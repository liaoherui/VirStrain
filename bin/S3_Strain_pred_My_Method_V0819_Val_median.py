#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

BASE_ORDER = ["A", "T", "G", "C"]
BASE_P = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="VirStrain strain prediction script"
    )
    parser.add_argument("-i", "--input_reads", required=True, help="SE reads or read1")
    parser.add_argument("-p", "--input_reads2", default="", help="PE reads read2")
    parser.add_argument("-k", "--kmer_size", type=int, default=25, help="k-mer size")
    parser.add_argument("-m", "--matrix_file", required=True, help="SNP-position matrix file")
    parser.add_argument("-s", "--snp_kmer_file", required=True, help="kmer -> snp_pos file")
    parser.add_argument("-f", "--snp_kmer_fa", required=True, help="kmer fasta")
    parser.add_argument("-c", "--cls_file", required=True, help="cluster info file")
    parser.add_argument("-b", "--sub_kmr_file", required=True, help="subcluster kmer file")
    parser.add_argument("-o", "--output", required=True, help="output text report")
    return parser.parse_args()


def run_cmd(cmd: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path is not None:
        with stdout_path.open("w", encoding="utf-8") as handle:
            subprocess.run(cmd, check=True, stdout=handle)
    else:
        subprocess.run(cmd, check=True)


def load_kmer_to_snp(snp_kmr_file: Path) -> tuple[dict[str, str], dict[str, int]]:
    dkps: dict[str, str] = {}
    dpsc: dict[str, int] = {}
    with snp_kmr_file.open("r", encoding="utf-8") as f1:
        for raw_line in f1:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split("\t")
            dkps[ele[0]] = ""
            ps = ele[1].split(",")
            for e in ps:
                dkps[ele[0]] = e
                dpsc[e] = 0
    return dkps, dpsc


def load_matrix_header(matrix_file: Path) -> tuple[list[str], list[str]]:
    with matrix_file.open("r", encoding="utf-8") as f3:
        header = f3.readline().strip()
    pos_snp = header.split("\t")
    carr: list[str] = []
    for p in pos_snp:
        c = p.split("-")[0]
        if c not in carr:
            carr.append(c)
    return pos_snp, carr


def count_kmers(
    file_dir: Path,
    snp_kmr_fa: Path,
    read_1: Path,
    read_2: Path | None,
    jf_file: Path,
    fa_file: Path,
    k: int,
) -> None:
    jellyfish_bin = file_dir / "jellyfish-linux"
    count_cmd = [
        str(jellyfish_bin),
        "count",
        "-m",
        str(k),
        "-s",
        "100M",
        "-t",
        "8",
        "--if",
        str(snp_kmr_fa),
        "-o",
        str(jf_file),
        str(read_1),
    ]
    if read_2 is not None:
        count_cmd.append(str(read_2))
    run_cmd(count_cmd)
    run_cmd([str(jellyfish_bin), "dump", "-c", str(jf_file)], stdout_path=fa_file)


def update_dpsc_from_dump(fa_file: Path, dkps: dict[str, str], dpsc: dict[str, int]) -> None:
    with fa_file.open("r", encoding="utf-8") as fnew:
        for raw_line in fnew:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split()
            if len(ele) < 2:
                continue
            kmer = ele[0]
            if kmer in dkps:
                dpsc[dkps[kmer]] += int(ele[1])


def safe_percentile_bounds(arr: np.ndarray, low: float, high: float) -> tuple[float, float]:
    if arr.size == 0:
        return 0.0, 0.0
    return tuple(np.percentile(arr, [low, high]))  # type: ignore[return-value]


def build_freq_arr(pos_snp: list[str], dpsc: dict[str, int]) -> np.ndarray:
    freq_arr = [dpsc.get(p, 0) for p in pos_snp]
    return np.array(freq_arr, dtype=float)


def compute_filtered_freq(freq_arr: np.ndarray) -> tuple[np.ndarray, float]:
    min_depth_percentile = 10
    max_depth_percentile = 90
    min_depth_absolute = 2
    min_depth_rate = 0.05

    keep = freq_arr != 0
    check_arr = freq_arr[keep]

    if check_arr.size == 0:
        return freq_arr.copy(), float(min_depth_absolute)

    min_depth, max_depth = safe_percentile_bounds(
        check_arr, min_depth_percentile, max_depth_percentile
    )
    keep2 = np.logical_and.reduce((check_arr >= min_depth, check_arr <= max_depth))
    check_arr2 = check_arr[keep2]

    if check_arr2.size == 0:
        min_depth_adf = float(min_depth_absolute)
    else:
        min_depth_adf = min_depth_rate * float(np.mean(check_arr2))
        if min_depth_adf < min_depth_absolute:
            min_depth_adf = float(min_depth_absolute)

    new_freq = freq_arr.copy()
    new_freq[new_freq <= min_depth_adf] = 0
    return new_freq, min_depth_adf


def load_strain_matrix(
    matrix_file: Path,
    freq_arr: np.ndarray,
    weighted_freq_arr: np.ndarray,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, float],
    dict[str, str],
    dict[str, float],
    dict[str, dict[int, float]],
    dict[str, dict[int, int]],
    np.ndarray,
]:
    ds_pos: dict[str, np.ndarray] = {}
    ds_freq: dict[str, np.ndarray] = {}
    dmap_rate: dict[str, float] = {}
    ds_num: dict[str, str] = {}
    dmr: dict[str, float] = {}
    dscf: dict[str, dict[int, float]] = {}
    dscl: dict[str, dict[int, int]] = {}
    all_ps: list[np.ndarray] = []

    with matrix_file.open("r", encoding="utf-8") as f3:
        _ = f3.readline()
        for raw_line in f3:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split("\t")
            tem = np.array([int(e) for e in ele[1:]], dtype=int)
            all_ps.append(tem)
            dscf[ele[0]] = {}
            dscl[ele[0]] = {}

            nt = freq_arr * tem
            raw_c = len(tem[tem == 1])
            map_c = len(nt[nt > 0])
            map_rate = float(np.sum(tem * weighted_freq_arr)) if weighted_freq_arr.size > 0 else 0.0

            dmap_rate[ele[0]] = map_rate
            ds_num[ele[0]] = f"{map_c}/{raw_c}" if raw_c > 0 else "0/0"
            dmr[ele[0]] = float(map_c) / float(raw_c) if raw_c > 0 else 0.0
            ds_pos[ele[0]] = tem
            ds_freq[ele[0]] = nt

    return ds_pos, ds_freq, dmap_rate, ds_num, dmr, dscf, dscl, np.array(all_ps)


def compute_top_map_strains(
    dmr: dict[str, float],
    dmap_rate: dict[str, float],
) -> tuple[list[str], list[tuple[str, float]]]:
    res = sorted(dmap_rate.items(), key=lambda d: d[1], reverse=True)
    top10_score_s = res[:10]
    if not res:
        return [], top10_score_s
    top_map_strain: list[str] = []
    for r in res:
        if r[1] == res[0][1]:
            top_map_strain.append(r[0])
        else:
            break
    return top_map_strain, top10_score_s


def resolve_top_map_ties(
    top_map_strain: list[str],
    ds_pos: dict[str, np.ndarray],
    freq_arr: np.ndarray,
    pos_snp: list[str],
    dmap_rate: dict[str, float],
    dmr: dict[str, float],
) -> list[str]:
    if not top_map_strain:
        return top_map_strain

    snp_arr = [ds_pos[s] for s in top_map_strain]
    pre_freq_arr = np.array([], dtype=float)

    for s in top_map_strain:
        pre_pa = ds_pos[s] * (-1)
        pre_pa = np.array(pre_pa, dtype=float)
        pre_pa[pre_pa == 0] = 1
        pre_freq_arr = freq_arr * pre_pa
        pre_freq_arr[pre_freq_arr < 0] = 0

    keep = pre_freq_arr != 0
    pre_freq_arr = pre_freq_arr[keep]
    pre_pos_snp = np.array(pos_snp)[keep]
    pre_ds_pos = {s: ds_pos[s][keep] for s in ds_pos}

    if pre_freq_arr.size == 0:
        return top_map_strain

    pre_wf_arr = pre_freq_arr / np.sum(pre_freq_arr)
    strain_num: dict[str, str] = {}
    sn = 0

    while True:
        if len(pre_freq_arr) == 0:
            break
        smr: dict[str, float] = {}
        for r in ds_pos:
            if r in top_map_strain:
                continue
            tt = pre_ds_pos[r]
            nt = tt * pre_wf_arr
            mr = float(np.sum(nt))
            smr[r] = mr

        if not smr:
            break

        res = sorted(smr.items(), key=lambda d: d[1], reverse=True)
        ts: list[str] = []
        for r in res:
            if r[1] == res[0][1]:
                ts.append(r[0])

        if not ts:
            break

        if len(ts) > 1:
            rmr = {s: dmap_rate[s] for s in ts}
            res2 = sorted(rmr.items(), key=lambda d: d[1], reverse=True)
            strain_num[res2[0][0]] = ""
            chosen = res2[0][0]
        else:
            strain_num[ts[0]] = ""
            chosen = ts[0]

        vm1 = len(pre_freq_arr)
        pre_pa = pre_ds_pos[chosen] * (-1)
        pre_pa[pre_pa == 0] = 1
        pre_freq_arr = pre_freq_arr * pre_pa
        pre_freq_arr[pre_freq_arr < 0] = 0
        keep = pre_freq_arr != 0
        pre_freq_arr = pre_freq_arr[keep]

        if np.sum(pre_freq_arr) != 0:
            pre_wf_arr = pre_freq_arr / np.sum(pre_freq_arr)
        pre_pos_snp = pre_pos_snp[keep]
        for s in pre_ds_pos:
            pre_ds_pos[s] = pre_ds_pos[s][keep]

        vm = vm1 - len(pre_freq_arr)
        if vm > 1:
            sn += 1

    if sn > 1:
        for s in top_map_strain:
            strain_num[s] = ""
        sscore: dict[str, float] = {}
        sna = np.array([ds_pos[s] for s in strain_num], dtype=float)
        if sna.size == 0:
            return top_map_strain
        ssum = sna.sum(axis=0)
        ssum[ssum == 0] = 1
        for s in strain_num:
            snt = ds_pos[s] / ssum
            ns = dmap_rate[s] * snt
            sscore[s] = float(ns.sum(axis=0))
        res = sorted(sscore.items(), key=lambda d: d[1], reverse=True)
        max_map = sorted(dmr.items(), key=lambda d: d[1], reverse=True)[0][1]
        tem_map_strain: list[str] = []
        for r in res:
            if dmr[r[0]] != max_map:
                continue
            tem_map_strain.append(r[0])
            break
        if len(tem_map_strain) > 0:
            top_map_strain = tem_map_strain

    return top_map_strain


def collect_unique_snps(
    top_map_strain: list[str],
    ds_pos: dict[str, np.ndarray],
    pos_snp: list[str],
    pos_freq_map: dict[str, float],
    min_depth_absolute: int = 2,
) -> tuple[dict[str, dict[str, float]], dict[str, int]]:
    if not top_map_strain:
        return {}, {}

    snp_arr = np.array([ds_pos[s] for s in top_map_strain], dtype=int)
    pos_sum = snp_arr.sum(axis=0)
    pos_sum[pos_sum > 1] = 0

    strain_unique: dict[str, dict[str, float]] = {}
    strain_unique_count: dict[str, int] = {}

    for i, p in enumerate(pos_sum):
        column = pos_snp[i]
        if p == 1:
            if pos_freq_map[column] <= min_depth_absolute:
                continue
            window = snp_arr[:, i]
            for i2, w in enumerate(window):
                if w == 1:
                    strain = top_map_strain[i2]
                    if strain not in strain_unique:
                        strain_unique[strain] = {column: pos_freq_map[column]}
                        strain_unique_count[strain] = 1
                    else:
                        strain_unique[strain][column] = pos_freq_map[column]
                        strain_unique_count[strain] += 1

    return strain_unique, strain_unique_count


def compute_depth_stats(
    mp_strain: list[str],
    ds_freq: dict[str, np.ndarray],
    ds_pos: dict[str, np.ndarray],
    freq_arr: np.ndarray,
    min_depth_percentile: int = 10,
    max_depth_percentile: int = 90,
) -> tuple[dict[str, float], dict[str, float], np.ndarray]:
    ds_avgd: dict[str, float] = {}
    ds_medd: dict[str, float] = {}
    new_freq_arr = freq_arr.copy()

    for m in mp_strain:
        keep = ds_freq[m] != 0
        vals = ds_freq[m][keep]
        if vals.size == 0:
            ds_avgd[m] = 0.0
            ds_medd[m] = 0.0
        else:
            min_depth, max_depth = safe_percentile_bounds(vals, min_depth_percentile, max_depth_percentile)
            keep2 = np.logical_and.reduce((ds_freq[m] >= min_depth, ds_freq[m] <= max_depth))
            filtered = ds_freq[m][keep2]
            if filtered.size == 0:
                ds_avgd[m] = 0.0
                ds_medd[m] = 0.0
            else:
                ds_avgd[m] = float(np.mean(filtered))
                ds_medd[m] = float(np.median(filtered))

        pos_arr = ds_pos[m] * (-1)
        pos_arr = np.array(pos_arr, dtype=float)
        pos_arr[pos_arr == 0] = 1
        new_freq_arr = new_freq_arr * pos_arr
        new_freq_arr[new_freq_arr < 0] = 0

    return ds_avgd, ds_medd, new_freq_arr


def compute_other_possible_strains(
    mp_strain: list[str],
    ds_pos: dict[str, np.ndarray],
    ds_freq: dict[str, np.ndarray],
    dmap_rate: dict[str, float],
    freq_arr: np.ndarray,
    pos_snp: list[str],
    ds_avgd: dict[str, float],
    ds_medd: dict[str, float],
    min_depth_percentile: int = 10,
    max_depth_percentile: int = 90,
) -> tuple[dict[str, list[str]], list[str], dict[str, int], dict[str, float], dict[str, float]]:
    keep = freq_arr != 0
    left_freq_arr = freq_arr[keep]
    pos_snp_np = np.array(pos_snp)
    left_pos_snp = pos_snp_np[keep]
    left_ds_pos = {s: ds_pos[s][keep] for s in ds_pos}
    left_ps_freq_map = dict(zip(left_pos_snp, left_freq_arr))

    if left_freq_arr.size == 0:
        return {}, [], {}, ds_avgd, ds_medd

    resl = sorted(left_ps_freq_map.items(), key=lambda d: d[1], reverse=True)
    os_strain: dict[str, list[str]] = {}
    os_arr: list[str] = []
    left_weighted_freq_arr = left_freq_arr / np.sum(left_freq_arr)
    vmap: dict[str, int] = {}

    max_iter_times = len(left_freq_arr)
    for _ in range(max_iter_times):
        if len(left_freq_arr) == 0:
            break

        strain_map_rate: dict[str, float] = {}
        for r in ds_pos:
            if r in mp_strain:
                continue
            tem = left_ds_pos[r]
            nt = tem * left_weighted_freq_arr
            map_rate = float(np.sum(nt))
            strain_map_rate[r] = map_rate

        if not strain_map_rate:
            break

        res = sorted(strain_map_rate.items(), key=lambda d: d[1], reverse=True)
        top_s: list[str] = []

        for r in res:
            if r[1] == res[0][1]:
                top_s.append(r[0])
                nt = left_ds_pos[r[0]] * left_freq_arr
                keep_nt = nt != 0
                nt = nt[keep_nt]
                if nt.size == 0:
                    ds_avgd[r[0]] = 0.0
                    ds_medd[r[0]] = 0.0
                else:
                    min_depth, max_depth = safe_percentile_bounds(
                        nt, min_depth_percentile, max_depth_percentile
                    )
                    keep_nt2 = np.logical_and.reduce((nt >= min_depth, nt <= max_depth))
                    if len(nt[keep_nt2]) != 0:
                        nt = nt[keep_nt2]
                    ds_avgd[r[0]] = float(np.mean(nt)) if nt.size else 0.0
                    ds_medd[r[0]] = float(np.median(nt)) if nt.size else 0.0

        if len(top_s) > 1:
            rank_map_rate = {s: dmap_rate[s] for s in top_s}
            res2 = sorted(rank_map_rate.items(), key=lambda d: d[1], reverse=True)
            top_s = [r[0] for r in res2]

        pre = [f"{r[0]}:{r[1]}" for r in resl]
        top_pos_snp = ",".join(pre) + "\t" + str(len(top_s))
        os_arr.append(top_pos_snp)
        os_strain[top_pos_snp] = []

        for s in top_s:
            os_strain[top_pos_snp].append(s)
            pos_arr = left_ds_pos[s] * (-1)
            pos_arr[pos_arr == 0] = 1
            left_freq_arr = left_freq_arr * pos_arr
            left_freq_arr[left_freq_arr < 0] = 0

        keep = left_freq_arr != 0
        valid_map = len(left_freq_arr)
        left_freq_arr = left_freq_arr[keep]
        valid_map = valid_map - len(left_freq_arr)
        vmap[top_pos_snp] = valid_map

        if np.sum(left_freq_arr) != 0:
            left_weighted_freq_arr = left_freq_arr / np.sum(left_freq_arr)

        left_pos_snp = left_pos_snp[keep]
        left_ps_freq_map = dict(zip(left_pos_snp, left_freq_arr))
        for s in left_ds_pos:
            left_ds_pos[s] = left_ds_pos[s][keep]
        resl = sorted(left_ps_freq_map.items(), key=lambda d: d[1], reverse=True)

    return os_strain, os_arr, vmap, ds_avgd, ds_medd


def load_cluster_info(
    cls_file: Path,
    all_s: list[str],
    sub_kmr_file: Path,
    read_1: Path,
    read_2: Path | None,
    file_dir: Path,
    min_depth_adf: float,
) -> tuple[dict[str, str], dict[str, str]]:
    s2cls: dict[str, str] = {}
    s2sub: dict[str, str] = {}
    candidate_cls: dict[str, str] = {}

    with cls_file.open("r", encoding="utf-8") as fcls:
        for raw_line in fcls:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split("\t")
            if ele[0] not in all_s:
                s2cls[ele[0]] = ele[2]
                continue
            if ele[2] == ele[3]:
                s2cls[ele[0]] = ele[2]
                s2sub[ele[0]] = "NA"
            else:
                s2cls[ele[0]] = ele[2]
                candidate_cls[ele[2]] = ""

    if len(candidate_cls) == 0:
        return s2cls, s2sub

    ksub: dict[str, dict[str, dict[str, str]]] = {}
    cls_sub: dict[str, dict[str, int]] = {}

    with sub_kmr_file.open("r", encoding="utf-8") as fsk:
        for raw_line in fsk:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split("\t")
            kmr = ele[0]
            if ele[1] in candidate_cls:
                if kmr not in ksub:
                    ksub[kmr] = {ele[1]: {}}
                if ele[1] not in ksub[kmr]:
                    ksub[kmr][ele[1]] = {}
                if ele[1] not in cls_sub:
                    cls_sub[ele[1]] = {}
                sub = ele[-1].split(",")
                for s in sub:
                    ksub[kmr][ele[1]][s] = ""
                    cls_sub[ele[1]][s] = 0

    temp_sub_fa = Path("Tem_Vs2Sub.fa")
    with temp_sub_fa.open("w", encoding="utf-8") as ok:
        co = 1
        for kmr in ksub:
            ok.write(f">{co}\n{kmr}\n")
            co += 1

    jellyfish_bin = file_dir / "jellyfish-linux"
    jf2 = Path("Tem_VS2.jf")
    fa2 = Path("Tem_Vs2.fa")

    count_cmd = [
        str(jellyfish_bin),
        "count",
        "-m",
        "25",
        "-s",
        "100M",
        "-t",
        "8",
        "--if",
        str(temp_sub_fa),
        "-o",
        str(jf2),
        str(read_1),
    ]
    if read_2 is not None:
        count_cmd.append(str(read_2))
    run_cmd(count_cmd)
    run_cmd([str(jellyfish_bin), "dump", "-c", str(jf2)], stdout_path=fa2)

    with fa2.open("r", encoding="utf-8") as ft2:
        for raw_line in ft2:
            line = raw_line.strip()
            if not line:
                continue
            ele = line.split()
            if len(ele) < 2:
                continue
            if int(ele[1]) >= min_depth_adf and ele[0] in ksub:
                for c in ksub[ele[0]]:
                    for c2 in ksub[ele[0]][c]:
                        cls_sub[c][c2] += int(ele[1])

    for s in all_s:
        if s in s2sub:
            continue
        if s2cls[s] not in cls_sub:
            s2sub[s] = "NA"
        else:
            res = sorted(cls_sub[s2cls[s]].items(), key=lambda d: d[1], reverse=True)
            if len(res) < 2:
                s2sub[s] = res[0][0] if res else "NA"
            elif res[0][1] == res[1][1]:
                s2sub[s] = "NA"
            else:
                s2sub[s] = res[0][0]

    return s2cls, s2sub


def write_text_report(
    output_file: Path,
    mp_strain: list[str],
    os_arr: list[str],
    os_strain: dict[str, list[str]],
    vmap: dict[str, int],
    s2cls: dict[str, str],
    s2sub: dict[str, str],
    dmap_rate: dict[str, float],
    ds_num: dict[str, str],
    ds_avgd: dict[str, float],
    ds_medd: dict[str, float],
    strain_unique: dict[str, dict[str, float]],
    dmr: dict[str, float],
    all_s: list[str],
    top10_score_s: list[tuple[str, float]],
) -> None:
    with output_file.open("w", encoding="utf-8") as o:
        o.write("\t\tStrain_ID\tCls_info\tSubCls_info\tMap_Score\tValid_Map_Rate\tTotal_Map_Rate\tStrain_Depth\tStrain_info\tUnique_SNP\n")
        o.write(">>Most possible strains:\n")

        for s in mp_strain:
            if s in strain_unique:
                o.write(
                    f"\t\t{s}\t{s2cls.get(s, 'NA')}\t{s2sub.get(s, 'NA')}\t{dmap_rate.get(s, 0)}\t"
                    f"{ds_num.get(s, '0/0')}\t{ds_num.get(s, '0/0')}\t"
                    f"{ds_avgd.get(s, 0)},{ds_medd.get(s, 0)}\t\t\t{strain_unique[s]}\n"
                )
            else:
                o.write(
                    f"\t\t{s}\t{s2cls.get(s, 'NA')}\t{s2sub.get(s, 'NA')}\t{dmap_rate.get(s, 0)}\t"
                    f"{ds_num.get(s, '0/0')}\t{ds_num.get(s, '0/0')}\t"
                    f"{ds_avgd.get(s, 0)},{ds_medd.get(s, 0)}\t\t\tNA\n"
                )

        o.write(">>Other possible strains:\n")
        if len(os_arr) > 0:
            for s in os_arr:
                ele = s.split("\t")
                if vmap[s] == 1:
                    o.write(f"\t>>(Could be FP){ele[0]},Genome_num: {ele[1]}\n")
                else:
                    o.write(f"\t>>{ele[0]},Genome_num: {ele[1]}\n")
                for n in os_strain[s]:
                    a = ds_num[n].split("/")[-1]
                    vm = f"{vmap[s]}/{a}"
                    if n in s2sub:
                        o.write(
                            f"\t\t{n}\t{s2cls.get(n, 'NA')}\t{s2sub.get(n, 'NA')}\t{dmap_rate.get(n, 0)}\t"
                            f"{vm}\t{ds_num.get(n, '0/0')}\t{ds_avgd.get(n, 0)},{ds_medd.get(n, 0)}\t\t\tNot_record\n"
                        )
                    else:
                        o.write(
                            f"\t\t{n}\t{s2cls.get(n, 'NA')}\tNot_record\t{dmap_rate.get(n, 0)}\t"
                            f"{vm}\t{ds_num.get(n, '0/0')}\t{ds_avgd.get(n, 0)},{ds_medd.get(n, 0)}\t\t\tNot_record\n"
                        )
        else:
            o.write("\tCan not detect other strains.\n")

        res = sorted(dmr.items(), key=lambda d: d[1], reverse=True)
        o.write("\n>>Highest_Map_Strains (Could be FP):\n")
        final: dict[str, float] = {}
        if res:
            for r in res:
                if r[1] == res[0][1]:
                    if r[0] not in all_s:
                        final[r[0]] = dmap_rate[r[0]]
        if len(final) != 0:
            res2 = sorted(final.items(), key=lambda d: d[1], reverse=True)
            for s in res2:
                if s[0] in s2sub:
                    o.write(
                        f"\t\t{s[0]}\t{s2cls.get(s[0], 'NA')}\t{s2sub.get(s[0], 'NA')}\t{dmap_rate.get(s[0], 0)}\t"
                        f"{ds_num.get(s[0], '0/0')}\t{ds_num.get(s[0], '0/0')}\tNA\tNA\n"
                    )
                else:
                    o.write(
                        f"\t\t{s[0]}\t{s2cls.get(s[0], 'NA')}\tNot_record\t{dmap_rate.get(s[0], 0)}\t"
                        f"{ds_num.get(s[0], '0/0')}\t{ds_num.get(s[0], '0/0')}\tNA\tNA\n"
                    )

        o.write(">>Top10_Score_Strains:\n")
        for t in top10_score_s:
            if t[0] in s2sub:
                o.write(
                    f"\t\t{t[0]}\t{s2cls.get(t[0], 'NA')}\t{s2sub.get(t[0], 'NA')}\t{t[1]}\t"
                    f"{ds_num.get(t[0], '0/0')}\t{ds_num.get(t[0], '0/0')}\tNA\tNA\n"
                )
            else:
                o.write(
                    f"\t\t{t[0]}\t{s2cls.get(t[0], 'NA')}\tNot_record\t{t[1]}\t"
                    f"{ds_num.get(t[0], '0/0')}\t{ds_num.get(t[0], '0/0')}\tNA\tNA\n"
                )


def populate_visualization_dicts(
    ds_freq: dict[str, np.ndarray],
    ds_pos: dict[str, np.ndarray],
    vs_sd: list[str],
    vs_so: list[str],
    pos_snp: np.ndarray,
    pos_label: dict[str, int],
    dscf: dict[str, dict[int, float]],
    dscl: dict[str, dict[int, int]],
) -> None:
    for s in ds_freq:
        check = 0
        if s in vs_sd:
            check = 1
        if s in vs_so:
            check = 2
        if check == 0:
            continue

        for i, c in enumerate(ds_freq[s]):
            if c == 0:
                continue
            current_ps = pos_snp[i]
            column = int(current_ps.split("-")[0])
            if column not in dscf[s]:
                dscf[s][column] = float(c)

        for i2, c in enumerate(ds_pos[s]):
            if c == 0:
                continue
            current_ps = pos_snp[i2]
            column = int(current_ps.split("-")[0])
            if column not in dscl[s]:
                dscl[s][column] = 0
            if c == 1:
                dscl[s][column] = pos_label[current_ps]


def write_depth_csvs(
    vs_sd: list[str],
    vs_so: list[str],
    pos_snp: np.ndarray,
    dscf: dict[str, dict[int, float]],
    dscl: dict[str, dict[int, int]],
) -> None:
    carr: list[int] = []
    for c in pos_snp:
        col = c.split("-")[0]
        if int(col) not in carr:
            carr.append(int(col))

    with Path("Mps_ps_depth_SCOV2.csv").open("w", encoding="utf-8") as ov1:
        ov1.write("ID,Column_ID")
        for s in vs_sd:
            ov1.write(f",{s}_Freq,{s}_LNum")
        ov1.write("\n")

        i = 1
        for c in carr:
            ov1.write(f"{i},{c}")
            for s in vs_sd:
                fval = dscf.get(s, {}).get(c, 0)
                lval = dscl.get(s, {}).get(c, 0)
                ov1.write(f",{fval},{lval}")
            ov1.write("\n")
            i += 1

    with Path("Ops_ps_depth_SCOV2.csv").open("w", encoding="utf-8") as ov2:
        ov2.write("ID,Column_ID")
        if len(vs_so) == 0:
            ov2.write(",None,None\n")
        else:
            for s in vs_so:
                ov2.write(f",{s}_Freq,{s}_LNum")
            ov2.write("\n")
            i = 1
            for c in carr:
                ov2.write(f"{i},{c}")
                for s in vs_so:
                    if s not in dscf or s not in dscl:
                        print(f"Warning: {s} not in final dict!")
                        continue
                    fval = dscf.get(s, {}).get(c, 0)
                    lval = dscl.get(s, {}).get(c, 0)
                    ov2.write(f",{fval},{lval}")
                ov2.write("\n")
                i += 1


def cleanup_temp_files() -> None:
    for pattern in ["Tem_Vs*", "Tem_VS*"]:
        for p in Path(".").glob(pattern):
            try:
                p.unlink()
            except IsADirectoryError:
                pass
            except FileNotFoundError:
                pass


def main() -> int:
    args = parse_args()

    read_1 = Path(args.input_reads).resolve()
    read_2 = Path(args.input_reads2).resolve() if args.input_reads2 else None
    snp_kmr_file = Path(args.snp_kmer_file).resolve()
    snp_kmr_fa = Path(args.snp_kmer_fa).resolve()
    matrix_file = Path(args.matrix_file).resolve()
    cls_file = Path(args.cls_file).resolve()
    sub_kmr_file = Path(args.sub_kmr_file).resolve()
    out_dir = Path(args.output).resolve()
    k = args.kmer_size

    file_dir = Path(__file__).resolve().parent

    dkps, dpsc = load_kmer_to_snp(snp_kmr_file)
    pos_snp, _ = load_matrix_header(matrix_file)

    jf_file = Path("Tem_VS.jf")
    fa_file = Path("Tem_Vs.fa")

    try:
        count_kmers(file_dir, snp_kmr_fa, read_1, read_2, jf_file, fa_file, k)
        update_dpsc_from_dump(fa_file, dkps, dpsc)

        freq_arr = build_freq_arr(pos_snp, dpsc)
        filtered_freq_arr, min_depth_adf = compute_filtered_freq(freq_arr)

        if np.sum(filtered_freq_arr) == 0:
            weighted_freq_arr = np.zeros_like(filtered_freq_arr, dtype=float)
        else:
            weighted_freq_arr = filtered_freq_arr / np.sum(filtered_freq_arr)

        pos_freq_map = dict(zip(pos_snp, filtered_freq_arr))

        (
            ds_pos,
            ds_freq,
            dmap_rate,
            ds_num,
            dmr,
            dscf,
            dscl,
            all_ps,
        ) = load_strain_matrix(matrix_file, filtered_freq_arr, weighted_freq_arr)

        all_sum = np.sum(all_ps, axis=0) if all_ps.size else np.array([0] * len(pos_snp))
        pos_label = dict(zip(pos_snp, list(all_sum)))

        top_map_strain, top10_score_s = compute_top_map_strains(dmr, dmap_rate)
        top_map_strain = resolve_top_map_ties(
            top_map_strain, ds_pos, filtered_freq_arr, pos_snp, dmap_rate, dmr
        )

        strain_unique, _strain_unique_count = collect_unique_snps(
            top_map_strain, ds_pos, pos_snp, pos_freq_map
        )

        mp_strain: list[str] = []
        if len(strain_unique) != 0:
            for s in strain_unique:
                mp_strain.append(s)
        else:
            for r in top_map_strain:
                mp_strain.append(r)

        ds_avgd, ds_medd, left_after_mp = compute_depth_stats(
            mp_strain, ds_freq, ds_pos, filtered_freq_arr
        )

        (
            os_strain,
            os_arr,
            vmap,
            ds_avgd,
            ds_medd,
        ) = compute_other_possible_strains(
            mp_strain,
            ds_pos,
            ds_freq,
            dmap_rate,
            left_after_mp,
            pos_snp,
            ds_avgd,
            ds_medd,
        )

        all_s: list[str] = []
        vs_sd: list[str] = []
        vs_so: list[str] = []

        for s in mp_strain:
            all_s.append(s)
            vs_sd.append(s)

        if len(os_arr) > 0:
            for s in os_arr:
                all_s.append(os_strain[s][0])
                if vmap[s] > 1:
                    vs_so.append(os_strain[s][0])

        s2cls, s2sub = load_cluster_info(
            cls_file, all_s, sub_kmr_file, read_1, read_2, file_dir, min_depth_adf
        )

        write_text_report(
            out_dir,
            mp_strain,
            os_arr,
            os_strain,
            vmap,
            s2cls,
            s2sub,
            dmap_rate,
            ds_num,
            ds_avgd,
            ds_medd,
            strain_unique,
            dmr,
            all_s,
            top10_score_s,
        )

        print("Txt report is done. Now will generate CSV report!")

        vs_so = vs_so[:5]
        populate_visualization_dicts(
            ds_freq,
            ds_pos,
            vs_sd,
            vs_so,
            np.array(pos_snp),
            pos_label,
            dscf,
            dscl,
        )
        write_depth_csvs(vs_sd, vs_so, np.array(pos_snp), dscf, dscl)

    finally:
        cleanup_temp_files()

    return 0


if __name__ == "__main__":
    sys.exit(main())
