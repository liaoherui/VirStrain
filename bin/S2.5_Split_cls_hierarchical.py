#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
from Bio import AlignIO
from Bio.Seq import Seq
from networkx.algorithms.components import connected_components


def estimate_shannon_entropy(dna_seq: str) -> tuple[float, dict[str, int]]:
    m = 0
    base_dict = {"a": "", "t": "", "g": "", "c": "", "A": "", "T": "", "G": "", "C": ""}
    bases_raw = dict(collections.Counter(dna_seq))

    # We do not consider dash or special characters when calculating entropy
    bases: dict[str, int] = {}
    for b, count in bases_raw.items():
        if b not in base_dict:
            continue
        bases[b] = count
        m += count

    if m == 0:
        return 0.0, bases_raw

    shannon_entropy_value = 0.0
    for base, n_i in bases.items():
        p_i = n_i / float(m)
        entropy_i = p_i * math.log(p_i)
        shannon_entropy_value += entropy_i

    return shannon_entropy_value * (-1), bases_raw


def estimate_shannon_entropy_spe(alignment, strain_set: dict[int, str], left_c: dict[int, str]):
    centropy: dict[int, float] = {}
    cbase_freq: dict[int, dict[str, int]] = {}

    for c in left_c:
        seq = []
        for i, n in enumerate(alignment[:, c]):
            if i not in strain_set:
                continue
            seq.append(n)

        entropy, base_count = estimate_shannon_entropy("".join(seq))
        cbase_freq[c] = base_count
        centropy[c] = 1000 if entropy == 0 else entropy

    return centropy, cbase_freq


def estimate_shannon_entropy_iterate(
    alignment,
    used_strains: dict[int, str],
    used_columns: dict[int, str],
    centropy: dict[int, float],
    strain_set: dict[int, str],
):
    new_centropy: dict[int, float] = {}
    new_cbase_freq: dict[int, dict[str, int]] = {}
    base_dict = {"a": "", "t": "", "g": "", "c": "", "A": "", "T": "", "G": "", "C": ""}

    for c in centropy:
        tem_bf: dict[str, int] = {}
        m = 0

        if c in used_columns:
            new_centropy[c] = 1000
            continue

        for i, s in enumerate(alignment[:, c]):
            if i not in strain_set:
                continue
            if i in used_strains:
                continue
            if s not in base_dict:
                continue

            tem_bf[s] = tem_bf.get(s, 0) + 1
            m += 1

        if m == 0:
            new_centropy[c] = 1000
            new_cbase_freq[c] = tem_bf
            continue

        shannon_entropy_value = 0.0
        for base, n_i in tem_bf.items():
            p_i = n_i / float(m)
            entropy_i = p_i * math.log(p_i)
            shannon_entropy_value += entropy_i

        new_centropy[c] = shannon_entropy_value * (-1)
        new_cbase_freq[c] = tem_bf

    return new_centropy, new_cbase_freq


def to_edges(items):
    it = iter(items)
    try:
        last = next(it)
    except StopIteration:
        return

    for current in it:
        yield last, current
        last = current


def to_graph(groups):
    graph = nx.Graph()
    for part in groups:
        graph.add_nodes_from(part)
        graph.add_edges_from(to_edges(part))
    return graph


def cls_with_entropy_hier(alignment, left_c: dict[int, str], strain_set: dict[int, str], k: int):
    used_strains: dict[int, str] = {}
    used_columns: dict[int, str] = {}
    selected_columns: dict[int, str] = {}
    base_dict = {"a": "", "t": "", "g": "", "c": "", "A": "", "T": "", "G": "", "C": ""}

    cqnum = len(left_c)
    snum = len(strain_set)
    sarr = sorted(strain_set.keys())

    centropy, cbase_freq = estimate_shannon_entropy_spe(alignment, strain_set, left_c)

    for c in centropy:
        if centropy[c] == 1000:
            used_columns[c] = ""

    ibreak = 0
    while True:
        if len(used_strains) == snum or len(used_columns) == cqnum:
            break

        print("Progress: C:", len(used_columns), "/", cqnum, "S:", len(used_strains), "/", snum)
        res = sorted(centropy.items(), key=lambda d: d[1])

        if res[0][1] == 1000 and ibreak == 0:
            return {}, {}, {}, []

        ibreak += 1
        if res[0][1] == 1000:
            break

        for r in res:
            if r[0] not in used_columns:
                current_c = r[0]
                break

        rawseq = ""
        sid2ni = {}
        for i, s in enumerate(sarr):
            rawseq += alignment[s, current_c]
            sid2ni[i] = s

        res2 = sorted(cbase_freq[current_c].items(), key=lambda d: d[1])
        round_strain = []

        for r in res2:
            if r[0] not in base_dict:
                continue
            if r[0] == res2[-1][0]:
                continue

            strain_id = [sid2ni[i] for i, letter in enumerate(rawseq) if letter == r[0]]
            round_strain.extend(strain_id)

            if centropy[current_c] != 1000:
                selected_columns[current_c] = ""

            cbase_freq[current_c][r[0]] = 0

        for s in round_strain:
            used_strains[s] = ""

        used_columns[current_c] = ""
        centropy, cbase_freq = estimate_shannon_entropy_iterate(
            alignment, used_strains, used_columns, centropy, strain_set
        )

    sub_info = {}
    kmr_sub = {}
    kmr_pos = {}
    strain_sub = {}
    sort_sub = []

    if len(selected_columns) == 0:
        return sub_info, kmr_sub, kmr_pos, sort_sub

    carr = sorted(selected_columns.keys())
    scseq = []

    for s in strain_set:
        cseq = "".join(alignment[s, c] for c in carr)
        scseq.append([cseq, f"{s}_cls"])

    graph = to_graph(scseq)
    count = 1

    for ele in connected_components(graph):
        ele = list(ele)
        pre = f"SubCls{count}_{len(ele) - 1}"
        sort_sub.append(pre)
        sub_info[pre] = []

        for e in ele:
            if "_cls" not in e:
                continue
            e = re.sub(r"_cls$", "", e)
            e = int(e)
            sub_info[pre].append(e)
            strain_sub[e] = pre

        count += 1

    check_dict = {"a": "A", "t": "T", "g": "G", "c": "C", "A": "A", "T": "T", "G": "G", "C": "C", "-": ""}
    dmap = {"a": "A", "t": "T", "g": "G", "c": "C", "A": "A", "T": "T", "G": "G", "C": "C"}
    half_k = (k - 1) // 2

    for c in carr:
        for s in strain_set:
            base = alignment[s, c].upper()
            left = str(alignment[s].seq[:c])[::-1]
            right = str(alignment[s].seq[c + 1 :])
            pos_base = f"{c}-{base}"

            if c - 12 < 0:
                block_seq = str(alignment[s].seq[0:c]) + base + str(alignment[s].seq[c + 1 : c + 13])
            else:
                block_seq = str(alignment[s].seq[c - 12 : c]) + base + str(alignment[s].seq[c + 1 : c + 13])

            if any(b not in check_dict for b in block_seq):
                continue

            lseq = ""
            rseq = ""

            for l in left:
                if len(lseq) == half_k:
                    break
                if l in dmap:
                    lseq += dmap[l]
                elif l != "-":
                    lseq += l.upper()

            for r in right:
                if len(rseq) == half_k:
                    break
                if r in dmap:
                    rseq += dmap[r]
                elif r != "-":
                    rseq += r.upper()

            lseq = lseq[::-1]
            kmr = lseq + base + rseq

            if len(kmr) < k:
                if len(lseq) < half_k:
                    rseq = ""
                    rl = k - 1 - len(lseq)
                    for r in right:
                        if len(rseq) == rl:
                            break
                        if r in dmap:
                            rseq += dmap[r]
                        elif r != "-":
                            rseq += r.upper()

                if len(rseq) < half_k:
                    lseq = ""
                    ll = k - 1 - len(rseq)
                    for l in left:
                        if len(lseq) == ll:
                            break
                        if l in dmap:
                            lseq += dmap[l]
                        elif l != "-":
                            lseq += l.upper()
                    lseq = lseq[::-1]

                kmr = lseq + base + rseq

            if kmr not in kmr_sub:
                kmr_sub[kmr] = {strain_sub[s]: ""}
            else:
                kmr_sub[kmr][strain_sub[s]] = ""

            if kmr not in kmr_pos:
                kmr_pos[kmr] = {pos_base: ""}
            else:
                kmr_pos[kmr][pos_base] = ""

    return sub_info, kmr_sub, kmr_pos, sort_sub


def parse_args():
    parser = argparse.ArgumentParser(description="Split clusters using entropy over MSA columns.")
    parser.add_argument("-i", "--input_msa", required=True, help="Input MSA FASTA file")
    parser.add_argument("-c", "--column_file", required=True, help="Column file")
    parser.add_argument("-r", "--cls_file", required=True, help="Cluster file")
    parser.add_argument("-k", "--kmer", type=int, default=25, help="K-mer size")
    parser.add_argument(
        "--dash_cutoff",
        type=float,
        default=0.01,
        help="Maximum allowed fraction of gaps in a column",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_msa = Path(args.input_msa)
    column_file = Path(args.column_file)
    cls_file = Path(args.cls_file)
    dash_cutoff = args.dash_cutoff
    k = args.kmer

    # Load used column info
    bused_c = {}
    with column_file.open("r", encoding="utf-8") as cfile:
        for raw_line in cfile:
            line = raw_line.strip()
            if not line:
                continue
            if "column" not in line:
                continue
            ele = line.split()
            bused_c[int(ele[1])] = ""

    # Load MSA
    alignment = AlignIO.read(str(input_msa), "fasta")
    sid = list(range(len(alignment[:, 0])))
    sname = []
    sname2seq = {}
    cnum = 0

    for record in alignment:
        name = record.id.strip()
        sname.append(">" + name)
        cnum = len(record.seq)
        sname2seq[name] = record.seq

    cid = list(range(cnum))
    centropy = {}
    cbase_freq = {}

    # Filter candidate columns
    for column in cid:
        seq = alignment[:, column]
        entropy, base_count = estimate_shannon_entropy(seq)

        if "-" not in base_count:
            if entropy != 0:
                centropy[column] = entropy
                cbase_freq[column] = base_count
        else:
            if base_count["-"] <= dash_cutoff * len(sid):
                if entropy != 0:
                    cbase_freq[column] = base_count
                    centropy[column] = entropy

    sid2name = dict(zip(sid, sname))
    sname2id = dict(zip(sname, sid))

    stat = 0
    left_c = {}
    for c in centropy:
        if c not in bused_c:
            left_c[c] = ""
            stat += 1

    print("Total columns:", cnum)
    print("Used column in first iterate:", len(bused_c))
    print("Left Column:", stat, "/", len(centropy))

    # Parse cluster file
    strain_info = {}
    cls_strain = {}
    cls_rep = {}
    strain_rep = {}

    with cls_file.open("r", encoding="utf-8") as fc:
        cls = None
        rep = None
        for raw_line in fc:
            line = raw_line.strip()
            if not line:
                continue

            if re.search("Cluster", line):
                cls = re.sub(">", "", line)
                cls_strain[cls] = {}
            else:
                ele = line.split("\t")
                if "*" in line:
                    rep = ele[0]
                    strain_rep[rep] = rep
                    strain_info[rep] = cls
                    cls_strain[cls][rep] = int(ele[-1])
                    cls_rep[cls] = rep
                else:
                    strain_rep[ele[0]] = rep
                    strain_info[ele[0]] = cls
                    cls_strain[cls][ele[0]] = int(ele[-1])

    raw_cls_size = []
    new_cls_size = []

    with open("Strain_cls_info.txt", "w", encoding="utf-8") as o1, \
         open("Rebuild_cls.clstr", "w", encoding="utf-8") as o2, \
         open("SubCls_kmer.txt", "w", encoding="utf-8") as o3:

        for c in cls_strain:
            raw_cls_size.append(len(cls_strain[c]))

            if len(cls_strain[c]) == 1:
                new_cls_size.append(len(cls_strain[c]))
                for s in cls_strain[c]:
                    o1.write(
                        f"{s}\t{strain_rep[s]}\t{strain_info[s]}\t{strain_info[s]}\t{sname2id[s]}\t1\n"
                    )
                    o2.write(f">{c}\n")
                    o2.write(f"\t{s}\t*\t{sname2id[s]}\n")
            else:
                strain_set = {cls_strain[c][s]: "" for s in cls_strain[c]}
                sub_info, kmr_sub, kmr_pos, sort_sub = cls_with_entropy_hier(
                    alignment, left_c, strain_set, k
                )
                print(sub_info, sort_sub)

                if len(sub_info) == 0 or len(sort_sub) == 1:
                    new_cls_size.append(len(cls_strain[c]))
                    o2.write(f">{c}\n")
                    for s in cls_strain[c]:
                        o1.write(
                            f"{s}\t{strain_rep[s]}\t{strain_info[s]}\t{strain_info[s]}\t{sname2id[s]}\t{len(cls_strain[c])}\n"
                        )
                        if s == cls_rep[c]:
                            o2.write(f"\t{s}\t*\t{sname2id[s]}\n")
                        else:
                            o2.write(f"\t{s}\t{sname2id[s]}\n")
                else:
                    o2.write(f">{c}\n")
                    for sub in sort_sub:
                        o2.write(f"\t>{sub}\n")
                        new_cls_size.append(len(sub_info[sub]))
                        for ss in sub_info[sub]:
                            ss_name = sid2name[ss]
                            o1.write(
                                f"{ss_name}\t{strain_rep[ss_name]}\t{strain_info[ss_name]}\t{sub}\t{sname2id[ss_name]}\t{len(sub_info[sub])}\n"
                            )
                            if ss_name == cls_rep[c]:
                                o2.write(f"\t\t{ss_name}\t*\t{sname2id[ss_name]}\n")
                            else:
                                o2.write(f"\t\t{ss_name}\t{sname2id[ss_name]}\n")

                    for kl in kmr_sub:
                        info1 = ",".join(kmr_pos[kl].keys())
                        info2 = ",".join(kmr_sub[kl].keys())
                        o3.write(f"{kl}\t{c}\t{info1}\t{info2}\n")
                        k2 = str(Seq(kl).reverse_complement())
                        o3.write(f"{k2}\t{c}\t{info1}\t{info2}\n")

    plt.hist(raw_cls_size, bins=100)
    plt.xlabel("Cluster size")
    plt.ylabel("Cluster Number")
    plt.savefig("Before_split.png")
    plt.close()

    plt.figure()
    plt.hist(new_cls_size, bins=100)
    plt.xlabel("Cluster size")
    plt.ylabel("Cluster Number")
    plt.savefig("After_split.png")
    plt.close()


if __name__ == "__main__":
    main()
