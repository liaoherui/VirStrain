#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from bin import prescan

__author__ = "Liao Herui"

USAGE = "VirStrain - An RNA virus strain-level identification tool for short reads."


def run_cmd(cmd: list[str], cwd: Path | None = None) -> None:
    """Run a command and raise an exception if it fails."""
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def split_contig(
    ingenome: Path,
    kmers: Path,
    report_handle,
    min_length: int = 100,
) -> tuple[dict[Path, str], Path]:
    """
    Split contigs that match the merged DB into separate FASTA files.

    Returns
    -------
    tuple
        (dictionary of contig fasta paths, temp directory path)
    """
    temp_index_prefix = Path(tempfile.mkdtemp(prefix="virstrain_bt2_")) / "bt2_index"
    temp_contig_dir = Path(tempfile.mkdtemp(prefix="virstrain_contigs_"))
    sam_file = temp_index_prefix.parent / "alignment.sam"

    # Build Bowtie2 index and align kmers
    run_cmd(["bowtie2-build", str(ingenome), str(temp_index_prefix)])
    run_cmd(
        [
            "bowtie2",
            "-x",
            str(temp_index_prefix),
            "-f",
            str(kmers),
            "--score-min",
            "C,0,0",
            "-S",
            str(sam_file),
        ]
    )

    matched_contigs: dict[str, str] = {}

    with sam_file.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) > 2 and fields[2] != "*":
                matched_contigs[fields[2]] = ""

    # Clean Bowtie2 intermediate files
    for file in temp_index_prefix.parent.glob("bt2_index*"):
        file.unlink(missing_ok=True)
    sam_file.unlink(missing_ok=True)

    sequences: dict[str, str] = {}
    current_header: str | None = None
    keep = False

    with ingenome.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_header = line
                contig_id = line.split()[0].removeprefix(">")
                if contig_id in matched_contigs:
                    keep = True
                    sequences[current_header] = ""
                else:
                    keep = False
                    report_handle.write(f"{contig_id}\tNA\tSkip\n")
            else:
                if keep and current_header is not None:
                    sequences[current_header] += line

    result: dict[Path, str] = {}

    for header, seq in sequences.items():
        if len(seq) < min_length:
            continue

        contig_id = header.removeprefix(">").split()[0]
        out_fasta = temp_contig_dir / f"{contig_id}.fasta"
        with out_fasta.open("w", encoding="utf-8") as out_handle:
            out_handle.write(f"{header}\n{seq}\n")
        result[out_fasta] = ""

    shutil.rmtree(temp_index_prefix.parent, ignore_errors=True)
    return result, temp_contig_dir


def read_third_data_line(report_file: Path) -> str:
    """
    Read the third line from a VirStrain report and convert tabs to pipes.
    """
    with report_file.open("r", encoding="utf-8") as handle:
        lines = [handle.readline().strip() for _ in range(3)]

    if len(lines) < 3 or not lines[2]:
        return "NA"

    return lines[2].replace("\t", "|")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="VirStrain.py",
        description=USAGE,
    )
    parser.add_argument(
        "-i",
        "--input_contigs",
        dest="input_contigs",
        type=str,
        required=True,
        help="Input contig FASTA data --- Required",
    )
    parser.add_argument(
        "-d",
        "--database_dir",
        dest="db_dir",
        type=str,
        required=True,
        help="Database dir --- Required",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="out_dir",
        type=str,
        default="VirStrain_Out",
        help="Output dir (default: current dir/VirStrain_Out)",
    )
    parser.add_argument(
        "-c",
        "--site_filter_cutoff",
        dest="sf_cutoff",
        type=float,
        default=0.05,
        help="The cutoff of filtering one site (default: 0.05)",
    )
    parser.add_argument(
        "-s",
        "--rank_by_sites",
        dest="rk_site",
        type=int,
        default=0,
        help="If set to 1, then VirStrain will sort the most possible strain by matches to the sites. (default: 0)",
    )
    parser.add_argument(
        "-m",
        "--high_mutation_virus",
        dest="hm_virus",
        default="1",
        nargs="?",
        help="If the virus has high mutation rate, use this option. (default: Not use)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    in_contigs = Path(args.input_contigs).resolve()
    db_dir = Path(args.db_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    rank_sites = int(args.rk_site)

    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir = out_dir / "Report"
    report_dir.mkdir(parents=True, exist_ok=True)

    merged_db_fa = db_dir / "merge_db.fa"
    if not in_contigs.exists():
        raise FileNotFoundError(f"Input contig FASTA not found: {in_contigs}")
    if not merged_db_fa.exists():
        raise FileNotFoundError(f"Merged DB FASTA not found: {merged_db_fa}")

    summary_report = out_dir / "VirStrain_contig_report.txt"

    with summary_report.open("w", encoding="utf-8") as out_handle:
        out_handle.write("Contigs_ID\tSpecies_info\tStrain_info\n")

        contigs, temp_contig_dir = split_contig(in_contigs, merged_db_fa, out_handle)

        try:
            for genome in contigs:
                contig_name = genome.stem
                out_handle.write(contig_name)

                target_sp, kmatch = prescan.scan(str(genome), str(merged_db_fa))

                if int(kmatch) < 10:
                    out_handle.write(f"\t{target_sp}:{kmatch}\tSkip\n")
                    continue

                nd_dir = db_dir / target_sp
                per_contig_report = Path(f"{contig_name}_VirStrain_report.txt")

                cmd = [
                    sys.executable,
                    "bin/S3_Strain_pred_My_Method_V0819_Val_contig.py",
                    "-i",
                    str(genome),
                    "-s",
                    str(nd_dir / "Pos-snp-kmer-all.txt"),
                    "-m",
                    str(nd_dir / "Strain_pos_snp_matrix_not_redundant_MM_Call.txt"),
                    "-f",
                    str(nd_dir / "Pos-snp-kmer-all.fa"),
                    "-c",
                    str(nd_dir / "Strain_cls_info.txt"),
                    "-b",
                    str(nd_dir / "SubCls_kmer.txt"),
                    "-d",
                    str(nd_dir),
                    "-o",
                    str(per_contig_report),
                ]

                if rank_sites == 1:
                    cmd.extend(["-r", "1"])

                run_cmd(cmd)

                third_line = read_third_data_line(per_contig_report)
                out_handle.write(f"\t{target_sp}:{kmatch}\t{third_line}\n")

                shutil.move(str(per_contig_report), str(report_dir / per_contig_report.name))

        finally:
            shutil.rmtree(temp_contig_dir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
