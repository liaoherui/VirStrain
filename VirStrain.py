#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

__author__ = "Liao Herui"

USAGE = "VirStrain - An RNA virus strain-level identification tool for short reads."


def run_cmd(cmd: list[str], cwd: Path | None = None) -> None:
    """Run a command and fail immediately if it returns a non-zero exit code."""
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def move_existing(file_names: list[str], out_dir: Path) -> None:
    """Move output files to the destination directory if they exist."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in file_names:
        src = Path(name)
        if src.exists():
            shutil.move(str(src), str(out_dir / src.name))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="VirStrain.py", description=USAGE)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s v1.17")
    parser.add_argument(
        "-i",
        "--input_reads",
        dest="input_reads",
        type=str,
        required=True,
        help="Input fastq data --- Required",
    )
    parser.add_argument(
        "-p",
        "--input_reads2",
        dest="input_reads2",
        type=str,
        default=None,
        help="Input fastq data for PE reads.",
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
        "-f",
        "--turn_off_figures",
        dest="close_fig",
        type=int,
        default=0,
        help="If set to 1, then VirStrain will not generate figures. (default: 0)",
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

    in_read1 = Path(args.input_reads).resolve()
    in_read2 = Path(args.input_reads2).resolve() if args.input_reads2 else None
    db_dir = Path(args.db_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    rks = int(args.rk_site)
    cf = int(args.close_fig)
    hm = args.hm_virus

    out_dir.mkdir(parents=True, exist_ok=True)

    # Note:
    # The original script parses sf_cutoff but never uses it.
    # Preserved here for compatibility, but intentionally unused.
    _sfc = float(args.sf_cutoff)

    if not in_read1.exists():
        raise FileNotFoundError(f"Input reads not found: {in_read1}")
    if in_read2 is not None and not in_read2.exists():
        raise FileNotFoundError(f"Paired-end reads file not found: {in_read2}")
    if not db_dir.exists():
        raise FileNotFoundError(f"Database directory not found: {db_dir}")

    # whether this is the virus with high mutation rate, like HIV
    # Preserved original logic:
    # if hm is falsy / absent, run HIV-specific script.
    if not hm:
        hiv_cmd = [
            sys.executable,
            "bin/S3_Strain_pred_My_Method_For_HIV.py",
            "-i",
            str(in_read1),
            "-s",
            str(db_dir / "Pos-snp-kmer-all.txt"),
            "-f",
            str(db_dir / "Pos-snp-kmer-all.fa"),
            "-c",
            str(db_dir / "Remove_redundant_matrix_MM_Call.clstr"),
            "-d",
            str(db_dir / "ID2Name.txt"),
            "-b",
            str(db_dir),
            "-o",
            "VirStrain_report.txt",
        ]

        if in_read2 is not None:
            hiv_cmd.extend(["-p", str(in_read2)])

        run_cmd(hiv_cmd)

        if Path("VirStrain_report.txt").exists():
            run_cmd([sys.executable, "bin/S4_Plot_strain_cov.py"])
            move_existing(
                [
                    "VirStrain_report.txt",
                    "VirStrain_report.html",
                    "Mps_ps_depth.csv",
                    "Ops_ps_depth.csv",
                ],
                out_dir,
            )
        return 0

    # Start to identify strains
    strain_cmd = [
        sys.executable,
        "bin/S3_Strain_pred_My_Method_V0819_Val.py",
        "-i",
        str(in_read1),
        "-s",
        str(db_dir / "Pos-snp-kmer-all.txt"),
        "-m",
        str(db_dir / "Strain_pos_snp_matrix_not_redundant_MM_Call.txt"),
        "-f",
        str(db_dir / "Pos-snp-kmer-all.fa"),
        "-c",
        str(db_dir / "Strain_cls_info.txt"),
        "-b",
        str(db_dir / "SubCls_kmer.txt"),
        "-d",
        str(db_dir),
        "-o",
        "VirStrain_report.txt",
    ]

    if in_read2 is not None:
        strain_cmd.extend(["-p", str(in_read2)])

    if rks != 0:
        strain_cmd.extend(["-r", "1"])

    run_cmd(strain_cmd)

    # Plot HTML page
    if Path("VirStrain_report.txt").exists():
        if cf == 0:
            run_cmd([sys.executable, "bin/S4_Plot_strain_cov.py"])
            move_existing(
                [
                    "VirStrain_report.txt",
                    "VirStrain_report.html",
                    "Mps_ps_depth.csv",
                    "Ops_ps_depth.csv",
                ],
                out_dir,
            )
        else:
            move_existing(
                [
                    "VirStrain_report.txt",
                    "Mps_ps_depth.csv",
                    "Ops_ps_depth.csv",
                ],
                out_dir,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
