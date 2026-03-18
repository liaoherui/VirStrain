#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

__author__ = "Liao Herui"
USAGE = "VirStrain - An RNA virus strain-level identification tool for short reads."


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="VirStrain_contigDB_merge.py",
        description=USAGE,
    )
    parser.add_argument(
        "-i",
        "--merge_db",
        dest="merge_db",
        type=str,
        required=True,
        help=(
            "Comma-separated directories of merged VirStrain databases. "
            "Example: -i VirStrain_DB/SCOV2 or "
            "-i VirStrain_DB/SCOV2,VirStrain_DB/H1N1"
        ),
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="out_dir",
        type=str,
        default="VirStrain_DB_merge",
        help="Output directory of merged database. Default: ./VirStrain_DB_merge",
    )
    return parser.parse_args()


def merge_databases(db_dirs: list[Path], out_dir: Path) -> None:
    """
    Merge Pos-snp-kmer-all.fa entries across VirStrain database directories.

    Original behavior preserved:
    - For each DB, read Pos-snp-kmer-all.fa
    - Append directory basename to FASTA header
    - Keep only sequences that occur exactly once across all DBs
    - Copy each source DB directory into the output directory
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_counts: dict[str, int] = {}
    seq_annotation: dict[str, str] = {}

    merged_fasta = out_dir / "merge_db.fa"

    with merged_fasta.open("w", encoding="utf-8") as out_handle:
        for db_dir in db_dirs:
            fasta_file = db_dir / "Pos-snp-kmer-all.fa"
            prefix = db_dir.name

            if not db_dir.exists():
                raise FileNotFoundError(f"Database directory not found: {db_dir}")

            if not fasta_file.exists():
                raise FileNotFoundError(f"Missing FASTA file: {fasta_file}")

            current_header: str | None = None

            with fasta_file.open("r", encoding="utf-8") as in_handle:
                for raw_line in in_handle:
                    line = raw_line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        current_header = f"{line}_{prefix}"
                    else:
                        if current_header is None:
                            raise ValueError(
                                f"Sequence found before FASTA header in file: {fasta_file}"
                            )

                        seq_counts[line] = seq_counts.get(line, 0) + 1
                        seq_annotation[line] = current_header

            # Copy DB directory into the output directory.
            dest_dir = out_dir / db_dir.name
            if dest_dir.exists():
                shutil.rmtree(dest_dir)
            shutil.copytree(db_dir, dest_dir)

        for seq, count in seq_counts.items():
            if count > 1:
                continue
            out_handle.write(f"{seq_annotation[seq]}\n{seq}\n")


def main() -> int:
    args = parse_args()

    db_dirs = [Path(x.strip()) for x in args.merge_db.split(",") if x.strip()]
    out_dir = Path(args.out_dir)

    if not db_dirs:
        raise ValueError("No valid database directories were provided in --merge_db")

    merge_databases(db_dirs, out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
