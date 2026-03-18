#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild annotated cluster file for VirStrain."
    )
    parser.add_argument(
        "-i",
        "--input_fa",
        required=True,
        help="Input FASTA file",
    )
    parser.add_argument(
        "-c",
        "--cls_file",
        required=True,
        help="Cluster file",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        default="VirStrain_DB",
        help="Output directory (default: VirStrain_DB)",
    )
    return parser.parse_args()


def load_headers(input_fa: Path) -> dict[str, str]:
    """
    Parse FASTA headers and extract annotation text.
    """
    dhead: dict[str, str] = {}

    with input_fa.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or not line.startswith(">"):
                continue

            if "|" in line:
                pre = line.split("|")[0].strip()
            else:
                pre = line.split()[0].strip()

            info = re.sub(r".*human", "", line)
            info = re.sub(r",.*", "", info)
            dhead[pre] = info

    return dhead


def rebuild_cluster_annotations(input_fa: Path, cls_file: Path, out_dir: Path) -> Path:
    """
    Rebuild the annotated cluster file.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    dhead = load_headers(input_fa)
    output_file = out_dir / "Anno_rebuild_cls.clstr"

    with cls_file.open("r", encoding="utf-8") as f_in, output_file.open(
        "w", encoding="utf-8"
    ) as f_out:
        for raw_line in f_in:
            line = raw_line.strip()
            if not line:
                continue

            if "Cluster" in line or "Cls" in line:
                if "Cluster" in line:
                    f_out.write(f"{line}\n")
                else:
                    f_out.write(f"\t{line}\n")
            else:
                fields = line.split("\t")
                key = fields[0]
                annotation = dhead.get(key, "NA")
                f_out.write(f"\t\t{line}\t{annotation}\n")

    return output_file


def main() -> int:
    args = parse_args()

    input_fa = Path(args.input_fa).resolve()
    cls_file = Path(args.cls_file).resolve()
    out_dir = Path(args.out_dir).resolve()

    if not input_fa.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fa}")
    if not cls_file.exists():
        raise FileNotFoundError(f"Cluster file not found: {cls_file}")

    rebuild_cluster_annotations(input_fa, cls_file, out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
