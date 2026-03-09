#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import math
import sys
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate per-column entropy from an MSA and iteratively cluster strains."
    )
    parser.add_argument(
        "-i",
        "--input_msa",
        required=True,
        help="Input MSA in FASTA format",
    )
    parser.add_argument(
        "--dash-cutoff",
        type=float,
        default=0.01,
        help="Maximum allowed dash fraction per column (default: 0.01)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output_test.txt",
        help="Output file (default: output_test.txt)",
    )
    return parser.parse_args()


def estimate_shannon_entropy(dna_seq: str) -> tuple[float, dict[str, int]]:
    """
    Compute Shannon entropy for a DNA sequence column.
    Only A/T/G/C (upper or lower case) contribute to the entropy calculation.
    """
    valid_bases = {"a", "t", "g", "c", "A", "T", "G", "C"}
    bases_raw = dict(collections.Counter(dna_seq))

    bases: dict[str, int] = {}
    m = 0
    for base, count in bases_raw.items():
        if base not in valid_bases:
            continue
        bases[base] = count
        m += count

    if m == 0:
        return 0.0, bases_raw

    shannon_entropy_value = 0.0
    for count in bases.values():
        p_i = count / float(m)
        shannon_entropy_value += p_i * math.log(p_i)

    return -shannon_entropy_value, bases_raw


def estimate_shannon_entropy_iterate(
    alignment: MultipleSeqAlignment,
    used_strains: set[int],
    used_columns: set[int],
    centropy: dict[int, float],
) -> tuple[dict[int, float], dict[int, dict[str, int]]]:
    """
    Recompute per-column entropies after excluding used strains and columns.
    """
    new_centropy: dict[int, float] = {}
    new_cbase_freq: dict[int, dict[str, int]] = {}
    valid_bases = {"a", "t", "g", "c", "A", "T", "G", "C"}

    for column_id in centropy:
        if column_id in used_columns:
            new_centropy[column_id] = 1000.0
            continue

        temp_base_freq: dict[str, int] = {}
        m = 0

        for strain_index, base in enumerate(alignment[:, column_id]):
            if strain_index in used_strains:
                continue
            if base not in valid_bases:
                continue

            temp_base_freq[base] = temp_base_freq.get(base, 0) + 1
            m += 1

        if m == 0:
            new_centropy[column_id] = 1000.0
            new_cbase_freq[column_id] = temp_base_freq
            continue

        shannon_entropy_value = 0.0
        for count in temp_base_freq.values():
            p_i = count / float(m)
            shannon_entropy_value += p_i * math.log(p_i)

        new_centropy[column_id] = -shannon_entropy_value
        new_cbase_freq[column_id] = temp_base_freq

    return new_centropy, new_cbase_freq


def cls_with_entropy(
    alignment: MultipleSeqAlignment,
    centropy: dict[int, float],
    cbase_freq: dict[int, dict[str, int]],
    sid2name: dict[int, str],
    sname2id: dict[str, int],
    output_path: Path,
) -> None:
    """
    Iteratively select low-entropy columns and group strains by minor alleles.
    """
    used_strains: set[int] = set()
    used_columns: set[int] = set()
    valid_bases = {"a", "t", "g", "c", "A", "T", "G", "C"}

    total_columns = len(centropy)
    total_strains = len(sname2id)

    with output_path.open("w", encoding="utf-8") as out_handle:
        out_handle.write(
            f"finish inputting {len(sname2id)} sequences from the alignment\n"
        )

        while True:
            if len(used_strains) == len(sname2id) or len(used_columns) == total_columns:
                break

            print(
                f"Progress: C: {len(used_columns)}/{total_columns} "
                f"S: {len(used_strains)}/{total_strains}"
            )

            ranked_columns = sorted(centropy.items(), key=lambda item: item[1])

            current_column: int | None = None
            for column_id, _entropy in ranked_columns:
                if column_id not in used_columns:
                    current_column = column_id
                    break

            if current_column is None:
                break

            rawseq = alignment[:, current_column]
            if current_column not in cbase_freq or not cbase_freq[current_column]:
                used_columns.add(current_column)
                centropy, cbase_freq = estimate_shannon_entropy_iterate(
                    alignment, used_strains, used_columns, centropy
                )
                continue

            ranked_bases = sorted(cbase_freq[current_column].items(), key=lambda item: item[1])

            round_strains: list[int] = []

            for base, _count in ranked_bases:
                if base not in valid_bases:
                    continue
                if base == ranked_bases[-1][0]:
                    continue

                strain_ids = [i for i, letter in enumerate(rawseq) if letter == base]
                round_strains.extend(strain_ids)

                out_handle.write(f"column {current_column}")
                for b in ["a", "t", "g", "c"]:
                    out_handle.write(f" {b}({cbase_freq[current_column].get(b, 0)})")
                out_handle.write("\n")

                for strain_id in strain_ids:
                    out_handle.write(f">{sid2name[strain_id]} ")
                out_handle.write(f"{current_column} {base}\n")

                cbase_freq[current_column][base] = 0

            for strain_id in round_strains:
                used_strains.add(strain_id)

            used_columns.add(current_column)

            centropy, cbase_freq = estimate_shannon_entropy_iterate(
                alignment, used_strains, used_columns, centropy
            )


def build_entropy_tables(
    alignment: MultipleSeqAlignment,
    dash_cutoff: float,
) -> tuple[dict[int, float], dict[int, dict[str, int]], dict[int, str]]:
    """
    Build per-column entropy and base frequency tables after dash filtering.
    """
    strain_ids = range(len(alignment[:, 0]))
    num_columns = len(alignment[0].seq)

    centropy: dict[int, float] = {}
    cbase_freq: dict[int, dict[str, int]] = {}
    cid_seq: dict[int, str] = {}

    for column in range(num_columns):
        seq = alignment[:, column]
        entropy, base_count = estimate_shannon_entropy(seq)

        if "-" not in base_count:
            if entropy != 0:
                centropy[column] = entropy
                cbase_freq[column] = base_count
                cid_seq[column] = seq
        else:
            if base_count["-"] <= dash_cutoff * len(list(strain_ids)):
                if entropy != 0:
                    centropy[column] = entropy
                    cbase_freq[column] = base_count
                    cid_seq[column] = seq

    return centropy, cbase_freq, cid_seq


def main() -> int:
    args = parse_args()

    input_msa = Path(args.input_msa).resolve()
    output_path = Path(args.output).resolve()
    dash_cutoff = args.dash_cutoff

    if not input_msa.exists():
        raise FileNotFoundError(f"Input MSA not found: {input_msa}")

    alignment = AlignIO.read(str(input_msa), "fasta")

    strain_ids = list(range(len(alignment[:, 0])))
    strain_names: list[str] = []
    sname2seq: dict[str, str] = {}

    for record in alignment:
        name = record.id.strip()
        strain_names.append(name)
        sname2seq[name] = str(record.seq)

    sid2name = dict(zip(strain_ids, strain_names))
    sname2id = dict(zip(strain_names, strain_ids))

    centropy, cbase_freq, _cid_seq = build_entropy_tables(alignment, dash_cutoff)

    cls_with_entropy(
        alignment=alignment,
        centropy=centropy,
        cbase_freq=cbase_freq,
        sid2name=sid2name,
        sname2id=sname2id,
        output_path=output_path,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
