#!/usr/bin/env python3

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path


def scan(ingenome: str | Path, db: str | Path) -> tuple[str, int]:
    """
    Scan an input genome against the VirStrain merged DB using jellyfish.

    Returns
    -------
    tuple[str, int]
        Top matching species/database label and its matched kmer count.
        Returns ('NA', 0) if no match is found.
    """
    file_dir = Path(__file__).resolve().parent
    ingenome = Path(ingenome).resolve()
    db = Path(db).resolve()

    if not ingenome.exists():
        raise FileNotFoundError(f"Input genome file not found: {ingenome}")
    if not db.exists():
        raise FileNotFoundError(f"Database file not found: {db}")

    # Map kmer sequence -> species/database label
    kmer_to_species: dict[str, str] = {}
    current_species: str | None = None

    with db.open("r", encoding="utf-8") as fd:
        for raw_line in fd:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                parts = line.split("_")
                current_species = "_".join(parts[1:])
            else:
                if current_species is not None:
                    kmer_to_species[line] = current_species

    jellyfish_bin = file_dir / "jellyfish-linux"

    with tempfile.TemporaryDirectory(prefix="virstrain_prescan_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        jf_file = tmpdir_path / "Tem_Vs.jf"
        fa_file = tmpdir_path / "Tem_Vs.fa"

        try:
            subprocess.run(
                [
                    str(jellyfish_bin),
                    "count",
                    "-m",
                    "25",
                    "-s",
                    "100M",
                    "-t",
                    "8",
                    "--if",
                    str(db),
                    "-o",
                    str(jf_file),
                    str(ingenome),
                ],
                check=True,
            )

            with fa_file.open("w", encoding="utf-8") as out_handle:
                subprocess.run(
                    [
                        str(jellyfish_bin),
                        "dump",
                        "-c",
                        str(jf_file),
                    ],
                    check=True,
                    stdout=out_handle,
                )

            species_counts: dict[str, int] = {}

            with fa_file.open("r", encoding="utf-8") as f:
                for raw_line in f:
                    line = raw_line.strip()
                    if not line:
                        continue

                    fields = line.split()
                    if len(fields) < 2:
                        continue

                    kmer = fields[0]
                    try:
                        count = int(fields[-1])
                    except ValueError:
                        continue

                    if count > 0 and kmer in kmer_to_species:
                        species = kmer_to_species[kmer]
                        species_counts[species] = species_counts.get(species, 0) + 1

            if not species_counts:
                return "NA", 0

            top_species, top_count = sorted(
                species_counts.items(),
                key=lambda item: item[1],
                reverse=True,
            )[0]

            return top_species, top_count

        finally:
            # TemporaryDirectory handles cleanup automatically,
            # but keeping the structure explicit here is fine.
            pass
