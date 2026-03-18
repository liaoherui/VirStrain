#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple
import networkx as nx
from networkx.algorithms.components import connected_components


INPUT_FASTA = Path("H1N1-16273-one-line.fasta")
OUTPUT_FASTA = Path("H1N1-16273-one-line-removeSame.fasta")
OUTPUT_CLUSTER = Path("H1N1-16273-one-line-removeSame.clstr")


def read_one_line_fasta(fasta_file: Path) -> Dict[str, str]:
    """
    Read a one-line FASTA file into a dictionary:
    {header: sequence}

    Assumes each sequence is on a single line immediately after its header.
    """
    sequences: Dict[str, str] = {}
    current_header: str | None = None

    with fasta_file.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_header = line
                sequences[current_header] = ""
            else:
                if current_header is None:
                    raise ValueError(
                        f"Found sequence line before any FASTA header in {fasta_file}"
                    )
                sequences[current_header] = line

    return sequences


def to_edges(items: Iterable[str]) -> Iterator[Tuple[str, str]]:
    """
    Convert an iterable like [a, b] or [a, b, c] into edges:
    (a, b), (b, c), ...
    """
    iterator = iter(items)
    try:
        last = next(iterator)
    except StopIteration:
        return

    for current in iterator:
        yield last, current
        last = current


def to_graph(groups: List[List[str]]) -> nx.Graph:
    """
    Build an undirected graph from grouped items.
    Each group connects its elements as a chain.
    """
    graph = nx.Graph()
    for group in groups:
        graph.add_nodes_from(group)
        graph.add_edges_from(to_edges(group))
    return graph


def main() -> None:
    dseq = read_one_line_fasta(INPUT_FASTA)

    # Original structure: [header, sequence]
    arr: List[List[str]] = [[header, seq] for header, seq in dseq.items()]

    graph = to_graph(arr)

    with OUTPUT_FASTA.open("w", encoding="utf-8") as fasta_out, \
         OUTPUT_CLUSTER.open("w", encoding="utf-8") as cluster_out:

        for component in connected_components(graph):
            component_list = list(component)

            # Keep only FASTA headers from the component
            headers = [item for item in component_list if item.startswith(">")]

            if not headers:
                continue

            representative = headers[0]
            fasta_out.write(f"{representative}\n{dseq[representative]}\n")

            cluster_entries = [
                f"{header}:{len(dseq[header])}"
                for header in headers
            ]
            cluster_out.write("\t".join(cluster_entries) + "\n")


if __name__ == "__main__":
    main()
