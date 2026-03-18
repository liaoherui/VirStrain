#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def extract_strain_prefixes(df: pd.DataFrame, skip_none: bool = False) -> list[str]:
    """
    Extract unique strain prefixes from dataframe columns.
    """
    prefixes: list[str] = []
    for col in df.columns:
        if "ID" in col:
            continue
        if skip_none and "None" in col:
            continue
        name = re.split(r"_", col)[0]
        if name not in prefixes:
            prefixes.append(name)
    return prefixes


def main() -> int:
    mps_file = Path("Mps_ps_depth_HIV.csv")
    ops_file = Path("Ops_ps_depth_HIV.csv")
    output_html = Path("Test.html")

    if not mps_file.exists():
        raise FileNotFoundError(f"Missing input file: {mps_file}")
    if not ops_file.exists():
        raise FileNotFoundError(f"Missing input file: {ops_file}")

    sfl = pd.read_csv(mps_file)
    sf2 = pd.read_csv(ops_file)

    if "ID" not in sfl.columns:
        raise ValueError("Column 'ID' not found in Mps_ps_depth_HIV.csv")
    if "ID" not in sf2.columns:
        raise ValueError("Column 'ID' not found in Ops_ps_depth_HIV.csv")

    s = extract_strain_prefixes(sfl)
    so = extract_strain_prefixes(sf2, skip_none=True)

    if not s:
        raise ValueError("No strain-related columns found in Mps_ps_depth_HIV.csv")

    # Keep the original behavior: only use the first detected strain prefix
    strain = s[0]
    freq_col = f"{strain}_Freq"
    lnum_col = f"{strain}_LNum"

    if freq_col not in sfl.columns:
        raise ValueError(f"Expected column not found: {freq_col}")
    if lnum_col not in sfl.columns:
        raise ValueError(f"Expected column not found: {lnum_col}")

    zero_mask = sfl[freq_col] == 0
    zx = sfl.loc[zero_mask, "ID"]
    zy = sfl.loc[zero_mask, freq_col]

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    fig.add_trace(
        go.Bar(
            x=sfl["ID"],
            y=sfl[freq_col],
            name="depth",
        ),
        secondary_y=False,
    )

    fig.add_trace(
        go.Scatter(
            x=sfl["ID"],
            y=(-1) * sfl[lnum_col],
            mode="lines+markers",
            opacity=0.5,
            marker=dict(size=5),
            name="(-1)*(strain number)",
        ),
        secondary_y=True,
    )

    fig.add_trace(
        go.Scatter(
            x=zx,
            y=zy,
            mode="markers",
            opacity=0.5,
            marker=dict(size=5, color="red"),
            name="Zero depth pos",
        ),
        secondary_y=False,
    )

    max_lnum = sfl[lnum_col].max()
    max_freq = sfl[freq_col].max()

    fig.update_layout(
        autosize=False,
        width=1200,
        height=300,
        title={
            "text": f"Most Possible Strain: {strain}",
            "xanchor": "center",
        },
    )

    fig.update_yaxes(range=[(-3) * max_lnum, 3 * max_lnum], secondary_y=True)
    fig.update_yaxes(range=[(-1) * max_freq, max_freq], secondary_y=False)

    output_html.write_text(
        fig.to_html(full_html=False, include_plotlyjs="cdn"),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
