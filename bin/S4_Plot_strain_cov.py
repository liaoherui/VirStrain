#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def extract_strain_names(df: pd.DataFrame, skip_none: bool = False) -> list[str]:
    """
    Extract unique strain names from columns like:
    <strain>_Freq and <strain>_LNum
    """
    strains: list[str] = []

    for col in df.columns:
        if "ID" in col:
            continue
        if skip_none and "None" in col:
            continue

        name = re.sub(r"_Freq$", "", col)
        name = re.sub(r"_LNum$", "", name)

        if name not in strains:
            strains.append(name)

    return strains


def zero_depth_points(df: pd.DataFrame, strain: str) -> tuple[pd.Series, pd.Series]:
    """
    Get x/y points where depth is zero for a given strain.
    """
    freq_col = f"{strain}_Freq"
    zero_mask = df[freq_col] == 0
    zx = df.loc[zero_mask, "ID"]
    zy = df.loc[zero_mask, freq_col]
    return zx, zy


def make_strain_figure(df: pd.DataFrame, strain: str, title: str) -> go.Figure:
    """
    Build a plotly figure for one strain.
    """
    freq_col = f"{strain}_Freq"
    lnum_col = f"{strain}_LNum"

    if freq_col not in df.columns:
        raise ValueError(f"Missing expected column: {freq_col}")
    if lnum_col not in df.columns:
        raise ValueError(f"Missing expected column: {lnum_col}")
    if "ID" not in df.columns:
        raise ValueError("Missing expected column: ID")

    zx, zy = zero_depth_points(df, strain)

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Bar(
            x=df["ID"],
            y=df[freq_col],
            name="depth",
            marker=dict(color="#1f77b4"),
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=df["ID"],
            y=(-1) * df[lnum_col],
            mode="lines+markers",
            opacity=0.7,
            marker=dict(size=5, color="#f0ad6d"),
            line=dict(color="#f0ad6d"),
            name="(-1)*(strain number)",
        ),
        secondary_y=True,
    )
    fig.add_trace(
        go.Scatter(
            x=zx,
            y=zy,
            mode="markers",
            opacity=0.7,
            marker=dict(size=5, color="red"),
            name="Zero depth pos",
        ),
        secondary_y=False,
    )

    max_lnum = df[lnum_col].max()
    max_freq = df[freq_col].max()

    fig.update_layout(
        autosize=False,
        width=1200,
        height=300,
        template="plotly_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        title={"text": title, "x": 0.5, "xanchor": "center"},
    )
    fig.update_yaxes(range=[(-3) * max_lnum, 3 * max_lnum], secondary_y=True)
    fig.update_yaxes(range=[(-1) * max_freq, max_freq], secondary_y=False)

    return fig


def main() -> int:
    mps_file = Path("Mps_ps_depth.csv")
    ops_file = Path("Ops_ps_depth.csv")
    output_html = Path("VirStrain_report.html")

    if not mps_file.exists():
        raise FileNotFoundError(f"Missing input file: {mps_file}")
    if not ops_file.exists():
        raise FileNotFoundError(f"Missing input file: {ops_file}")

    sfl = pd.read_csv(mps_file)
    sf2 = pd.read_csv(ops_file)

    if "ID" not in sfl.columns:
        raise ValueError("Column 'ID' not found in Mps_ps_depth.csv")
    if "ID" not in sf2.columns:
        raise ValueError("Column 'ID' not found in Ops_ps_depth.csv")

    main_strains = extract_strain_names(sfl)
    other_strains = extract_strain_names(sf2, skip_none=True)

    if not main_strains:
        raise ValueError("No strain columns found in Mps_ps_depth.csv")

    with output_html.open("w", encoding="utf-8") as ft:
        main_strain = main_strains[0]
        fig = make_strain_figure(
            sfl,
            main_strain,
            f"Most Possible Strain: {main_strain}",
        )
        ft.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))

        if other_strains:
            rank = 1
            for strain in other_strains:
                fig = make_strain_figure(
                    sf2,
                    strain,
                    f"Other Possible Strain: {strain} (Rank:{rank})",
                )
                ft.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
                rank += 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
