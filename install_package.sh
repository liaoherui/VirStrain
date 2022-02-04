#!/usr/bin/env bash
# Use conda to install latest version of plotly to avoid error
conda install -y -c plotly plotly &&\
# plotly
pip install plotly==3.10.0  &&\
# Networkx
conda install -y networkx &&\
# Biopython
conda install -y -c conda-forge biopython &&\
# Numpy
conda install -y numpy &&\
# Pandas
conda install -y pandas


