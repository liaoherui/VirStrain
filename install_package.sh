# Use conda to install latest version of plotly to avoid error
conda install -c ployly plotly &&\
# plotly
pip install plotly==3.10.0  &&\
# Networkx
conda install networkx &&\
# Biopython
conda install -c conda-forge biopython &&\
# Numpy
conda install numpy &&\
# Pandas
conda install pandas


