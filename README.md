# CLIPNET paper

 Scripts and notebooks to reproduce figures in the CLIPNET paper (preprint [here](https://www.biorxiv.org/content/10.1101/2024.03.13.583868)). Code for generating new predictions and feature interpretions with CLIPNET is available in a [separate repo](https://github.com/Danko-Lab/clipnet/).

## Install dependencies

Dependencies used to produce plots are listed in `requirements.txt`. To install them, run:

```bash
# Either use conda/mamba or python venv to isolate installation
mamba create -n clipnet-paper -c conda-forge -c bioconda python=3.9
mamba activate clipnet-paper
pip install -r requirements.txt
```

## Code/Notebooks

Jupyter notebooks to reproduce figures used in the paper are organized by analysis type. While we provide some scripts/instructions to calculate predictions and feature interpretations (we haven't written many of these up, so please just raise an issue if you want clarification on something), many of these calculations will be very time/resource-intensive. As a result, for the purposes of reproducing the figures in the paper, we assume precalculated predictions and feature interpretations, which we have archived on Zenodo.

## Downloads

Training data and data to reproduce figuress are available at [10.5281/zenodo.10597358](https://zenodo.org/doi/10.5281/zenodo.10597358). To preserve directory structure, we packaged the data into tar files, divided roughly by figure/analysis. For a description of these files, see the `DOWNLOADS_README.md` file.
