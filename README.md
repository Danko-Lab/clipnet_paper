# CLIPNET paper

 Scripts and notebooks to reproduce figures in the CLIPNET paper. Code for generating new predictions and feature interpretions with CLIPNET is available in a [separate repo](https://github.com/Danko-Lab/clipnet/).

## Install dependencies

Dependencies used to produce plots are listed in `requirements.txt`. To install them, run:

```bash
# Either use conda/mamba or python venv to isolate installation
mamba create -n clipnet-paper -c conda-forge -c bioconda python=3.9
mamba activate clipnet-paper
pip install -r requirements.txt
```

## Code/Notebooks

Jupyter notebooks to reproduce figures used in the paper are organized by figure. Please see READMEs in the subdirectories for more details. While we provide some scripts/instructions to calculate predictions and feature interpretations, many of these calculations will be very time/resource-intensive. As a result, for the purposes of reproducing the figures in the paper, we assume precalculated predictions and feature interpretations, which we provide as specified below.

## Downloads

Data to reproduce figures and training data are available at [10.5281/zenodo.10597358](https://zenodo.org/doi/10.5281/zenodo.10597358). To preserve directory structure, we packaged the data into tar files, divided roughly by figure/analysis. Here is a brief description of the contents of each tar file:

- `training_data.tar.gz`: Contains processed data used to train the CLIPNET models.
  - `individual_pints_peaks/`: Contains the PINTS peaks for each individual PRO-cap library.
  - `individual_jittered_windows/`: Contains the jittered (uniformly random, +/- 250bp around center of each peak) 1 kb windows for each individual PRO-cap library.
  - `processed_data/`: Contains the processed data used to train the models, including the individualized sequences and PRO-cap signal (RPM normalized). Packaged as npz arrays. Data were concatenated across libraries, then split into the data folds described in `processed_data/data_fold_assignments.csv.gz`. We note that the PRO-cap data are structured as N x 2000 arrays (1000 bp pl strand, 1000 bp mn strand). The sequence data are structured as N x 1000 x 4 arrays (N = number of sequences, 1000 = sequence length, 4 = two-hot encoding of sequences).
- `evaluation_data.tar.gz`: Contains data used to evaluate the performance of the CLIPNET models.
  - `fixed_uniq_windows.bed.gz`: A fixed set of 48,058 1 kb windows used to evaluate the models. We selected PRO-cap peaks that were present in at least 60 of the 67 libraries, then selected 1 kb windows around each of them (with 250 bp jittering).
  - `processed_data/`: Contains the processed data used to evaluate the models.
    - `procap/`: Contains the processed PRO-cap signal (csv) for each data fold. Contains both concatenated and individualized data.
    - `sequences/`: Contains the sequences (fasta) for each data fold. Contains both concatenated and individualized data.
  - `predictions/`: Contains the predictions of the CLIPNET models on the fixed set of 1 kb windows.
    - `ensemble_test`: Contains the predictions and evaluation metrics of the ensemble model and of model folds on the complete hold out data set (fold 0).
    - `individual_test`: Contains the predictions and evaluation metrics of the model folds on the individual model hold out folds (model 1 used fold 1 as a hold out, model 2 used fold 2, etc).
- `deepshap_scores.tar.gz`: Contains DeepSHAP contribution scores.
  - `merged_windows_all.bed.gz`: A nonredundant set of 212,777 windows around PRO-cap peaks (union across all libraries) used for calculating DeepSHAP scores.
  - `all_tss_windows_reference_seq.fna.gz`: The reference (hg38) sequence for the windows in `merged_windows_all.bed.gz`.
  - `all_seqs_onehot.npz`: A one-hot encoded version of the reference sequence. This and the score arrays are structured as N x 4 x 1000 arrays for compatibility with TF-MoDISco.
  - `mean_across_folds_all_profile.npz`: The profile contribution scores (mean across model folds).
  - `mean_across_folds_all_quantity.npz`: The quantity contribution scores (mean across model folds).
- `tfmodisco_results.tar.gz`: Contains TF-MoDISco results.
  - `mean_across_folds_all_profile_modisco.h5`: The TF-MoDISco results for the profile contribution scores.
  - `mean_across_folds_all_quantity_modisco.h5`: The TF-MoDISco results for the quantity contribution scores.
  - `mean_across_folds_all_profile_modisco/`: A report of the TF-MoDISco results for the profile contribution scores.
  - `mean_across_folds_all_quantity_modisco/`: A report of the TF-MoDISco results for the quantity contribution scores.
  - `mean_across_folds_all_modisco_positions.h5`: Distribution of TF-MoDISco motif positions around the max TSS for each window.
- `qtl_analysis`: Contains analysis of both tiQTLs and diQTLs.
  - `tiqtl/`: Contains the tiQTL analysis.
    - `tiqtl_results/`: Contains the tiQTL results.
    - `tiqtl_windows.bed.gz`: The windows used for the tiQTL analysis.
  - `diqtl/`: Contains the diQTL analysis.
    - `diqtl_results/`: Contains the diQTL results.
    - `diqtl_windows.bed.gz`: The windows used for the diQTL analysis.
