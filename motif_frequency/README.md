# DeepSHAP and TF-MoDISco calculation

This directory contains the code and notebooks to calculate DeepSHAP and TF-MoDISco scores for the CLIPNET paper.

For the purposes of motif discovery with TF-MoDISco, we used a set of ~200k sequence windows (the union across all PRO-cap libraries). We provide the DeepSHAP scores and TF-MoDISco results for these windows in the `deepshap_scores.tar.gz` and `tfmodisco_results.tar.gz` files, respectively. The `deepshap_scores.tar.gz` file contains the DeepSHAP scores for the sequence windows, and the `tfmodisco_results.tar.gz` file contains the TF-MoDISco results for the DeepSHAP scores. The reference sequence for the windows is also provided in the `deepshap_scores.tar.gz` file.

To recalculate DeepSHAP scores, first download and unpack the DeepSHAP tar:

```bash
wget https://zenodo.org/record/10597358/files/deepshap_scores.tar.gz -P /path/to/scratch/
tar -xvf /path/to/scratch/deepshap_scores.tar.gz -C /path/to/scratch/
```

Then, run the `calculate_deepshap_scores.py` on script from our [CLIPNET](https://github.com/Danko-Lab/clipnet) package. This script will likely take a long time, as calculating DeepSHAP scores on so many examples is very compute intensive, even on a GPU (please do not run on CPU). If GPU parallelization is possible (i.e., you have multiple GPUs available), you can parallelize the calculation by calculating DeepSHAP scores for each model fold in parallel (specify the model file path with `--model_fp`), then taking the mean across folds. This is the approach we used to calculate the DeepSHAP scores for the CLIPNET paper. Here are some example scripts (see the documentation in the CLIPNET repo for more details).

```bash
python calculate_deepshap_scores.py \
    /path/to/scratch/deepshap_scores/all_tss_windows_reference_seq.fna.gz \
    /path/to/scratch/deepshap_scores/all_tss_windows_deepshap_profile_1.npz \
    /path/to/scratch/all_tss_windows_reference_seq_onehot.npz \
    --model_fp /path/to/models/fold_1.h5 \
    --mode profile \
    --gpu

python calculate_deepshap_scores.py \
    /path/to/scratch/deepshap_scores/all_tss_windows_reference_seq.fna.gz \
    /path/to/scratch/deepshap_scores/all_tss_windows_deepshap_quantity_1.npz \
    /path/to/scratch/all_tss_windows_reference_seq_onehot.npz \
    --model_fp /path/to/models/fold_1.h5 \
    --mode quantity \
    --gpu
```

TF-MoDISco can then be used to cluster profile and quantity seqlets into motifs. We used the [lite implementation](https://github.com/jmschrei/tfmodisco-lite/) of TF-MoDISco. Please note that tfmodisco-lite is not installed by default either in the requirements for this repo or for the main CLIPNET repo. For dependency conflict reasons, I installed it in a separate conda environment. Here is a sketch of how to run TF-MoDISco (with the parameters we used):

```bash
output=profile # or quantity

# Calculate TF-MoDISco motifs, up to 1000,000 seqlets, Leiden clustering = 50
modisco motifs \
    -s /path/to/scratch/all_tss_windows_reference_seq_onehot.npz \
    -a /path/to/scratch/deepshap_scores/mean_across_folds_all_${output}.npz \
    -n 1000000 -l 50 -v \
    -o /path/to/scratch/tfmodisco_results/mean_across_folds_all_${output}_modisco.h5

# Download JASPAR motif database (used for motif annotation)
wget https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt \
    -O /path/to/scratch/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
    
# Report results
modisco report \
    -i /path/to/scratch/tfmodisco_results/mean_across_folds_all_${output}_modisco.h5 \
    -o /path/to/scratch/tfmodisco_results/mean_across_folds_all_${output}_modisco/ \
    -t -m /path/to/scratch/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
