from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
import tqdm
from tangermeme.annotate import annotate_seqlets
from tangermeme.seqlet import recursive_seqlets


def load_motif_cwms(modisco_results_path, include=None, trim_threshold=0.2):
    with h5py.File(modisco_results_path, "r") as f:
        cwm_dict = defaultdict(lambda: dict())
        cutoff_dict = defaultdict(lambda: dict())
        for patterns_group_name in ["pos_patterns", "neg_patterns"]:
            if (include is not None) and (patterns_group_name not in include.keys()):
                continue
            # if the results include pos/neg patterns...
            if patterns_group_name in f.keys():
                new_patterns_grp = f[patterns_group_name]
                # if there are any patterns for this metacluster...
                if len(new_patterns_grp.keys()) > 0:
                    pattern_names = list(new_patterns_grp.keys())
                    pattern_names = sorted(
                        pattern_names, key=lambda name: int(name.split("_")[1])
                    )
                    # for each hit...
                    for pattern in pattern_names:
                        if include is not None:
                            if (
                                int(pattern.replace("pattern_", ""))
                                not in include[patterns_group_name]
                            ):
                                continue
                        pattern_grp = new_patterns_grp[pattern]
                        cwm = pattern_grp["contrib_scores"][:].T
                        score = np.sum(np.abs(cwm), axis=0)
                        trim_thresh = (
                            np.max(score) * trim_threshold
                        )  # Cut off anything less than 20% of max score
                        pass_inds = np.where(score >= trim_thresh)[0]
                        trimmed = cwm[:, np.min(pass_inds) : np.max(pass_inds) + 1]
                        cwm_dict[patterns_group_name][pattern] = trimmed / trimmed.sum()
                        cutoff_dict[patterns_group_name][pattern] = np.percentile(
                            np.abs(
                                pattern_grp["seqlets"]["contrib_scores"][:].sum(
                                    axis=(1, 2)
                                )
                            ),
                            1,
                        )
    return cwm_dict, cutoff_dict


def load_motif_ppms(
    modisco_results_path, include=None, pseudocount=1, trim_threshold=0.2
):
    with h5py.File(modisco_results_path, "r") as f:
        cwm_dict = defaultdict(lambda: dict())
        for patterns_group_name in ["pos_patterns", "neg_patterns"]:
            if (include is not None) and (patterns_group_name not in include.keys()):
                continue
            # if the results include pos/neg patterns...
            if patterns_group_name in f.keys():
                new_patterns_grp = f[patterns_group_name]
                # if there are any patterns for this metacluster...
                if len(new_patterns_grp.keys()) > 0:
                    pattern_names = list(new_patterns_grp.keys())
                    pattern_names = sorted(
                        pattern_names, key=lambda name: int(name.split("_")[1])
                    )
                    # for each hit...
                    for pattern in pattern_names:
                        if include is not None:
                            if (
                                int(pattern.replace("pattern_", ""))
                                not in include[patterns_group_name]
                            ):
                                continue
                        pattern_grp = new_patterns_grp[pattern]
                        n = pattern_grp["seqlets"]["n_seqlets"][:][0]
                        ppm = (pattern_grp["sequence"][:].T * n + pseudocount) / (
                            n + pseudocount
                        )
                        cwm = pattern_grp["contrib_scores"][:].T
                        score = np.sum(np.abs(cwm), axis=0)
                        trim_thresh = (
                            np.max(score) * trim_threshold
                        )  # Cut off anything less than 20% of max score
                        pass_inds = np.where(score >= trim_thresh)[0]
                        trimmed = ppm[:, np.min(pass_inds) : np.max(pass_inds) + 1]
                        cwm_dict[patterns_group_name][pattern] = trimmed
    return cwm_dict


cwms, cutoffs = load_motif_cwms(
    "all_tss_windows_reference_seq_quantity_shap.modisco.h5"
)
# ppms = load_motif_ppms("all_tss_windows_reference_seq_quantity_shap.modisco.h5")

motif_idx = {
    "ETS-0": 0,
    "SP/KLF-0": 1,
    "NFY": 2,
    "NRF1-0": 3,
    "YY1-0": 4,
    "CREB-1": 5,
    "CA-Inr-0": 6,
    "IRF": 8,
    "CA-Inr-1": 9,
    "CTCF-0": 12,
    "TFEC-0": 14,
    "THAP-0": 15,
    "TBP-0": 16,
    "ZBTB33": 19,
    "POU": 22,
    "TBP-1": 24,
    "IRF-1": 26,
    "ZNF76": 28,
    "SP/KLF-1": 30,
    "IRF-2": 31,
    "ZNF76-0": 33,
    "TBP-2": 34,
    "ETS-1": 35,
    "CA-Inr-2": 36,
    "SP/KLF::ETS-0": 37,
    "ETS-2": 39,
    "ETS::CREB": 40,
    "THAP-1": 42,
    "NRF1-1": 44,
    "NFKB2": 45,
    "IRF-3": 46,
    "SP/KLF::ETS-1": 47,
    "CREB": 49,
    "ETS-3": 50,
    "IRF-4": 51,
    "CA-Inr-3": 52,
    "TFEC-1": 53,
    "CA-Inr-4": 54,
    "YY1-1": 55,
    "TBP-3": 56,
    "NRF1-2": 58,
    "ETS-4": 59,
    "THAP-2": 60,
    "SP/KLF": 61,
    "CA-Inr-5": 62,
    "CTCF-1": 63,
}

motif_idx_main = {
    "ETS-0": 0,
    "SP/KLF": 1,
    "NFY": 2,
    "NRF1-0": 3,
    "YY1-0": 4,
    "CREB-0": 5,
    "CA-Inr-0": 6,
    "IRF": 8,
    "CTCF-0": 12,
    "TBP-0": 16,
    "TBP-1": 24,
    "IRF-1": 26,
    "IRF-2": 31,
    "TBP-2": 34,
}

motif_idx_main = {
    "ETS": 0,
    "SP/KLF": 1,
    "NFY": 2,
    "NRF1": 3,
    "YY1": 4,
    "CREB": 5,
    "CA-Inr": 6,
    "IRF": 8,
    "CTCF": 12,
    "TBP": 16,
}

motifs = {
    k: cwms["pos_patterns"][f"pattern_{motif_idx[k]}"]
    for k in motif_idx.keys()
}

motif_cutoffs = {
    k: cutoffs["pos_patterns"][f"pattern_{motif_idx[k]}"]
    for k in motif_idx.keys()
}

attr = np.load("all_tss_windows_reference_seq_quantity_shap.npz")["arr_0"]
seqlets = recursive_seqlets(attr.sum(axis=1), threshold=0.05)

X = np.load("all_tss_windows_reference_seq_ohe.npz")["arr_0"]
ann = annotate_seqlets(attr, seqlets, motifs, n_nearest=3, n_jobs=32)

seqlets["motif_idx"] = ann[0][:, 0]
seqlets["motif_p"] = ann[1][:, 0]
seqlets["motif"] = [list(motifs.keys())[k] for k in seqlets["motif_idx"]]

seqlets.columns = [
    "example_idx",
    "start",
    "end",
    "attribution",
    "seqlet_p",
    "motif_idx",
    "motif_p",
    "motif",
]
sequences = []
for i, row in tqdm.tqdm(seqlets.iterrows()):
    alphabet_idx = np.abs(X[row["example_idx"], :, row["start"] : row["end"]]).argmax(
        axis=0
    )
    sequences.append("".join(["A", "C", "G", "T"][a] for a in alphabet_idx))

seqlets["sequence"] = sequences

filtered_seqlets = [
    seqlets[
        (seqlets.motif_p <= 0.05)
        & (seqlets.attribution >= motif_cutoffs[name] * 0.5)
        & seqlets.motif.str.startswith(name.replace("-0", ""))
    ]
    for name in tqdm.tqdm(list(motif_idx.keys()))
]

filtered_seqlets = pd.concat(filtered_seqlets).reset_index(drop=True)

filtered_seqlets.to_csv("quantity_annotated_seqlets_full.csv.gz", index=False)

windows = pd.read_csv(
    "../../data/lcl/merged_windows_all.bed.gz",
    sep="\t",
    header=None,
    names=["chrom", "start", "end"],
)
filtered_seqlets_bed = pd.DataFrame(
    {
        "chrom": [windows.chrom[i] for i in filtered_seqlets.example_idx],
        "start": [windows.start[i] for i in filtered_seqlets.example_idx]
        + filtered_seqlets.start,
        "end": [windows.start[i] for i in filtered_seqlets.example_idx]
        + filtered_seqlets.end,
        "example_idx": filtered_seqlets.example_idx,
        "attribution": filtered_seqlets.attribution,
        "seqlet_p": filtered_seqlets.seqlet_p,
        "motif": filtered_seqlets.motif,
        "tomtom_p": filtered_seqlets.motif_p,
        "sequence": filtered_seqlets.sequence,
    }
)
filtered_seqlets_bed.to_csv(
    "quantity_annotated_seqlets_full_0.04.bed.gz", sep="\t", index=False, header=None
)
