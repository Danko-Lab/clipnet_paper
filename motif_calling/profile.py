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


def load_motif_ppms(modisco_results_path, include=None, pseudocount=1):
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
                        cwm_dict[patterns_group_name][pattern] = ppm
    return cwm_dict


cwms, cutoffs = load_motif_cwms("all_tss_windows_reference_seq_profile_shap.modisco.h5")
ppms = load_motif_ppms("all_tss_windows_reference_seq_profile_shap.modisco.h5")


motif_idx = {
    "CA-Inr-0": 0,
    "SP/KLF-0": 1,
    "TA-Inr-0": 2,
    "ETS-0": 3,
    "TBP-0": 4,
    "DPR-0": 5,
    "NFY": 6,
    "DPR-1": 7,
    "DPR-2": 8,
    "DPR-3": 9,
    "DPR-4": 10,
    "TBP-1": 11,
    "DPR-5": 12,
    "DPR-6": 13,
    "DPR-7": 14,
    "YY1-0": 15,
    "DPR-8": 16,
    "IRF-0": 17,
    "NRF1-0": 18,
    "CA-Inr-1": 19,
    "POU-0": 21,
    "IRF-1": 22,
    "USF2": 23,
    "CA-Inr-2": 24,
    "CREB": 25,
    "CTCF": 26,
    "ATF-0": 27,
    "ATF-1": 28,
    "GA-Inr-0": 29,
    "TBP-2": 30,
    "ZBTB33": 32,
    "ETS-1": 35,
    "TA-Inr-1": 37,
    "ETS-2": 38,
    "TBP-3": 39,
    "GA-Inr-1": 40,
    "RUNX-0": 41,
    "CA-Inr-3": 42,
    "RUNX-1": 43,
    "ETS-3": 44,
    "TCT-Inr": 45,
    "XBP": 46,
    "POU-1": 47,
    "CA-Inr-4": 48,
    "GA-Inr-2": 49,
    "TBP-4": 50,
    "TBP-5": 51,
    "ETS-4": 52,
    "YY1-1": 53,
    "CA-Inr-5": 54,
    "CA-Inr-6": 55,
    "CA-Inr-7": 56,
    "JUN": 57,
    "CA-Inr-8": 58,
    "CA-Inr-9": 59,
    "NRF1-1": 60,
    "SP/KLF-1": 61,
    "NRF1-2": 62,
    "TA-Inr-2": 64,
}

motif_idx_main = {
    "CA-Inr-0": 0,
    "SP/KLF-0": 1,
    "TA-Inr-0": 2,
    "ETS-0": 3,
    "TBP-0": 4,
    "DPR-0": 5,
    "NFY": 6,
    "DPR-1": 7,
    "DPR-2": 8,
    "DPR-3": 9,
    "DPR-4": 10,
    "TBP-1": 11,
    "DPR-5": 12,
    "DPR-6": 13,
    "DPR-7": 14,
    "YY1-0": 15,
    "DPR-8": 16,
    "IRF-0": 17,
    "NRF1-0": 18,
    "CA-Inr-1": 19,
    "IRF-1": 22,
    "TBP-2": 30,
}

motif_idx_main = {
    "CA-Inr": 0,
    "SP/KLF": 1,
    "TA-Inr": 2,
    "ETS": 3,
    "TBP": 4,
    "NFY": 6,
    "YY1": 15,
    "IRF": 17,
    "NRF1": 18,
}

motifs = {k: cwms["pos_patterns"][f"pattern_{motif_idx[k]}"] for k in motif_idx.keys()}
motif_cutoffs = {
    k: cutoffs["pos_patterns"][f"pattern_{motif_idx[k]}"] for k in motif_idx.keys()
}

attr_raw = np.load("all_tss_windows_reference_seq_profile_shap.npz")["arr_0"]
attr = attr_raw / attr_raw.sum(axis=(1, 2), keepdims=True)
# X_scale = np.sign(X) * (np.abs(X)) ** (1 / 3)
seqlets = recursive_seqlets(attr.sum(axis=1), threshold=0.05)

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

seqlet_attr_raw = []
for i, row in tqdm.tqdm(seqlets.iterrows()):
    seqlet_attr_raw.append(attr_raw[row["example_idx"], :, row["start"] : row["end"]].sum())

seqlets["seqlet_attr_raw"] = seqlet_attr_raw

filtered_seqlets = [
    seqlets[
        (seqlets.motif_p <= 0.05)
        & (seqlets.seqlet_attr_raw >= motif_cutoffs[name] * 0.5)
        & (seqlets.motif == name)
    ]
    for name in tqdm.tqdm(list(motif_idx.keys()))
]

filtered_seqlets = pd.concat(filtered_seqlets).reset_index(drop=True)

X = np.load("all_tss_windows_reference_seq_ohe.npz")["arr_0"]
sequences = []
for i, row in tqdm.tqdm(filtered_seqlets.iterrows()):
    alphabet_idx = np.abs(X[row["example_idx"], :, row["start"] : row["end"]]).argmax(
        axis=0
    )
    sequences.append("".join(["A", "C", "G", "T"][a] for a in alphabet_idx))

filtered_seqlets["sequence"] = sequences
