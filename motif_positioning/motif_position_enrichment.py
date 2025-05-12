import numpy as np
import polars as pl
import tqdm

# LOAD DATA

# Read in windows around PRO-cap peaks in LCLs
all_bed = pl.read_csv(
    "merged_windows_all.bed.gz",
    separator="\t",
    has_header=False,
    new_columns=["chrom", "start", "end"],
)
# Load the experimental data:
signals = pl.read_csv("all_tss_windows_procap_mean.csv.gz", has_header=False).to_numpy()
signals = signals.to_numpy()[:, np.r_[251:751, 1251:1751]]


# Load the motif calls
profile_motifs = pl.read_csv(
    "profile_annotated_seqlets.bed.gz",
    separator="\t",
    has_header=False,
    new_columns=["chrom", "start", "end", "peak_idx", "motif"],
)

# Print motif counts
profile_motifs = profile_motifs.with_columns(
    (pl.col("motif").str.split("-").list.get(0)).alias("motif_prefix")
)
print(profile_motifs["motif_prefix"].value_counts())

# Load the motif calls
quantity_motifs = pl.read_csv(
    "quantity_annotated_seqlets.bed.gz",
    separator="\t",
    has_header=False,
    new_columns=["chrom", "start", "end", "peak_idx", "motif"],
)

# Print motif counts
quantity_motifs = quantity_motifs.with_columns(
    (pl.col("motif").str.split("-").list.get(0)).alias("motif_prefix")
)
print(quantity_motifs["motif_prefix"].value_counts())

profile_motifs_ = profile_motifs[["chrom", "start", "end", "peak_idx", "motif_prefix"]]
quantity_motifs_ = quantity_motifs[
    ["chrom", "start", "end", "peak_idx", "motif_prefix"]
]
all_motifs = pl.concat([profile_motifs_, quantity_motifs_]).unique()

# Get TSS coordinates
tss = signals.argmax(axis=1)
is_mn = tss >= 500
tss_df = pl.DataFrame({"tss": tss, "is_mn": is_mn})
# Subtract 500 from tss if is_mn
tss_df = tss_df.with_columns(
    pl.when(pl.col("tss") >= 500)
    .then(pl.col("tss") - 500)
    .otherwise(pl.col("tss"))
    .alias("tss")
    + 250
)

# Join relevant info from tss_df and all_bed onto all_motifs using peak_idx
motif_positions = (
    all_motifs.join(
        tss_df.with_row_index("peak_idx").cast({"peak_idx": pl.Int64}), on="peak_idx"
    )
    .join(
        all_bed.with_row_index("peak_idx")
        .select(["peak_idx", "start"])
        .rename({"start": "window_start"})
        .cast(pl.Int64),
        on="peak_idx",
    )
    .with_columns(
        [
            ((pl.col("end") + pl.col("start")) // 2).alias("motif_center"),
            (pl.col("tss") + pl.col("window_start")).alias("tss_coord"),
            pl.when(pl.col("is_mn"))
            .then(
                (pl.col("tss") + pl.col("window_start"))
                - ((pl.col("end") + pl.col("start")) // 2)
            )
            .otherwise(
                ((pl.col("end") + pl.col("start")) // 2)
                - (pl.col("tss") + pl.col("window_start"))
            )
            .alias("motif_position"),
        ]
    )
)

# Extract final result as a list or array
motif_positions = motif_positions["motif_position"].to_numpy()

all_motifs = all_motifs.with_columns(pl.Series("distance_from_tss", motif_positions))

# Generate permutations
permutations = np.array(
    [np.random.permutation(motif_positions) for _ in tqdm.trange(10_000)]
)

# Precompute total counts per motif_prefix
total_counts = all_motifs.group_by("motif_prefix").len().rename({"len": "total_count"})

center_freq_perms = [
    all_motifs.filter((p <= -40) & (p >= -100))
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
    for p in tqdm.tqdm(permutations)
]

tata_freq_perms = [
    all_motifs.filter((p <= -20) & (p >= -30))
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
    for p in tqdm.tqdm(permutations)
]

tss_freq_perms = [
    all_motifs.filter((p <= 5) & (p >= -5))
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
    for p in tqdm.tqdm(permutations)
]

dpe_freq_perms = [
    all_motifs.filter((p <= 30) & (p >= 20))
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
    for p in tqdm.tqdm(permutations)
]

pause_freq_perms = [
    all_motifs.filter((p <= 50) & (p >= 40))
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
    for p in tqdm.tqdm(permutations)
]

# ~~~~~~~~~~~~~~~

center_freq = (
    all_motifs.filter(
        (all_motifs["distance_from_tss"] <= -40)
        & (all_motifs["distance_from_tss"] >= -100)
    )
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
)

tata_freq = (
    all_motifs.filter(
        (all_motifs["distance_from_tss"] <= -20)
        & (all_motifs["distance_from_tss"] >= -30)
    )
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
)

tss_freq = (
    all_motifs.filter(
        (all_motifs["distance_from_tss"] <= 5) & (all_motifs["distance_from_tss"] >= -5)
    )
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
)

dpe_freq = (
    all_motifs.filter(
        (all_motifs["distance_from_tss"] <= 30)
        & (all_motifs["distance_from_tss"] >= 20)
    )
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
)

pause_freq = (
    all_motifs.filter(
        (all_motifs["distance_from_tss"] <= 50)
        & (all_motifs["distance_from_tss"] >= 40)
    )
    .group_by("motif_prefix")
    .len()
    .rename({"len": "filtered_count"})
    .join(total_counts, on="motif_prefix")
    .with_columns((pl.col("filtered_count") / pl.col("total_count")).alias("freq"))
    .select(["motif_prefix", "freq"])
)

# ~~~~~~~~~~~~~~~

center_p = {
    prefix: np.sum(
        [
            perm.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            > center_freq.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            for perm in center_freq_perms
        ]
    )
    for prefix in tqdm.tqdm(center_freq["motif_prefix"])
}

tata_p = {
    prefix: np.sum(
        [
            perm.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            > tata_freq.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            for perm in tata_freq_perms
        ]
    )
    for prefix in tqdm.tqdm(tata_freq["motif_prefix"])
}

tss_p = {
    prefix: np.sum(
        [
            perm.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            > tss_freq.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            for perm in tss_freq_perms
        ]
    )
    for prefix in tqdm.tqdm(tss_freq["motif_prefix"])
}

dpe_p = {
    prefix: np.sum(
        [
            perm.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            > dpe_freq.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            for perm in dpe_freq_perms
        ]
    )
    for prefix in tqdm.tqdm(dpe_freq["motif_prefix"])
}

pause_p = {
    prefix: np.sum(
        [
            perm.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            > pause_freq.filter(pl.col("motif_prefix") == prefix)["freq"][0]
            for perm in pause_freq_perms
        ]
    )
    for prefix in tqdm.tqdm(pause_freq["motif_prefix"])
}


pd.DataFrame({
    "center_p": center_p,
    "tata_p": tata_p,
    "tss_p": tss_p,
    "dpe_p": dpe_p,
    "pause_p": pause_p,
})
