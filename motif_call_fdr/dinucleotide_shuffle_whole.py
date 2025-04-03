import pyfastx
import tqdm
from clipnet.clipnet import shuffle

seqs = pyfastx.Fasta("attr/lcl/quantity_annotated_seqlets_SPKLF_1kb.fa")


shuffled_seq = [
    shuffle.kshuffle(rec.seq, random_seed=0)[0]
    for rec in tqdm.tqdm(seqs, desc="Shuffling sequences")
]
shuffled_center = [
    rec.seq[: (len(rec) - 30) // 2]
    + shuffled_seq[i][42:72]
    + rec.seq[(len(rec) + 30) // 2 :]
    for i, rec in enumerate(seqs)
]

with open("attr/lcl/quantity_annotated_seqlets_SPKLF_1kb_shuf.fa", "a") as f:
    for i, seq in enumerate(shuffled_center):
        f.write(f">{seqs[i].name}\n")
        f.write(f"{seq}\n")
    

