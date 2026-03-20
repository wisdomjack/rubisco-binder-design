from Bio import AlignIO
from Bio.Align import substitution_matrices
from collections import Counter
import math

aas = list("ACDEFGHIKLMNPQRSTVWY")

chem_groups = {
"A":"hydrophobic","V":"hydrophobic","I":"hydrophobic","L":"hydrophobic","M":"hydrophobic",
"F":"aromatic","W":"aromatic","Y":"aromatic",
"S":"polar","T":"polar","N":"polar","Q":"polar",
"K":"positive","R":"positive","H":"positive",
"D":"negative","E":"negative",
"G":"special","P":"special","C":"special"
}

blosum_matrix = substitution_matrices.load("BLOSUM62")

GAP_THRESHOLD = 0.7


def shannon_entropy(counts):

    total = sum(counts.values())
    H = 0

    for c in counts.values():
        p = c/total
        H -= p * math.log2(p)

    return H


def process_alignment(file, prefix):

    alignment = AlignIO.read(file,"fasta")

    n_seq = len(alignment)
    length = alignment.get_alignment_length()

    raw_cons = ""
    chem_cons = ""
    blosum_cons = ""

    meta = open(prefix+"_metadata.txt","w")

    meta.write(f"Alignment: {file}\n")
    meta.write(f"Sequences: {n_seq}\n")
    meta.write(f"Length: {length}\n")

    pos_out = 0

    for i in range(length):

        column_all = [rec.seq[i] for rec in alignment]

        gap_count = column_all.count("-")
        gap_fraction = gap_count / len(column_all)

        # skip low coverage columns
        if gap_fraction > GAP_THRESHOLD:

            meta.write(
                f"Pos {i+1:4d} | gap {gap_fraction:.3f} | skipped\n"
            )

            continue

        column = [aa for aa in column_all if aa not in ["-","X"]]

        if len(column) == 0:
            continue

        counts = Counter(column)

        identity = max(counts.values()) / len(column)
        entropy = shannon_entropy(counts)

        # RAW consensus
        raw = counts.most_common(1)[0][0]

        # CHEMICAL consensus
        group_counts = Counter()

        for aa,count in counts.items():
            group_counts[chem_groups.get(aa,"other")] += count

        best_group = group_counts.most_common(1)[0][0]

        group_res = [aa for aa in column if chem_groups.get(aa,"other")==best_group]

        chem = Counter(group_res).most_common(1)[0][0]

        # BLOSUM consensus
        best_score = -1e9
        best_aa = None

        for candidate in aas:

            score = 0

            for aa,count in counts.items():
                score += blosum_matrix[candidate,aa] * count

            if score > best_score:
                best_score = score
                best_aa = candidate

        blosum = best_aa

        raw_cons += raw
        chem_cons += chem
        blosum_cons += blosum

        pos_out += 1

        meta.write(
            f"Pos {i+1:4d} | identity {identity:.3f} | gap {gap_fraction:.3f} | entropy {entropy:.3f} | raw={raw} | chem={chem} | blosum={blosum} | {dict(counts)}\n"
        )

    meta.close()

    with open(prefix+"_consensus_raw.fasta","w") as f:
        f.write(f">{prefix}_raw_consensus\n{raw_cons}\n")

    with open(prefix+"_consensus_chemical.fasta","w") as f:
        f.write(f">{prefix}_chemical_consensus\n{chem_cons}\n")

    with open(prefix+"_consensus_blosum.fasta","w") as f:
        f.write(f">{prefix}_blosum_consensus\n{blosum_cons}\n")


process_alignment("rbcl_aligned.fasta","rbcl")
process_alignment("rbcs_aligned.fasta","rbcs")

print("Done.")
