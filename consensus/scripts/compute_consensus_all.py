from collections import Counter
from Bio import AlignIO
from Bio.Align import substitution_matrices

# amino acids
aas = list("ACDEFGHIKLMNPQRSTVWY")

# chemical groups
chem_groups = {
"A":"hydrophobic","V":"hydrophobic","I":"hydrophobic","L":"hydrophobic","M":"hydrophobic",
"F":"aromatic","W":"aromatic","Y":"aromatic",
"S":"polar","T":"polar","N":"polar","Q":"polar",
"K":"positive","R":"positive","H":"positive",
"D":"negative","E":"negative",
"G":"special","P":"special","C":"special"
}

blosum = substitution_matrices.load("BLOSUM62")


def compute_consensus(alignment):

    length = alignment.get_alignment_length()

    raw_consensus = ""
    chem_consensus = ""
    blosum_consensus = ""

    for i in range(length):

        column = [record.seq[i] for record in alignment if record.seq[i] != "-"]

        counts = Counter(column)

        # ---------- raw ----------
        raw_res = counts.most_common(1)[0][0]
        raw_consensus += raw_res

        # ---------- chemical ----------
        group_counts = Counter()

        for aa, count in counts.items():
            if aa in chem_groups:
                group_counts[chem_groups[aa]] += count

        best_group = group_counts.most_common(1)[0][0]

        group_res = [aa for aa in column if chem_groups.get(aa) == best_group]
        chem_res = Counter(group_res).most_common(1)[0][0]

        chem_consensus += chem_res

        # ---------- BLOSUM ----------
        best_score = -1e9
        best_aa = None

        for candidate in aas:

            score = 0
            for aa, count in counts.items():
                score += blosum[candidate, aa] * count

            if score > best_score:
                best_score = score
                best_aa = candidate

        blosum_consensus += best_aa

    return raw_consensus, chem_consensus, blosum_consensus


def process_alignment(file, prefix):

    print("Processing:", file)

    alignment = AlignIO.read(file, "fasta")

    raw, chem, blosum = compute_consensus(alignment)

    with open(f"{prefix}_raw_consensus.fasta","w") as f:
        f.write(f">{prefix}_raw_consensus\n{raw}\n")

    with open(f"{prefix}_chemical_consensus.fasta","w") as f:
        f.write(f">{prefix}_chemical_consensus\n{chem}\n")

    with open(f"{prefix}_blosum_consensus.fasta","w") as f:
        f.write(f">{prefix}_blosum_consensus\n{blosum}\n")


# ----------------------------
# run both datasets
# ----------------------------

process_alignment("rbcs_aligned.fasta", "rbcs")
process_alignment("rbcl_aligned.fasta", "rbcl")

print("Done.")
