from Bio import AlignIO
from collections import Counter
import math
import csv


def compute_entropy(aln_file):

    aln = AlignIO.read(aln_file, "fasta")

    n_seq = len(aln)
    n_pos = aln.get_alignment_length()

    print("\nAlignment:", aln_file)
    print("Sequences:", n_seq)
    print("Length:", n_pos)

    results = []

    for i in range(n_pos):

        col = aln[:, i]
        counts = Counter(col)

        # remove gaps if other residues exist
        if "-" in counts and len(counts) > 1:
            del counts["-"]

        total = sum(counts.values())

        # identity
        top_aa, top_count = counts.most_common(1)[0]
        identity = top_count / total

        # Shannon entropy
        H = 0
        for aa, c in counts.items():
            p = c / total
            H -= p * math.log2(p)

        results.append((i+1, identity, H, dict(counts)))

        print(
            f"Pos {i+1:4d} | identity {identity:0.3f} | entropy {H:0.3f} | {dict(counts)}"
        )

    return results


def write_csv(results, outfile):

    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow(["position", "identity", "entropy", "distribution"])

        for r in results:
            writer.writerow(r)


# run for RbcL
rbcl = compute_entropy("rbcl_aligned.fasta")
write_csv(rbcl, "rbcl_entropy_profile.csv")

# run for RbcS
rbcs = compute_entropy("rbcs_aligned.fasta")
write_csv(rbcs, "rbcs_entropy_profile.csv")

print("\nCSV files written:")
print("rbcl_entropy_profile.csv")
print("rbcs_entropy_profile.csv")
