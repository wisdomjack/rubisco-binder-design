import csv

input_file = "consensus_expanded_data_set.tsv"

rbcl_out = open("rbcl.fasta", "w")
rbcs_out = open("rbcs.fasta", "w")

with open(input_file) as f:
    reader = csv.reader(f, delimiter="\t")

    header = next(reader)  # skip header

    for row in reader:

        if len(row) < 6:
            continue

        species = row[0].strip()
        common = row[1].strip()

        rbcl_uniprot = row[2].strip()
        rbcl_seq = row[3].strip()

        rbcs_uniprot = row[4].strip()
        rbcs_seq = row[5].strip()

        if rbcl_seq:
            rbcl_out.write(f">{species}|{common}|{rbcl_uniprot}\n")
            rbcl_out.write(f"{rbcl_seq}\n")

        if rbcs_seq:
            rbcs_out.write(f">{species}|{common}|{rbcs_uniprot}\n")
            rbcs_out.write(f"{rbcs_seq}\n")

rbcl_out.close()
rbcs_out.close()

print("Done! Wrote rbcl.fasta and rbcs.fasta")
