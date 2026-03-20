#!/usr/bin/env python3

import sys

input_file = "pocket348_vert.pqr"
output_file = "pocket348_gly.pdb"

radius_cutoff = 4.5   # keep only large spheres (optional)

atoms = []

with open(input_file) as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue

        fields = line.split()

        x = float(fields[5])
        y = float(fields[6])
        z = float(fields[7])
        radius = float(fields[9])

        if radius < radius_cutoff:
            continue

        atoms.append((x,y,z))

with open(output_file,"w") as out:

    for i,(x,y,z) in enumerate(atoms,1):

        pdb_line = (
            f"ATOM  {i:5d}  CA  GLY A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00           C\n"
        )

        out.write(pdb_line)

    out.write("END\n")

print(f"Wrote {len(atoms)} glycine CA atoms → {output_file}")
