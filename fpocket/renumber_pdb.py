#!/usr/bin/env python3

import sys
import os
import string

if len(sys.argv) != 2:
    print("Usage: renumber_pdb.py structure.pdb")
    sys.exit(1)

input_pdb = sys.argv[1]
base = os.path.splitext(input_pdb)[0]
output_pdb = base + "_renumbered.pdb"

chain_letters = list(string.ascii_uppercase)

current_chain = 0
res_counter = 1
last_res_key = None

with open(input_pdb) as f, open(output_pdb, "w") as out:

    for line in f:

        if line.startswith(("ATOM", "HETATM")):

            # unique residue identifier (original chain + residue + insertion)
            res_key = line[21:27]

            if res_key != last_res_key:
                last_res_key = res_key
                new_res = res_counter
                res_counter += 1

            new_chain = chain_letters[current_chain]

            # rebuild line while clearing segID (columns 73-76)
            new_line = (
                line[:21] +                # up to chain column
                new_chain +                # new chain ID
                f"{new_res:4d}" +          # new residue number
                line[26:72] +              # rest of line
                "    " +                   # blank segID
                line[76:]                  # element / charge
            )

            out.write(new_line)

        elif line.startswith("TER"):

            out.write(line)
            current_chain += 1

        else:
            out.write(line)

print(f"Wrote cleaned file: {output_pdb}")
