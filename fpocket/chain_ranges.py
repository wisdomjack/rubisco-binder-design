#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    print("Usage: chain_ranges.py structure.pdb")
    sys.exit(1)

pdb = sys.argv[1]

chains = {}

with open(pdb) as f:
    for line in f:

        if not line.startswith(("ATOM", "HETATM")):
            continue

        chain = line[21]
        resnum = int(line[22:26])

        if chain not in chains:
            chains[chain] = set()

        chains[chain].add(resnum)

print()

for chain in sorted(chains):

    residues = sorted(chains[chain])

    start = residues[0]
    end = residues[-1]
    length = len(residues)

    print(f"Chain {chain} : {chain}{start}-{chain}{end} ({length} residues)")

print()
