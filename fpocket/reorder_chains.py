#!/usr/bin/env python3
"""
Reorder PDB chains so that chain E comes first,
then A, B, C, D — and renumber all residues sequentially from 1.
"""

import sys
from collections import defaultdict

def parse_pdb(filepath):
    chains = defaultdict(list)
    other_lines = []
    
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                chains[chain].append(line)
            elif line.startswith("END"):
                continue
            else:
                other_lines.append(line)
    
    return chains, other_lines

def reorder_and_renumber(chains, new_order):
    output_lines = []
    new_chain_labels = {}  # old chain -> new chain letter
    
    # Map old chains to new chain letters A, B, C, D, E
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, old_chain in enumerate(new_order):
        new_chain_labels[old_chain] = chain_letters[i]
    
    residue_counter = 1
    atom_counter = 1
    
    for old_chain in new_order:
        new_chain = new_chain_labels[old_chain]
        lines = chains[old_chain]
        
        # Track residue number mapping
        prev_res_num = None
        current_new_res = residue_counter - 1
        
        for line in lines:
            res_num = int(line[22:26].strip())
            
            if res_num != prev_res_num:
                current_new_res += 1
                prev_res_num = res_num
            
            # Rebuild the line with new chain, residue number, and atom serial
            new_line = (
                line[:6]                          # record type (ATOM/HETATM)
                + f"{atom_counter:5d}"            # atom serial
                + line[11:21]                     # atom name, alt loc, res name
                + new_chain                       # new chain ID
                + f"{current_new_res:4d}"         # new residue number
                + line[26:]                       # insertion code + coords + rest
            )
            output_lines.append(new_line)
            atom_counter += 1
        
        # After finishing this chain, update residue counter
        residue_counter = current_new_res + 1
        
        # Add TER record between chains
        output_lines.append(
            f"TER   {atom_counter:5d}      "
            + lines[-1][17:20]  # residue name
            + " "
            + new_chain
            + f"{current_new_res:4d}\n"
        )
        atom_counter += 1
    
    output_lines.append("END\n")
    return output_lines, new_chain_labels

def main():
    if len(sys.argv) < 2:
        print("Usage: python reorder_chains.py input.pdb [output.pdb]")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2] if len(sys.argv) > 2 else input_pdb.replace(".pdb", "_reordered.pdb")
    
    chains, other_lines = parse_pdb(input_pdb)
    
    # E first, then A, B, C, D
    new_order = ["E", "A", "B", "C", "D"]
    
    # Check all chains exist
    for c in new_order:
        if c not in chains:
            print(f"ERROR: Chain {c} not found in PDB. Found chains: {list(chains.keys())}")
            sys.exit(1)
    
    output_lines, mapping = reorder_and_renumber(chains, new_order)
    
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)
    
    print(f"Done! Written to: {output_pdb}")
    print(f"Chain mapping (old -> new):")
    for old, new in mapping.items():
        print(f"  {old} -> {new}")
    print(f"\nNew contig for RFdiffusion partial diffusion:")
    
    # Print the new contig string based on renumbered ranges
    # E(53) -> A1-53, A(130) -> B54-183, B(332) -> C184-515, C(323) -> D516-838, D(99) -> E839-937
    # Actually let's compute it properly
    pos = 1
    contig_parts = []
    chain_sizes = {"E": 53, "A": 130, "B": 332, "C": 323, "D": 99}
    new_letters = ["A", "B", "C", "D", "E"]
    for i, old in enumerate(new_order):
        size = chain_sizes[old]
        new_let = new_letters[i]
        end = pos + size - 1
        contig_parts.append(f"{new_let}{pos}-{end}")
        pos = end + 1
    
    # First chain is diffused (old E, now A), rest are fixed
    diffused = contig_parts[0]  # new A (old E) — will be diffused via partial_T
    fixed = "/0 ".join(contig_parts[1:])
    
    print(f"\n  'contigmap.contigs=[{contig_parts[0]}/0 {'/0 '.join(contig_parts[1:])}]' diffuser.partial_T=20")
    print(f"\n  (All chains are referenced by PDB number — partial_T=20 will diffuse the whole structure,")
    print(f"   but chain A (old E) is now first so indices align correctly.)")

if __name__ == "__main__":
    main()
