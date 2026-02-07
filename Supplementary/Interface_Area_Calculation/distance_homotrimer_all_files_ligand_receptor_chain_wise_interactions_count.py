import os
import glob
from Bio import PDB
import numpy as np
import pandas as pd

pdb_parser = PDB.PDBParser(QUIET=True)

def calculate_interactions(structure, clash_cutoff=3.5):
    """Compute atomic interactions between chains A, B, C and D, E, F."""
    group1 = {"A", "B", "C"}
    group2 = {"D", "E", "F"}

    atoms_group1 = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id in group1]
    atoms_group2 = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id in group2]

    coords_group2 = np.array([atom.coord for atom in atoms_group2], dtype="d")
    kdt_group2 = PDB.kdtrees.KDTree(coords_group2)

    clashes = []
    interaction_counts = {}

    for atom_1 in atoms_group1:
        kdt_search = kdt_group2.search(np.array(atom_1.coord, dtype="d"), clash_cutoff)
        for kdt_atom in kdt_search:
            j, atom_distance = kdt_atom.index, kdt_atom.radius
            atom_2 = atoms_group2[j]

            chain_1 = atom_1.get_parent().get_parent().id
            chain_2 = atom_2.get_parent().get_parent().id

            pair_key = f"{chain_1}-{chain_2}"
            interaction_counts[pair_key] = interaction_counts.get(pair_key, 0) + 1

            clashes.append([
                atom_distance,
                chain_1, atom_1.get_parent().id[1], atom_1.get_parent().get_resname(), atom_1.name,
                chain_2, atom_2.get_parent().id[1], atom_2.get_parent().get_resname(), atom_2.name
            ])

    return clashes, interaction_counts

# Process all PDB files in the current directory
pdb_files = glob.glob("*.pdb")
if not pdb_files:
    print("No PDB files found in the current directory.")
else:
    for pdb_file in pdb_files:
        structure_id = os.path.splitext(os.path.basename(pdb_file))[0]
        structure = pdb_parser.get_structure(structure_id, pdb_file)
        clashes, interaction_counts = calculate_interactions(structure, clash_cutoff=3.5)

        df = pd.DataFrame(clashes, columns=[
            'Distance', 
            'Atom1_chain', 'Atom1_resi_num', 'Atom1_resi_name', 'Atom1_atom',
            'Atom2_chain', 'Atom2_resi_num', 'Atom2_resi_name', 'Atom2_atom'
        ])

        # Create a summary from interaction_counts
        summary_rows = []
        for pair_key, count in interaction_counts.items():
            summary_rows.append([f"Total {pair_key} Interactions", count] + [""] * 7)

        summary_df = pd.DataFrame(summary_rows, columns=df.columns)
        df = pd.concat([df, summary_df], ignore_index=True)

        output_file = f"{structure_id}_ABC_DEF_interactions.csv"
        df.to_csv(output_file, index=False)

        print(f"\nProcessed: {pdb_file}")
        print(f"Saved output as: {output_file}")
        print(f"Total inter-group interactions: {len(df) - len(summary_rows)}")

