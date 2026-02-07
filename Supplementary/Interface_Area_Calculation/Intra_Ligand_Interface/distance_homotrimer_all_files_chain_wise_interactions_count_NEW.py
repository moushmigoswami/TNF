import os
import glob
from Bio import PDB
import numpy as np
import pandas as pd

pdb_parser = PDB.PDBParser(QUIET=True)

def calculate_interactions(structure, clash_cutoff=3.5):
    """Compute atomic interactions within Chain A, B, and C of a PDB structure."""
    chains_of_interest = {"A", "B", "C"}
    atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id in chains_of_interest]
    
    coords = np.array([atom.coord for atom in atoms], dtype="d")
    kdt = PDB.kdtrees.KDTree(coords)
    clashes = []
    interaction_counts = {"A_B": 0, "B_C": 0, "A_C": 0}

    for i, atom_1 in enumerate(atoms):
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), clash_cutoff)
        for kdt_atom in kdt_search:
            j, atom_distance = kdt_atom.index, kdt_atom.radius
            atom_2 = atoms[j]
            if atom_1.parent.id == atom_2.parent.id:
                continue
            
            chain_1 = atom_1.get_parent().get_parent().id
            chain_2 = atom_2.get_parent().get_parent().id
            
            if chain_1 != chain_2:
                clashes.append([
                    atom_distance,
                    chain_1, atom_1.get_parent().id[1], atom_1.get_parent().get_resname(), atom_1.name,
                    chain_2, atom_2.get_parent().id[1], atom_2.get_parent().get_resname(), atom_2.name
                ])
                if (chain_1, chain_2) in [("A", "B"), ("B", "A")]:
                    interaction_counts["A_B"] += 1
                elif (chain_1, chain_2) in [("B", "C"), ("C", "B")]:
                    interaction_counts["B_C"] += 1
                elif (chain_1, chain_2) in [("A", "C"), ("C", "A")]:
                    interaction_counts["A_C"] += 1

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

        summary_df = pd.DataFrame([
            ["Total A-B Interactions", interaction_counts["A_B"], "", "", "", "", "", "", ""],
            ["Total B-C Interactions", interaction_counts["B_C"], "", "", "", "", "", "", ""],
            ["Total A-C Interactions", interaction_counts["A_C"], "", "", "", "", "", "", ""]
        ], columns=df.columns)

        df = pd.concat([df, summary_df], ignore_index=True)

        output_file = f"{structure_id}_A_B_C_interactions.csv"
        df.to_csv(output_file, index=False)

        print(f"\nProcessed: {pdb_file}")
        print(f"Saved output as: {output_file}")
        print(f"Total inter-chain interactions: {len(df) - 3}")
        print(f"A-B: {interaction_counts['A_B']}, B-C: {interaction_counts['B_C']}, A-C: {interaction_counts['A_C']}")

