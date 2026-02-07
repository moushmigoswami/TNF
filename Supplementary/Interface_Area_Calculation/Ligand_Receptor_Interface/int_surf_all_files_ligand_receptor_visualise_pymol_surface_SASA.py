import pandas as pd
import freesasa
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.vectors import Vector
import numpy as np
import os

ligand_chains = {"A", "B", "C"}
receptor_chains = {"D", "E", "F"}

available_colors = [
    [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.5, 0.0, 0.5], [1.0, 1.0, 0.0],
    [1.0, 0.5, 0.0], [1.0, 0.75, 0.8], [0.0, 1.0, 1.0], [0.6, 0.2, 0.8], [0.5, 0.5, 0.5]
]
color_names = ["red", "blue", "purple", "yellow", "orange", "pink", "cyan", "violet", "gray"]
interaction_colors = {}

class ChainSelect(Select):
    def __init__(self, chains): self.chains = chains
    def accept_chain(self, chain): return chain.id in self.chains

def save_subset(structure, chains, filename):
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, select=ChainSelect(chains))

def compute_centroid(residue_coords):
    if not residue_coords:
        return [0.0, 0.0, 0.0]
    x = sum(coord[0] for coord in residue_coords) / len(residue_coords)
    y = sum(coord[1] for coord in residue_coords) / len(residue_coords)
    z = sum(coord[2] for coord in residue_coords) / len(residue_coords)
    return [round(x, 3), round(y, 3), round(z, 3)]

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm else v

def compute_offset_centroid(original, reference, distance=28.0):
    direction = np.array(original) - np.array(reference)
    direction_unit = normalize(direction)
    offset_centroid = np.array(original) + distance * direction_unit
    return [round(c, 3) for c in offset_centroid]

# Collect all results for final CSV
all_results = []

csv_files = [f for f in os.listdir('.') if f.endswith('grouped_all_residue_interactions.csv')]

if not csv_files:
    print("No matching CSV files found in the current directory.")
    exit(1)

for interaction_csv_file in csv_files:
    print(f"Processing interaction CSV: {interaction_csv_file}...")

    df = pd.read_csv(interaction_csv_file).dropna(subset=["Atom1_resi_num", "Atom2_resi_num"])
    pdb_id = os.path.splitext(interaction_csv_file)[0].split('_ABC_DEF_grouped_all_residue_interactions')[0]
    pdb_file = f"{pdb_id}.pdb"

    if not os.path.isfile(pdb_file):
        print(f"Skipping {interaction_csv_file}: No corresponding PDB file {pdb_file} found.")
        continue

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)

    interacting_residues = {}
    interaction_areas = {}
    color_idx = 0

    for _, row in df.iterrows():
        try:
            chain1, resi1 = row["Atom1_chain"], int(row["Atom1_resi_num"])
            chain2, resi2 = row["Atom2_chain"], int(row["Atom2_resi_num"])
        except ValueError:
            continue

        if chain1 in ligand_chains and chain2 in receptor_chains:
            key = f"{chain1}-{chain2}"
        elif chain2 in ligand_chains and chain1 in receptor_chains:
            key = f"{chain2}-{chain1}"
        else:
            continue

        if key not in interacting_residues:
            interacting_residues[key] = set()
            color = available_colors[color_idx % len(available_colors)]
            interaction_colors[key] = color_names[color_idx % len(color_names)]
            color_idx += 1

        interacting_residues[key].add((chain1, resi1))
        interacting_residues[key].add((chain2, resi2))

    for key in interacting_residues:
        c1, c2 = key.split("-")
        save_subset(structure, [c1], f"{pdb_id}_{c1}.pdb")
        save_subset(structure, [c2], f"{pdb_id}_{c2}.pdb")
        save_subset(structure, [c1, c2], f"{pdb_id}_{c1}_{c2}_combined.pdb")

        area1 = freesasa.calc(freesasa.Structure(f"{pdb_id}_{c1}.pdb")).totalArea()
        area2 = freesasa.calc(freesasa.Structure(f"{pdb_id}_{c2}.pdb")).totalArea()
        combined = freesasa.calc(freesasa.Structure(f"{pdb_id}_{c1}_{c2}_combined.pdb")).totalArea()
        interface = round(area1 + area2 - combined, 1)
        interaction_areas[key] = interface

        all_results.append({
            "pdb_id": pdb_id,
            "interaction_pair": key,
            "interface_area": interface
        })

    all_coords = [
        atom.coord for model in structure
        for chain in model
        for res in chain
        if "CA" in res
        for atom in [res["CA"]]
    ]
    global_center = compute_centroid(all_coords)

    pymol_script_name = f"{pdb_id}_ligand_receptor_interactions.pml"
    with open(pymol_script_name, "w") as f:
        f.write(f"reinitialize\nload {pdb_file}\n\n")
        f.write("hide everything\nshow cartoon\n")

        for key, colorname in interaction_colors.items():
            rgb = available_colors[color_names.index(colorname)]
            f.write(f"set_color {colorname} = [{rgb[0]}, {rgb[1]}, {rgb[2]}]\n")

        f.write("\n")

        for key, residues in interacting_residues.items():
            color = interaction_colors[key]
            selection_name = f"interacting_{key.replace('-', '_')}"
            selection = " or ".join([f"(chain {ch} and resi {res})" for ch, res in residues])
            f.write(f"select {selection_name}, {selection}\n")
            f.write(f"color {color}, {selection_name}\n")
            f.write(f"show surface, {selection_name}\n\n")

        f.write("set transparency, 0.2\nzoom all\n\n")

        for key, area in interaction_areas.items():
            residues = interacting_residues[key]
            coords = []
            for chain_id, res_id in residues:
                try:
                    res = structure[0][chain_id][res_id]
                    if "CA" in res:
                        coords.append(res["CA"].coord)
                except:
                    continue

            centroid = compute_centroid(coords)
            offset_centroid = compute_offset_centroid(centroid, global_center, distance=28.0)

            color_name = interaction_colors[key]
            label_obj = f"label_{key.replace('-', '_')}"
            label_text = f"{key} ({color_name}) area: {area} Å²"

            f.write(f"pseudoatom {label_obj}, pos={offset_centroid}, label=\"{label_text}\"\n")
            f.write(f"set label_size, 18, {label_obj}\n")
            f.write(f"set label_color, black, {label_obj}\n")
            f.write(f"color {color_name}, {label_obj}\n\n")

    print(f"PyMOL script '{pymol_script_name}' created.\n")

# Save all surface area data to a CSV
results_df = pd.DataFrame(all_results)
results_df.to_csv("all_interface_areas.csv", index=False)
print("All interface areas saved to 'all_interface_areas.csv'")

