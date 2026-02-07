import pandas as pd
import freesasa
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
import glob
import os

# === STEP 1: Find all grouped CSV files and corresponding PDBs ===
csv_files = glob.glob("*_grouped_all_residue_interactions.csv")
if not csv_files:
    raise FileNotFoundError("No '_grouped_all_residue_interactions.csv' files found.")

# Initialize a list to collect data for CSV output
all_surface_data = []

# === STEP 2: Process each CSV file ===
for csv_file in csv_files:
    prefix = csv_file.replace("_grouped_all_residue_interactions.csv", "")
    pdb_file = f"{prefix}.pdb"

    if not os.path.exists(pdb_file):
        print(f"Warning: Corresponding PDB file {pdb_file} not found for {csv_file}")
        continue

    # === STEP 3: Load interaction data and structure ===
    df = pd.read_csv(csv_file).dropna(subset=["Atom1_resi_num", "Atom2_resi_num"])
    structure = PDBParser(QUIET=True).get_structure(prefix, pdb_file)

    # === STEP 4: Interaction color mapping ===
    interaction_colors = {
        "A-B": "custom_red",
        "B-C": "custom_blue",
        "A-C": "custom_yellow"
    }

    # === STEP 5: Extract interacting residues ===
    interacting_residues = {key: set() for key in interaction_colors}
    for _, row in df.iterrows():
        try:
            chain1, resi1 = row["Atom1_chain"], int(row["Atom1_resi_num"])
            chain2, resi2 = row["Atom2_chain"], int(row["Atom2_resi_num"])
        except ValueError:
            continue
        interaction_type = row["Interaction_Group"]
        if interaction_type in interaction_colors:
            interacting_residues[interaction_type].add((chain1, resi1))
            interacting_residues[interaction_type].add((chain2, resi2))

    # === STEP 6: Helper functions ===
    class ChainSelect(Select):
        def __init__(self, chains): self.chains = chains
        def accept_chain(self, chain): return chain.id in self.chains

    def save_subset(structure, chains, filename):
        io = PDBIO()
        io.set_structure(structure)
        io.save(filename, select=ChainSelect(chains))

    def compute_centroid(coords):
        return [round(float(np.mean(axis)), 3) for axis in zip(*coords)] if coords else [0.0, 0.0, 0.0]

    def normalize(v):
        norm = np.linalg.norm(v)
        return v / norm if norm else v

    def compute_offset_centroid(original, reference, distance=28.0):
        direction_unit = normalize(np.array(original) - np.array(reference))
        offset = np.array(original) + distance * direction_unit
        return [round(c, 3) for c in offset]

    # === STEP 7: Interface surface area calculation ===
    interaction_areas = {}
    for key in interaction_colors:
        c1, c2 = key.split("-")
        save_subset(structure, [c1], f"{prefix}_{c1}.pdb")
        save_subset(structure, [c2], f"{prefix}_{c2}.pdb")
        save_subset(structure, [c1, c2], f"{prefix}_{c1}_{c2}_combined.pdb")

        area1 = freesasa.calc(freesasa.Structure(f"{prefix}_{c1}.pdb")).totalArea()
        area2 = freesasa.calc(freesasa.Structure(f"{prefix}_{c2}.pdb")).totalArea()
        combined = freesasa.calc(freesasa.Structure(f"{prefix}_{c1}_{c2}_combined.pdb")).totalArea()
        interface = area1 + area2 - combined
        interaction_areas[key] = round(interface, 1)

        # Collect for CSV
        all_surface_data.append({
            "Prefix": prefix,
            "Interaction": key,
            "Surface_Area": round(interface, 1)
        })

    # === STEP 8: Global CA centroid ===
    all_coords = [
        atom.coord for model in structure
        for chain in model
        for res in chain
        if "CA" in res
        for atom in [res["CA"]]
    ]
    global_center = compute_centroid(all_coords)

    # === STEP 9: Generate PyMOL script ===
    pml_filename = f"{prefix}_interactions_surface_SASA.pml"
    with open(pml_filename, "w") as f:
        f.write("reinitialize\n")
        f.write(f"load {pdb_file}\n\n")
        f.write("hide everything\nshow cartoon\n")

        # Define custom colors
        for color_name in interaction_colors.values():
            rgb = {
                "custom_red": "[1.0, 0.0, 0.0]",
                "custom_blue": "[0.0, 0.0, 1.0]",
                "custom_yellow": "[1.0, 1.0, 0.0]",
            }[color_name]
            f.write(f"set_color {color_name} = {rgb}\n")
        f.write("\n")

        # Highlight residues
        for key, residues in interacting_residues.items():
            color = interaction_colors[key]
            selection = " or ".join([f"(chain {ch} and resi {res})" for ch, res in residues])
            selection_name = f"interacting_{key.replace('-', '_')}"
            f.write(f"select {selection_name}, {selection}\n")
            f.write(f"color {color}, {selection_name}\n")
            f.write(f"show surface, {selection_name}\n\n")

        f.write("set transparency, 0.2\nzoom all\n\n")

        # Add offset pseudoatoms
        for key, area in interaction_areas.items():
            coords = [
                structure[0][chain][res]["CA"].coord
                for chain, res in interacting_residues[key]
                if res in structure[0][chain] and "CA" in structure[0][chain][res]
            ]
            centroid = compute_centroid(coords)
            offset_centroid = compute_offset_centroid(centroid, global_center, 28.0)
            label = f"{key} ({interaction_colors[key].replace('custom_', '')}) area: {area} Å²"
            label_obj = f"label_{key.replace('-', '_')}"
            f.write(f"pseudoatom {label_obj}, pos={offset_centroid}, label=\"{label}\"\n")
            f.write(f"set label_size, 18, {label_obj}\n")
            f.write(f"set label_color, black, {label_obj}\n")
            f.write(f"color {interaction_colors[key]}, {label_obj}\n\n")

    print(f"PyMOL script '{pml_filename}' created for {prefix}.")

# === STEP 10: Saving all surface area results to a CSV file ===
surface_area_df = pd.DataFrame(all_surface_data)
surface_area_df.to_csv("all_surface_areas.csv", index=False)
print("All surface areas from all files have been saved to 'all_surface_areas.csv'.")
