import pandas as pd
import glob
import os

# Get all CSV files ending with _A_B_C_interactions.csv in the current directory
csv_files = glob.glob("*_A_B_C_interactions.csv")

for csv_file in csv_files:
    # Derive output file name based on input file prefix
    base_name = csv_file.replace("_A_B_C_interactions.csv", "")
    output_file = f"{base_name}_grouped_all_residue_interactions.csv"

    # Load CSV
    df = pd.read_csv(csv_file)

    # Ensure residue numbers are treated as integers
    df["Atom1_resi_num"] = pd.to_numeric(df["Atom1_resi_num"], errors="coerce").fillna(0).astype(int)
    df["Atom2_resi_num"] = pd.to_numeric(df["Atom2_resi_num"], errors="coerce").fillna(0).astype(int)

    # Initialize dictionary to store interactions
    grouped_rows = {"A-B": [], "B-C": [], "A-C": []}

    # Categorize and store all residue-residue interactions
    for _, row in df.iterrows():
        chain_1, chain_2 = str(row["Atom1_chain"]), str(row["Atom2_chain"])

        # Determine interaction type
        if {chain_1, chain_2} == {"A", "B"}:
            interaction_type = "A-B"
        elif {chain_1, chain_2} == {"B", "C"}:
            interaction_type = "B-C"
        elif {chain_1, chain_2} == {"A", "C"}:
            interaction_type = "A-C"
        else:
            continue  # Skip other interactions

        # Store the interaction
        row["Interaction_Group"] = interaction_type
        grouped_rows[interaction_type].append(row)

    # Assign numbering and compile final data
    final_rows = []
    interaction_counts = {}

    for interaction_type, rows in grouped_rows.items():
        interaction_counts[interaction_type] = len(rows)
        for i, row in enumerate(rows, start=1):
            row["Interaction_Number"] = i
            final_rows.append(row)

    filtered_df = pd.DataFrame(final_rows)

    # Create summary rows
    summary_rows = pd.DataFrame([
        {"Interaction_Group": "Summary", "Interaction_Number": "",
         "Atom1_chain": "", "Atom1_resi_num": "",
         "Atom2_chain": "", "Atom2_resi_num": "",
         "Message": f"Total {group}: {count} interactions"}
        for group, count in interaction_counts.items()
    ])

    # Append summary rows
    filtered_df = pd.concat([filtered_df, summary_rows], ignore_index=True)

    # Save output
    filtered_df.to_csv(output_file, index=False)

    # Print summary
    print(f"Processed: {csv_file}")
    for group, count in interaction_counts.items():
        print(f"  Total {group} interactions: {count}")
    print(f"  Output saved to: {output_file}\n")

