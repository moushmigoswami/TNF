import pandas as pd
import glob

# Define valid interaction pairs
valid_pairs = {
    ("A", "D"), ("A", "E"), ("A", "F"),
    ("B", "D"), ("B", "E"), ("B", "F"),
    ("C", "D"), ("C", "E"), ("C", "F")
}

# Get all CSV files ending with _interactions.csv in the current directory
csv_files = glob.glob("*ABC_DEF_interactions.csv")

for csv_file in csv_files:
    base_name = csv_file.replace("_interactions.csv", "")
    output_file = f"{base_name}_grouped_all_residue_interactions.csv"

    df = pd.read_csv(csv_file)

    df["Atom1_resi_num"] = pd.to_numeric(df["Atom1_resi_num"], errors="coerce").fillna(0).astype(int)
    df["Atom2_resi_num"] = pd.to_numeric(df["Atom2_resi_num"], errors="coerce").fillna(0).astype(int)

    grouped_rows = {
        "A-D": [], "A-E": [], "A-F": [],
        "B-D": [], "B-E": [], "B-F": [],
        "C-D": [], "C-E": [], "C-F": []
    }

    for _, row in df.iterrows():
        c1, c2 = str(row["Atom1_chain"]), str(row["Atom2_chain"])
        pair = (c1, c2)
        reverse_pair = (c2, c1)

        for valid_pair in valid_pairs:
            if pair == valid_pair or reverse_pair == valid_pair:
                group_name = f"{valid_pair[0]}-{valid_pair[1]}"
                row["Interaction_Group"] = group_name
                grouped_rows[group_name].append(row)
                break  # Found a matching pair; skip to next row

    # Build final rows with interaction numbers
    final_rows = []
    interaction_counts = {}

    for group, rows in grouped_rows.items():
        interaction_counts[group] = len(rows)
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
        if count > 0
    ])

    # Append summary
    filtered_df = pd.concat([filtered_df, summary_rows], ignore_index=True)
    filtered_df.to_csv(output_file, index=False)

    # Print summary
    print(f"Processed: {csv_file}")
    for group, count in interaction_counts.items():
        if count > 0:
            print(f"  Total {group} interactions: {count}")
    print(f"  Output saved to: {output_file}\n")



