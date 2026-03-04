import pandas as pd
import glob
import os
from collections import defaultdict

# Get all Excel file paths
file_paths = glob.glob("/Volumes/Extreme SSD/NGS8/30-1181097629/00_fastq/10m/*.xlsx")

# Dictionary mapping value → set of filenames it appears in
value_to_files = defaultdict(set)

# Go through each file
for path in file_paths:
    filename = os.path.basename(path)
    try:
        # Only read third column (Column C)
        df = pd.read_excel(path, engine="openpyxl", usecols=[2])  # Column C is index 2
        col_c = df.iloc[:, 0].dropna()

        for val in col_c:
            value_to_files[val].add(filename)

    except Exception as e:
        print(f"Error reading {filename}: {e}")

# Filter to values that appear in more than one file
shared_values = {val: files for val, files in value_to_files.items() if len(files) > 1}

# Print to terminal
print("\nShared values in Column C (found in multiple files):")
for val, files in shared_values.items():
    print(f"Value: {val} — Found in: {', '.join(files)}")

# Save to CSV
output_data = [
    {"Value": val, "Files": ", ".join(sorted(files))}
    for val, files in shared_values.items()
]
output_df = pd.DataFrame(output_data)
output_df.to_csv("shared_values_column_c.csv", index=False)

print("\n✅ Results saved to 'shared_values_column_c.csv'")
