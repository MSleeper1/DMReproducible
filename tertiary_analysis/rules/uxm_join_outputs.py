import pandas as pd
import glob

# Example: Get list of all CSV files to combine, adjust path and pattern as needed
csv_files = glob.glob("uxm_output*.csv")

# List to hold each dataframe
dfs = []

for file in csv_files:
    df = pd.read_csv(file)
    df.set_index('CellType', inplace=True)  # Set CellType as index for easy concat
    dfs.append(df)

# Concatenate along columns (axis=1) to combine samples side-by-side
combined_df = pd.concat(dfs, axis=1)

# Optional: If you want to remove duplicate columns (e.g., same sample in multiple files)
# combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

# Save combined DataFrame as CSV
combined_df.to_csv("uxm_output.all.csv")
