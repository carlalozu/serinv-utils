import pandas as pd

# Load the CSV files
def load_csv(file_path):
    return pd.read_csv(file_path)

# Compute the difference between two DataFrames
def compute_difference(df1, df2):
    diff_df = df1.copy()
    
    # Iterate over unique combinations of bandwidth and n_offdiags
    for bandwidth_value in df1['bandwidth'].unique():
        for n_offdiags_value in df1['n_offdiags'].unique():
            # Filter the rows matching the combination in both DataFrames
            mask1 = (df1['bandwidth'] == bandwidth_value) & (df1['n_offdiags'] == n_offdiags_value)
            mask2 = (df2['bandwidth'] == bandwidth_value) & (df2['n_offdiags'] == n_offdiags_value)

            # Compute mean values from the first file
            mean_time = df2.loc[mask2, 'time'].mean()
            mean_flops = df2.loc[mask2, 'FLOPS'].mean()

            # Subtract the mean from the second file for time and FLOPS
            diff_df.loc[mask1, 'time'] = df1.loc[mask1, 'time'] - mean_time
            diff_df.loc[mask1, 'FLOPS'] = df1.loc[mask1, 'FLOPS'] - mean_flops

    return diff_df

# Main function to load files, compute differences, and save the result
def main(file1, file2, output_file):
    df1 = load_csv(file1)
    df2 = load_csv(file2)

    diff_df = compute_difference(df1, df2)
    diff_df.to_csv(output_file, index=False)
    print(f"Difference file saved as {output_file}")

# Example usage
if __name__ == "__main__":

    cluster = 'fritz'
    arrowhead = 64
    file1 = f"../jobs/{cluster}/results/pobtasi_{arrowhead}.txt"
    file2 = f"../jobs/{cluster}/results/pobtaf_{arrowhead}.txt"
    output_file = f"../jobs/{cluster}/results/pobtasi_d_{arrowhead}.txt"
    main(file1, file2, output_file)
