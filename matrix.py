import os
import pandas as pd
import numpy as np

def custom_aggregate(series):
    # Check if the series contains numeric data
    if pd.api.types.is_numeric_dtype(series):
        # Use max for numeric data, ignoring NaNs
        return series.max(skipna=True)
    else:
        # For non-numeric data, return the first non-null, non-zero value
        non_null_values = series[series.notna() & (series != 0)]
        if not non_null_values.empty:
            return non_null_values.iloc[0]
        return 0

def normalize_process_order(df):
    # Convert 'Process A' and 'Process B' to lowercase
    df['Process A'] = df['Process A'].str.lower()
    df['Process B'] = df['Process B'].str.lower()

    # Ensure 'Process A' and 'Process B' are in a consistent order
    df[['Process A', 'Process B']] = pd.DataFrame(np.sort(df[['Process A', 'Process B']], axis=1), index=df.index)

    # Remove rows where 'Process A' or 'Process B' are 'Process A', 'Process B', or 'Annotation'
    invalid_entries = ['process a', 'process b', 'annotation']
    df = df[~df['Process A'].isin(invalid_entries) & ~df['Process B'].isin(invalid_entries)]

    return df

def consolidate_csv_files(folder_path):
    # Initialize an empty DataFrame for consolidation
    consolidated_df = pd.DataFrame()

    # Traverse the folder to find CSV files
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                
                # Read the CSV file
                df = pd.read_csv(file_path)

                # Normalize the order of 'Process A' and 'Process B' and remove invalid rows
                df = normalize_process_order(df)

                # Identify the unique similarity score column (third column)
                similarity_column = df.columns[2]

                # Rename this column to reflect the similarity score type
                df.rename(columns={similarity_column: similarity_column}, inplace=True)

                # Append to the consolidated DataFrame
                consolidated_df = pd.concat([consolidated_df, df])

    # Group by 'Process A' and 'Process B', and apply custom aggregation
    consolidated_df = consolidated_df.groupby(['Process A', 'Process B']).agg(custom_aggregate).reset_index()

    # Fill missing values with 0
    consolidated_df.fillna(0, inplace=True)

    return consolidated_df

# Replace 'path_to_folder' with the path to your folder containing CSV files
folder_path = 'Metrics'
consolidated_df = consolidate_csv_files(folder_path)

# Save the consolidated DataFrame to a CSV file named 'BIG.csv'
output_path = os.path.join(folder_path, 'BIG.csv')
consolidated_df.to_csv(output_path, index=False)

print(f"Consolidated CSV created at: {output_path}")
