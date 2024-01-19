import os
import pandas as pd

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

                # Identify the unique similarity score column (third column)
                similarity_column = df.columns[2]

                # Rename this column to reflect the similarity score type
                df.rename(columns={similarity_column: similarity_column}, inplace=True)

                # If the consolidated DataFrame is empty, initialize it with the first DataFrame
                if consolidated_df.empty:
                    consolidated_df = df
                else:
                    # Merge with the consolidated DataFrame
                    consolidated_df = pd.merge(consolidated_df, df, on=['Process A', 'Process B'], how='outer')

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
