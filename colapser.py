import os
import pandas as pd

def generate_unique_comparisons(folder_path):
    unique_comparisons = set()

    # Traverse the folder to find CSV files
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                df = pd.read_csv(file_path)

                # Iterate through each row to extract 'Process A' and 'Process B'
                for index, row in df.iterrows():
                    process_a = row['Process A']
                    process_b = row['Process B']

                    # Add the comparison in a sorted manner to ensure uniqueness
                    comparison = tuple(sorted([process_a, process_b]))
                    unique_comparisons.add(comparison)

    return unique_comparisons

# Replace 'path_to_folder' with the path to your folder containing CSV files
folder_path = 'Metrics'
unique_comparisons = generate_unique_comparisons(folder_path)

# Convert the set to a DataFrame
comparison_df = pd.DataFrame(unique_comparisons, columns=['Process A', 'Process B'])

# Save the DataFrame as a CSV file
comparison_df.to_csv('unique_comparisons.csv', index=False)

print("Unique comparisons saved as 'unique_comparisons.csv'")
