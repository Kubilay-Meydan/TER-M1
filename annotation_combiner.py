import pandas as pd
from itertools import combinations as comb
import os

def extract_process_data(df):
    """
    Extract process names and values from a dataframe.
    """
    # Assuming first column is unnamed and contains process names
    processes = df.columns[1:]
    values = {}
    for i, row in df.iterrows():
        process_row = row[0].lower().replace('.py', '')
        for process_col in processes:
            process_col_cleaned = process_col.lower().replace('.py', '')
            values[(process_row, process_col_cleaned)] = row[process_col]
    return [process.lower().replace('.py', '') for process in processes], values

def process_csv_files(root_folder_path):
    """
    Process all CSV files in the given folder and its subfolders.
    """
    all_processes = []
    process_values = {}
    for root, dirs, files in os.walk(root_folder_path):
        for file_name in files:
            if file_name.endswith('.csv'):
                file_path = os.path.join(root, file_name)
                df = pd.read_csv(file_path)
                processes, values = extract_process_data(df)
                all_processes.extend([(process, file_path) for process in processes])
                process_values.update(values)

    # Generating all unique combinations of processes
    combinations = list(comb(all_processes, 2))

    # Preparing data for the new dataframe
    data = []
    for (process_a, file_a), (process_b, file_b) in combinations:
        key = (process_a, process_b) if file_a == file_b else ('D', 'D')
        annotation = process_values.get(key, 'D')  # Default to 'D' if not found

        # Skip comparisons where the value is 'I'
        if annotation == 'I':
            continue

        data.append({
            'Process A': process_a,
            'Process B': process_b,
            'Annotation': annotation
        })

    # Creating the final dataframe
    final_df = pd.DataFrame(data)

    # Saving the final dataframe to a new CSV file
    output_file_path = os.path.join(root_folder_path, 'combined_processes.csv')
    final_df.to_csv(output_file_path, index=False)

    return output_file_path

# Example usage
folder_path = 'Annotation_Kubilay'
output_file = process_csv_files(folder_path)
print(f"Combined processes CSV file created at: {output_file}")
