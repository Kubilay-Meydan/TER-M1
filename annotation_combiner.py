import pandas as pd
from itertools import product
import os

def extract_process_data(df):
    """
    Extract process names and values from a dataframe.
    """
    processes = df.columns[1:]  # Assuming first column is unnamed and contains process names
    values = df.set_index(df.columns[0]).to_dict()
    return processes, values

def process_csv_files(root_folder_path):
    """
    Process all CSV files in the given folder and its subfolders.
    """
    all_processes = []
    file_process_mapping = {}
    for root, dirs, files in os.walk(root_folder_path):
        for file_name in files:
            if file_name.endswith('.csv'):
                file_path = os.path.join(root, file_name)
                df = pd.read_csv(file_path)
                processes, values = extract_process_data(df)
                all_processes.extend([(process, file_path) for process in processes])
                file_process_mapping[file_path] = values

    # Generating all combinations of processes
    combinations = list(product(all_processes, repeat=2))

    # Preparing data for the new dataframe
    data = []
    for (process_a, file_a), (process_b, file_b) in combinations:
        # Check if both processes are from the same file
        if file_a == file_b:
            annotation = file_process_mapping[file_a][process_b].get(process_a, 'N/A')
        else:
            annotation = 'D'  # Different files

        # Skip comparisons where the value is 'I'
        if annotation == 'I':
            continue

        # Remove '.py' extension from process names
        process_a = process_a.replace('.py', '')
        process_b = process_b.replace('.py', '')

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
