import os

rna_seq_tools = ["STAR", "Cufflinks", "Bedtools", "STAR-Fusion", "Picard", "RSeQC", "DESeq2",
                 "3D RNA-seq", "BEAVR", "rMATS", "RaNA-Seq", "RAP", "kallisto", 
                 "RNASeqPower", "SeqFold", "EDASeq", "IUTA", "StringTie", "HISAT2", 
                 "Cufflinks", "RSeQC", "iDEP", "scBatch", "PLncDB", "Scotty", 
                 "RnaSeqSampleSize", "MISO", "samtools", "samtools index", "bedtools", 
                 "bedtools genomecov", "cutadapt", "CutAdapt", "fastp", "fastq", "fastqc", 
                 "FastQC", "hisat2", "picard", "MultiQC", "multiqc", "salmon", 
                 "salmon index", "salmon quant", "quant", "Salmon", "samtools flagstat", 
                 "samtools sort"]

def extract_rna_seq_tools_from_shell(file_path):
    found_tools = set()
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            for tool in rna_seq_tools:
                if tool.lower() in line.lower():
                    found_tools.add(tool)
    return found_tools

def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def save_found_tools(extracted_tools_dir, file_name, tools):
    ensure_directory_exists(extracted_tools_dir)
    with open(os.path.join(extracted_tools_dir, f'{file_name}_tools_in_shell.txt'), 'w') as f:
        for tool in tools:
            f.write(f'{tool}\n')

def main(repo_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    extracted_tools_dir = os.path.join(script_dir, 'Extracted_RNA_Seq_Tools')

    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                tools = extract_rna_seq_tools_from_shell(file_path)
                if tools:
                    save_found_tools(extracted_tools_dir, os.path.splitext(file)[0], tools)

    print(f'RNA-seq tools have been saved into the "{extracted_tools_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the directory: ')
    main(repo_path)
