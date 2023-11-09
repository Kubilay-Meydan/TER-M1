# TER-M1


# Snakemake rules Extractor

This Python script extracts module names from Snakemake workflow files. It searches through all the files in a specified local repository and extracts words that follow the word "rule" at the beginning of a line, provided that the word is followed by a colon `:` or ` :`. This is the typical snakemake sythax for defining rules or modules

## How to Use

1. Download the script to your local machine.
2. Open your terminal or command prompt and navigate to the location where you downloaded the script.
3. Run the script
4. When prompted, enter the path to your local repository containing the Snakemake workflow files.

The script will then process the files in the specified repository, workflow by workflow and save a list of txt files for each one containing a list of module names. It will save those txt files in a directory called `Extracted_Rules` it creates in the Workflow directory. So to access your text files containg the rules, you will have to go to `Workflows/Extracted_Rules`

Searches through '.py', '.java', '.c', '.cpp', '.js', '.ts', '.html', '.css', '.smk' file types `OR` files that have no type but are named "Snakefile". (The capitalisation doesn't matter)

## Requirements

- Python 3.x

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
