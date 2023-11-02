# TER-M1


# Snakemake rules Extractor

This Python script extracts module names from Snakemake workflow files. It searches through all the files in a specified local repository and extracts words that follow the word "rule" at the beginning of a line, provided that the word is followed by a colon `:`.

## How to Use

1. Download the script to your local machine.
2. Open your terminal or command prompt and navigate to the location where you downloaded the script.
3. Run the script
4. When prompted, enter the path to your local repository containing the Snakemake workflow files.

The script will then process the files in the specified repository and print out a list of module names sorted by file, followed by a list of all module names found in all the files and the number of modules found.

Searches through .py', '.java', '.c', '.cpp', '.js', '.ts', '.html', '.css', '.smk' file types

## Requirements

- Python 3.x

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
