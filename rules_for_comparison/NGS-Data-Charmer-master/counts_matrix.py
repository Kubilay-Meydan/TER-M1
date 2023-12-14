
# Compile counts for RNA-seq
rule counts_matrix:
    input:
        gene_counts = expand("output/counts/{sample}.gene.counts.txt", sample=SAMPLES),
        exon_counts = expand("output/counts/{sample}.exon.counts.txt", sample=SAMPLES)
    output:
        gene_matrix = "output/counts_genic_matrix.txt",
        exon_matrix = "output/counts_exonic_matrix.txt"
    run:
        import pandas as pd
        import platform

        ## Process the genic counts
        dict_of_counts = {}
        for file in input.gene_counts:
            sample = file.split(".")[0]
            if platform.system() != 'Windows':
                sample = sample.split("/")[2]
            else:
                sample = sample.split("\\")[2]
            
            dict_of_counts[sample] = {}
            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    dict_of_counts[sample][lines[0]] = int(float(lines[-1]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[0], sep='\t')

        ## Process the exonic counts
        dict_of_counts = {}
        for file in input.exon_counts:
            sample = file.split(".")[0]
            if platform.system() != 'Windows':
                sample = sample.split("/")[2]
            else:
                sample = sample.split("\\")[2]
            
            dict_of_counts[sample] = {}
            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    dict_of_counts[sample][lines[0]] = int(float(lines[-1]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[1], sep='\t')