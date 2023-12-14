rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = merged_outdir + "{name}_1.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
        shell:
            """
            cat {params.one} > {output.out_one}
            """
