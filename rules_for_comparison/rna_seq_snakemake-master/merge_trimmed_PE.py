if config['end_type'] == "pe":
    rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0],
            two = lambda wildcards: get_trimmed(wildcards.name)[1]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = merged_outdir + "{name}_1.merged.fastq.gz",
            out_two = merged_outdir + "{name}_2.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
            two = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[1])
        shell:
            """
            cat {params.one} > {output.out_one}
            cat {params.two} > {output.out_two}
            """