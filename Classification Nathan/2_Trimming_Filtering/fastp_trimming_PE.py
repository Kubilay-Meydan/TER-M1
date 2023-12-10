
if config['end_type'] == "pe":
    rule fastp_trimming:
        input:
        #get the value in the fast1 column
            fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
            fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False)
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES),
            unit="|".join(UNITS)
        conda:
            "../env/align.yaml"
        output:
            out_fastqc = fastp_outdir + "{unit}_{name}_1.trimmed.fastq.gz",
            out_fastqc2 = fastp_outdir + "{unit}_{name}_2.trimmed.fastq.gz",
            fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
            fastphtml = fastp_outdir + "{unit}_{name}_fastp.html",
        params:
            fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
            fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
            fastphtml = fastp_outdir + "{unit}_{name}_fastp.html"

        log:
            os.path.join(config['project_top_level'], "logs", "{unit}_{name}.fastp_stdout.log")

        shell:
            """
            #free -h
            fastp \
            --in1 {input.fastq_file} \
            --in2 {input.fastq_file2} \
            --out1 {output.out_fastqc} \
            --out2 {output.out_fastqc2} \
            --json {output.fastpjson} \
            --html {output.fastphtml} \
            {params.fastp_parameters} \
            2> {log}
            """