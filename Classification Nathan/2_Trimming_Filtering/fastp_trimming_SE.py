rule fastp_trimming:
            input:
            #get the value in the fast1 column
                fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True)
            wildcard_constraints:
                name="|".join(SAMPLE_NAMES),
                unit="|".join(UNITS)
            conda:
                "../env/align.yaml"
            output:
                out_fastqc = fastp_outdir + "{unit}_{name}_1.trimmed.fastq.gz",
                fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
                fastphtml = fastp_outdir + "{unit}_{name}_fastp.html",
            params:
                fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
                fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False),
            #out_fastqc2 = lambda wildcards: return_fastq2_name(wildcards.name,wildcards.unit),
                fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
                fastphtml = fastp_outdir + "{unit}_{name}_fastp.html"

            log:
                os.path.join(config['project_top_level'], "logs", "{unit}_{name}.fastp_stdout.log")

            shell:
                """
                fastp -i {input.fastq_file} \
                -o {output.out_fastqc} \
                --json {output.fastpjson} \
                --html {output.fastphtml} \
                {params.fastp_parameters} \
                2> {log}
                """