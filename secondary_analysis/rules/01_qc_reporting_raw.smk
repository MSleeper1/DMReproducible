### Snakemake rules for quality control reporting of raw sequence files (01) ###
# fastqc and fastq_screen are used to generate quality control reports for raw sequence files

### FASTQC ###
# FastQC is a quality control tool for high throughput sequence data. It reads in sequence data in a variety of formats and will create an HTML report to view the results.
# input: raw sequence files (fastq)
# output: fastqc report (html, zip)

# Rule to run fastqc on single-end sequence files
rule fastqc_se:
    input:
        expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])
    
    log:
        "logs/secondary_rules/01_fastqc_se/01_fastqc_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
    
    conda:
        "../../environment_files/fastqc.yaml"
    
    params:
        output_dir = expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    wildcard_constraints:
        layout = "se"
    
    shell:
        """
        mkdir -p {params.output_dir}
        echo "running fastqc on {input}"
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        echo "done"
        """

# Rule to run fastqc on paired-end sequence files
rule fastqc_pe:
    input: 
        r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir = config["data_dir"]),
        r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir = config["data_dir"])
    
    output:
        r1 = expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        r2 = expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])
    
    log:
        "logs/secondary_rules/01_fastqc_pe/01_fastqc_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    params:
        output_dir = expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    wildcard_constraints:
        layout = "pe"
    
    shell: 
        """
        mkdir -p {params.output_dir}
        echo "running fastqc on {input.r1} and {input.r2}"
        fastqc -o {params.output_dir} {input.r1} {input.r2} > {log} 2>&1
        echo "done"
        """

### FASTQ SCREEN ###
# FastQ Screen is a tool to screen a set of sequences in FastQ format against a set of sequence databases so you can see if the composition of the library matches with what you expect.
# input: raw sequence files (fastq)
# output: fastq_screen report (txt, html)

# Rule to run fastq_screen on single-end sequence files
rule fastq_screen_se:
    input:
        conf = expand("{root}/{genomes_dir}/{genome}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        fq_file = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"])

    log:
        "logs/secondary_rules/01_fastq_screen_se/01_fastq_screen_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastq-screen.yaml"

    params:
        out_dir = expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"])
    
    # shadow: 
    #     "shallow"
    
    threads: 6
    
    wildcard_constraints:
        layout = "se"
    
    shell:
        """
        mkdir -p {params.out_dir}
        echo "running fastq_screen on {input.fq_file}"
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.fq_file}
        echo "done"
        """

# Rule to run fastq_screen on paired-end sequence files
rule fastq_screen_pe:
    input:
        conf = expand("{root}/{genomes_dir}/{genome}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir = config["data_dir"]),
        r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir = config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"]),
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_screen.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["txt", "html"])

    log:
        "logs/secondary_rules/01_fastq_screen_pe/01_fastq_screen_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
   
    conda:
        "../../environment_files/fastq-screen.yaml"

    params:
        out_dir = expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"])
    
    threads: 6
    
    # shadow: 
    #     "shallow"
    
    wildcard_constraints:
        layout = "pe"
    
    shell:
        """
        mkdir -p {params.out_dir}
        echo "running fastq_screen on {input.r1} and {input.r2}"
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.r1}
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.r2}
        echo "done"
        """

# ## compile reports for raw reads
# # Note that multiqc is not working properly in this rule. The shell command is not running as expected and failing when multiqc is called. The commands work when run interactively in the terminal..
# # 01 raw read reports
# rule multiqc_compile_reports_01:
#     input:
#         expand("{root}/{rep_dir}/01_fastqc/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
#         expand("{root}/{rep_dir}/01_fastqc/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
#         expand("{root}/{rep_dir}/01_fastq_screen/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", root = config["root"], se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
#         expand("{root}/{rep_dir}/01_fastq_screen/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", root = config["root"], pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")

#     output:
#         input_list = expand("{root}/{rep_dir}/summary/{sample.ref}/01_qc_report_list.txt", root = config["root"], rep_dir=config["reports_dir"]),
#         multiqc_command = expand("{root}/{rep_dir}/summary/{sample.ref}/01_multiqc_command.sh", root = config["root"], rep_dir=config["reports_dir"]),
#         file = expand("{root}/{rep_dir}/summary/{sample.ref}/01_raw_multiqc.html", root = config["root"], rep_dir=config["reports_dir"])
#     log:
#         "logs/secondary_rules/multiqc/01_multiqc.log"
#     conda:
#         "../../environment_files/multiqc.yaml"
#     params:
#         output_filename = "01_raw_multiqc", # Required: do not change without adjusting the rule output
#         outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
#         comment = "Multiqc report for raw reads produced by fastqc and fastq_screen. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
#         report_title = "01 Raw Sequence QC Reports", # Optional: title for multiqc report (can be changed freely)
#         extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
#     shell:
#         """
#         touch {log}
#         echo "starting multiqc rule" > {log}
#         mkdir -p {params.outdir}
#         echo "removing old input list file if it exists" >> {log}
#         if [ -f {output.input_list} ]; then 
#             rm {output.input_list}
#         fi
#         echo "creating new input list file" >> {log}
#         touch {output.input_list}
#         for i in {input}; do 
#             echo $i >> {output.input_list}
#         done
#         echo "creating multiqc command file for reference" >> {log}
#         echo 'multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title "{params.report_title}" --comment "{params.comment}" {params.extra}' > {output.multiqc_command}
#         echo "running multiqc" >> {log}
#         multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title "{params.report_title}" --comment "{params.comment}" {params.extra} >> {log} 2>&1
#         sleep 5
#         echo "done" >> {log}
#         """