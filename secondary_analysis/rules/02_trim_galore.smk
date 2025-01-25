
### trim_galore rules ###
# trim_galore is a wrapper around cutadapt and fastqc that trims adapters and low quality bases from fastq files and generates fastqc reports for the trimmed files
# input: raw fastq files
# output: trimmed fastq files and fastqc reports for the trimmed files

# trim_galore_se: rule to trim single-end fastq files
rule trim_galore_se:
    input:
        expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", root = config["root"], data_dir = config["data_dir"])
    
    output:
        trimmed_fq = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", root = config["root"], data_dir=config["data_dir"]),
        fastqc_reports = expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["html", "zip"]),
        trim_reports = expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", root = config["root"], rep_dir=config["reports_dir"])
        
    log:
        stdout = "logs/secondary_rules/02_trim_galore_se/02_trim_galore_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "logs/secondary_rules/02_trim_galore_se/02_trim_galore_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    # shadow: 
    #     "shallow"

    conda:
        "../../environment_files/trim_galore.yaml"

    params:
        output_dir = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),  # trimmed files and reports will be saved in this directory
        temp_fastqc_reports = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["html", "zip"]),  # temporary fastqc reports will be saved in this directory by trim galore but moved to a report directory after
        fastqc_rep_dir = expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        temp_trim_reports = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", root = config["root"], data_dir=config["data_dir"]),  # temporary trim reports will be saved in this directory by trim galore but moved to a report directory after
        trim_rep_dir = expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),  # trim reports will be moved to this directory
        user_args = config["prep_args"]["trim_galore_se"]  # user args can be adjusted in the config file

    wildcard_constraints:
        layout="se"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        echo "trimming with: trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input}" >> {log.stdout} 2>> {log.stderr}
        trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        echo "moving fastqc reports to {params.fastqc_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.fastqc_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.fastqc_rep_dir} {params.temp_fastqc_reports} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        echo "moving trim reports to {params.trim_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.trim_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.trim_rep_dir} {params.temp_trim_reports} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        """

# trim_galore_pe: rule to trim paired-end fastq files
rule trim_galore_pe:
    input:
        r1 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", root = config["root"], data_dir=config["data_dir"]),
        r2 = expand("{root}/{data_dir}/01_raw_sequence_files/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", root = config["root"], data_dir=config["data_dir"])

    output:
        trimmed_fq = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_val_{read}.fq", root = config["root"], data_dir=config["data_dir"], read=["1", "2"]),
        fastqc_reports = expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_val_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], read=["1", "2"], suf=["html", "zip"]),
        trim_reports = expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}.fastq_trimming_report.txt", root = config["root"], rep_dir=config["reports_dir"], read=["1", "2"])
        
    log:
        stdout = "logs/secondary_rules/02_trim_galore_pe/02_trim_galore_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "logs/secondary_rules/02_trim_galore_pe/02_trim_galore_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    # shadow: 
    #     "shallow"

    conda:
        "../../environment_files/trim_galore.yaml"

    params:
        output_dir = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),  # trimmed files and reports will be saved in this directory
        temp_fastqc_reports = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_val_{read}_fastqc.{suf}", root = config["root"], data_dir=config["data_dir"], read=["1", "2"], suf=["html", "zip"]),  # temporary fastqc reports will be saved in this directory by trim galore but moved to a report directory after
        fastqc_rep_dir = expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        temp_trim_reports = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}.fastq_trimming_report.txt", root = config["root"], data_dir=config["data_dir"], read=["1", "2"]),  # temporary trim reports will be saved in this directory by trim galore but moved to a report directory after
        trim_rep_dir = expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),  # trim reports will be moved to this directory
        user_args = config["prep_args"]["trim_galore_se"]  # user args can be adjusted in the config file

    wildcard_constraints:
        layout="pe"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        echo "trimming with: trim_galore {params.user_args} --paired --fastqc --output_dir {params.output_dir} {input.r1} {input.r2}" >> {log.stdout} 2>> {log.stderr}
        trim_galore {params.user_args} --paired --fastqc --output_dir {params.output_dir} {input.r1} {input.r2} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        echo "moving fastqc reports to {params.fastqc_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.fastqc_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.fastqc_rep_dir} {params.temp_fastqc_reports} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        echo "moving trim reports to {params.trim_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.trim_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.trim_rep_dir} {params.temp_trim_reports} >> {log.stdout} 2>> {log.stderr}
        echo "done" >> {log.stdout} 2>> {log.stderr}
        """

## compile reports for trimmed reads
# 02 trimmed read reports
rule multiqc_compile_reports_02:
    input:
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), 
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_val_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", root = config["root"], se=sample_info_se.itertuples(), rep_dir=config["reports_dir"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", root = config["root"], pe=sample_info_pe.itertuples(), rep_dir=config["reports_dir"], read=["1", "2"])

    output:
        file = expand("{root}/{rep_dir}/summary/02_trimmed_multiqc.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/02_qc_report_list.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/02_multiqc_command.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/02_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    params:
        output_filename = "02_trimmed_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed reads that includes trimgalore and fastqc report. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "02 Trimmed Sequence QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        touch {log}
        echo "starting multiqc rule" > {log}
        mkdir -p {params.outdir}
        echo "removing old input list file if it exists" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi
        echo "creating new input list file" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done
        echo "creating multiqc command file for reference" >> {log}
        echo 'multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title "{params.report_title}" --comment "{params.comment}" {params.extra}' > {output.multiqc_command}
        echo "running multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title "{params.report_title}" --comment "{params.comment}" {params.extra} >> {log} 2>&1
        sleep 5
        echo "done" >> {log}
        """