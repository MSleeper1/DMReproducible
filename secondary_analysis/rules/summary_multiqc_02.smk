#--------------------------------------------------------------------------------
# Rules for compiling multiqc reports for raw and trimmed fastq files
#--------------------------------------------------------------------------------

# 01 raw read reports: compile reports for raw reads to look for quality issues from sequencing
rule multiqc_compile_reports_01:
    input:
        expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", root = config["root"], se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", root = config["root"], pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]) # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
    output:
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--01_qc_report_list.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--01_multiqc_command.sh", root = config["root"], rep_dir=config["reports_dir"]),
        file = expand("{root}/{rep_dir}/summary/{{ref}}--01_raw_multiqc.html", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--01_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--01_raw_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for raw reads produced by fastqc and fastq_screen. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "01 Raw Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 2"  # Optional: extra parameters for multiqc (can be changed freely)
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

# 02 trimmed read reports: compile reports for trimmed reads to look for quality issues still present after trimming
rule multiqc_compile_reports_02:
    input:
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), 
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_val_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", root = config["root"], se=sample_info_se.itertuples(), rep_dir=config["reports_dir"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", root = config["root"], pe=sample_info_pe.itertuples(), rep_dir=config["reports_dir"], read=["1", "2"])
    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--02_trimmed_multiqc.html", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--02_qc_report_list.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--02_multiqc_command.sh", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--02_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--02_trimmed_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed reads that includes trimgalore and fastqc report. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "02 Trimmed Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 2"  # Optional: extra parameters for multiqc (can be changed freely)
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

# 01 - 02 summary multiqc report: compile a summary multiqc report for comparison of raw and trimmed reads
rule multiqc_compile_reports_01_02:
    input:
        expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{root}/{rep_dir}/01_fastqc/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", root = config["root"], se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{root}/{rep_dir}/01_fastq_screen/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", root = config["root"], pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), 
        expand("{root}/{rep_dir}/02_fastqc_post_trim/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_val_{read}_fastqc.{suf}", root = config["root"], rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", root = config["root"], se=sample_info_se.itertuples(), rep_dir=config["reports_dir"]), 
        expand("{root}/{rep_dir}/02_trim_galore/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", root = config["root"], pe=sample_info_pe.itertuples(), rep_dir=config["reports_dir"], read=["1", "2"])
    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--01-02_trimmed_multiqc.html", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--01-02_qc_report_list.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--01-02_multiqc_command.sh", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--02_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--01-02_trimmed_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed reads that includes trimgalore and fastqc reports for raw and trimmed reads. This report is useful for comparing the raw and trimmed fastq files. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "01-02 Raw and Trimmed Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 2"  # Optional: extra parameters for multiqc (can be changed freely)
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

