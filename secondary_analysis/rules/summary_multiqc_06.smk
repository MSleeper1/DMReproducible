# This rule compiles all the reports from the 06 post merge step into a single multiqc report
# The reports include featureCounts, fastqc, samtools, mosdepth, and qualimap reports

#------------------------------------------------
# MultiQC rule for 06_post_merge reports
#------------------------------------------------
rule multiqc_compile_bwa_reports_06:
    input:
        expand("{root}/{rep_dir}/06_feature_counts_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bwa
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bwa
        expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),        
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/06_qualimap_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bwa
    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--06_post_merge_multiqc_bwa.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--06_qc_report_list_bwa.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--06_multiqc_command_bwa.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--06_multiqc_bwa.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--06_post_merge_multiqc_bwa", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed and deduplicated reads after merging. Includes reports from fastqc, samtools, mosdepth featurecounts... Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "06 Post-merge Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose"  # Optional: extra parameters for multiqc (can be changed freely)
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

rule multiqc_compile_bis_reports_06:
    input:
        expand("{root}/{rep_dir}/06_feature_counts_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bis
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bis
        expand("{root}/{rep_dir}/06_samtools_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats",  root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bis
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_qualimap_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bis
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis per-base report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis summary report
        expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])
    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--06_post_merge_multiqc_bis.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--06_qc_report_list_bis.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--06_multiqc_command_bis.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--06_multiqc_bis.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--06_post_merge_multiqc_bis", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed and deduplicated reads after merging. Includes reports from fastqc, samtools, mosdepth featurecounts... Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "06 Post-merge Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose"  # Optional: extra parameters for multiqc (can be changed freely)
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

# "scp 'farm:{params.files}' ."
# scp 'farm:{/home/msleeper/scratch/reports/summary/896--06_post_merge_multiqc_bwa.html,/home/msleeper/scratch/reports/summary/896--01-02_trimmed_multiqc.html,/home/msleeper/scratch/reports/summary/896--02_trimmed_multiqc.html,/home/msleeper/scratch/reports/summary/896--07_final_multiqc_bwa.html,/home/msleeper/scratch/reports/summary/896--04_pre_merge_multiqc_bwa.html,/home/msleeper/scratch/reports/summary/896--01_raw_multiqc.html}' .