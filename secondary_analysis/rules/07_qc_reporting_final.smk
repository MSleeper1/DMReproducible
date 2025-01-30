
#---------------------------------------------------------------
# MultiQC rule for 06-07 post merge and methylation reports
#---------------------------------------------------------------
# Compliling a list of reports based on the alignment pathways in the config file
# Bismark and Bwameth pathways both used
if config["bismark"] == True and config["bwameth"]==True:
    post_merge_report_list = [
        expand("{root}/{rep_dir}/06_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bwa
        expand("{root}/{rep_dir}/06_feature_counts_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bis
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bwa
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bis
        expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
        expand("{root}/{rep_dir}/06_samtools_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats",  root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bis
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        # expand("{root}/{rep_dir}/06_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bwa
        # expand("{root}/{rep_dir}/06_qualimap_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bis
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis per-base report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis summary report
        expand("{root}/{rep_dir}/06_bismark_summary/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_merged_bismark_deduplication/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.deduplication_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.log", root = config["root"], sample = sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.cytosine_context_summary.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"])
    ]

# Bismark pathway only
# only add bismark pathway reports to the list
if config["bismark"] == True and config["bwameth"]==False:
    post_merge_report_list = [
        expand("{root}/{rep_dir}/06_feature_counts_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bis
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bis
        expand("{root}/{rep_dir}/06_samtools_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats",  root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bis
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        # expand("{root}/{rep_dir}/06_qualimap_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bis
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis per-base report
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # mosdepth bis summary report
        expand("{root}/{rep_dir}/06_bismark_summary/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_merged_bismark_deduplication/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.deduplication_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"]),
        expand("{root}/{data_dir}/07_bismark_methyl_extractor/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.cytosine_context_summary.txt", sample = sample_info.itertuples(), root = config["root"], data_dir=config["data_dir"])    
    ]

# Bwameth pathway only
# only add bwa pathway reports to the list
if config["bismark"]==False and config["bwameth"]==True:
    post_merge_report_list = [
        expand("{root}/{rep_dir}/06_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample = sample_info.itertuples()), # featureCounts bwa
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"], sample = sample_info.itertuples()), # fastqc report bwa
        expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam.stats", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # samtools stats report bwa
        expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.flagstat.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),        
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/06_merged_deduped_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.log", root = config["root"], sample = sample_info.itertuples(), rep_dir=config["reports_dir"]),
        # expand("{root}/{rep_dir}/06_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}", root = config["root"], rep_dir=config["reports_dir"], sample = sample_info.itertuples()), # qualimap report bwa
        ]

# multiqc rule for 06-07 post merge and methylation reports
rule multiqc_compile_reports_07:
    input:
        post_merge_report_list
    output:
        file = expand("{root}/{rep_dir}/summary/07_final_multiqc.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/07_final_qc_report_list.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/07_multiqc_command.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/07_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    params:
        output_filename = "07_final_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for files post merge and methylation extraction. Includes reports from fastqc, samtools, mosdepth featurecounts... Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "07 Final QC Reports for Secondary Analysis", # Optional: title for multiqc report (can be changed freely)
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

# rule bismark_sample_report:
#     input:
#         bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bam", root = config["root"], data_dir=config["data_dir"]),
#         mbias_report = expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
#         splitting_report = expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
#         nucleotide_report = expand("{rep_dir}/06_bismark_summary/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.nucleotide_stats.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
#         dedup_report = expand("{rep_dir}/06_merged_bismark_deduplication/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.deduplication_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
#         cytosine_report = expand("{data_dir}/07_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.cytosine_context_summary.txt", sample=sample_info.itertuples(), data_dir=config["data_dir"])
#     output:
 