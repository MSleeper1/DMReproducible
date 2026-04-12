# This rule compiles reports from the 03 alignment and 04 deduplication rules for multiqc.
# Inputs are determined based on the alignment pathways specified in the config file.


# # Report directories for qualimap output
# # Bismark and Bwameth pathways both used
# if config["bismark"] == True and config["bwameth"]==True:
#     qualimap_list = [
#         directory(expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])),
#         directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]))
# ]
# # Bismark pathway only
# # only add bismark pathway reports to the list
# if config["bismark"] == True and config["bwameth"]==False:
#     qualimap_list = [
#         directory(expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]))
# ]
# # Bwameth pathway only
# # only add bwa pathway reports to the list
# if config["bismark"]==False and config["bwameth"]==True:
#     qualimap_list = [
#         directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]))
#     ]

#------------------------------------------------
# MultiQC rule for 04_pre_merge reports
#------------------------------------------------
rule multiqc_compile_bwa_reports_04:
    input:
        expand("{root}/{rep_dir}/03_bwameth/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # bwameth alignment report for PE and SE
        expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        # expand("{root}/{rep_dir}/04_fastqc_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        expand("{root}/{rep_dir}/04_feature_counts_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])

    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--04_pre_merge_multiqc_bwa.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--04_qc_report_list_bwa.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--04_multiqc_command_bwa.sh", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "logs/secondary_rules/multiqc/{ref}--04_multiqc_bwa.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--04_pre_merge_multiqc_bwa", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed and deduplicated reads before merging. Includes reports from fastqc, mosdepth, featurecounts, samtools, deduplication, and alignment. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Pre-merge Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
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

rule multiqc_compile_bis_reports_04:
    input:
        expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", pe = sample_info_pe.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # PE bis alignment report
        expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", se = sample_info_se.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # SE bis alignment report
        expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplication_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        expand("{root}/{rep_dir}/04_feature_counts_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}__bismark_deduplicated.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])
    output:
        file = expand("{root}/{rep_dir}/summary/{{ref}}--04_pre_merge_multiqc_bis.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/{{ref}}--04_qc_report_list_bis.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/{{ref}}--04_multiqc_command_bis.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/{ref}--04_multiqc_bis.log"
    conda:
        "../../environment_files/multiqc.yaml"
    wildcard_constraints:
        ref = "|".join(sample_info["ref"].tolist())
    params:
        output_filename = "{ref}--04_pre_merge_multiqc_bis", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed and deduplicated reads before merging. Includes reports from fastqc, mosdepth, featurecounts, samtools, deduplication, and alignment. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Pre-merge Sequence QC Reports for dataset ref:{ref}", # Optional: title for multiqc report (can be changed freely)
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

