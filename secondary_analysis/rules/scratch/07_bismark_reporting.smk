# rule move_bis_reports_se:
#     input:
        
#     output:
#         # Reports
#         mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
#         splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylone_CpG_cov = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),
#         ucsc_bedgraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bedGraph_UCSC.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),
#         # Primary output files: methylation status at each read cytosine position: (extremely large)
#         read_base_meth_state_cpg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CpG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         cpg_ot = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CpG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         cpg_ob = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CpG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         # * You could merge CHG, CHH using: --merge_non_CpG
#         read_base_meth_state_chg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         read_base_meth_state_chh = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHH_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"])
#         chh_ot = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHH_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         chh_ob = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHH_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         chg_ot = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"]),
#         chg_ob = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.txt.gz", root = config["root"], data_dir=config["data_dir"])

#         mbias_r1 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], data_dir=config["data_dir"]),
#         mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
#         splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylone_CpG_cov = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"])
    
#     output:
#         mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
#         mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
#         splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
#         methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
#         methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"])
    
#     log:
#         "logs/secondary_rules/06_move_bis_reports/06_move_bis_reports_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
#     # shadow:
#     #     "shallow"    
    
#     wildcard_constraints:
#         layout="se"
    
#     params:
#         output_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor", root = config["root"], rep_dir=config["reports_dir"])

#     shell:
#         """
#         mkdir -p {params.output_dir}
#         mv {input.mbias_r1} {output.mbias_r1}
#         mv {input.mbias_report} {output.mbias_report}
#         mv {input.splitting_report} {output.splitting_report}
#         mv {input.methylone_CpG_cov} {output.methylone_CpG_cov}
#         mv {input.methylome_CpG_mlevel_bedGraph} {output.methylome_CpG_mlevel_bedGraph}
#         """

    
# rule move_bis_reports_pe:
#     input:
#         mbias_r1 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], data_dir=config["data_dir"]),
#         # Only for PE BAMS:
#         mbias_r2 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R2.png", root = config["root"], data_dir=config["data_dir"]),
#         mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
#         splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylone_CpG_cov = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"])
    
#     output:
#         mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
#         mbias_r2 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R2.png", root = config["root"], rep_dir=config["reports_dir"]),
#         mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
#         splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
#         methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
#         methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"])

#     log:
#         "logs/secondary_rules/06_move_bis_reports/06_move_bis_reports_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

#     wildcard_constraints:
#         layout="pe"

#     params:
#         output_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor", root = config["root"], rep_dir=config["reports_dir"])
    
#     shell:
#         """
#         mkdir -p {params.output_dir}
#         mv {input.mbias_r1} {output.mbias_r1}
#         mv {input.mbias_r2} {output.mbias_r2}
#         mv {input.mbias_report} {output.mbias_report}
#         mv {input.splitting_report} {output.splitting_report}
#         mv {input.methylone_CpG_cov} {output.methylone_CpG_cov}
#         mv {input.methylome_CpG_mlevel_bedGraph} {output.methylome_CpG_mlevel_bedGraph}
#         """
        # alignment_report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        # dedup_report = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        
rule bismark_merged_sample_reporting_se:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        alignment_report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        dedup_report = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
        splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_splitting_report.txt", root = config["root"], data_dir=config["data_dir"])

    output:
        bis_sample_report = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_SE_report.html", root = config["root"], rep_dir=config["reports_dir"]),
        bis_report_data = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        nucleotide_freq_report = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.nucleotide_stats.txt", root = config["root"], rep_dir=config["reports_dir"]),
        mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
        splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "logs/secondary_rules/06_bismark_reporting_se/06_bismark_reporting_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    
    wildcard_constraints:
        layout="se"
    
    params:
        work_dir = expand("{root}/{data_dir}/06_bismark_methyl_extractor/", root = config["root"], data_dir=config["data_dir"]),
        report_dir = expand("{root}/{rep_dir}/06_bismark_summary/", root = config["root"], rep_dir=config["reports_dir"]),
        meth_report_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/", root = config["root"], rep_dir=config["reports_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"])

    shell:
        """
        echo "Creating symbolic link to {input.bam} in {params.work_dir}" > {log}
        ln -s {input.bam} {params.work_dir}
        echo "Creating symbolic links for alignment report and deduplication report" >> {log}
        ln -s {input.alignment_report} {params.work_dir}
        ln -s {input.dedup_report} {params.work_dir}
        echo " Running bam2nuc to create nucleotide frequency report" >> {log}
        bam2nuc --dir {params.work_dir} --genome_folder {params.genome} {input.bam}
        echo "Running bismark2report to create sample report" >> {log}
        bismark2report --dir {params.output_dir}
        echo "Moving methylatiuon extractor reports to report directory {params.meth_report_dir}" >> {log}
        mkdir -p {params.meth_report_dir}
        mv -f -v --target-directory={params.meth_report_dir} {input.mbias_report} {input.splitting_report}
        echo "Moving summary reports to report directory {params.report_dir}" >> {log}
        mkdir -p {params.report_dir}
        mv -f -v --target-directory={params.report_dir} {params.work_dir}/*_merged_SE_report.html {params.work_dir}/*_merged_SE_report.txt
        echo "Done" >> {log}
        """
    

rule bismark_premerged_sample_reporting_se:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        alignment_report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        dedup_report = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"]),

    output:
        bis_sample_report = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_SE_report.html", root = config["root"], rep_dir=config["reports_dir"]),
        bis_report_data = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        nucleotide_freq_report = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.nucleotide_stats.txt", root = config["root"], rep_dir=config["reports_dir"]),

    log:
        "logs/secondary_rules/06_bismark_reporting_se/06_bismark_reporting_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    
    wildcard_constraints:
        layout="se"
    
    params:
        work_dir = expand("{root}/{data_dir}/06_bismark_methyl_extractor/", root = config["root"], data_dir=config["data_dir"]),
        report_dir = expand("{root}/{rep_dir}/06_bismark_summary/", root = config["root"], rep_dir=config["reports_dir"]),
        meth_report_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/", root = config["root"], rep_dir=config["reports_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"])

    shell:
        """
        echo "Creating symbolic link to {input.bam} in {params.work_dir}" > {log}
        ln -s {input.bam} {params.work_dir}
        echo "Creating symbolic links for alignment report and deduplication report" >> {log}
        ln -s {input.alignment_report} {params.work_dir}
        ln -s {input.dedup_report} {params.work_dir}
        echo " Running bam2nuc to create nucleotide frequency report" >> {log}
        bam2nuc --dir {params.work_dir} --genome_folder {params.genome} {input.bam}
        echo "Running bismark2report to create sample report" >> {log}
        bismark2report --dir {params.output_dir}
        echo "Moving methylatiuon extractor reports to report directory {params.meth_report_dir}" >> {log}
        mkdir -p {params.meth_report_dir}
        mv -f -v --target-directory={params.meth_report_dir} {input.mbias_report} {input.splitting_report}
        echo "Moving summary reports to report directory {params.report_dir}" >> {log}
        mkdir -p {params.report_dir}
        mv -f -v --target-directory={params.report_dir} {params.work_dir}/*_merged_SE_report.html {params.work_dir}/*_merged_SE_report.txt
        echo "Done" >> {log}
        """
    