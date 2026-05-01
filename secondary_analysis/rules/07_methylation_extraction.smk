### Methylation extraction rules for bwa and bismark pathways
#--------------------------------------------------------------------------------
# bismark methylation extractor rules for bismark pathway
#--------------------------------------------------------------------------------

# # bismark_methylation_extractor --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder /path_to_genome_folder/ sample_bismark_bt2.bam
# # bismark_methylation_extractor --gzip --single-end --output_dir [dir] --split_by_chromosome --cytosine_report --bedGraph --ucsc --genome_folder <path> 
# bismark_methylation_extractor --gzip --paired-end --output_dir [dir] --cytosine_report  --bedGraph --ucsc
# bismark_methylation_extractor --gzip --single-end --output_dir {params.output_dir} --cytosine_report --bedGraph --ucsc --genome_folder {input.genome} {input.bam}        
# bismark_methylation_extractor --gzip --single-end --output_dir /home/msleeper/scratch/data/07_bismark_methyl_extractor --split_by_chromosome --cytosine_report --bedGraph --ucsc --genome_folder /home/msleeper/scratch/genomes/hg38/bismark/Bisulfite_Genome/ /home/msleeper/scratch/data/05_sambamba_bis/402--241-Cancer_M-SRX17589484-se.bam

rule bismark_methylation_extractor_se:
    input: 
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"])
    
    output:
        # Reports
        mbias_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
        splitting_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
        
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),
        ucsc_bedgraph = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph_UCSC.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        # read_base_meth_state_cpg = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        cpg_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        cpg_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),

        # * You could merge CHG, CHH using: --merge_non_CpG
        # read_base_meth_state_chg = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        # read_base_meth_state_chh = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chh_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chh_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chg_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chg_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),

        cph_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.CpG_report.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        c_summary = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.cytosine_context_summary.txt", root = config["root"], data_dir=config["data_dir"])
        
    log:
        "logs/secondary_rules/07_bismark_methylation_extractor_se/07_bismark_methylation_extractor_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/bismark.yaml"
    
    wildcard_constraints:
        layout="se"
    
    params:
        output_dir = expand("{root}/{data_dir}/07_bismark_methyl_extractor", root = config["root"], data_dir=config["data_dir"])
    
    threads: 2

    resources:
        mem_mb=4000

    shell:
        '''
        echo "Running bismark_methylation_extractor for {input.bam}" > {log}
        bismark_methylation_extractor --gzip --single-end --output_dir {params.output_dir} --cytosine_report --bedGraph --ucsc --genome_folder {input.genome} {input.bam}
        echo "done" >> {log}
        '''

rule bismark_methylation_extractor_pe:
    input: 
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"])
    
    output:
        # Reports
        mbias_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
        splitting_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
        
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),
        ucsc_bedgraph = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph_UCSC.bedGraph.gz", root = config["root"], data_dir=config["data_dir"]),

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        # read_base_meth_state_cpg = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        cpg_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        cpg_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CpG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),

        # * You could merge CHG, CHH using: --merge_non_CpG
        # read_base_meth_state_chg = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        # read_base_meth_state_chh = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chh_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chh_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHH_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chg_ot = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_OT_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        chg_ob = expand("{root}/{data_dir}/07_bismark_methyl_extractor/CHG_OB_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        
        cph_report = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.CpG_report.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        c_summary = expand("{root}/{data_dir}/07_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.cytosine_context_summary.txt", root = config["root"], data_dir=config["data_dir"])


    log:
        "logs/secondary_rules/07_bismark_methylation_extractor_pe/07_bismark_methylation_extractor_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/bismark.yaml"
    
    wildcard_constraints:
        layout="pe"
    
    params:
        output_dir = expand("{root}/{data_dir}/07_bismark_methyl_extractor", root = config["root"], data_dir=config["data_dir"])
    
    threads: 2

    resources:
        mem_mb=4000

    shell:
        '''
        echo "Running bismark_methylation_extractor for {input.bam}" > {log}
        bismark_methylation_extractor --gzip --paired-end --output_dir {params.output_dir} --cytosine_report --bedGraph --ucsc --genome_folder {input.genome} {input.bam}
        echo "done" >> {log}
        '''

#--------------------------------------------------------------------------------
# wgbstools methylation extractor rules for bwa and bismark pathways
#--------------------------------------------------------------------------------

# note: --mbias is only compatible with pe samples. If I run into issues with mbias plotting error messages I may need to split se and pe samples at this step

rule wgbstools_convert_bam_to_beta_bwa:
    ''' convert bams into pat and beta files'''
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        ref = expand("{root}/{wgbstools_ref_dir}/{fasta}", root = config["root"], wgbstools_ref_dir = config["ref"]["wgbstools_idx_dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{root}/{data_dir}/07_wgbstools_betas_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"])
    
    log:
        "logs/secondary_rules/07_wgbstools_convert_bam_to_beta_bwa/07_wgbstools_convert_bam_to_beta_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    params:
        outdir = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"]),
        tempdir = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/temp", root = config["root"], data_dir=config["data_dir"]),
        genome_name = config["ref"]["fasta"]

    threads: 2

    resources:
        mem_mb=4000

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        echo "Converting bam files {input.bam} to pat and beta files" > {log}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.genome_name} --temp_dir {params.tempdir} {input.bam} >> {log} 2>&1
        echo "done" >> {log}
        """

rule wgbstools_convert_bam_to_beta_bis:
    ''' convert bams into pat and beta files'''
    input:
        bam = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        ref = expand("{root}/{wgbstools_ref_dir}/{fasta}", root = config["root"], wgbstools_ref_dir = config["ref"]["wgbstools_idx_dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{root}/{data_dir}/07_wgbstools_betas_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"])
    
    log:
        "logs/secondary_rules/07_wgbstools_convert_bam_to_beta_bis/07_wgbstools_convert_bam_to_beta_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    threads: 2

    resources:
        mem_mb=4000

    params:
        outdir = expand("{root}/{data_dir}/07_wgbstools_betas_bis/", root = config["root"], data_dir=config["data_dir"]),
        tempdir = expand("{root}/{data_dir}/07_wgbstools_betas_bis/temp", root = config["root"], data_dir=config["data_dir"]),
        genome_name = config["ref"]["fasta"]

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        echo "Converting bam files {input.bam} to pat and beta files" > {log}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.genome_name} --temp_dir {params.tempdir} {input.bam} >> {log} 2>&1
        echo "done" >> {log}
        """

# bismark_methylation_extractor --gzip --single-end --output_dir /home/msleeper/scratch/data/07_bis_test/ --cytosine_report --bedGraph --genome_folder /home/msleeper/scratch/genomes/hg38/bismark/ /home/msleeper/scratch/data/06_merged_deduped_bis/215--N37-Colon-Normal-SRX1631736-pe.bam
# bismark_methylation_extractor --gzip --single-end --output_dir /home/msleeper/scratch/data/07_bis_test/ --cytosine_report --bedGraph --genome_folder /home/msleeper/scratch/genomes/hg38/bismark/ /home/msleeper/scratch/data/06_merged_deduped_bwa/215--N37-Colon-Normal-SRX1631736-pe.bam
# bismark_methylation_extractor --gzip --paired-end --output_dir /home/msleeper/scratch/data/07_bis_test/ --comprehensive --genome_folder /home/msleeper/scratch/genomes/hg38/bismark/ /home/msleeper/scratch/data/06_merged_deduped_bwa/215--N37-Colon-Normal-SRX1631736-pe.bam

rule methyldackel_bwa:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        ref = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"], fasta=config["ref"]["fasta"])

    output:
        cpg_context = expand("{root}/{data_dir}/07_methyldackel_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_CpG.bedGraph", root = config["root"], data_dir=config["data_dir"]),
        merged_cpg_context = expand("{root}/{data_dir}/07_methyldackel_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_mergedContext_CpG.bedGraph", root = config["root"], data_dir=config["data_dir"]),
        methylkit_in = expand("{root}/{data_dir}/07_methyldackel_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_CpG.methylKit", root = config["root"], data_dir=config["data_dir"]),
        cytosine_report = expand("{root}/{data_dir}/07_methyldackel_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.cytosine_report.txt", root = config["root"], data_dir=config["data_dir"])

    log:
        "logs/secondary_rules/07_methyldackel_bwa/07_methyldackel_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/methyldackel.yaml"

    threads: 2

    resources:
        mem_mb=4000

    params:
        out_dir = expand("{root}/{data_dir}/07_methyldackel_bwa/", root = config["root"], data_dir=config["data_dir"]),
        out_prefix = "{root}/{data_dir}/07_methyldackel_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}".format(root = config["root"], data_dir=config["data_dir"])

    shell:
        '''
        echo "Creating output directory {params.out_dir}" > {log}
        mkdir -p {params.out_dir}
        echo "Running methyldackel for {input.bam}" > {log}
        MethylDackel extract -o {params.out_prefix} {input.ref} {input.bam}
        MethylDackel extract --mergeContext -o {params.out_prefix}_mergedContext {input.ref} {input.bam}
        MethylDackel extract --cytosine_report -o {params.out_prefix} {input.ref} {input.bam}
        MethylDackel extract --methylKit -o {params.out_prefix} {input.ref} {input.bam}
        echo "done" >> {log}
        '''