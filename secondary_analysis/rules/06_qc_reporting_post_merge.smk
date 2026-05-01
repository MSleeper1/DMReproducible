### Snakemake rules for quality control of merged bam files (05) ###
# qualimap, samtools stats, fastqc, and featureCounts are used to generate quality control reports for the merged bam files.
# reports are saved to the reports directory. 
# directories will have 05 prefix to indicate that they are part of the 05_quality_control section of the pipeline reporting on merged bam files.

### QUALIMAP ###
# Qualimap is used to generate quality control reports for the merged bam files.
# Qualimaps flag --outdir is a relative path. To avoid issues with this, the ouput will be moved to the reports directory after qualimap is finished.

# qualimap_post_merge: qualimap is run on the merged bam files. The output is moved to the reports directory after qualimap is finished.
# input: merged bam files and gtf file (must be unzipped)
# output: qualimap reports

rule qualimap_post_merge_bwa:
    input:
        bwa_bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
    
    output:
        directory(expand("{root}/{rep_dir}/06_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]))
        
    log:
        "logs/secondary_rules/06_qualimap_post_merge_bwa/06_qualimap_post_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    threads: 4

    resources:
        mem_mb=8000
    
    params:
        out_dir = expand("{root}/{rep_dir}/06_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        extra = "--java-mem-size=35G"
    shell:
        """
        echo making output directory {params.out_dir} > {log}
        mkdir -p {params.out_dir}
        echo "Running qualimap on {input.bwa_bam}" > {log}
        qualimap bamqc -bam {input.bwa_bam} -c -sd -os -gd hg38 -gff {input.gtf} {params.extra} --outdir {params.out_dir} >> {log} 2>&1
        echo "Done" >> {log}
        """

rule qualimap_post_merge_bis:
    input:
        bis_bam = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])

    output:
        directory(expand("{root}/{rep_dir}/06_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]))
        
    log:
        "logs/secondary_rules/06_qualimap_post_merge_bis/06_qualimap_post_merge_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    threads: 4

    resources:
        mem_mb=8000
    
    params:
        out_dir = expand("{root}/{rep_dir}/06_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        extra = "--java-mem-size=35G"
    shell:
        """
        echo making output directory {params.out_dir} > {log}
        mkdir -p {params.out_dir}
        echo "Running qualimap on {input.bis_bam}" > {log}
        qualimap bamqc -bam {input.bis_bam} -c -sd -os -gd hg38 -gff {input.gtf} {params.extra} --outdir {params.out_dir} >> {log} 2>&1
        echo "Done" >> {log}
        """


### FEATURE COUNTS ###
# FeatureCounts is a program that counts the number of reads that map to each feature in a GTF file
# input.sam can be a bam file but must be called input.sam to function with wrapper

rule feature_counts_post_merge_bwa_se:
    input:
        sam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/06_feature_counts_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/06_feature_counts_bwa_se/06_feature_counts_bwa_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads: 2

    resources:
        mem_mb=4000

    conda:
        "../../environment_files/feature_counts.yaml"

    wildcard_constraints:
        layout = "se"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule feature_counts_post_merge_bwa_pe:
    input:
        sam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/06_feature_counts_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/06_feature_counts_bwa_pe/06_feature_counts_bwa_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads: 2

    resources:
        mem_mb=4000

    conda:
        "../../environment_files/feature_counts.yaml"

    wildcard_constraints:
        layout = "pe"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f -p"

    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule feature_counts_post_merge_bis_se:
    input:
        sam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/06_feature_counts_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/05_feature_counts_bis_se/05_feature_counts_bis_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads: 2

    resources:
        mem_mb=4000

    conda:
        "../../environment_files/feature_counts.yaml"

    wildcard_constraints:
        layout = "se"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule feature_counts_post_merge_bis_pe:
    input:
        sam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/06_feature_counts_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/05_feature_counts_bis_pe/05_feature_counts_bis_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads: 2

    resources:
        mem_mb=4000

    conda:
        "../../environment_files/feature_counts.yaml"

    wildcard_constraints:
        layout = "pe"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f -p"

    wrapper:
        "0.72.0/bio/subread/featurecounts"   


### FASTQC ###
# added sleep 5 due to latency issues


rule fastqc_post_merge_bwa:
    input:
        expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),

    output:
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/06_fastqc_post_merge_bwa/06_fastqc_post_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000

    params:
        output_dir = expand("{root}/{rep_dir}/06_fastqc_post_merge_bwa", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
        mkdir -p {params.output_dir}
        echo "Running fastqc on {input}" > {log}
        fastqc -o {params.output_dir} {input} >> {log} 2>&1
        echo "Done" >> {log}
        """

rule fastqc_post_merge_bis:
    input:
        expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])  

    output:
        expand("{root}/{rep_dir}/06_fastqc_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/06_fastqc_post_merge_bis/06_fastqc_post_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000
    
    params:
        output_dir = expand("{root}/{rep_dir}/06_fastqc_post_merge_bis", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
        mkdir -p {params.output_dir}
        echo "Running fastqc on {input}" > {log}
        fastqc -o {params.output_dir} {input} >> {log} 2>&1
        echo "Done" >> {log}
        """


### SAMTOOLS STATS ###

rule samtools_stats_post_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
 
    output:
        report = expand("{root}/{rep_dir}/06_samtools_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/06_samtools_post_merge_bwa/06_samtools_post_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000
    
    shell:
        """
        echo "Running samtools stats on {input}" > {log}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {log}
        """

rule samtools_stats_post_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])  

    output:
        report = expand("{root}/{rep_dir}/06_samtools_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.stats",  root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/06_samtools_post_merge_bis/06_samtools_post_merge_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000

    shell:
        """
        echo "Running samtools stats on {input}" > {log}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {log}
        """

### Samtools flagstat rule ###
# samtools flagstat is a program that counts the number of reads that pass various filters
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: samtools flagstat report for deduplicated sequence files

rule samtools_flagstat_post_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),

    output:
        report = expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.flagstat.txt", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000

    shell:
        """
        echo "Running samtools flagstat on {input.bam}" > {output.report}
        samtools flagstat {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

rule samtools_flagstat_post_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        report = expand("{root}/{rep_dir}/06_samtools_flagstat_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.flagstat.txt", root = config["root"], rep_dir=config["reports_dir"])
    
    conda:
        "../../environment_files/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000

    shell:
        """
        echo "Running samtools flagstat on {input.bam}" > {output.report}
        samtools flagstat {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

### MOSDEPTH RULE ###
# mosdepth is a program that calculates the depth of coverage for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: mosdepth report for deduplicated sequence files (global distribution text file, per-base bed file, and summary text file)

# ADD mosdepth rule for bwa

rule mosdepth_post_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/06_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.bai", root = config["root"], data_dir=config["data_dir"])

    output:
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.per-base.bed.gz.csi", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing

    log:
        "logs/secondary_rules/05_mosdepth_bis/05_mosdepth_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads: 2

    resources:
        mem_mb=4000

    conda:
        "../../environment_files/mosdepth.yaml"

    params:
        extra="--fast-mode",  # optional
        mapping_quality = 10,
        out_prefix = expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        outdir = expand("{root}/{rep_dir}/06_mosdepth_post_merge_bis", root = config["root"], rep_dir=config["reports_dir"])
    
    shell:
        '''
        echo "making output directory" > {log}
        mkdir -p {params.outdir}
        echo "Running mosdepth on {input.bam}" >> {log}
        mosdepth -x -Q {params.mapping_quality} {params.out_prefix} {input.bam}
        echo "Done" >> {log}
        '''

rule mosdepth_post_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/06_merged_deduped_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam.bai", root = config["root"], data_dir=config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.per-base.bed.gz.csi", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing

    log:
        "logs/secondary_rules/06_mosdepth_bwa/06_mosdepth_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/mosdepth.yaml"

    threads: 2

    resources:
        mem_mb=4000

    params:
        extra="--fast-mode",  # optional
        mapping_quality = 10,
        out_prefix = expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        outdir = expand("{root}/{rep_dir}/06_mosdepth_post_merge_bwa", root = config["root"], rep_dir=config["reports_dir"])
    
    shell:
        '''
        echo "making output directory" > {log}
        mkdir -p {params.outdir}
        echo "Running mosdepth on {input.bam}" >> {log}
        mosdepth -x -Q {params.mapping_quality} {params.out_prefix} {input.bam}
        echo "Done" >> {log}
        '''

# # Get the depth for each sample
# rule mosdepth:
#     input:
#         '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai',
#         bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
#     output:
#         '6_mosdepth/{sample}.sorted.markdupes.mosdepth.global.dist.txt',
#         '6_mosdepth/{sample}.sorted.markdupes.mosdepth.summary.txt',
#         '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz',
#         '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz.csi'
#     threads:
#         config['mosdepth']['threads']
#     params:
#         mapping_quality = config['mosdepth']['mapping_quality'],
#         mosdepth_path = config['paths']['mosdepth_path'],
#         out_prefix = '6_mosdepth/{sample}.sorted.markdupes'
#     shell:
#         '''
#         {params.mosdepth_path} \
#         -x \
#         -t {threads} \
#         -Q {params.mapping_quality} \
#         {params.out_prefix} \
#         {input.bam}
#         '''

# # Calculate the coverage from the mosdepth output
# rule calc_coverage:
#     input:
#         bed = '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz'
#     output:
#         '6_mosdepth/{sample}.sorted.markdupes.coverage.txt'
#     params:
#         genome = REFERENCE_GENOME
#     shell:
#         '''
#         scripts/mosdepth_to_x_coverage.py \
#         -f {params.genome} \
#         -m {input.bed} \
#         > {output}
#         '''

# Bismark bam2nuc rule
# Calculate nucleotide frequency report for deduplicated sequence files
rule bismark_post_merge_nuc_freq:
    input:
        bam = expand("{root}/{data_dir}/06_merged_deduped_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bam", root = config["root"], data_dir=config["data_dir"])  

    output:
        nucleotide_freq_report = expand("{root}/{rep_dir}/06_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.nucleotide_stats.txt", root = config["root"], rep_dir=config["reports_dir"])
 
    log:
        "logs/secondary_rules/06_bismark_post_merge_nuc_freq/06_bismark_post_merge_nuc_freq-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/bismark.yaml"
    
    threads: 2

    resources:
        mem_mb=4000

    params:
        out_dir = expand("{root}/{rep_dir}/06_bismark_summary/", root = config["root"], rep_dir=config["reports_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"])

    shell:
        """
        echo "Making output directory" > {log}
        mkdir -p {params.out_dir}
        echo "Removing old nucleotide frequency report if it exists" >> {log}
        if [ -f {output.nucleotide_freq_report} ]; then 
            rm {output.nucleotide_freq_report}
        fi
        echo "Creating output report file" >> {log}
        touch {output.nucleotide_freq_report}
        echo " Running bam2nuc to create nucleotide frequency report" >> {log}
        bam2nuc --dir {params.out_dir} --genome_folder {params.genome} {input.bam}
        echo "Done" >> {log}
        """