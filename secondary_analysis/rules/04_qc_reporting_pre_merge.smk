### Snakemake rules for quality control reporting of trimmed, aligned, and deduplicated sequence files (04) ###
# fastqc and samtools stats are run on the deduplicated sequence files

### FASTQC RULES ###
# input: trimmed, alignned, and deduplicated sequence files (bam)
# output: fastqc reports for deduplicated sequence files (html and zip)

# Rule to run fastqc on bam files after deduplication
rule fastqc_pre_merge_bwa:
    input: 
        bwa_bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir = config["data_dir"]),
    
    output:
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/04_fastqc_pre_merge_bwa/04_fastqc_pre_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    params:
        output_dir = expand("{root}/{rep_dir}/04_fastqc_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
		mkdir -p {params.output_dir}
		echo "Running fastqc on {input.bwa_bam}" > {log}
		fastqc -o {params.output_dir} {input.bwa_bam} >> {log} 2>&1
		echo "Done" >> {log}
		"""

rule fastqc_pre_merge_bis:
    input: 
        bis_bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"])
        
    output:
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/04_fastqc_pre_merge_bis/04_fastqc_pre_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    params:
        output_dir = expand("{root}/{rep_dir}/04_fastqc_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
		mkdir -p {params.output_dir}
		echo "Running fastqc on {input.bis_bam}" > {log}
		fastqc -o {params.output_dir} {input.bis_bam} >> {log} 2>&1
		echo "Done" >> {log}
		"""

### SAMTOOLS STATS RULE ###
# samtools stats is a program that generates general statistics for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: samtools stats report for deduplicated sequence files (.stats text file)
rule samtools_stats_pre_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"]), # sambamba output
        ref = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])
    
    output:
        report = expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input.bam}" > {output.report}
        samtools stats -p -d -r {input.ref} {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

rule samtools_stats_pre_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"]),
        ref = expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"])
    
    output:
        report = expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input.bam}" > {output.report}
        samtools stats -p -d -r {input.ref} {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

### Samtools Flagstat Rule ###
# samtools flagstat is a program that generates general statistics for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: samtools flagstat report for deduplicated sequence files (.flagstat text file)
rule samtools_flagstat_pre_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        report = expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.flagstat", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools flagstat on {input.bam}" > {output.report}
        samtools flagstat {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

rule samtools_flagstat_pre_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        report = expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam.flagstat", root = config["root"], rep_dir=config["reports_dir"])
    
    conda:
        "../../environment_files/samtools.yaml"
    
    shell:
        """
        echo "Running samtools flagstat on {input.bam}" > {output.report}
        samtools flagstat {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

### QUALIMAP RULES ###
# qualimap is a program that helps identify contamination and other issues in sequence files
rule qualimap_pre_merge_bwa:
    input:
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        bwa_bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir = config["data_dir"]),

    output:
        directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge", root = config["root"], rep_dir=config["reports_dir"])),
        # expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge/genome_results.txt", root = config["root"], rep_dir=config["reports_dir"]),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge/raw_data_qualimapReport", root = config["root"], rep_dir=config["reports_dir"]))

    log:
        "logs/secondary_rules/04_qualimap_pre_merge_bwa/04_qualimap_pre_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    params:
        out_dir = expand("{root}/{rep_dir}/04_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge", root = config["root"], rep_dir=config["reports_dir"]),
        extra = "--java-mem-size=14G"

    shell:
        """
        echo "making output directory" > {log}
        mkdir -p {params.out_dir}
        echo "Running qualimap on {input.bwa_bam}" > {log}
        qualimap bamqc -bam {input.bwa_bam} -c -sd -os -gd hg38 -gff {input.gtf} {params.extra} --outdir {params.out_dir} >> {log} 2>&1
        echo "Done" >> {log}
        """

rule qualimap_pre_merge_bis:
    input:
        bis_bam = expand("{root}/{data_dir}/04_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/04_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.bam.bai", root = config["root"], data_dir=config["data_dir"]),
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
    
    output:
        directory(expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge", root = config["root"], rep_dir=config["reports_dir"])),
        # expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge/genome_results.txt", root = config["root"], rep_dir=config["reports_dir"]),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge/raw_data_qualimapReport", root = config["root"], rep_dir=config["reports_dir"]))

    log:
        "logs/secondary_rules/04_qualimap_pre_merge_bis/04_qualimap_pre_merge_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    params:
        out_dir = expand("{root}/{rep_dir}/04_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_pre_merge", root = config["root"], rep_dir=config["reports_dir"]),
        extra = "--java-mem-size=14G"

    shell:
        """
        echo "making output directory" > {log}
        mkdir -p {params.out_dir}
        echo "Running qualimap on {input.bis_bam}" > {log}
        qualimap bamqc -bam {input.bis_bam} -c -sd -os -gd hg38 -gff {input.gtf} {params.extra} --outdir {params.out_dir} >> {log} 2>&1
        echo "Done" >> {log}
        """


### FEATURE COUNTS ###
# FeatureCounts is a program that counts the number of reads that map to each feature in a GTF file
# input.sam can be a bam file but must be called input.sam to function with wrapper

rule feature_counts_pre_merge_bwa_se:
    input:
        sam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir = config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/04_feature_counts_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/04_feature_counts_bwa/04_feature_counts_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    threads:
        3

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

rule feature_counts_pre_merge_bis_se:
    input:
        sam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag

    output:
        expand("{root}/{rep_dir}/04_feature_counts_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}__bismark_deduplicated.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/04_feature_counts_bis/04_feature_counts_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    threads:
        3

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


rule feature_counts_pre_merge_bwa_pe:
    input:
        sam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir = config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/04_feature_counts_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/04_feature_counts_bwa/04_feature_counts_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    threads:
        3

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

rule feature_counts_pre_merge_bis_pe:
    input:
        sam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"]),
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        fasta=expand("{root}/{genomes_dir}/{genome}/{fasta}.fa", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], fasta = config["ref"]["fasta"]) # implicitly sets the -G flag

    output:
        expand("{root}/{rep_dir}/04_feature_counts_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}__bismark_deduplicated.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/04_feature_counts_bis/04_feature_counts_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    threads:
        3

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

### MOSDEPTH RULE ###
# mosdepth is a program that calculates the depth of coverage for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: mosdepth report for deduplicated sequence files (global distribution text file, per-base bed file, and summary text file)
rule mosdepth_pre_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/04_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/04_bismark_sorted_by_coordinate/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.bam.bai", root = config["root"], data_dir=config["data_dir"])

    output:
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.per-base.bed.gz.csi", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing

    log:
        "logs/secondary_rules/04_mosdepth_bis/04_mosdepth_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/mosdepth.yaml"

    params:
        extra="--fast-mode",  # optional
        mapping_quality = 10,
        out_prefix = expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_sorted", root = config["root"], rep_dir=config["reports_dir"]),
        outdir = expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis", root = config["root"], rep_dir=config["reports_dir"])
    
    shell:
        '''
        echo "making output directory" > {log}
        mkdir -p {params.outdir}
        echo "Running mosdepth on {input.bam}" >> {log}
        mosdepth -x -Q {params.mapping_quality} {params.out_prefix} {input.bam}
        echo "Done" >> {log}
        '''

rule mosdepth_pre_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.bai", root = config["root"], data_dir=config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.per-base.bed.gz.csi", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing

    log:
        "logs/secondary_rules/04_mosdepth_bwa/04_mosdepth_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/mosdepth.yaml"

    params:
        extra="--fast-mode",  # optional
        mapping_quality = 10,
        out_prefix = expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup", root = config["root"], rep_dir=config["reports_dir"]),
        outdir = expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa", root = config["root"], rep_dir=config["reports_dir"])
    
    shell:
        '''
        echo "making output directory" > {log}
        mkdir -p {params.outdir}
        echo "Running mosdepth on {input.bam}" >> {log}
        mosdepth -x -Q {params.mapping_quality} {params.out_prefix} {input.bam}
        echo "Done" >> {log}
        '''

# Bismark bam2nuc rule
# Calculate nucleotide frequency report for deduplicated sequence files
rule bismark_pre_merge_nuc_freq:
    input:
        bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        nucleotide_freq_report = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.nucleotide_stats.txt", root = config["root"], rep_dir=config["reports_dir"]),
 
    log:
        "logs/secondary_rules/04_bismark_pre_merge_nuc_freq/04_bismark_pre_merge_nuc_freq-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/bismark.yaml"
    
    params:
        out_dir = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/", root = config["root"], rep_dir=config["reports_dir"]),
        genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"]),

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

#------------------------------------------------
# MultiQC rule for 04_pre_merge reports
#------------------------------------------------
# Compliling a list of reports based on the alignment pathways in the config file
# Bismark and Bwameth pathways both used
if config["bismark"] == True and config["bwameth"]==True:
    pre_merge_report_list = [
        expand("{root}/{rep_dir}/03_bismark_bwt2/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", pe = sample_info_pe.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # PE bis alignment report
        expand("{root}/{rep_dir}/03_bismark_bwt2/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", se = sample_info_se.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # SE bis alignment report
        expand("{root}/{rep_dir}/03_bwameth/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # bwameth alignment report for PE and SE
        expand("{root}/{rep_dir}/04_bismark_deduplication/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplication_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])),
        expand("{root}/{rep_dir}/04_feature_counts_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}__bismark_deduplicated.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/04_bismark_summary/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    ]

# Bismark pathway only
# only add bismark pathway reports to the list
if config["bismark"] == True and config["bwameth"]==False:
    pre_merge_report_list = [
        expand("{root}/{rep_dir}/03_bismark_bwt2/{pe.ref}--{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", pe = sample_info_pe.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # PE bis alignment report
        expand("{root}/{rep_dir}/03_bismark_bwt2/{se.ref}--{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", se = sample_info_se.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # SE bis alignment report
        expand("{root}/{rep_dir}/04_bismark_deduplication/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplication_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])),
        expand("{root}/{rep_dir}/04_feature_counts_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}__bismark_deduplicated.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bis/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_sorted.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # this named output is required for prefix parsing
        expand("{root}/{rep_dir}/04_bismark_summary/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.nucleotide_stats.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
    ]

# Bwameth pathway only
# only add bwa pathway reports to the list
if config["bismark"]==False and config["bwameth"]==True:
    pre_merge_report_list = [
        expand("{root}/{rep_dir}/03_bwameth/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # bwameth alignment report for PE and SE
        expand("{root}/{rep_dir}/04_sambamba_bwameth_dedup/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_fastqc_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"]),
        # directory(expand("{root}/{rep_dir}/04_qualimap_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_pre_merge", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"])),
        expand("{root}/{rep_dir}/04_feature_counts_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.featureCounts{suf}", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_samtools_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.flagstat", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.global.dist.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.per-base.bed.gz.csi", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        expand("{root}/{rep_dir}/04_mosdepth_pre_merge_bwa/{sample.ref}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.mosdepth.summary.txt", sample = sample_info.itertuples(), root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing

    ]

# multiqc rule for 04 pre-merge reports
rule multiqc_compile_reports_04:
    input:
        pre_merge_report_list
    output:
        file = expand("{root}/{rep_dir}/summary/04_pre_merge_multiqc.html", root = config["root"], rep_dir=config["reports_dir"]),
        input_list = expand("{root}/{rep_dir}/summary/04_qc_report_list.txt", root = config["root"], rep_dir=config["reports_dir"]),
        multiqc_command = expand("{root}/{rep_dir}/summary/04_multiqc_command.sh", root = config["root"], rep_dir=config["reports_dir"])
    log:
        "logs/secondary_rules/multiqc/04_multiqc.log"
    conda:
        "../../environment_files/multiqc.yaml"
    params:
        output_filename = "04_pre_merge_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["root"] + "/" + config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed and deduplicated reads before merging. Includes reports from fastqc, mosdepth, featurecounts, samtools, deduplication, and alignment. Reports compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Pre-merge Sequence QC Reports", # Optional: title for multiqc report (can be changed freely)
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


# #ADD Bismarl report 

# # rule bismark_report_se:
# #     input:
# #         bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"]) # sambamba output
# #         alignment_report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
# #         dedup_report = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"]),

# #     output:
# #         bis_sample_report = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup_SE_report.html", root = config["root"], rep_dir=config["reports_dir"]),
# #         bis_report_data = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
# #         nucleotide_freq_report = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.nucleotide_stats.txt", root = config["root"], rep_dir=config["reports_dir"]),

# #     log:
# #         "logs/secondary_rules/04_bismark_report_se/04_bismark_report_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    
# #     wildcard_constraints:
# #         layout="se"
    
# #     params:
# #         work_dir = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/", root = config["root"], data_dir=config["data_dir"]),
# #         output_dir = expand("{root}/{rep_dir}/04_bismark_summary/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/", root = config["root"], rep_dir=config["reports_dir"]),
# #         genome = expand("{root}/{genomes_dir}/{genome}/bismark/", root = config["root"], genomes_dir=config["genomes_dir"], genome=config["ref"]["genome"]),

# #     shell:
# #         """
# #         echo " Running bam2nuc to create nucleotide frequency report" >> {log}
# #         bam2nuc --dir {params.work_dir} --genome_folder {params.genome} {input.bam}

# #         echo "Running bismark2report to create sample report" >> {log}
# #         bismark2report --dir {params.work_dir}

# #         cho "Moving methylatiuon extractor reports to report directory {params.meth_report_dir}" >> {log}
# #         mkdir -p {params.meth_report_dier}
# #         mv -f -v --target-directory={params.meth_report_dir} {input.mbias_report} {input.splitting_report}
# #         echo "Moving summary reports to report directory {params.report_dir}" >> {log}
# #         mkdir -p {params.report_dir}
# #         mv -f -v --target-directory={params.report_dir} {params.work_dir}/*_merged_SE_report.html {params.work_dir}/*_merged_SE_report.txt
# #         echo "Done" >> {log}
# #         """

# # bismark2report --dir test_bis_report --dedup_report ~/home/msleeper/scratch/reports/04_bismark_deduplication/318/1-crc-SRX8409040-se/SRR11858959_bismark.deduplication_report.txt --alignment_report ~/home/msleeper/scratch/reports/03_bismark_bwt2/318/1-crc-SRX8409040-se/SRR11858959_trimmed_bismark_bt2_SE_report.txt
