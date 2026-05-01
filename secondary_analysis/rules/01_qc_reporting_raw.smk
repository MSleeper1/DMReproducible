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
    
    threads: 4

    resources:
        mem_mb=8000

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

    threads: 4

    resources:
        mem_mb=8000

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
    
    threads: 4

    resources:
        mem_mb=8000
    
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
    
    threads: 4

    resources:
        mem_mb=8000
    
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
