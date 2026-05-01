### Genome preparation for fastq_screen ###

# ref_prep_fastqscreen: rule to download the fastq_screen genomes that have been already been prepared with bismark genome preparation
# This rule will download the fastq_screen genomes that have been already been prepared with bismark genome preparation
# Fastq_screen will create the directory and download the genomes to the directory and create a fastq_screen.conf file
# The fastq_screen.conf file will be used to specify the genomes to be used in the fastq_screen rule
# default genomes included by fastqc are human, mouse, rat, drosophilia, worm, yeast, arabdopsis, ecoli, and PhiX
# the default genomes can be changed by editing the fastq_screen.conf file
# An alternative method is to download the genomes of interest from ncbi and prepare them with bismark_genome_preparation

rule ref_prep_fastqscreen:
    output:
        directory(expand("{root}/{genomes_dir}/{genome}/FastQ_Screen_Genomes_Bisulfite/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"])),
        expand("{root}/{genomes_dir}/{genome}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"])

    log:
        "logs/initialize_rules/ref_prep_fastqscreen.err"

    conda:
        "../../environment_files/fastq-screen.yaml"

    threads: 4
    
    resources:
        mem_mb=8000

    params:
        genomes_dir = expand("{root}/{genomes_dir}/{genome}", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"])
    
    shell:
        """
        mkdir -p {params.genomes_dir}
        cd {params.genomes_dir}
        fastq_screen --get_genomes --bisulfite --quiet --threads {threads}
        """

# if error occuring with fastq_screen get genomes command
# fastq_screen --bisulfite --get_genomes 
# work around by downloading them directlt with wget
# wget -r -np -nH --cut-dirs=3 https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/Genome_Data/FastQ_Screen_Genomes_Bisulfite/
# conf file needs to be adjusted if you use this work around       