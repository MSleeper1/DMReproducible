#! /bin/bash -login
#SBATCH -p low                                  # partition requested
#SBATCH -J blocks                               # job name
#SBATCH -t 3-10:00:00                           # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=40gb                              # memory pool to all cores
#SBATCH -e slurm.crossover.j%j.err              # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.crossover.j%j.out              # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS

# OBSERVED RESOURCE USAGE:

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
#
# sbatch beta_prep_for_plots.sh
#
#-----------------------------------------------------------------------------------#

# activate conda
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment
mamba activate dmr

# make things fail on errors
set -o nounset  # same as 'set -u' which treats unset variables as an error when substituting
set -o errexit  # same as 'set -e' which exits immediately if a command exits with a non-zero status
set -x          # same as set '-o xtrace' which prints commands and their arguments as they are executed

# create output and tmp dirs
mkdir -p /home/msleeper/scratch/data/07_wgbstools_betas_bwa/blocks
mkdir -p /home/msleeper/scratch/data/07_wgbstools_betas_bwa/avg_meth_tables

# navigate to directory 
cd /home/msleeper/scratch/data/07_wgbstools_betas_bwa/

# # random already run as sample (CRC tissues plus cfdna)
# # wgbstools segment -o blocks.271.438.896.bed --betas 438cf_betas/*beta 896_betas/*beta 271_betas/*.beta
# # mv blocks.271.438.896.bed blocks/blocks.271.438.896.bed
# # used 18 / 80 G mem
wgbstools index blocks/blocks.271.438.896.bed
wgbstools beta_to_table blocks/blocks.271.438.896.bed.gz --betas 438cf_betas/*beta 896_betas/*beta 271_betas/*.beta --output avg_meth_tables/avg_meth.271.438.896.tsv

# # # all refs
# # # 171 215 271 318 438 535 644 783 896
# # wgbstools segment -o blocks.171.215.271.318.438.535.644.783.896.bed --betas 171_betas/*.beta 215_betas/*.beta 271_betas/*.beta 318_betas/*.beta 438cf_betas/*.beta 535_betas/*.beta 644_betas/*.beta 783_betas/*.beta 896_betas/*beta 
# # mv blocks.171.215.271.318.438.535.644.783.896.bed blocks/blocks.171.215.271.318.438.535.644.783.896.bed
# wgbstools index blocks/blocks.171.215.271.318.438.535.644.783.896.bed 
# wgbstools beta_to_table blocks/blocks.171.215.271.318.438.535.644.783.896.bed.gz --betas 171_betas/*.beta 215_betas/*.beta 271_betas/*.beta 318_betas/*.beta 438cf_betas/*.beta 535_betas/*.beta 644_betas/*.beta 783_betas/*.beta 896_betas/*beta --output avg_meth_tables/avg_meth.all_studies.171.215.271.318.438.535.644.783.896.tsv

# # # all colon tissue (including norm and cancer)
# # # 171 271 318 535 644 783 896
# # wgbstools segment -o blocks.colon_tissue_studies.171.271.318.535.644.783.896.bed --betas 171_betas/*.beta 271_betas/*.beta 318_betas/*.beta 535_betas/*.beta 644_betas/*.beta 783_betas/*.beta 896_betas/*beta 
# # mv blocks.colon_tissue_studies.171.271.318.535.644.783.896.bed blocks/blocks.colon_tissue_studies.171.271.318.535.644.783.896.bed
# wgbstools index blocks/blocks.colon_tissue_studies.171.271.318.535.644.783.896.bed
# wgbstools beta_to_table blocks/blocks.colon_tissue_studies.171.271.318.535.644.783.896.bed.gz --betas 171_betas/*.beta 271_betas/*.beta 318_betas/*.beta 535_betas/*.beta 644_betas/*.beta 783_betas/*.beta 896_betas/*beta --output avg_meth_tables/avg_meth.colon_tissue_studies.171.271.318.535.644.783.896.tsv

# # different tissues and cfdna (does cfdna with cancer group with specific tissues?)
# # 215 438
# # wgbstools segment -o blocks.tissue_origin.215.438.bed --betas 215_betas/*beta 438cf_betas/*beta
# # mv blocks.tissue_origin.215.438.bed blocks/blocks.tissue_origin.215.438.bed
# wgbstools index blocks/blocks.tissue_origin.215.438.bed
# wgbstools beta_to_table blocks/blocks.tissue_origin.215.438.bed.gz --betas 215_betas/*beta 438cf_betas/*beta --output avg_meth_tables/avg_meth.tissue_origin.215.438.tsv

# # # old vs young colon normal
# # # 535 783
# # wgbstools segment -o blocks.old_v_young_colon.535.783.bed --betas 535_betas/*beta 783_betas/*beta
# # mv blocks.old_v_young_colon.535.783.bed blocks/blocks.old_v_young_colon.535.783.bed
# wgbstools index blocks/blocks.old_v_young_colon.535.783.bed 
# wgbstools beta_to_table blocks/blocks.old_v_young_colon.535.783.bed.gz --betas 535_betas/*beta 783_betas/*beta --output avg_meth_tables/avg_meth.old_v_young_colon.535.783.tsv

# # # norm v cancer tissue (colon)
# # # 271 318 644 896
# # wgbstools segment -o blocks.norm_v_crc.271.318.644.896.bed --betas 271_betas/*.beta 318_betas/*.beta 644_betas/*.beta 896_betas/*beta 
# # mv blocks.norm_v_crc.271.318.644.896.bed blocks/blocks.norm_v_crc.271.318.644.896.bed
# wgbstools index blocks/blocks.norm_v_crc.271.318.644.896.bed 
# wgbstools beta_to_table blocks/blocks.norm_v_crc.271.318.644.896.bed.gz --betas 271_betas/*.beta 318_betas/*.beta 644_betas/*.beta 896_betas/*beta --output avg_meth_tables/avg_meth.norm_v_crc.271.318.644.896.tsv

# # # pre crc and old v young
# # # 535 783 896
# # wgbstools segment -o blocks.old_v_young_precancer_v_cancer.535.783.896.bed  --betas 535_betas/*beta 783_betas/*beta 896_betas/*beta
# # mv blocks.old_v_young_precancer_v_cancer.535.783.896.bed blocks/blocks.old_v_young_precancer_v_cancer.535.783.896.bed 
# wgbstools index blocks/blocks.old_v_young_precancer_v_cancer.535.783.896.bed 
# wgbstools beta_to_table blocks/blocks.old_v_young_precancer_v_cancer.535.783.896.bed.gz --betas 535_betas/*beta 783_betas/*beta 896_betas/*beta --output avg_meth_tables/avg_meth.old_v_young_precancer_v_cancer.535.783.896.tsv

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM
