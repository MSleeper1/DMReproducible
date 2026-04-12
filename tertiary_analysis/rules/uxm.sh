#! /bin/bash -login
#SBATCH -p high
#SBATCH -J uxm_build
#SBATCH -t 3-0:00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem=30gb


#-----------------------------------------------------------------------------------#
# uxm.sh
# Created: 1/31/26 MS
#
# accepts indexed pat files and creates an atlas based on groups for deconvolution
#
# Required inputs:
#   -o indicates name of output atlas
#   -w indicates where to navigate to for access to file
#   -m markers.bed file to use
#   -p indexed pat files to use (.pat.gz)
#
# sbatch bwameth-workflow.sh -w [directory/path] -o [output-atlas-name] -m [markers.bed] -g [groups_file] [pat.gz_files ...]
#
#-----------------------------------------------------------------------------------#


# activate conda
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment
conda activate dmr

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:m:g:' flag
do
    case "${flag}" in
        o) output=${OPTARG}
            ;;
        w) where=${OPTARG}
            ;;
        m) markers_bed=${OPTARG}
            ;;
        g) groups_file=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# creating an array to keep track of marked bam files that are created
FILES=()

# iterate over files to add them to array
for i; do

    # add completed accession number to FILES array
    FILES+=($i)

done

echo "array of files has been created: ${FILES[@]}"

# navigate to directory specified by -w flag for where
echo "Going to directory: $where"
cd $where

echo "running uxm_build"
# Example usage
# uxm build --where /home/msleeper/scratch/data/07_wgbstools_betas_bwa --markers 215/Markers.25.all.bed --groups groups.215.uxm.csv --output 215/Atlas.25.all_tissues --pats 215*.pat.gz
# sbatch uxm.sh -w /home/msleeper/scratch/data/07_wgbstools_betas_bwa -m 215/Markers.25.all.bed -g groups.215.uxm.csv -o 215/Atlas.25.all_tissues -p 215*.pat.gz

uxm build  --threads 12 --markers $markers_bed --output $output --groups $groups_file --pats ${FILES[@]}
echo "done with uxm_build"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM