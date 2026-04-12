#! /bin/bash -login
#SBATCH -p low                                  # partition requested
#SBATCH -J crossmap                             # job name
#SBATCH -t 3-10:00:00                           # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=20gb                              # memory pool to all cores
#SBATCH -e slurm.crossover.j%j.err              # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.crossover.j%j.out              # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 
 
# OBSERVED RESOURCE USAGE:

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# crossmap_alignment_hg38Tohg19.sh
#
# input bam files
# outputs crossmapped bam files in crossover directory within the bam file directory
# 
# options:
#   -w indicates where to navigate to for access to bam files
#   non-option args should be the bam files being processed 
#
# sbatch crossmap_alignment_hg38Tohg19.sh -w [directory/path] [bam-file-1] [bam-file-2] bam-file-3]...
# 
#-----------------------------------------------------------------------------------#
 
# activate conda 
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"
 
# activate a specific conda environment
mamba activate crossmap
 
# make things fail on errors
set -o nounset  # same as 'set -u' which treats unset variables as an error when substituting
set -o errexit  # same as 'set -e' which exits immediately if a command exits with a non-zero status
set -x          # same as set '-o xtrace' which prints commands and their arguments as they are executed

# defining flags for arguments passed in
while getopts ':w:' flag

do
    case "${flag}" in
        w) where=${OPTARG}
           ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
           ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
           ;;
    esac
done

# navigate to directory specified by -w flag for where
echo "Going to directory: $where"
cd $where
mkdir -p crossmap_hg38Tohg19

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# creating an array to keep track of marked bam files that are created
FILES=()

# iterate over files to convert, sort, and mark duplicates
for i; do
    
    # crossmap each bam file
    echo "Running CrossMap for $i to output $i.hg38Tohg19.bam"
    CrossMap bam -a ~/scratch/genomes/liftover/hg38ToHg19.over.chain $i $i.hg38Tohg19
    FILES+=($i'.hg38Tohg19.bam')
    mv $i.hg38Tohg19.* crossmap_hg38Tohg19/

done

# letting the user know where to find the sorted and marked files
echo "all files have been crossmapped"
echo "${FILES[@]} can be found in $where/crossmap_hg38Tohg19"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM
