### liftover_betas.sh
# this script coverts betas to bed files and uses lift over to convert genome references
# outputs beta files after lifover conversion
# only accepts one beta file at a time
# must already have downloaded the chain file of choice

# example golden paths and file names for chain files / how to download and unzip
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# gunzip hg38ToHg19.over.chain.gz
# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
# gunzip hg19ToHg38.over.chain.gz

# my chain file dir: /home/msleeper/scratch/genomes/liftover

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':r:n:c:f:w:' flag
do
    case "${flag}" in
        r) ref_in=${OPTARG}
            ;;
        n) new_ref=${OPTARG}
            ;;
        c) chain_file=${OPTARG}
            ;;
        f) beta_in=${OPTARG}
            ;;
        w) where=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# navigate to directory with beta files specified by -w flag 
echo "Going to directory: $where"
cd $where

# activate wgbstools env
conda activate dmr

# convert beta to bed files
echo "running wgbstools beta2bed for $beta_in"
wgbstools beta2bed --genome $ref_in --outpath $beta_in.bed $beta_in 

# activate liftover env
conda activate liftover

# liftover bed file
echo "running liftover on $beta_in.bed"
liftOver $beta_in.bed $chain_file $beta_in.lift2$new_ref.bed $beta_in.UNMAPPED.bed

# activate wgbstools env
conda activate dmr

# convert bed to beta files
echo "running wgbstools bed2beta for $beta_in.lift2$new_ref.bed"
wgbstools bed2beta --genome $new_ref $beta_in.lift2$new_ref.bed
