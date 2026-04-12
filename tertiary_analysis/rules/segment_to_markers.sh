
#
# this script will:
#   -time alignment of fastq files to a reference 
#   -time conversion of resulting sam file to bam file
#   -time sorting bam file
#   -time marking duplicates in bam file
#
# input accession numbers for aligned bam files
#   - this should be the file name without extension output from sratools-workflow.sh
# outputs:
# alignment files in sam and sorted + marked bam files along with a final merged bam
#
# options:
#   -o indicates outfolder that will be created at location -w
#      output bam and sam files will be found here
#   -w indicates where to navigate to for access to file
#   -r indicates the reference genome to use
#   non-option args should be the accession numbers for fastq files being processed
#
# sbatch time-align2marked.sh -w [directory/path] -o [folder-name] -r [reference-genome] -f [accession-num-1] [accession-num-2] [accession-num-3]...
#
#-----------------------------------------------------------------------------------#


# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':r:f:w:b:m:' flag
do
    case "${flag}" in
        r) ref_num=${OPTARG}
            ;;
        f) files_beta=${OPTARG}
            ;;
        w) where=${OPTARG}
            ;;
        b) block_file=${OPTARG}
            ;;
        m) meth_table=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# example inputs
        # beta_dir = expand("{root}/{data_dir}/07_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"], sample = sample_info.itertuples()),
        # beta_files = expand("{{ref}}--{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.{suf}", suf=["beta"], sample = sample_info.itertuples()),
        # block_file = "blocks.{{ref}}.bed",
        # meth_table = "avg_meth.{{ref}}.tsv"

# navigate to directory with beta files specified by -w flag 
echo "Going to directory: $where"
cd $where

echo "running wgbstools segment to find_markers on beta files for study $ref_num"

# segmenting into homogenously methylated chunks of CpG sites according to the tutorial settings
echo "Segmenting beta files: $files_beta"
wgbstools segment --betas $files_beta --min_cpg 3 --max_bp 2000 -o $block_file  

# compress bed file and generate corresponding index file for visualization steps
echo " indexing blocks output by segment"
wgbstools index $block_file

# collapse beta files to homogenously methylated chunks found in segmentation step and writes to csv
echo "extracting average methylation per block with beta to table"
wgbstools beta_to_table $block_file.gz --betas $files_beta | column -t >> $meth_table

# NEED TO WRITE SCRIPT TO MAKE GROUPS FILE AND PASS INTO THIS SCRIPT

# # find DMRs for all 3 samples (lung, pancreas, and colon) according to tutorial settings
# # note that pval is set to 1 and change this after testing phase
wgbstools find_markers --blocks_path blocks.bed.gz --betas *beta --groups_file groups.csv --delta_quants .3 --pval 1
wgbstools find_markers --blocks_path blocks.[ref].bed.gz --betas [ref]*beta --groups_file groups.[ref].csv --delta_quants .3 --pval 1
wgbstools find_markers --blocks_path blocks.318.bed.gz --betas 318*beta --groups_file groups.318.csv --delta_quants .3 --pval 1



# # writes markers to csv file
# cat Markers.*.bed >> DMR_markers.csv


