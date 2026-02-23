# for beta file run wgbstools segment, index, beta to tables, 

# create ref array
ref_array=("318" "215" "171" "644" "271" "535" "783" "438cf" "896")

# loop through ref array
for ref in "${ref_array[@]}"; do
    # check if the file exists
    if [ -f "${ref}*beta" ]; then
        echo "segmenting ${ref} betas"
        wgbstools segment --betas ${ref}*beta --min_cpg 3 --max_bp 2000 -o blocks.${ref}.bed
    else
        echo "File ${ref} betas do not exist, skipping."
    fi
done

# loop through ref array to run wgbstools index
# wgbstools index blocks.[ref].bed
for ref in "${ref_array[@]}"; do
    # check if the file exists
    if [ -f "blocks.${ref}.bed" ]; then
        echo "indexing ${ref} blocks"
        wgbstools index blocks.${ref}.bed
    else
        echo "File blocks.${ref}.bed does not exist, skipping."
    fi
done

# loop through ref array to run beta_to_table
# wgbstools beta_to_table blocks.[ref].bed --betas *beta | column -t > avg_meth.[ref].tsv
for ref in "${ref_array[@]}"; do
    # check if the file exists
    if [ -f "blocks.${ref}.bed.gz" ]; then
        echo "running beta_to_table on ${ref} blocks"
        wgbstools beta_to_table blocks.${ref}.bed.gz --betas ${ref}*beta | column -t > avg_meth.${ref}.tsv
    else
        echo "File blocks.${ref}.bed does not exist, skipping."
    fi
done

### running everything above

# running wgbstools homog
# wgbstools homog [ref]*pat.gz -b blocks.[ref].bed -o homog_out --thresholds 0.25,0.75
# loop through ref array to run homog
for ref in "${ref_array[@]}"; do
    # check if the file exists
    if ls ${ref}*pat.gz 1> /dev/null 2>&1; then
        echo "running homog on ${ref} blocks"
        wgbstools homog ${ref}*pat.gz -b blocks.${ref}.bed -o homog_out --thresholds 0.25,0.75
    else
        echo "File ${ref}*pat.gz does not exist, skipping."
    fi
done
