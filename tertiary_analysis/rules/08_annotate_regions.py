# Python script for annotating Markers output from wgbstools
# Uses a gtf file to search for genes that overlap with DMRs in Markers file
# Outputs a tsv with genes all listed in one column (genes are sepoarated by commas within the column)

# Use with the conda environment annotate

import pandas as pd
import pybedtools
import re
# import pyranges as pr

# --- Input files ---
input_file = "Markers.CRC.bed"      # One region per line: chr#:start-end
gtf_file = "hg38.knownGene.gtf"        # GENCODE GTF annotation
output_file = "regions_with_genes_combined.tsv"

# --- Promoter definition (in bp) --- (cite a source)
PROMOTER_UPSTREAM = 2000   # 2kb upstream
PROMOTER_DOWNSTREAM = 500  # 0.5kb downstream

# --- Step 1: Load regions into DataFrame ---
df = pd.read_csv(input_file, sep='\t')
# Rename columns in-place
df.rename(columns={'#chrom': 'chrom'}, inplace=True)

# --- Step 2: Convert regions to BED format ---
bed = pybedtools.BedTool.from_dataframe(df)

# --- Step 3: Load gene annotation and filter for "gene" features ---
genes = pybedtools.BedTool(gtf_file)
genes = genes.filter(lambda x: x[2] == "gene")

# --- Step 4: Intersect regions with genes ---
intersect = bed.intersect(genes, wa=True, wb=True)

# --- Step 5: Extract gene attributes (name, ID, type) ---
gene_results = []
for interval in intersect:
    chrom, start, end = interval[0], int(interval[1]), int(interval[2])
    region_str = f"{chrom}:{start}-{end}"

    gtf_attributes = interval[17]
    gene_name_match = re.search(r'gene_name "([^"]+)"', gtf_attributes)
    gene_id_match = re.search(r'gene_id "([^"]+)"', gtf_attributes)
    gene_type_match = re.search(r'gene_type "([^"]+)"', gtf_attributes)

    gene_name = gene_name_match.group(1) if gene_name_match else "NA"
    gene_id = gene_id_match.group(1) if gene_id_match else "NA"
    gene_type = gene_type_match.group(1) if gene_type_match else "NA"

    gene_results.append((region_str, gene_name, gene_id, gene_type))

# Convert intersections to DataFrame
intersect_df = pd.DataFrame(gene_results, columns=['region', 'gene_name', 'gene_id', 'gene_type'])

# --- Step 6: Aggregate per region ---
agg_df = (
    intersect_df
    .groupby('region')
    .agg({
        'gene_name': lambda x: ','.join(sorted(set(x))),
        'gene_id': lambda x: ','.join(sorted(set(x))),
        'gene_type': lambda x: ','.join(sorted(set(x)))
    })
    .reset_index()
)

# --- Step 7: Build promoter intervals ---
promoters = []
for feature in genes:
    chrom = feature.chrom
    strand = feature.strand
    gene_start = int(feature.start)
    gene_end = int(feature.end)
    gtf_attributes = feature[8]

    gene_name_match = re.search(r'gene_name "([^"]+)"', gtf_attributes)
    gene_name = gene_name_match.group(1) if gene_name_match else "NA"

    # Calculate TSS
    tss = gene_start if strand == "+" else gene_end
    promoter_start = max(0, tss - PROMOTER_UPSTREAM)
    promoter_end = tss + PROMOTER_DOWNSTREAM

    promoters.append([chrom, promoter_start, promoter_end, gene_name])

# Convert to BedTool
promoter_bed = pybedtools.BedTool(promoters)

# Intersect regions with promoters
promoter_intersect = bed.intersect(promoter_bed, wa=True, wb=True)

# Parse promoter results
promoter_results = []
for interval in promoter_intersect:
    chrom, start, end = interval[0], int(interval[1]), int(interval[2])
    region_str = f"{chrom}:{start}-{end}"
    promoter_gene = interval[6]
    promoter_results.append((region_str, promoter_gene))

# Aggregate promoter info
promoter_df = pd.DataFrame(promoter_results, columns=['region', 'promoter_gene'])
promoter_agg = promoter_df.groupby('region')['promoter_gene'].apply(lambda x: ','.join(sorted(set(x)))).reset_index()

# --- Step 8: Merge all annotations ---
final_df = (
    df.merge(agg_df, on='region', how='left')
      .merge(promoter_agg, on='region', how='left')
      .merge(df, on='region', how='left')
      .fillna('NA')
)

# Fill NA for regions with no genes
final_df['gene_name'] = final_df['gene_name'].fillna('NA')
final_df['gene_id'] = final_df['gene_id'].fillna('NA')

# --- Step 9: Save output ---
final_df.to_csv(output_file, sep='\t', index=False)

print(f"✅ Annotated file saved to {output_file}")
