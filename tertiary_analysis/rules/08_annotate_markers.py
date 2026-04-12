# Conda environment to use: annot

# References for sources used to write this script:
#   - Gene annotation methods: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2
#   - Gencode annotation file format: https://www.gencodegenes.org/pages/data_format.html

# Information about interpretting the chromosome labels
# https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=


#-----------------#
# Import packages #
#-----------------#
import pandas as pd
import numpy as np
import pysam

#-------------------#
# Assign varieables #
#-------------------#
# Assign directories and files to use for annotation references
annotation_file = "gencode.v44.annotation.all.tsv.sorted.formatted.bed.gz"
annotation_file_genes = "gencode.v44.annotation.genes.tsv.sorted.formatted.bed.gz"
annotation_file_regulatory = "ensembl.GRch38.p14.regulatory.tsv.sorted.formatted.bed.gz"

# DMR regions to annotate (tab-delimited bed formatted file)
DMR_file = "Markers.CRC.bed"

# Output file names
output_file_per_gene = "colon_cancer_markers_annotated_per_gene_new.tsv"
output_file = "colon_cancer_markers_annotated_new.tsv"
output_file_per_reg = "colon_cancer_markers_annotated_per_reg_new.tsv"


#-------------------#
#  Define Functions #
#-------------------#

### Functions to Read in reference and DMR files ###

def gencode_genes_ref_to_df(gencode_genes_bed):
    # Reference for genes only: 
    # Read reference annotation file with tabix
    gencode_v44_genes = pysam.TabixFile(gencode_genes_bed)
    # read in the annotation file as a pandas dataframe
    annot_df_gene = pd.read_csv(gencode_genes_bed, sep='\t', header=None, names=['chr', 'start', 'end', 'attribute', 
                                                                                 'feature', 'strand', 'source', 'gene_name', 
                                                                                 'gene_id', 'gene_type', 'gene_level'])
    # drop unwanted columns from annotation dfs
    # genes only:
    annot_df_gene = annot_df_gene.drop(['attribute'], axis=1)
    # return cleaned dfs for gene annotation
    return gencode_v44_genes, annot_df_gene

def regulatory_ref_to_df(ensembl_regulatory_bed):
    # Reference for regulatory features from ensembl regulatory build:
    # Read reference annotation file with tabix
    regulatory_build = pysam.TabixFile(ensembl_regulatory_bed)
    annot_df_regulatory = pd.read_csv(ensembl_regulatory_bed, 
                                    sep='\t', 
                                    header=0, 
                                    low_memory=False, 
                                    names=['chr', 'start', 'end',
                                            'feature_type', 'feature_type_description'])
    
    return regulatory_build, annot_df_regulatory

def gencode_all_to_ref(gencode_all_annotation_file):
    # Reference for all features: 
    # Read reference annotation file with tabix
    gencode_v44 = pysam.TabixFile(gencode_all_annotation_file)
    # read in the annotation file as a pandas dataframe
    annot_df = pd.read_csv(gencode_all_annotation_file, sep='\t', header=None, names=['chr', 'start', 'end', 'attribute', 
                                                                                      'feature', 'strand', 'source', 'gene_name', 
                                                                                      'gene_id', 'gene_type', 'gene_level', 'exon_id', 
                                                                                      'exon_number', 'transcript_name','transcript_type'])
    # drop unwanted columns
    annot_df = annot_df.drop(['attribute'], axis=1)
    # return cleaned dfs for annotation of all features
    return gencode_v44, annot_df


def dmr_to_df(DMR_file):
    ### Read the DMR sample file in to a pandas dataframe
    df = pd.read_table(DMR_file, low_memory=False, header=0)

    # Clean DMR file
    # remove CpGs from lenCpG values dataframe
    df['lenCpG'] = df['lenCpG'].apply(lambda x: x.split("CpGs")[0])
    # remove bp from bp values in dataframe
    df['bp'] = df['bp'].apply(lambda x: x.split("bp")[0])
    # rename chromosome column
    df.rename(columns = {'#chr':'chr'}, inplace = True)
    # uncomment below line if you want to drop unwanted columns
    #df = df.drop(['startCpG', 'endCpG', 'len_bp', 'region', 'tg_mean', 'bg_mean', 'delta_quants', 'delta_maxmin'], axis=1)
    
    return df

### Identify overlap ###
    # Determine if there is overlap between genes and DMRs
    # query and reference
def overlap(q_st, q_end, ref_st, ref_end):
    overlap  = min(q_end, ref_end)-max(q_st, ref_st)
    return overlap

### Gencode annotate dmrs ###
    # annotate dmrs using gencode_annotate_loci function and gencode_v44 tabix annotation file
    # returns genes, feature, gene_type, gene_level, exon_id, exon_number, transcript_name, transcript_type that overlap with DMRs
    # returns amount of overlap as a percentage under gene column as gene-name(overlap %)(overlap %)
    # # loci is a row in a dataframe `a` with columns chr, start, end ['#chr', 'start', 'end']
    # # `tb` tabix is the indexed annotation file as a pysam.TabixFile object
def gencode_annotate_loci(a, tb):
    genes = []
    gene_type = []
    gene_level = []
    gene_id = []
    strand = []

    try:
        # reads the tabix file
        # columns in tabix file: ['chr', 'start', 'end', 'attribute', 'feature', 
                                # 'strand', 'source', 'gene_name', 'gene_id','
                                # gene_type', 'gene_level'])
        
        for region in tb.fetch(a['chr'], int(a['start']), int(a['end'])):
            if region:
                r = region.split('\t')

                ### Calculate length of the overlap between the query interval and the gencode gene
                overlap_len = overlap(int(a['start']), int(a['end']), int(r[1]), int(r[2]))

                ### first value is the percentage of the query interval that overlap with the gencode gene
                ### second value is the percentage of the gencode gene that overlaps with the query interval
                ret_val = '{}({})({})'.format(r[7], 
                                              np.round(overlap_len/float(int(a['end'])-int(a['start']))*100, 2), 
                                              np.round(overlap_len/float(int(r[2])-int(r[1]))*100, 2)
                                              ) 

                genes.append(ret_val) 
                gene_type.append(r[9])
                gene_level.append(r[10])
                gene_id.append(r[8])
                strand.append(r[5])
    
        if len(genes)>0:
            return [";".join(genes), ";".join(gene_type), ";".join(gene_level), ";".join(gene_id), ";".join(strand)]
        else:
            return ["NA(0)(0)", "NA", "NA", "NA", "NA"]
    
    except ValueError:
        return ["NA(0)(0)", "NA", "NA", "NA", "NA"]
    
### Regulatory Annotation ###
    # annotate dmrs using regulatory_annotate_loci function and regulatory build tabix annotation file
    # returns regulatory feature types that overlap with DMRs
    # returns amount of overlap as a percentage under regulatory_feature column as feature-name(overlap %)(overlap %)
    # # loci is a row in a dataframe `a` with columns chr, start, end ['#chr', 'start', 'end']
    # # `tb` tabix is the indexed annotation file as a pysam.TabixFile object
def regulatory_annotate_loci(a, tb):
    reg_features = []
    reg_feature_description = []

    try:
        # reads the tabix file
        # columns in tabix file: ['chr', 'start', 'end', 'feature_type', 'feature_type_description']
        
        for region in tb.fetch(a['chr'], int(a['start']), int(a['end'])):
            if region:
                r = region.split('\t')

                ### Calculate length of the overlap between the query interval and the gencode gene
                overlap_len = overlap(int(a['start']), int(a['end']), int(r[1]), int(r[2]))

                ### first value is the percentage of the query interval that overlap with the gencode gene
                ### second value is the percentage of the gencode gene that overlaps with the query interval
                ret_val = '{}({})({})'.format(r[3], 
                                              np.round(overlap_len/float(int(a['end'])-int(a['start']))*100, 2), 
                                              np.round(overlap_len/float(int(r[2])-int(r[1]))*100, 2)
                                              ) 

                reg_features.append(ret_val)
                reg_feature_description.append(r[4])
    
        if len(reg_features)>0:
            return [";".join(reg_features), ";".join(reg_feature_description)]
        else:
            return ["NA(0)(0)", "NA"]
    except ValueError:
        return ["NA(0)(0)", "NA"]
    
# Function to count the number of genes found for each region 
# by counting the number of semicolons in the gene_name or regulatory_feature column
def count_features(gene):
    if gene == 'NA(0)(0)':
        return 0
    else:
        return gene.count(";") + 1
    
# Function to return if region is intergenic or intragenic
def genic(gene):
    if gene == 'NA(0)(0)':
        return "intergenic"
    else:
        return "intragenic"

# Function to return if region has regulatory features
def regulatory(regulatory_feature):
    if regulatory_feature == 'NA(0)(0)':
        return "no"
    else:
        return "yes"
    
# write a function to split the deliminates gene columns into lists
def split_gene(gene_name, gene_type, gene_level, gene_id, strand):
    gene_list = gene_name.split(";")
    gt_list = str(gene_type).split(";") 
    gl_list = str(gene_level).split(";")
    gid_list = str(gene_id).split(";")
    strand_list = str(strand).split(";")
    return gene_list, gt_list, gl_list, gid_list, strand_list

def split_regulatory(regulatory_feature, reg_feature_description):
    regulatory_list = regulatory_feature.split(";")
    reg_feature_description_list = reg_feature_description.split(";")
    return regulatory_list, reg_feature_description_list

def count_genes(df):
    # count the number of genes found for each region and add to g_nums list
    g_nums = []
    for i,r in df.iterrows():
        g_num = count_features(str(r['gene_name']))
        g_nums.append(g_num)

    # add gene number column to dataframe
    df['gene_num'] = g_nums
    print("Gene counts added to dataframe.")

    # return max number of genes found in one region and preview of regions with max number of genes
    max_g = max(g_nums)
    print("Max number of genes in one region is {}.".format(max_g))
    print("Regions with max number of genes:")
    df.loc[(df['gene_num']==max_g)].head()
    return df

def count_regulatory_features(df):
    # count the number of regulatory features found for each region and add to r_nums list
    reg_nums = []
    for i,r in df.iterrows():
        reg_num = count_features(str(r['regulatory_feature']))
        reg_nums.append(reg_num)

    # add gene number column to dataframe
    df['reg_num'] = reg_nums
    print("Regulatory feature counts added to dataframe.")

    # return max number of genes found in one region and preview of regions with max number of genes
    max_r = max(reg_nums)
    print("Max number of regulatory features in one region is {}.".format(max_r))
    print("Regions with max number of regulatory features:")
    df.loc[(df['reg_num']==max_r)].head()
    return df

### split out the genes to multiple rows (regions are no longer a unique identifier)
def dup_regions_to_split_genes(df):
    # split values in the columns gene, gene_type, gene_level by delimiter ; and add to a list of lists called new_rows
    new_rows = []
    for i,r in df.iterrows():
        g_list, gt_list, gl_list, gid_list, strand_list = split_gene(r['gene_name'], r['gene_type'], r['gene_level'], r['gene_id'], r['strand'])

        for g in range(len(g_list)):
            new_rows.append(np.append(r[['chr', 'start', 'end', 'startCpG', 
                                        'endCpG', 'target', 'region', 'lenCpG', 
                                        'bp', 'tg_mean', 'bg_mean', 'delta_means', 
                                        'delta_quants', 'delta_maxmin', 'ttest', 
                                        'direction', 'gene_name', 'gene_type',
                                        'gene_level', 'gene_id', 'strand', 'genic', 'gene_num']
                                        ].values, [g_list[g], gt_list[g], gl_list[g], 
                                                    gid_list[g], strand_list[g]
                                                    ]))
    # # Check df against new rows to confirm proper splitting
    # df.head()
    # print(len(new_rows))
    # print(new_rows[0:5])
    # Add new rows to add to the dataframe

    df_perGene = pd.DataFrame(new_rows, columns=['chr', 'start', 'end', 'startCpG',
                                                'endCpG', 'target', 'region', 'lenCpG',
                                                'bp', 'tg_mean', 'bg_mean', 'delta_means',
                                                'delta_quants', 'delta_maxmin', 'ttest', 
                                                'direction', 'gene_name', 'gene_type',
                                                'gene_level', 'gene_id', 'strand', 
                                                'genic', 'gene_num', 'gene_name_split', 
                                                'gene_type_split', 'gene_level_split',
                                                'gene_id_split', 'strand_split'])
    
    # Split gene_name column into gene_name, region_coverage, gene_coverage
    df_perGene['gene_name_cleaned'] = df_perGene['gene_name_split'].apply(lambda x: x.split("(")[0].replace(")", ""))
    df_perGene['region_coverage'] = df_perGene['gene_name_split'].apply(lambda x: x.split("(")[1].replace(")", ""))
    df_perGene['gene_coverage'] = df_perGene['gene_name_split'].apply(lambda x: x.split("(")[2].replace(")", ""))

    # once I have confirmed proper annotation, I can clean up the df
    # columns in df_per_gene: ['chr', 'start', 'end', 'startCpG', 'endCpG', 'target', 'region',
                        #    'lenCpG', 'bp', 'tg_mean', 'bg_mean', 'delta_means', 'delta_quants',
                        #    'delta_maxmin', 'ttest', 'direction', 'gene_name', 'gene_type',
                        #    'gene_level', 'gene_id', 'strand', 'gene_num', 'gene_name_split',
                        #    'gene_type_split', 'gene_level_split', 'gene_id_split', 'strand_split',
                        #    'gene_name_cleaned', 'region_coverage', 'gene_coverage']

    ## drop the genes column
    df_perGene = df_perGene.drop(['gene_name', 'gene_type','gene_level', 
                                'gene_id', 'strand', 'gene_name_split'], axis=1)
    # rename columns
    df_perGene.rename(columns = {'gene_name_cleaned':'gene_name', 
                                'gene_id_split':'gene_ID',
                                'strand_split':'strand',
                                'gene_type_split':'gene_type',
                                'gene_level_split':'gene_level',
                                'gene_num':'num_genes',
                                'gene_coverage':'perc_gene_cov_region',
                                'region_coverage':'perc_region_cov_gene'
                                }, inplace = True)
    # reorder columns
    df_perGene = df_perGene.loc[:,['chr', 'start', 'end', 'startCpG', 
                                    'endCpG', 'target', 'region', 'lenCpG', 
                                    'bp', 'tg_mean', 'bg_mean', 'delta_means', 
                                    'delta_quants', 'delta_maxmin', 'ttest', 
                                    'direction', 'num_genes', 'gene_name', 
                                    'perc_gene_cov_region', 'perc_region_cov_gene',
                                    'gene_type', 'gene_level', 'gene_ID', 'strand', 
                                    'genic']]
    return df_perGene


### split out the regulatory features to multiple rows (regions are no longer a unique identifier)
def dup_regions_to_split_regulatory_features(df):
    # split values in the columns gene, gene_type, gene_level by delimiter ; and add to a list of lists called new_rows
    new_rows = []
    for i,r in df.iterrows():
        reg_list, reg_desc_list = split_regulatory(r['regulatory_feature'], r['regulatory_feature_description'])

        for reg in range(len(reg_list)):
            new_rows.append(np.append(r[['chr', 'start', 'end', 'startCpG', 
                                        'endCpG', 'target', 'region', 'lenCpG', 
                                        'bp', 'tg_mean', 'bg_mean', 'delta_means', 
                                        'delta_quants', 'delta_maxmin', 'ttest', 
                                        'direction', 'gene_name', 'gene_type',
                                        'gene_level', 'gene_id', 'strand', 'genic', 
                                        'gene_num','regulatory_feature', 'regulatory_feature_description', 'reg_num', 
                                        'regulatory']
                                        ].values, [reg_list[reg], reg_desc_list[reg]
                                                   ]))
    # # compare the df against new rows to confirm proper splitting
    # df.head()
    # print(len(new_rows))
    # print(new_rows[0:10])

    # Add new rows to add to the dataframe (df_perReg)
    df_perReg = pd.DataFrame(new_rows, columns=['chr', 'start', 'end', 'startCpG', 
                                                'endCpG', 'target', 'region', 'lenCpG', 
                                                'bp', 'tg_mean', 'bg_mean', 'delta_means', 
                                                'delta_quants', 'delta_maxmin', 'ttest', 
                                                'direction', 'gene_name', 'gene_type',
                                                'gene_level', 'gene_id', 'strand', 'genic', 
                                                'gene_num','regulatory_feature', 'regulatory_feature_description', 'reg_num', 
                                                'regulatory', 'reg_split', 'reg_desc_split'])

    # Split regulatory_feature column into reg feature type, region_coverage, feature_coverage
    df_perReg['reg_name_cleaned'] = df_perReg['reg_split'].apply(lambda x: x.split("(")[0].replace(")", ""))
    df_perReg['perc_reg_cov_region'] = df_perReg['reg_split'].apply(lambda x: x.split("(")[1].replace(")", ""))
    df_perReg['perc_region_cov_reg'] = df_perReg['reg_split'].apply(lambda x: x.split("(")[2].replace(")", ""))
    
    # Check columns and clean
    # df_perReg.columns
    # once I have confirmed proper annotation, I can clean up the df
    # columns in df_perReg:  ['chr', 'start', 'end', 'startCpG', 'endCpG', 'target', 'region',
                        #    'lenCpG', 'bp', 'tg_mean', 'bg_mean', 'delta_means', 'delta_quants',
                        #    'delta_maxmin', 'ttest', 'direction', 'gene_name', 'gene_type',
                        #    'gene_level', 'gene_id', 'strand', 'genic', 'gene_num',
                        #    'regulatory_feature', 'regulatory_feature_description', 'reg_num', 'regulatory', 'reg_split', 'reg_desc_split',
                        #    'reg_name_cleaned', 'perc_reg_cov_region', 'perc_region_cov_reg']

    ## drop the genes column
    df_perReg = df_perReg.drop(['regulatory_feature', 'regulatory_feature_description'], axis=1)

    # rename columns
    df_perReg.rename(columns = {'reg_name_cleaned':'regulatory_feature',
                                'reg_desc_split':'regulatory_description',
                                'reg_num':'num_regs',
                                'gene_num':'num_genes'
                                }, inplace = True)

    # reorder columns
    df_perReg = df_perReg.loc[:,['chr', 'start', 'end', 'startCpG', 
                                'endCpG', 'target', 'region', 'lenCpG', 
                                'bp', 'tg_mean', 'bg_mean', 'delta_means',
                                'delta_quants', 'delta_maxmin', 'ttest',
                                'direction', 'num_regs', 'regulatory_feature', 
                                'perc_reg_cov_region', 'perc_region_cov_reg', 
                                'num_genes', 'gene_name', 'gene_type', 
                                'gene_level', 'gene_id', 'strand', 'genic', 
                                'regulatory', 'regulatory_description']]
    return df_perReg







### MAIN ###

if __name__ == '__main__':
    
    # Read in and clean gencode genes only reference into df
    gencode_v44_genes, annot_df_gene = gencode_genes_ref_to_df(annotation_file_genes)

    # Read in and clean ensembl regulatory features into df
    regulatory_build, annot_df_regulatory = regulatory_ref_to_df(annotation_file_regulatory)

    ### Uncomment if you want to utilize more features in the gencode gtf
    # # Read in and clean gencode all features reference into df
    # gencode_v44, annot_df = gencode_all_to_ref(annotation_file)

    # Read in and clean DMR file
    df = dmr_to_df(DMR_file)

    # gene annotation with gencode v44
    df[['gene_name', 'gene_type', 'gene_level', 'gene_id', 'strand']] = df.apply(lambda x: gencode_annotate_loci(x[['chr', 'start', 'end']], gencode_v44_genes), axis=1, result_type='expand')
    # replace empty strings with NaN
    df = df.replace(r'^\s*$', np.nan, regex=True)

    # count number of genes in each region and add column 'gene_num'
    df = count_genes(df)

    # apply genic function to df to create a variable for intergenic vs intragenic
    df['genic'] = df['gene_name'].apply(genic)

    # annotate dmrs using regulatory_annotate_loci function and regulatory build tabix annotation file
    df[['regulatory_feature', 'regulatory_feature_description']] = df.apply(lambda x: regulatory_annotate_loci(x[['chr', 'start', 'end']], regulatory_build), axis=1, result_type='expand')
    # replace empty strings with NaN
    df = df.replace(r'^\s*$', np.nan, regex=True)

    # count the number of regulatory features in each region and add column 'reg_num'
    df = count_regulatory_features(df)

    # create a variable indicating if there is a regulatory feature in the region
    df['regulatory'] = df['regulatory_feature'].apply(regulatory)

    # duplicating regions to have one gene per row (regions are no longer a unique identifier)
    df_perGene = dup_regions_to_split_genes(df)

     # duplicating regions to have one regulatory feature per row (regions are no longer a unique identifier)
    df_perReg = dup_regions_to_split_regulatory_features(df)

    df.to_csv(output_file, index=False, header=True, sep="\t")
    df_perReg.to_csv(output_file_per_reg, index=False, header=True, sep="\t")
    df_perGene.to_csv(output_file_per_gene, index=False, header=True, sep="\t")

    print(f"✅ Annotated files saved to {output_file}, {output_file_per_reg}, and {output_file_per_gene}")






    ### Troublshooting and checking for accuracy of annotation from gtfs
    # checking proper annotation (uncomment below lines to check different things)

    ## Observing DMR annotations and comparing to annotation files
    # df.head()
    # annot_df_gene.head()
    # annot_df_regulatory.head()

    # # checking by chromosome
    # display(df.loc[(df['chr']=='chr4')])
    # display(annot_df_gene.loc[(annot_df_gene['chr']=='chr4')])
    # display(annot_df_regulatory.loc[(annot_df_gene['chr']=='chr4')])

    # # checking for specific genes
    # display(df.loc[(annot_df['gene_name']=='SEPTIN9') & (df['chr']=='chr17')])
    # display(annot_df_gene.loc[(annot_df_gene['gene_name']=='SEPTIN9') & (annot_df_gene['chr']=='chr17')])
    # display(annot_df_regulatory.loc[annot_df_regulatory['chr']=='chr17'])

    # # checking by chromosome and start/end
    # df_check = display(df.loc[(dmr['chr']=='chr17') & (df['start'].between(77280569,77500596))])
    # annot_check = display(annot_df_gene.loc[(annot_df_gene['chr']=='chr17') & (annot_df_gene['start'].between(77280569,77500596))])
    # annot_check_reg = display(annot_df_regulatory.loc[(annot_df_regulatory['chr']=='chr17') & (annot_df_regulatory['start'].between(77280569,77500596))])







