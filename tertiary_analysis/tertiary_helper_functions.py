import pandas as pd

#### Python functions for snakemake pipeline ####

''' Read sample info table into pandas dataframe and confirm column names that are needed for the pipeline'''
def get_sample_info_df(input_tsv):
    sample_info = pd.read_table(input_tsv, dtype=str)
    # Verify column names
    if not {'ref', 'patient_id', 'srx_id', 'accession', 'layout', 'group', 'age', 'sex', 'ethnicity'}.issubset(sample_info.columns.values):
        raise KeyError("The sample info file must contain the following named columns: 'ref', 'patient_id', 'srx_id', 'accession', 'layout', 'group', 'age', 'sex' 'ethnicity'")
    sample_info = sample_info.set_index(["srx_id"], drop=False)
    return sample_info

# ''' Get list of accessions for a given srx_id and determine if a merge is needed for bam files associated with that srx_id'''
# def accession_list_by_srx_id(sample_info_df, srx_id):
#     srx_samples = sample_info_df[sample_info_df['srx_id'] == srx_id]
#     accessions = srx_samples['accession'].tolist()
#     # if len(accessions) > 1:
#     #     print("More than 1 run provided for srx_id = {}: \n Merge needed for bam files of accessions = {}".format(srx_id, accessions))
#     # if len(accessions) == 1:
#     #     print("Only 1 run provided for srx_id = {}. \n Merge not needed for accession = {}".format(srx_id, accessions))
#     if len(accessions) == 0:
#         raise KeyError("ERROR: No accessions provided for runs from srx_id = {}".format(srx_id))
#     return accessions

# ''' Make a dictionary of srx_id and associated accessions'''
# def make_srx_accession_dict(sample_info_df):
#     srx_id_list = sample_info_df['srx_id'].unique().tolist()
#     srx_acc_dict = {}
#     for srx_id in srx_id_list:
#         accessions = accession_list_by_srx_id(sample_info_df, srx_id)
#         srx_acc_dict[srx_id] = (accessions)
#     return srx_acc_dict

# ''' Make a dictionary of srx_id and associated pandas dataframe rows'''
# def make_srx_df_dict(sample_info_df):
#     srx_id_list = sample_info_df['srx_id'].unique().tolist()
#     srx_df_dict = {}
#     for srx_id in srx_id_list:
#         srx_df_dict[srx_id] = sample_info_df[sample_info_df['srx_id'] == srx_id]
#     return srx_df_dict

# ''' drop the accessions column from the dataframe and merge the rows into a single row'''
# def merge_rows_by_srx_id(srx_df_dict):
#     merged_df_dict = {}
#     for srx_id, df in srx_df_dict.items():
#         df = df.drop(columns=['accession'])
#         df = df.drop_duplicates()
#         merged_df_dict[srx_id] = df
#         if len(df) > 1:
#             raise ValueError("ERROR: All df samlpe details for srx_id = {} should be the same except accession because srx_id represents an experiment and is unique to one sample. This allows for sample information to be collapsed by srx_id when accessions are removed. Please check the sample info file for discrepancies in columns with srx_id={}.".format(srx_id, srx_id))
#     return merged_df_dict

# ''' Return a reduced df for for sample info by srx_id'''
# def get_sample_info_by_srx_id(sample_info_df):
#     sample_info_by_srx = sample_info_df.drop(columns=['accession'])
#     sample_info_by_srx = sample_info_by_srx.drop_duplicates()
#     return sample_info_by_srx