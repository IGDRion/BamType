import pandas as pd
import numpy as np
import os
import csv

# Final file name
final_file = "/groups/dog/stage/victor/bamslam/csv/sample_features/all_samples.csv"

# List all csv files in the directory & create a new empty data frame
file_path = "/groups/dog/stage/victor/bamslam/csv/sample_features/"
file_list = os.listdir(file_path)

# If the final file already exist, get it off the list 
if "all_samples.csv" in file_list:
    file_list.remove("all_samples.csv")

df = pd.DataFrame()

# Merge all of the csv files of the directory
for file in file_list:
    df_temp = pd.read_csv(str(file_path + file))
    df = pd.concat((df, df_temp), ignore_index=True)

# Replace the sample column with juste the sample name (not the full path)
for i, row in df.iterrows():
    new_name = df.at[i,'sample'].rsplit("/",1)[1]
    df.at[i,'sample'] = new_name


# Adding a column for group (cancer type)
groups = map(lambda x: x.rsplit("_")[0], df['sample']) #map apply a fonction on a itterable object and return the result of the function in another itterable
df['group'] = list(groups)


# Adding a column for subgroup (resistant vs sensitive)
def subgroup_maker(sample:str):
    if sample.endswith(("R_R1_new_tr", "R1_R1_new_tr", "R2_R1_new_tr", "R3_R1_new_tr")):
        return "resistant"
    if sample.endswith(("S_R1_new_tr", "S1_R1_new_tr", "S2_R1_new_tr", "S3_R1_new_tr")):
        return "sensitive"   
    
subgroups = map(subgroup_maker, df['sample']) 
df['subgroup'] = list(subgroups)


# Merging with count_output.csv to get number of unmapped reads
df_unmapped = pd.read_csv("/groups/dog/stage/victor/count_output.csv",names=['sample1', 'nbr_unmapped_reads'])
df_unmapped = df_unmapped.loc[df_unmapped['sample1'].str.endswith("unmapped")]
df_unmapped['sample1'] = df_unmapped['sample1'].str.replace('-unmapped','')

df = pd.merge(df, df_unmapped, left_on='sample', right_on='sample1', how='left').drop('sample1', axis=1)


# Merging to get qscore of unampped reads
file_path = "/scratch/vlebars/unmapped_list/"
file_list = os.listdir(file_path)
sample_list = []
qscore_median_list = []


for file in file_list:
    df_unmapped_list = pd.read_csv(str(file_path + file))
    df_unmapped_list.columns = ['read']

    sample = file.split("_R1_new")[0]

    df_qscore = pd.read_csv(str("/groups/dog/nanopore/lncrna_resist_cgo/primary/human/" + sample + "/cdna_SQK-DCS109/flowcell/guppy_6.0.0/sequencing_summary.txt"), usecols=["read_id","mean_qscore_template"], sep="\t")
    df_qscore_unmapped = pd.merge(df_unmapped_list, df_qscore, left_on='read', right_on='read_id', how='left').drop('read_id', axis=1)
    
    sample_list.append(file.split("-unmapped")[0])
    qscore_median_list.append(df_qscore_unmapped['mean_qscore_template'].median())
    
df_qscore_unmapped_final = pd.DataFrame(list(zip(sample_list, qscore_median_list)), columns=['sample01', 'unmapped_qscore'])

df = pd.merge(df, df_qscore_unmapped_final, left_on='sample', right_on='sample01', how='left').drop('sample01', axis=1)

# Create new csv file
df.to_csv(final_file, sep=",", index=False)