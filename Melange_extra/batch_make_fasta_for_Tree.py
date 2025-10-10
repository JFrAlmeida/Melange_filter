import os
import itertools as itool
import pandas as pd
import re
from Bio import SeqIO
import time
from warnings import simplefilter



time_start = time.time()


# COG categories
COGS_GH18_endochitinases = ["COG3325", "COG3469"]

COGS_GH19_endochitinases = ["COG3179"]

COGS_exochi = ["COG4724", "COG3525", "COG3979"]

COGS_LPMO = ["COG3397"]

COGS_deacet = ["COG3394", "COG0726", "COG2861"]

#Paths storage
path_Pfam_table = "Outputs/PFAM_grouped_table.csv"

#Check if output folder exists, if not make one
if not os.path.exists("Outputs/Fasta_files"):
    os.makedirs("Outputs/Fasta_files")

###Import PFAM table and make one dataframe for each COG category
pfam_table = pd.read_csv(path_Pfam_table)

#Cut dataframe down into name and amino acid sequence of all enzymes
names_selected = list(pfam_table)
del names_selected[0] #eliminates the sequence name from the list
del names_selected[0] #eliminates COG column  from the list
del names_selected[1] #eliminates aa sequence column  from the list

#eliminate unwanted columns from pfam dataframe
pfam_table = pfam_table.drop(names_selected, axis=1)

#Cut down dataframe by COG category and cut it down into name, sequence;
GH18_endochitinases_df = pfam_table.loc[pfam_table["COG"].isin(COGS_GH18_endochitinases),:].copy()
GH19_endochitinases_df = pfam_table.loc[pfam_table["COG"].isin(COGS_GH19_endochitinases),:].copy()
COGS_exochi_df = pfam_table.loc[pfam_table["COG"].isin(COGS_exochi),:].copy()
COGS_LPMO_df = pfam_table.loc[pfam_table["COG"].isin(COGS_LPMO),:].copy()
COGS_deacet_df = pfam_table.loc[pfam_table["COG"].isin(COGS_deacet),:].copy()

#lists for announcements to be changed below!
no_cogs_found_list = [] #to store cogs with no enzymes found for an announcement at the end
cogs_found_list = [] #which cogs had results
#print me some fastas now boy!

#GH18_endochitinases
if len(GH18_endochitinases_df.index) > 0: # seq_list is ready to be exported into fasta
    Gh18_name_list = GH18_endochitinases_df["genome_prokka_feats"].to_list()
    Gh18_sequence_list = GH18_endochitinases_df["aa_sequence"].to_list()
    seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(Gh18_name_list, Gh18_sequence_list)]

    # write the fasta file with each item from short_seq_list
    with open("Outputs/Fasta_files/GH18_endochitinases.fasta", "w") as file:
        for item in seq_list:
            file.write(item)  # Add a newline after each item

    #empty out seq_list
    seq_list.clear()
    cogs_found_list = cogs_found_list + COGS_GH18_endochitinases

else:
    no_cogs_found_list = no_cogs_found_list + COGS_GH18_endochitinases

#GH19_endochitinases
if len(GH19_endochitinases_df.index) > 0: # seq_list is ready to be exported into fasta
    Gh19_name_list = GH19_endochitinases_df["genome_prokka_feats"].to_list()
    Gh19_sequence_list = GH19_endochitinases_df["aa_sequence"].to_list()
    seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(Gh19_name_list, Gh19_sequence_list)]

    # write the fasta file with each item from short_seq_list
    with open("Outputs/Fasta_files/GH19_endochitinases.fasta", "w") as file:
        for item in seq_list:
            file.write(item)  # Add a newline after each item
    seq_list.clear()
    cogs_found_list = cogs_found_list + COGS_GH19_endochitinases

else:
    no_cogs_found_list = no_cogs_found_list + COGS_GH19_endochitinases


#COGS_exochitinases
if len(COGS_exochi_df.index) > 0: # seq_list is ready to be exported into fasta
    COGS_exochi_name_list = COGS_exochi_df["genome_prokka_feats"].to_list()
    COGS_exochi_sequence_list = COGS_exochi_df["aa_sequence"].to_list()
    seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(COGS_exochi_name_list, COGS_exochi_sequence_list)]

    # write the fasta file with each item from short_seq_list
    with open("Outputs/Fasta_files/COGS_exochitinases.fasta", "w") as file:
        for item in seq_list:
            file.write(item)  # Add a newline after each item
    seq_list.clear()
    cogs_found_list = cogs_found_list + COGS_exochi

else:
    no_cogs_found_list = no_cogs_found_list + COGS_exochi


#COGS_LPMO
if len(COGS_LPMO_df.index) > 0: # seq_list is ready to be exported into fasta
    COGS_LPMO_name_list = COGS_LPMO_df["genome_prokka_feats"].to_list()
    COGS_LPMO_sequence_list = COGS_LPMO_df["aa_sequence"].to_list()
    seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(COGS_LPMO_name_list, COGS_LPMO_sequence_list)]

    # write the fasta file with each item from short_seq_list
    with open("Outputs/Fasta_files/COGS_LPMOs.fasta", "w") as file:
        for item in seq_list:
            file.write(item)  # Add a newline after each item
    seq_list.clear()
    cogs_found_list = cogs_found_list + COGS_LPMO

else:
    no_cogs_found_list = no_cogs_found_list + COGS_LPMO


#COGS_deacet
if len(COGS_deacet_df.index) > 0: # seq_list is ready to be exported into fasta
    COGS_deacet_name_list = COGS_deacet_df["genome_prokka_feats"].to_list()
    COGS_deacet_sequence_list = COGS_deacet_df["aa_sequence"].to_list()
    seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(COGS_deacet_name_list, COGS_deacet_sequence_list)]

    # write the fasta file with each item from short_seq_list
    with open("Outputs/Fasta_files/COGS_deacetylases.fasta", "w") as file:
        for item in seq_list:
            file.write(item)  # Add a newline after each item
    seq_list.clear()
    cogs_found_list = cogs_found_list + COGS_deacet

else:
    no_cogs_found_list = no_cogs_found_list + COGS_deacet

#Cogs found announcement
if len(cogs_found_list) > 0 and len(no_cogs_found_list) > 0:
    print("The following cogs had sequences in your dataset :", cogs_found_list)

#No cogs found announcement
if len(no_cogs_found_list) > 0:
    print("No genes were found for the the following cogs :", no_cogs_found_list)

#All cogs were found in your dataset
if len(cogs_found_list) > 0 and len(no_cogs_found_list) == 0:
    print("All the cogs existed in your dataset")
    
#Runtime calculation
time_finished = time.time()

print(os.path.basename(__file__), " took ", (time_finished - time_start), " seconds to run")

print("Thank you for using this script - João Fragoso de Almeida")