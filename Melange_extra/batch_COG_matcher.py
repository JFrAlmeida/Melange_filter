import os
import itertools as itool
import pandas as pd
import re
from Bio import SeqIO
import time
from warnings import simplefilter
from config import all_cats

#This Script takes in the Melange annotation (https://sandragodinhosilva.github.io/melange/) of any number of genomes,
# selects only specific COGs, adds to it the dna and amino acid sequences corresponding to those COGs, and produces two tables,
# one with only metadata and the sequences, another with PFAM counts for statistics

############################################# HOW TO USE THE SCRIPT:

# Put it in the folder where your Melange results are, after moving them to a new folder of your choice (for example :
# a folder called "melange_results" in your desktop)

# That folder should have in it the "Annotation" and "Annotation_results" folders only, to try and keep it clean

# open this script in an IDE of your choice, and make sure you have the python packages in the very top of the script
# installed (if you dont the script will not run)

#############################################

time_start = time.time()


# default path to ORFS, if you want to change it go ahead:
# ORFs_path = os.path.normpath("Annotation_results/Orfs_per_genome_dummy/") #Remove the Dummy part at the end
ORFs_path = os.path.normpath("Annotation_results/Orfs_per_genome/")


#Get COGs from config file, and eliminate duplicates
all_COGs = list(set(all_cats))
for element in all_COGs:
    if not re.match(r"^COG\d+$", element):
        if element == "":
            print('One or more of your categories has the characters "" instead of a COG number, please eliminate'
                  ' them/correct it and run again!')
            quit()
        if element == "COG":
            print("One or more of your categories is missing the number after the COG letters, please correct it and"
                  " try again!")
            quit()
        if re.match(r"^\d+$", element):
            print("One of your categories is composed only of numbers, please correct it and try again!")
            quit()
        else:
            print("Something strange is wrong with your categories, please verify them and try again!")
            quit()

#location of Outputs and related tables, if it doesnt exist, make it
if not os.path.exists("Outputs"):
    os.makedirs("Outputs")

if not os.path.exists("Outputs/All_Features_per_Genome"):
    os.makedirs("Outputs/All_Features_per_Genome")

if not os.path.exists("Outputs/Genomes_with_PFAM_list"):
    os.makedirs("Outputs/Genomes_with_PFAM_list")

if not os.path.exists("Outputs/Transposed"):
    os.makedirs("Outputs/Transposed")

if not os.path.exists("Outputs/Selected"):
    os.makedirs("Outputs/Selected")



#Location of the fna and aa files
aa_fna_location = "Annotation/"

#Function to extract sequences from FASTA files
def extract_sequences(fasta_file, ids):
    seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return [seq_dict.get(seq_id) for seq_id in ids]  # Return sequence if ID exists

### Generate a description of the PFAM universe specific to each genome set

#Location of PFAM description file
pfam_desc_path = "Annotation_results/Pfam_description.csv"

#Make a dict with all entries of PFAMs
pfam_desc_full = pd.read_csv(pfam_desc_path)
pfam_index = dict(zip(pfam_desc_full.PFAM_ACC,pfam_desc_full.PFAM_desc))
pfam_list = [] #open an empty list to use later

#List of all genomes to process, to pass to loops
genome_path_list = os.listdir(ORFs_path)

for genome_name in genome_path_list:
    genome_path = os.path.join(ORFs_path, genome_name)
    genome = pd.read_csv(genome_path)
    genome["PFAM"] = genome["PFAM"].astype(str).str.replace(r"\.\d+", "", regex=True) #delete the ".number" that exists
    # after every PFAM accession number

    temp = list(genome.loc[:,"PFAM"].dropna().str.split("+")) #unfold the sequential PFAM numbers in the PFAM_ACC column
    pfam_temp = list(itool.chain.from_iterable(temp)) #make them into a list
    pfam_list = pfam_list + pfam_temp

#make the final list by removing duplicates between the different genomes
pfam_list = list(dict.fromkeys(pfam_list)) #dict.fromkeys will preserve order and eliminate repeats
pfam_list_df = pd.DataFrame(data= pfam_list) #and into a dataframe

#rename columns and remove NaN
pfam_list_df = pfam_list_df.rename(axis=1, mapper = {0:"PFAM_ACC"})
pfam_list_df = pfam_list_df.dropna()

#Pull hits of "PFAM_ACC" on the pfam index to a new column "description"
pfam_list_df["description"] = pfam_list_df["PFAM_ACC"].apply(lambda x: pfam_index.get(x))
pfam_list_df = pfam_list_df.sort_values(by= ["PFAM_ACC"], ascending= True)

#make a dictionary of Pfam number :  concatenation of the PFAM+descriptions, to substitue further down
pfam_list_intermediate = pfam_list_df
pfam_list_intermediate["pfam_desc"] = pfam_list_intermediate["PFAM_ACC"].astype(str) + ": " + pfam_list_intermediate["description"].astype(str)
pfam_list_desc = pfam_list_intermediate.loc[:,["PFAM_ACC", "pfam_desc"]]
pfam_list_dict = dict(zip(pfam_list_desc["PFAM_ACC"], pfam_list_desc["pfam_desc"]))

#Start processing genomes
#Make some essential DFs, that need to be outside the loops
grouped_features_table = pd.DataFrame()
merged_df = pd.DataFrame()

#Counter for genome progress management
counter = 0
counter_percent = 0


#Ignoring performance warning from Pandas, which does not appear do be relevant, but spit out anyways past some loops
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)



#Open genome and build the metainfo-table
for genome_name in genome_path_list:

    genome_path = os.path.join(ORFs_path, genome_name)

    #paths to fna and aa files are prepared
    aa_name = genome_name.replace("_all_features.csv", ".faa")
    aa_joined = os.path.join(aa_fna_location, aa_name)
    aa_path = os.path.normpath(aa_joined)

    #path to fna
    ffn_name = genome_name.replace("_all_features.csv", ".ffn")
    ffn_joined = os.path.join(aa_fna_location, ffn_name)
    ffn_path = os.path.normpath(ffn_joined)

    #Now ready to start building the table
    genome = pd.read_csv(genome_path)

    #Make a mask stating yes where a COG is in our list, and false when it isnt, use it to select lines of interest
    # create the empty bool series first
    filter_bool = pd.Series(dtype=bool)

    # This loop makes a filter_bool boolean series that adds a true for only those lines containing any of the COGs listed above
    for item in all_COGs:
        finder = genome["COG"].str.contains(item, case=False, regex=True, na=False)

        # This allows the attribution of the first entry, then addition of the next ones to the same thing
        if filter_bool.size == 0:
            filter_bool = finder
        else:
            filter_bool = (filter_bool + finder)

    #This is when the mask is used to select the COGs, behaves similarly to match in R (I think)
    genome_selected = genome.loc[filter_bool]
    genome_selected = genome_selected.rename(mapper={"row_0":"prokka_features"}, axis=1)
    genome_selected = genome_selected.reset_index()

    #In case a dataframe is empty after the mask is applied, this skips the rest of the iteration
    if genome_selected.index.tolist() == []:
        continue


    #Build a first column containing the genome name plus the prokka feature considered on that line
    #unfucking the genome name
    cond = r"GCA_\d+\.\d+_(.*?)_all_features\.csv"
    simpler_genome_name = re.search(cond, genome_name).group(1) #leaves the genome name plus strain ID

    #Make col with the unique name per genome per prokka feature
    genome_selected["genome_prokka_feats"] = simpler_genome_name + "_" + genome_selected["prokka_features"]

    #Reorder dataframe so ["genome_prokka_feats"]
    cols = genome_selected.columns.tolist()
    cols2 = cols[-1:] + cols[:-1]
    genome_selected = genome_selected[cols2]

    # add two new columns, one for nucleotide, one for amino acid sequences into genome dataframe; This also calls
    # the function that parses the fna and aa files. I think .loc[:, "prokka_features"] passes an iterator
    # to the function, processing one line of the dataframe and one feature (of those remaining) at a time
    genome_selected.loc[:, "dna_sequence"] = extract_sequences(ffn_path, genome_selected.loc[:, "prokka_features"])
    genome_selected.loc[:, "aa_sequence"] = extract_sequences(aa_path, genome_selected.loc[:, "prokka_features"])

    #eliminate an Index column that remains leftover
    genome_selected = genome_selected.drop(["index"], axis=1)

    out_genome_path = os.path.join("Outputs/All_Features_per_Genome", genome_name)
    genome_selected.to_csv(out_genome_path, index=False)


    #making a dataframe that aggregates all the general features of all genomes into one file
    if grouped_features_table.empty == True:
        grouped_features_table = genome_selected
    else:
        genome_selected_no_head = genome_selected.drop([0])
        grouped_features_table = pd.concat([grouped_features_table, genome_selected], ignore_index=True, )

    #PFAM table generation
    #At this point genome_selected was already outputed and joined to the grouped table, and is free to me modified
    # further by the PFAM table code below
    genome_selected_PFAM = genome_selected.copy()
    genome_selected_PFAM["PFAM"] = genome_selected_PFAM["PFAM"].astype(str).str.replace(r"\.\d+", "", regex=True) #delete the
    # .number that exists after every PFAM number
    temp = list(genome_selected_PFAM.loc[:,"PFAM"].dropna().str.split("+")) #unfold the sequential PFAM numbers in the PFAM_ACC column
    pfam_temp = list(itool.chain.from_iterable(temp)) #make them into a list
    pfam_temp_list = list(dict.fromkeys(pfam_temp))  # dict.fromkeys will preserve order and eliminate repeats

    for item in pfam_temp_list:
        genome_selected_PFAM[item] = ""

    #objective is counting how many PFAMS exist in each row
    temp_pfam_table = genome_selected_PFAM.loc[:, genome_selected_PFAM.filter(regex="PF\d+$", axis=1).columns]
    temp_pfam_table = temp_pfam_table.columns.values.tolist()


    #Count how many of each pfam there is in each COG
    for pfam in temp_pfam_table:
        pfam_series = genome_selected_PFAM["PFAM"] #now its a series and I can use .count() on it
        genome_selected_PFAM[pfam] = pfam_series.str.count(pfam)

    # Get rid of any NaN columns that for some reason keep appearing from the PFAM lists above #
    genome_selected_PFAM = genome_selected_PFAM.drop(labels= "nan", axis= 1, errors = "ignore")

    # Now just rename the column with pfam_list_dict, which has this structure: "PFAM number" : "pfam number : description"
    #genome_selected_PFAM = genome_selected_PFAM.rename(mapper=pfam_list_dict, axis=1)

    #For some reason I am getting a column "P" in some genomes, and I cant figure out why, so I am removing it here
    col_temp_pfam = genome_selected_PFAM.columns.tolist()
    pfam_temp_set = set(col_temp_pfam)
    if "P" in pfam_temp_set:
        genome_selected_PFAM = genome_selected_PFAM.drop("P", axis=1)

    #export this genome with its PFAM list to:
    genome_selected_PFAM.to_csv("Outputs/Genomes_with_PFAM_list/" + genome_name.replace("_all_features.csv",
                                                                                   "_PFAM_list.csv"),index=False)

    #Below code makes both the intermediate transposed genome tables and a concatenated table with all PFAMs over
    # all genomes

    #Start by dropping unnecessary columns
    genome_selected_PFAM_all = genome_selected_PFAM.copy()
    genome_selected_PFAM_all = genome_selected_PFAM.set_index("genome_prokka_feats")
    genome_selected_PFAM_all = genome_selected_PFAM_all.drop(["A", "CAZymes", "KO", "MEROPS", "prokka_features", "PFAM"], axis=1, errors="ignore")
    genome_selected_PFAM_all_T = genome_selected_PFAM_all.T #Transposing now
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.reset_index() #Fixing Index after transposing
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.rename(columns = {"index": "a_genome_prokka_feats"}) #Still
    # fixing Index, leaving the name with an a before so that the sorting of the table in the end is simpler

    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.set_index("a_genome_prokka_feats") #Finally index is "a_a_genome_prokka_feats"

    #changing row names of dna and aa sequence for the sorting at the end of the script places columns as intended
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.rename(mapper={"dna_sequence": "Dna_sequence",
                                                                            "aa_sequence": "Eaa_sequence"}, axis= 0)

    #Start merging them all iteratively into the same table
    if merged_df.empty:
        merged_df = genome_selected_PFAM_all_T
    else:
        merged_df = merged_df.merge(genome_selected_PFAM_all_T, on="a_genome_prokka_feats", how="outer")

    genome_selected_PFAM_all_T.to_csv("Outputs/Transposed/" + genome_name.replace("_all_features.csv",
                                                                        "_PFAM_list_transposed.csv"), index=True)
    counter = counter + 1
    counter_percent = (counter/len(genome_path_list))*100
    counter_percent = round(counter_percent, 2)
    counter_announcement = ("Finished genome " + str(counter) + " of " + str(len(genome_path_list)) +
                            " (" + str(counter_percent) + "%)")
    print(counter_announcement)

#sort the final PFAM table, and handle index
merged_df = merged_df.sort_values(by="a_genome_prokka_feats", ascending=True)
merged_df_untransposed = merged_df.copy().T.reset_index()

#make the now column names make sense again
merged_df_untransposed = merged_df_untransposed.rename(columns={"Dna_sequence": "dna_sequence",
                                                                "Eaa_sequence": "aa_sequence"})
#and set an index
merged_df_untransposed = merged_df_untransposed.set_index("genome_prokka_feats").fillna(value= 0)

#Now add summary statistics to it
merged_df_untransposed.loc["Total sum of each PFAM"] = merged_df_untransposed.sum(
    numeric_only= True, axis= 0) #col totals
merged_df_untransposed.loc[: , "Total PFAMs per COG"] = merged_df_untransposed.iloc[:, 3:].sum(
    numeric_only= True, axis= 1) #row totals


#Export general tables
grouped_features_table.to_csv("Outputs/All_Genomes_grouped_features.csv", index=False)
merged_df_untransposed.to_csv("Outputs/PFAM_grouped_table.csv", index=True)



#Runtime calculation
time_finished = time.time()

print(os.path.basename(__file__), " took ", (time_finished - time_start), " seconds to run")