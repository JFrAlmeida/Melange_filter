import os
import itertools as itool
import pandas as pd
import re
from Bio import SeqIO
import time
import warnings
from warnings import simplefilter
import shutil
import logging

#My module imports
from config import name_change_condition, annotation_columns_target
from category_checks import cats_with_stuff, cats_with_types
from config import run_dbcan_parser




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

#################### Put functions here ####################

def increment_counter ():
    global counter
    global counter_percent
    counter = counter + 1
    counter_percent = (counter / len(genome_path_list)) * 100
    counter_percent = round(counter_percent, 2)

#Function to extract sequences from FASTA files
def extract_sequences(fasta_file, seqs):
    seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return [seq_dict.get(seq_id) for seq_id in seqs]  # Return sequence if ID exists

############################################################


################### logging basic config ###################

if os.path.exists("Logs/"):
    shutil.rmtree("Logs")
    os.makedirs("Logs")
    open("Logs/empty_genomes.log", "w")
else:
    os.makedirs("Logs")
    open("Logs/empty_genomes.log", "w")

logging.basicConfig(
     filename="Logs/empty_genomes.log",
     filemode="a",
     format="{asctime} - {message}",
     style="{",
     datefmt="%Y-%m-%d %H:%M",
    level=logging.DEBUG
)

logger = logging.getLogger(__name__)
####################  hardset variables ####################

#Counter for genome progress management
counter = 0
counter_percent = 0
empty_genome_counter = 0
empty_genome_exists = False
no_change = False

############################################################

#Paths if rundbcan_parser has been run, vs not
if run_dbcan_parser:
    ORFs_path = os.path.normpath("Annotation_results/Orfs_per_genome_CAZy/")
else:
    ORFs_path = os.path.normpath("Annotation_results/Orfs_per_genome/")

#location of Outputs and related tables, if it doesnt exist, make it
if not os.path.exists("Outputs"):
    os.makedirs("Outputs")

if os.path.exists("Outputs/All_Features_per_Genome"): #cleans up folder if it has a bunch of trash in it, same for the rest
    shutil.rmtree("Outputs/All_Features_per_Genome")
    os.makedirs("Outputs/All_Features_per_Genome")
else:
    os.makedirs("Outputs/All_Features_per_Genome")

if os.path.exists("Outputs/Genomes_with_PFAM_list"):
    shutil.rmtree("Outputs/Genomes_with_PFAM_list")
    os.makedirs("Outputs/Genomes_with_PFAM_list")
else:
    os.makedirs("Outputs/Genomes_with_PFAM_list")

if os.path.exists("Outputs/Transposed"):
    shutil.rmtree("Outputs/Transposed")
    os.makedirs("Outputs/Transposed")
else:
    os.makedirs("Outputs/Transposed")

if os.path.exists("Outputs/Fasta_files"):
    shutil.rmtree("Outputs/Fasta_files")
    os.makedirs("Outputs/Fasta_files")
else:
    os.makedirs("Outputs/Fasta_files")




#Location of the fna and aa files
aa_fna_location = "Annotation/"

#condition for arranging the name of the genome into the first column

#Controls the condition to name each genome_name + PROKKA_feat concatenation, in the column "genome_prokka_feats"
try:
    if name_change_condition == "no_change":
        no_change = True # Name will not be changed
    elif name_change_condition != r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv":
        cond = name_change_condition
        print(f"your name_change_condition is : {name_change_condition}")
    else:
        cond = r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv"  # condition for name change will be this
        print(f"name_change_condition was not set as any valid option in config.py, running as default...")

except:
    print("something went wrong processing name_change_condition variable, set in config, please ensure it is set \n"
          "quitting...")
    quit()



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




    # This loop makes a filter_bool pandas series of type boolean with true for values matching any of the annotations
    # specified in any category in config.
    for item in annotation_columns_target: #item is the target column ("COG" "KO" etc.)
        for type in list(cats_with_types.keys()): #type here is the key to the type of annotation in config ("COG" "KO"
            # etc.), its value needs to match item to proceed
            print(f"catswith types: {cats_with_types[type]}, and item: {item}")
            if cats_with_types[type] == item:
                for category in list(cats_with_stuff.keys()): #category is cat_1, cat_2 etc.
                    for anno in cats_with_stuff[category]: #anno is the annotation targets, "COG132" "COG312" etc.
                        # print(f"looking for {anno} in {item}")
                        finder = genome[item].str.contains(anno, case=False, na=False) # True where match is found
                        if sum(finder) != 0: # to avoid operations when no information is added
                            if not filter_bool.size == 0:
                                filter_bool = filter_bool + finder # add finder to the boolean mask of other loops
# inside the same genome. Builds a large boolean mask aggregating all True hits over every cat_
                            else:
                                filter_bool = finder # if its not made yet, make it


    #This is when the boolean filter is used to return only rows for which the filter is TRUE
    if sum(filter_bool) == 0:
        #If filter_bool empty make a log into list of empty genomes
        empty_genome_counter = empty_genome_counter + 1
        logger.debug(genome_name)
        continue

    genome_selected = genome.loc[filter_bool]
    genome_selected = genome_selected.rename(mapper={"row_0":"prokka_features"}, axis=1) #QoL change for the col name
    genome_selected = genome_selected.reset_index() #Integer index reset

    #In case a dataframe is empty after the mask is applied, this skips the rest of the iteration
    if genome_selected.index.tolist() == []:
        increment_counter()
        counter_announcement = ("Finished genome " + str(counter) + " of " + str(len(genome_path_list)) +
                                " (" + str(counter_percent) + "%) -> " + genome_name + " did not have any of the required cogs")
        print(counter_announcement)
        continue

    #Build a first column containing the genome name plus the prokka feature considered on that line
    # Also unfucking the genome name a bit; cond is a regex pattern settable in config
    if no_change == False:
        simpler_genome_name = re.search(cond, genome_name).group(1) #leaves the genome name plus strain ID
        genome_selected["genome_prokka_feats"] = simpler_genome_name + "_" + genome_selected["prokka_features"]
    else:
        genome_selected["genome_prokka_feats"] = genome_name + "_" + genome_selected["prokka_features"]

    #Reorder dataframe so ["genome_prokka_feats"] is the first column
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
    genome_selected.to_csv(os.path.normpath(out_genome_path), index=False)


    #making a dataframe that aggregates all the general features of all genomes into one file
    if grouped_features_table.empty == False:
        genome_selected_no_head = genome_selected.drop([0])
        grouped_features_table = pd.concat([grouped_features_table, genome_selected], ignore_index=True, )
    else:
        grouped_features_table = genome_selected

    #PFAM table generation
    #The .copy() means I finally learned how python actually uses pointers :)
    genome_selected_PFAM = genome_selected.copy()
    genome_selected_PFAM["PFAM"] = genome_selected_PFAM["PFAM"].astype(str).str.replace(r"\.\d+", "", regex=True) #delete the
    # .number that exists after every PFAM number
    temp = list(genome_selected_PFAM.loc[:,"PFAM"].dropna().str.split("+")) #unfold the sequential PFAM numbers in the PFAM_ACC column
    pfam_temp = list(itool.chain.from_iterable(temp)) #make them into a list
    pfam_temp_list = list(dict.fromkeys(pfam_temp))  # dict.fromkeys will preserve order and eliminate duplicates

    for item in pfam_temp_list: # Opens an empty column for every PFAM entry
        genome_selected_PFAM[item] = ""

    #Creates a list of columns names that start with PF followed by any numbers, so every PFAM found in the genome
    temp_pfam_table = genome_selected_PFAM.loc[:, genome_selected_PFAM.filter(regex="PF\d+$", axis=1).columns]
    temp_pfam_table = temp_pfam_table.columns.values.tolist()


    #Count how many of each pfam there is in each COG
    for pfam in temp_pfam_table:
        pfam_series = genome_selected_PFAM["PFAM"] #now its a series and I can use .count() on it
        genome_selected_PFAM[pfam] = pfam_series.str.count(pfam)

    # Get rid of any NaN columns that for some reason keep appearing from the PFAM lists above #
    genome_selected_PFAM = genome_selected_PFAM.drop(labels= "nan", axis= 1, errors = "ignore")

    #export this genome with its PFAM list to:
    genome_selected_PFAM.to_csv("Outputs/Genomes_with_PFAM_list/" +
                                genome_name.replace("_all_features.csv",
                                                    "_PFAM_list.csv"),index=False)

    #Below code makes both the intermediate transposed genome tables and a concatenated table with all PFAMs over
    # all genomes

    #Start by dropping unnecessary columns
    genome_selected_PFAM_all = genome_selected_PFAM.copy()
    genome_selected_PFAM_all = genome_selected_PFAM.set_index("genome_prokka_feats")
    genome_selected_PFAM_all = genome_selected_PFAM_all.drop(["A", "prokka_features", "PFAM"], axis=1, errors="ignore")
    genome_selected_PFAM_all_T = genome_selected_PFAM_all.T #Transposing now
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.reset_index() #Fixing Index after transposing
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.rename(columns = {"index": "a_genome_prokka_feats"}) #Still
    # fixing Index, leaving the name with an a before so that the sorting of the table in the end is simpler

    #Finally index is "a_a_genome_prokka_feats"
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.set_index("a_genome_prokka_feats")

    #changing row names of dna and aa sequence for the sorting at the end of the script places columns as intended
    genome_selected_PFAM_all_T = genome_selected_PFAM_all_T.rename(mapper={"dna_sequence": "Dna_sequence",
                                                                            "aa_sequence": "Eaa_sequence"}, axis= 0)

    #Start merging them all iteratively into the same table
    if merged_df.empty:
        merged_df = genome_selected_PFAM_all_T
    else:
        merged_df = merged_df.merge(genome_selected_PFAM_all_T, on="a_genome_prokka_feats", how="outer")

    #writing PFAM table to csv
    genome_selected_PFAM_all_T.to_csv("Outputs/Transposed/" + genome_name.replace("_all_features.csv",
                                                                        "_PFAM_list_transposed.csv"), index=True)

    #Icrement counter and announce genome finished!
    increment_counter()
    counter_announcement = (f"Finished genome {counter} of {len(genome_path_list)} ({counter_percent}%)")
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

#Emptu genomes warning
if empty_genome_counter != 0:
    warnings.warn(f"Script ran successfully! But {empty_genome_counter} genomes had none of the annotations you requested..."
      f"\n a list of their names is in ./Logs/empty_genomes.log")


#nice goodbye message
print("Thank you for using my scripts, hope it helped (｡◕‿◕｡) -- JFA")

#Runtime calculation
time_finished = time.time()

print(f"{os.path.basename(__file__)}, took , {time_finished - time_start}, seconds to run")