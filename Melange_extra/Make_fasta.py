import os
import shutil
import time
import pandas as pd
import textwrap

# Import variables from modules:
from config import name_change_condition, annotation_columns_target, make_fasta_groups, select_by_dna
from category_checks import cats_with_stuff, cats_with_types, cats_with_names


time_start = time.time()

#make folders
if os.path.exists("Outputs/Fasta_files"): #cleans up folder if it has a bunch of trash in it, same for the rest
    shutil.rmtree("Outputs/Fasta_files")
    os.makedirs("Outputs/Fasta_files")
else:
    os.makedirs("Outputs/Fasta_files")

#paths
path_to_outputs = "Outputs/Fasta_files/"



#import the table!
df = pd.read_csv("Outputs/All_Genomes_grouped_features.csv")

#aggregate annotation columns, in hindsight making the "which column is it in" problem very easy to solve... oh well
df["anno_sum"] = df[["COG", "KO", "MEROPS", "CAZymes", "PFAM"]].astype(str).agg('+'.join, axis=1)




#start selecting by
for category in make_fasta_groups: #selects groups of categories

    #MISSING IMPLEMENTATION OF THIS FOR ONE CATEGORY OR ALL CATEGORIES!!!

    name_list = []
    seq_list = []
    #make sure to set the name only once per cycle
    name_set = False

    if "!" in category:

        #now have a list with all cats that go in one fasta file
        divided = category.split("!")

        # Set the name for the fasta file as the name of the first category
        fasta_name = cats_with_names[("name_" + divided[0])]

        #In case several cats end up with the same name
        while fasta_name in os.listdir(path_to_outputs):
            fasta_name = fasta_name + "_1"

        #Now make the file path for the fasta
        fasta_name = fasta_name + ".fasta"
        path_to_write = os.path.join(path_to_outputs,fasta_name)
        # print(f"path_to_write: {path_to_write}")
        name_set = True

        #to hold all the annnotations
        all_annos = []

        #gets all annotations of each category into all_annos
        for cat in divided:
            all_annos = all_annos + cats_with_stuff[cat] #Running as intended

    else:
        # Set the name for the fasta file
        fasta_name = cats_with_names[("name_" + category)]

        #In case several cats end up with the same name
        while fasta_name in os.listdir(path_to_outputs):
            fasta_name = fasta_name + "_1"

        #Now make the file path for the fasta
        fasta_name = fasta_name + ".fasta"
        path_to_write = os.path.join(path_to_outputs, fasta_name)
        # print(f"path_to_write: {path_to_write}")
        name_set = True

        # to hold all the annnotations
        all_annos = []
        all_annos = all_annos + cats_with_stuff[cat]  # Running as intended


    #now build the bool mask to select which sequences go into the fasta file
    filter_bool = pd.Series(dtype=bool)

    for anno in all_annos: #iterates though each annotation
        finder = df["anno_sum"].str.contains(anno, case=False, na=False)  # True where match is found
        if sum(finder) != 0:  # to avoid operations when no information is added
            if not filter_bool.size == 0:
                filter_bool = filter_bool + finder  # add finder to the boolean mask of other loops
            # inside the same genome. Builds a large boolean mask aggregating all True hits over every cat_n
            else:
                filter_bool = finder  # if its not made yet, make it

        # This is when the boolean filter is used to return only rows for which the filter is TRUE
        if sum(filter_bool) == 0:
            # If filter_bool empty skip this iteration
            continue

        df_selected = df.loc[filter_bool].copy()
        df_selected = df_selected.reset_index()  # Integer index reset

        # In case a dataframe is empty after the mask is applied, this skips the rest of the iteration
        if df_selected.index.tolist() == []:
            continue

        # makes two lists, one with the names of each sequence, another with the dna or aa sequences, depening on select_by_dna
        name_list = df_selected["genome_prokka_feats"].tolist()
        if select_by_dna: #Set in config, dna selected if True, amino acids if False
            seq_list= df_selected["dna_sequence"].tolist()
        else:
            seq_list= df_selected["aa_sequence"].tolist()

        #Wrap seq_list so it is at most 60 characters per line
        wrapped_seqs = [textwrap.fill(seq, width=60) for seq in seq_list]

        #builds a list where each entry is one full fasta block to write to file:
        fasta_list = [">" + name + "\n" + seq + "\n" for name, seq in
                    zip(name_list, wrapped_seqs)]  # Already made in fasta format, just needs writting

        with open(path_to_write, "w") as file:
            for item in fasta_list:
                file.write(item)

time_finished = time.time()

print(f"{os.path.basename(__file__)}, took , {time_finished - time_start}, seconds to run")