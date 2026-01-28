import os
import shutil
import time
import pandas as pd

#import modules
from ..config import annotation_columns_target, presence_counter_groups, path_main_melange
from ..categories.category_checks import cats_with_stuff, cats_with_types, cats_with_names

#This script counts the hits per category (presence_counter_groups) in all genomes, uses
# "All_Genomes_grouped_features.csv" as a base file

time_start = time.time()
print(f"starting {os.path.basename(__file__)}...")

#Paths
input = os.path.normpath(os.path.join(path_main_melange,"Outputs/All_Genomes_grouped_features.csv"))
output = os.path.normpath(os.path.join(path_main_melange,"Outputs/Statistics/"))
file_name = "presence_counter.csv"

#make folders
if os.path.exists(output): #cleans up folder if it has a bunch of trash in it, same for the rest
    # Check if the file exists
    files = [file for file in os.listdir(output) if "presence_counter.csv" in file]

    if files:  # If list not empty basically
        for item in files:
            os.remove(os.path.normpath(os.path.join(output, item)))  # removes old files
else:
    os.makedirs(output)

#For use in loop to make a list of the names of things
name_list = []
count_series = pd.Series(dtype="int")
series_list =[]


######################################################## Code ########################################################

#import the table, then combine the annotation columns
df = pd.read_csv(input, index_col= "genome_prokka_feats")
anno_serie = df[annotation_columns_target].astype(str).agg('+'.join, axis=1)

#make another df with only the index
top_index = df.index.str.replace(r"_PROKKA_\d*$", "", regex=True)

#Make a set containing genome names only
genome_set = set(top_index.to_list())

#Now make a loop and count for each "group" how many counts it has
for group in presence_counter_groups: #selects groups of categories

    #Controls atributing names
    name_set = False

    if "!" in group:
        divided = group.split("!")
        # Set the name for the group as the name of the first item in the group
        if not name_set:
            group_name = cats_with_names[("name_" + divided[0])]
    else:
        # Set the name for the fasta file
        if not name_set:
            group_name = cats_with_names[("name_" + group)]
            divided = [group]
        name_set = True

    #Now get te values in each group from the cats with stuff dict
    all_annos = []
    for item in divided:
        all_annos = all_annos + cats_with_stuff[item]

    #Make a list that ends up with all the names of all the different groups
    name_list.append(group_name)

    #make an empty series BUT WITH THE RIGHT INDEX!!! to store the counts; also name already done
    count_series = pd.Series(data= 0, index=anno_serie.index)

    #goes through every annotation group by group, and adds them together!
    for anno in all_annos:
        count_series = count_series.add(anno_serie.str.count(anno))

    #add the name now after all the operations
    count_series.name = group_name

    #then append it to the storage list
    series_list.append(count_series)

#Merge the counts series back into the df!
df = pd.concat(series_list, axis=1).sort_index(axis="index").fillna(value=0)

#Make multi_level index for the df
df.index = pd.MultiIndex.from_arrays(
    [top_index, df.index],
    names=["genome", "genome_prokka_feats"])

#Sum up the counts, since they are the only numerical columns, the others get dumped by pandas
df_toplevel = df.groupby(level=0).sum(numeric_only=True)

if "__" or "-" in df_toplevel.index:
    index_real = df_toplevel.index.to_list()
    index_new = []
    for item in index_real:
        item = item.replace("__","_")
        item = item.replace("-","_")
        index_new.append(item)

df_toplevel.index = index_new

#Export the table to output and filename
df_toplevel.to_csv(os.path.normpath(os.path.join(output,file_name)),index=True, sep=",")

time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")
