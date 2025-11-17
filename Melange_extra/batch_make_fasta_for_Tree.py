import os
import time
import pandas as pd

# Import variables from config now:
from config import cat_1, cat_2, cat_3, cat_4, cat_5, cat_6, cat_7, cat_8, cat_9, cat_10
from config import cat_1_name, cat_2_name, cat_3_name, cat_4_name, cat_5_name, cat_6_name, cat_7_name, cat_8_name, \
    cat_9_name, cat_10_name

#For future reference, you can probably just wrap the search part in a function. Anyways a better way to do it would
# be to do the search once and then divide the results after I think. But really why bother, this is fast enough



time_start = time.time()




#Paths storage
path_Pfam_table = "Outputs/PFAM_grouped_table.csv"
path_to_outputs = "Outputs/Fasta_files/"

#open lists for announcements
no_cogs_found_list = [] #to store cogs with no hits found for an announcement at the end
cogs_found_list = [] #which cogs had positive results

#Check if output folder exists, if not make one
if not os.path.exists("Outputs/Fasta_files"):
    os.makedirs("Outputs/Fasta_files")

###Import PFAM table and make one dataframe for each COG category
pfam_table = pd.read_csv(path_Pfam_table)

#Make list of columns to KEEP for processing later (cogs, prokka names, amino acid sequence)
names_selected = list(pfam_table)
del names_selected[0] #eliminates the sequence name from the list
del names_selected[0] #eliminates COG column  from the list
del names_selected[1] #eliminates aa sequence column  from the list

#eliminate unwanted columns from pfam dataframe
pfam_table = pfam_table.drop(names_selected, axis=1)


#Category 1 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_1_df = pfam_table.loc[pfam_table["COG"].isin(cat_1),:].copy()

if not cat_1 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_1_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_1_list_of_name = cat_1_df["genome_prokka_feats"].to_list()
        cat_1_list_of_sequence = cat_1_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_1_list_of_name, cat_1_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_1_output = path_to_outputs + cat_1_name + ".fasta"
        with open(path_cat_1_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_1_cogs_check_list = cat_1_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_1:
            if element in cat_1_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_1


#Category 2 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_2_df = pfam_table.loc[pfam_table["COG"].isin(cat_2),:].copy()

if not cat_2 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_2_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_2_list_of_name = cat_2_df["genome_prokka_feats"].to_list()
        cat_2_list_of_sequence = cat_2_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_2_list_of_name, cat_2_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_2_output = path_to_outputs + cat_2_name + ".fasta"
        with open(path_cat_2_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_2_cogs_check_list = cat_2_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_2:
            if element in cat_2_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_2

#Category 3 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_3_df = pfam_table.loc[pfam_table["COG"].isin(cat_3),:].copy()

if not cat_3 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_3_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_3_list_of_name = cat_3_df["genome_prokka_feats"].to_list()
        cat_3_list_of_sequence = cat_3_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_3_list_of_name, cat_3_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_3_output = path_to_outputs + cat_3_name + ".fasta"
        with open(path_cat_3_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_3_cogs_check_list = cat_3_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_3:
            if element in cat_3_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_3

#Category 4 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_4_df = pfam_table.loc[pfam_table["COG"].isin(cat_4),:].copy()

if not cat_4 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_4_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_4_list_of_name = cat_4_df["genome_prokka_feats"].to_list()
        cat_4_list_of_sequence = cat_4_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_4_list_of_name, cat_4_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_4_output = path_to_outputs + cat_4_name + ".fasta"
        with open(path_cat_4_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_4_cogs_check_list = cat_4_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_4:
            if element in cat_4_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_4

#Category 5 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_5_df = pfam_table.loc[pfam_table["COG"].isin(cat_5),:].copy()

if not cat_5 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_5_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_5_list_of_name = cat_5_df["genome_prokka_feats"].to_list()
        cat_5_list_of_sequence = cat_5_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_5_list_of_name, cat_5_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_5_output = path_to_outputs + cat_5_name + ".fasta"
        with open(path_cat_5_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_5_cogs_check_list = cat_5_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_5:
            if element in cat_5_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_5

#Category 6 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_6_df = pfam_table.loc[pfam_table["COG"].isin(cat_6),:].copy()

if not cat_6 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_6_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_6_list_of_name = cat_6_df["genome_prokka_feats"].to_list()
        cat_6_list_of_sequence = cat_6_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_6_list_of_name, cat_6_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_6_output = path_to_outputs + cat_6_name + ".fasta"
        with open(path_cat_6_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_6_cogs_check_list = cat_6_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_6:
            if element in cat_6_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_6

#Category 7 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_7_df = pfam_table.loc[pfam_table["COG"].isin(cat_7),:].copy()

if not cat_7 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_7_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_7_list_of_name = cat_7_df["genome_prokka_feats"].to_list()
        cat_7_list_of_sequence = cat_7_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_7_list_of_name, cat_7_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_7_output = path_to_outputs + cat_7_name + ".fasta"
        with open(path_cat_7_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_7_cogs_check_list = cat_7_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_7:
            if element in cat_7_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_7

#Category 8 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_8_df = pfam_table.loc[pfam_table["COG"].isin(cat_8),:].copy()

if not cat_8 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_8_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_8_list_of_name = cat_8_df["genome_prokka_feats"].to_list()
        cat_8_list_of_sequence = cat_8_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_8_list_of_name, cat_8_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_8_output = path_to_outputs + cat_8_name + ".fasta"
        with open(path_cat_8_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_8_cogs_check_list = cat_8_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_8:
            if element in cat_8_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_8

#Category 9 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_9_df = pfam_table.loc[pfam_table["COG"].isin(cat_9),:].copy()

if not cat_9 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_9_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_9_list_of_name = cat_9_df["genome_prokka_feats"].to_list()
        cat_9_list_of_sequence = cat_9_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_9_list_of_name, cat_9_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_9_output = path_to_outputs + cat_9_name + ".fasta"
        with open(path_cat_9_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_9_cogs_check_list = cat_9_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_9:
            if element in cat_9_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_9

#Category 10 starts here
#Cut down dataframe by COG category and cut it down into name, sequence;
cat_10_df = pfam_table.loc[pfam_table["COG"].isin(cat_10),:].copy()

if not cat_10 == [""]: #Make sure this shit actually works, but cmon shouldnt be thaat hard no?
    if len(cat_10_df.index) > 0: # seq_list is ready to be exported into fasta
        cat_10_list_of_name = cat_10_df["genome_prokka_feats"].to_list()
        cat_10_list_of_sequence = cat_10_df["aa_sequence"].to_list()
        seq_list = [">" + name + "\n" + seq + "\n" for name, seq in zip(cat_10_list_of_name, cat_10_list_of_sequence)] #Already made
        # in fasta format!

        # write the fasta file with each item from short_seq_list
        path_cat_10_output = path_to_outputs + cat_10_name + ".fasta"
        with open(path_cat_10_output, "w") as file:
            for item in seq_list:
                file.write(item)

        #Find out which COGs are actually present in the dataset and sort them into the (no_)cogs_found variables
        cat_10_cogs_check_list = cat_10_df.COG.to_list() #Also can be a function, we dont care about this variable after this
        for element in cat_10:
            if element in cat_10_cogs_check_list:
                cogs_found_list.append(element)
            else:
                no_cogs_found_list.append(str(element))

        # empty out seq_list
        seq_list.clear()

    else:
        no_cogs_found_list = no_cogs_found_list + cat_10


#Cogs found announcement
if len(cogs_found_list) > 0 and len(no_cogs_found_list) > 0:
    print("The following cogs had sequences in your dataset :", cogs_found_list)

#No cogs found announcement
if len(no_cogs_found_list) > 0:
    print("No genes were found for the the following cogs :", no_cogs_found_list)

#All cogs were found in your dataset
if len(cogs_found_list) > 0 and len(no_cogs_found_list) == 0:
    print("All the cogs existed in your dataset")

#nice goodbye message
print("Thank you for using my scripts, hope it helped (｡◕‿◕｡) -- JFA")

#Runtime calculation
time_finished = time.time()

print(os.path.basename(__file__), " took ", (time_finished - time_start), " seconds to run")

print("Thank you for using this script - João Fragoso de Almeida")
