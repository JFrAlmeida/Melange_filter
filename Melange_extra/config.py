#Objective is to hold variables to be share bewteen the melange processing scripts here
#no need to add file extensions (e.g. .fasta) to the names

#edit each category

########################## Edit this part of the script!!! ################################################

#If True, runs Process_dbcan_results and parse_dbcan_to_melange_ORF, adding the results of running dbcan (previously
# by the user) to the melange_Orfs_per_genome_table
run_dbcan_parser = True

#This specifies the target column of MELANGE annotation to select for, in s list format -> ["item", "item2", "item3"]
# Valid options are COG", "KO", "MEROPS", "PFAM" or "CAZy", in any order, and any number and combination of
# them, as long as inside "" and separated by a comma; only COGS has been tested, but all exept PFAM should work
annotation_columns_target = ["COG","CAZymes"]

#This controls the outcome of the "genome_prokka_feats" column. It concatenates what the regex expression matches, with
# the PROKKA feature nr of each annotation match
#By default the space between the GCA/GCF number of a genome and the _allfeatures.csv termination Melange attributes it
#valid options are "no_change" or your own condition within "". Defaults to r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv"
name_change_condition = r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv"

#Below is a condition that matches the entire genome name before _all_features\.csv. Use it in case your genomes don't have a pattern
# name_change_condition = r"(.*)_all_features\.csv"


########################################### fasta maker options ###########################################
#If this is True, at the end of batch_COG_matcher.py, Make_fasta.py will be run, if False it won't
run_make_fasta = True

#Whether to run Make_fasta.py -> makes a fasta file for each category, groups of categories, or all together.
#Valid choices are "cat" (one for each category) "all" (all together) or a list with the cat groupings, example below
make_fasta_groups = ["cat_1!cat_6","cat_2!cat_7","cat_3!cat_8","cat_4","cat_5"]
#Note: the name_cat_n of the first item in each run_fasta element will be used to name the fasta file (e.g in
# "cat_1!cat_6", the name will be the one defined in name_cat_1)

#this lets you select if dna or amino acids sequences are used to build the fasta file
# Valid options are False -> amino acids; True -> dna
select_by_dna = False

#example of grouping many cats ["cat_1!cat_2!cat_3!cat_4","cat_5!cat_6"]
#make a list, each element separated by "," and inside each element of the list put the cats you want to be in the same
# fasta file, with an "!" as separator

########################################### Annotaion features to search for ###########################################
#IMPORTANT! -> each cat can only take one type of annotation, only COGs, or only Cazy, or only MEROPS etc.
# But different cats can be of different types (cat_1 CAZYmes and cat_2 MEROPS)
    #valid choices for cat_n is any COG, KO, MEROPS, PFAM or CAZY number
    #Valid choices for name_cat_n is anything, this determines the name of the fasta file with amino acids/DNA
    #Valid choices for type_n is "COG", "KO", "MEROPS", "PFAM" or "CAZymes"

#Note that you can delete or add more categories at your will, as long as they keep this format
cat_1 = ["COG3325", "COG3469"]
name_cat_1 = "GH18_endochitinases"
type_1 = "COG" #case sensitive! double check!

cat_2 = ["COG3179"]
name_cat_2 = "GH19_endochitinases"
type_2 = "COG" #case sensitive! double check!

cat_3 = ["COG4724", "COG3525", "COG3979"]
name_cat_3 = "Exochitiniases"
type_3 = "COG" #case sensitive! double check!

cat_4 = ["COG3397"]
name_cat_4 = "LPMOs"
type_4 = "COG" #case sensitive! double check!

cat_5 = ["COG3394", "COG0726", "COG2861"]
name_cat_5 = "Deacetylases"
type_5 = "COG" #case sensitive! double check!

cat_6 = ["GH18"]
name_cat_6 = "GH18_endo"
type_6 = "CAZymes" #case sensitive! double check!

cat_7 = ["GH19"]
name_cat_7 = "GH19_endo"
type_7 = "CAZymes" #case sensitive! double check!

cat_8 = ["GH20"]
name_cat_8 = "CAZYmes"
type_8 = "GH20 exochitinases" #case sensitive! double check!



####################################### End of editable part of the script #######################################

#Please don't touch or delete this! it is necessary for other scripts to see the options in this one
config_locals = locals().items()


########################## Example of a config file ################################################


#Example of usage for chitinase COGs

#cat_1 = ["COG3325", "COG3469"] # GH18 endochitinases
#name_cat_1 = "GH18_endochitinases" #will be the name of the file outputed
#
#cat_2 = ["COG3179"] # GH19 endochitinases
#name_cat_2 = "GH19_endochitinases" #will be the name of the file outputed
#
#cat_3 = ["COG4724", "COG3525", "COG3979"] # exochitiniases
#name_cat_3 = "Exochitiniases" #will be the name of the file outputed
#
#cat_4 = ["COG3397"] # LPMOs
#name_cat_4 = "LPMOs"
#
#cat_5 = ["COG3397"] #
#name_cat_5 = "LPMOs"
#
#cat_6 = ["COG3394", "COG0726", "COG2861"] # Deacitylases
#name_cat_6 = "Deacitylases" #will be the name of the file outputed
#
#cat_7 = [] # empty
#name_cat_7 = "empty" #will be the name of the file outputed
#
#cat_8 = [] # empty
#name_cat_8 = "empty" #will be the name of the file outputed
#
#cat_9 = [] # empty
#name_cat_9 = "empty" #will be the name of the file outputed
#
#cat_10 = [] # empty
#name_cat_10 = "empty" #will be the name of the file outputed






























