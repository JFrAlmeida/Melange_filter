import os
########################################### Annotation finder options ###########################################

#This specifies the target column of MELANGE annotation to select for, in s list format -> ["item", "item2", "item3"]
# Valid options are COG", "KO", "MEROPS", "PFAM" or "CAZy", in any order, and any number and combination of
# them, as long as inside "" and separated by a comma; only COGS has been tested, but all exept PFAM should work
annotation_columns_target = ["CAZymes"]

#This controls the outcome of the "genome_prokka_feats" column. It concatenates what the regex expression matches, with
# the PROKKA feature nr of each annotation match
#By default the space between the GCA/GCF number of a genome and the _allfeatures.csv termination Melange attributes it
#valid options are "no_change" or your own condition within "". Defaults to r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv"
name_change_condition = r"GC[AF]_\d+\.\d+_(.*?)_all_features\.csv"

#Below is a condition that matches the entire genome name before _all_features\.csv. Use it in case your genomes don't
# have an easily recognizable pattern
# name_change_condition = r"(.*)_all_features\.csv"

########################################### dbcan parser options ###########################################

#If True, runs Process_dbcan_results and parse_dbcan_to_melange_ORF, adding the results of running dbcan (previously
# by the user) to the melange_Orfs_per_genome_table
run_dbcan_parser = True

#set condition for the name changing of the dbcan results run_dbcan_parser.py and process_dbcan_results.py use. Will
# change the name of the genomes in the CAZYmes_counts table.
#Valid options are "no change" or your own condition within "". Defaults to "^(?:[^_]*_){2}(.*)$"
dbcan_name_change = "^(?:[^_]*_){2}(.*)$"

########################################### fasta maker options ###########################################
#If this is True, at the end of batch_COG_matcher.py, Make_fasta.py will be run, if False it won't
run_make_fasta = False

#list with the category groups to make fasta files for. For each element ["element","thing_in_element!another_thing"], a
# fasta file will be output to Outputs/Fasta_files; Example at the bottom of fasta maker options
make_fasta_groups = []
#Note: the name_cat_n of the first item in each run_fasta element will be used to name the fasta file (e.g in
# "cat_1!cat_6", the name will be the one defined in name_cat_1)

#this lets you select if dna or amino acids sequences are used to build the fasta file
# Valid options are False -> amino acids; True -> dna
select_by_dna = False

#example of grouping many cats ["cat_1!cat_2!cat_3!cat_4","cat_5!cat_6"]
#make a list, each element separated by "," and inside each element of the list put the categories you want to be in
# the same fasta file, with an "!" as separator


########################################### Statistics options ###########################################

statistics_prefix = "Algae" # To be added as prefix to several statistics outputs (stacked bar, presence counter)
statistics_multiple_runs = True # To save files used for statistics (presence counter) of different runs into another directory
save_directory = "Statistics_general"

            ####################### presence counter options #######################

# whether to run presence_counter:
run_presence_counter = True #True for yes run it, False for don't run it
save_to_perm = True # if true saves the result to Outputs/Stats_perm, which is not cleaned each melange_extra run

#Category groups for presence counter, same syntax and valid options as for as make_fasta_groups above,
# will often be the exactly the same as make_fasta_groups
presence_counter_groups = ["cat_1","cat_2","cat_3","cat_4","cat_5","cat_6","cat_7","cat_8","cat_9","cat_10"]

            ####################### stacked bar plot options #######################

run_stacked_bar_plot = True #If true makes a stacked bar plot from the results of presence_finder

#!!!!Put this file inside /home/your/path/Melange_extra/statistics!!!!
stacked_order_color_file_name = "tree_order.csv" #file with two columns "Tree_order", with the intended order of the samples in the figure, and "Color", with the intended color of the LABELS in the figure;
# the figure, and "Color", with the intended color of the LABELS in the figure;
stacked_graph_format = "svg" #valid terminations are "png", "pdf", "svg" etc... check matplotlib plt.savefig for more
stacked_dpi = 300 #Of the figure
stacked_width = 0.4 #Of the bars
stacked_graph_title = "Algae polyssacharide targeting features" # Title of the graph
stacked_graph_name_prefix = "Algae" # would show as Algae_stacked_bar.svg, just changes the prefix

#stacked bar plot label options
stacked_fontsize = 10 #size of the label text
stacked_fontstyle = "italic" #can be bold, italic, normal
stacked_ha = "right" #set_horizontalalignment from matplotlib
stacked_rotation = 45 #rotation of the labels

            ####################### PCOA maker options #######################


"Not implemented yet, but the options will be more or less this"

run_pcoa_maker = True #Creates PCOAs from the counts.csv files in Annotation_results 
groups_file = "PCOA_group.csv"
colors_in_groupfile = True #If your groups file has a column with the colors for each group in RGB format 

pcoa_graph_format = "svg" #valid terminations are "png", "pdf", "svg" etc... check matplotlib plt.savefig for more
pcoa_dpi = 300 #Of the figure
pcoa_distance_matrix_metric = "braycurtis"
names_in_graph = False

########################################### PFAM description options ###########################################

pfam_topn_list = 10 #Make a small table with pfams that have at least n hits; set it to the minimum number you want, or leave empty if you don't want this table

########################################### Annotaion features to search for ###########################################
#IMPORTANT! -> each cat can only take one type of annotation, only COGs, or only Cazy, or only MEROPS etc.
# But different cats can be of different types (cat_1 CAZYmes and cat_2 MEROPS)
    #valid choices for cat_n is any COG, KO, MEROPS, PFAM or CAZY number
    #Valid choices for name_cat_n is anything, this determines the name of the fasta file with amino acids/DNA, and of
        # the presence_counter files
    #Valid choices for type_n is "COG", "KO", "MEROPS", "PFAM" or "CAZymes"

last_filter = [] # To select, from your categories below only those that are present here too leave blank -> [] if you dont want another filter
last_filter_columns = ["PFAM"] #same syntax as annotation_columns_target above

#Note that you can delete or add more categories at your will, as long as they keep this format
cat_1 = ["GH117", "GH50"]
name_cat_1 = "Agar"
type_1 = "CAZymes" #case sensitive! double check!

cat_2 = ["GH150","GH167","GH82"]
name_cat_2 = "Carrageenan"
type_2 = "CAZymes" #case sensitive! double check!

cat_3 = ["PL14", "PL15", "PL17", "PL18", "PL31", "PL34", "PL36", "PL38", "PL39", "PL41", "PL44", "PL6","PL7","PL5"]
name_cat_3 = "Alginate"
type_3 = "CAZymes" #case sensitive! double check!

cat_4 = ["GH3","GH74","PL3","PL1"]
name_cat_4 = "Porphyran"
type_4 = "CAZymes" #case sensitive! double check!

cat_5 = ["AA15","AA16","GT2","GH9"]
name_cat_5 = "Cellulose"
type_5 = "CAZymes" #case sensitive! double check!

cat_6 = ["GH55","GH17"]
name_cat_6 = "Laminarin"
type_6 = "CAZymes" #case sensitive! double check!

cat_7 = ["AA14","CE12","CE15","CE2","CE20","CE3","CE5","CE6","CE7","GH11","GH115","GH120","GH141","GH5","GH54","GH62","GH67","GT43","GT47","GT61","CE1","GH39","GT8"]
name_cat_7 = "Xylan"
type_7 = "CAZymes" #case sensitive! double check!

cat_8 = ["GH174","GH107","GH187","PL43"]
name_cat_8 = "Fucoidan"
type_8 = "CAZymes" #case sensitive! double check!

cat_9 = ["GH105","PL25","PL28","PL37","PL40","GH78","GH88","PL24"]
name_cat_9 = "Ulvan"
type_9 = "CAZymes" #case sensitive! double check!

cat_10 = ["GH119","GH126","GH15","GH176","GH31","GH49","GH57","GH97","AA13","GH13","GH14","GH144","GT35","GT5"]
name_cat_10 = "Starch"
type_10 = "CAZymes" #case sensitive! double check!



####################################### End of editable part of the script #######################################

#Please don't touch or delete this! it is necessary for other scripts to see the options in this one
config_locals = locals().items()

#common paths
path_main_melange = os.path.dirname(os.getcwd())
