import os
import time
import pandas as pd
import numpy as np
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from statistics import SIMPER_func as smp

#my imports
from config import path_main_melange, groups_file, colors_in_groupfile, pcoa_graph_format, pcoa_dpi
from config import pcoa_distance_matrix_metric, names_in_graph


#Warning, this file assumes you only have one colour annotation column, does not take more than that


time_start = time.time()

print(f"starting {os.path.basename(__file__)}...")

#paths
path_to_groupfile = os.path.normpath(os.path.join(path_main_melange, "Melange_extra/statistics/PCOA_group.csv"))
counts_files_dir_path = os.path.normpath(os.path.join(path_main_melange, "Annotation_results"))
counts_files_identifier = "_counts.csv" #termination of the counts files
output = os.path.normpath(os.path.join(path_main_melange, "Outputs/Statistics"))

#options for the graph
PCOA_format = "svg"
dpi=300

#Hardset_variables
distance_matrix_metric = "braycurtis" #Metric for the distance matrix, from skbio.diversity.beta_diversity
names_in_graph = False #Whether to represent names of the genomes in the graph itself


##################################################### Code below ######################################################

#get the csv with the groups
df_groups = pd.read_csv(path_to_groupfile, header=0, names=["Index","Groups","Color"], index_col="Index")


#Make output if it doesn't exist
if os.path.exists(output):
    #Check if there are files with PCOA.svg termination
    pcoa_files = [file for file in os.listdir(output) if f"PCOA.{PCOA_format}" in file]

    if pcoa_files: #If list not empty basically
        for item in pcoa_files:
            os.remove(os.path.normpath(os.path.join(output,item))) #removes old files
else:
    os.makedirs(output) #if folder dont exist, make it

#iterate through all the counts files
for file in os.listdir(counts_files_dir_path):
    if counts_files_identifier in file:
        filepath = os.path.join(counts_files_dir_path, file) #build the path to the file
        df = pd.read_csv(filepath, index_col=0)  # working fine

        # Hellinger transformation
        df_hellinger = np.sqrt(
            df.div(
                other=df.sum(axis=0),  # From how data is prepared, features in rows and samples (genomes) in columns,
                # I need to do the summing by columns, vertically
                axis=1))  # Then here divide horizontally

        # Need to transpose the table to the genomes are considered as index, which the distmatrix considers the samples
        df_hellinger = df_hellinger.T

        temp_index = df_hellinger.index.to_series()
        temp_index = temp_index.str[16:]
        temp_index = temp_index.replace("__","_")
        df_hellinger.index = temp_index



        if not sorted(df_hellinger.index.to_list()) == sorted(df_groups.index.to_list()):
            print(f"The indexes for your groups file and counts file are not the same; Ensure they are and run again"
                  f"\n quitting...")
            quit()

        # calculate the bray curtis distance matrix
        dist_mat = beta_diversity(metric="braycurtis",
                      counts=df_hellinger,  # uses a numpy array, which is supported by skbio
                      ids=df_hellinger.index)

        # calculate the PCOA
        pcoa_result = pcoa(dist_mat, number_of_dimensions=2)  # Getting a warning that the positive eigenvalues arent
        # significantly larger than the negative ones. weird. people online are ignoring it, so will I

        # Make the coordinates from the PCOA
        coordinates = pcoa_result.samples

        # arrange a df with them
        df_pcoa = coordinates[['PC1', 'PC2']]  # results of PCOA for each PC
        df_pcoa["index"] = df_hellinger.index.to_numpy()  # index
        df_pcoa = df_pcoa.set_index("index")  # setting the index
        df_pcoa = df_pcoa.join(df_groups)  # Adds the groups to the df

        # Make a color map for them groups
        codes = df_pcoa["Groups"].astype(
            "category").cat.codes  # gets the groups into a category and then into sequential integers
        # aka into (1,2,3,4) etc...

        #Which color map to use, in relation to how many groups are in groups file ["Groups"]
        if len(set(codes)) < 9:
            cmap = plt.get_cmap("tab10")
        elif len(set(codes)) <= 20:
            cmap = plt.get_cmap("tab20")
        else:
            print("You have more than 20 groups in your PCOA, this script is built for max 20, please reduce the number!"
                  "\n quitting...")
            quit()
        colors = cmap(codes)  # makes a "tab10" color map for codes (above), now you can pass it to c= in ax.scatter

"""
Todo: 

Make the colors a list and pass it to scatter like this:

# Example: Mapping 3 points
colors = ['magento', '#4287f5', 'teal']
ax.scatter(df_pcoa["PC1"], df_pcoa["PC2"], c=colors)

make a map for the legend?

use pandas.DataFrame.drop_duplicates
to get unique values for the groups colors

# Define your palette
color_map = {
    'Control': 'skyblue',
    'Treatment_A': '#FF5733',
    'Treatment_B': 'forestgreen'
}

# Map the categories in your dataframe to the list of colors
colors = df_pcoa['Group'].map(color_map)

"""




        # Now a legend
        categories = df_pcoa["Groups"].astype("category").cat.categories
        legend_plot = [
            Line2D(
                [0], [0],
                marker="o",
                color="w",
                label=cat,
                markerfacecolor=cmap(i),
                markersize=8
            )
            for i, cat in enumerate(categories)
        ]

        # make the titles for the axes, combining the proportion of variance explained by each PC
        proportion_series = pcoa_result.proportion_explained * 100  # Now in %

        # labels for the axis
        label_PC1 = ("PC1 (" + str(round(proportion_series["PC1"], 2)) + " %)")
        label_PC2 = ("PC2 (" + str(round(proportion_series["PC2"], 2)) + " %)")

        # Make graph
        fig, ax = plt.subplots()
        ax.scatter(df_pcoa["PC1"], df_pcoa["PC2"], c=colors)
        plt.title(file.replace("_counts.csv",""))
        ax.legend(handles=legend_plot, title="Groups")

        # Annotates graph based on the names in df, check how iterrows actually works
        if names_in_graph:
            droped_Groups = df_pcoa[["PC1","PC2"]]
            for x, y in droped_Groups.iterrows():
                ax.annotate(x, y)

        ax.set_xlabel(label_PC1)
        ax.set_ylabel(label_PC2)

        hellinger_np = df_hellinger
        groups_list = df_groups["Groups"].to_list()

        # simp is a DF with a multi level index of ((group vs group), feature), e.g. ((abiotic - algae), CBM6), with
        # columns sp_mean (average dissimilaty)	sp_sd (standard deviation)	ratio	sp_pct	cumulative
        simp = smp.simper(hellinger_np, groups_list)

        #Calculate the Mean bray curtis dissimilarity BY PAIR, and map it back into simp
        summed_sp_means = simp.groupby(level=0).sum()["sp_mean"]
        simp["mean_BC"] = simp.index.get_level_values(0).map(summed_sp_means)

        #simp index level=1 is "index", which are features, level=0 are simper pairs
        # print("groups :",simp.groupby('index').groups)
        # print(simp.groupby("index").AA10)

        weighted_avg = {}
        #Calculate weighted (by the BC mean) mean of each feature across all groups
        for feature, group in simp.groupby(level=1): #groups by feature
            values = group["sp_mean"].to_numpy()
            weights = group["mean_BC"].to_numpy()
            weighted_avg[feature] = np.average(values, weights=weights)

        weighted_avg_df = (pd.DataFrame.from_dict(weighted_avg,
                                                  orient='index',
                                                  columns=['weighted_sp_mean']).sort_values(by='weighted_sp_mean',
                                                                                            ascending=False))

        #This selects the top 10 most dissimilar features
        top10_weighted = weighted_avg_df.iloc[0:10,:].copy()

        #Ensure identical sample order
        df_pcs = df_pcoa.loc[:,("PC1","PC2")].copy()
        hell_select_PCOA = df_hellinger.loc[df_pcs.index]

        # computes position of features in ordination, equivalent to wascores() in R::vegan
        feature_position_in_pcoa = (df_pcs.T @ hell_select_PCOA / hell_select_PCOA.sum(axis=0)).T

        #Select the top 10 features from top10_weighted in the feature_position_in_pcoa
        top10_positions = feature_position_in_pcoa.loc[top10_weighted.index,:]

        # try an actual arrowpatch
        for item in top10_positions.index:
            position = (top10_positions.loc[item,"PC1"],top10_positions.loc[item,"PC2"])
            arrow = patches.FancyArrowPatch(
                posA= (0.0,0.0),
                posB=position,
                arrowstyle="Simple,head_length=10,head_width=5,tail_width=2",
                color="black",
                label=item,
            )
            ax.add_patch(arrow)
            ts = 1.2 #to make the text go a bit further off from the arrow tip
            ax.text(
                x=position[0]*ts, y= position[1]*ts, s=item,
                ha='center',
                va='center',
                fontsize = 10,
                color = "black"

            )


        export_name = os.path.normpath(os.path.join(output, file.replace("counts.csv",f"PCOA.{PCOA_format}")))

        #export graphs to output
        plt.savefig(fname= export_name,dpi=dpi, format= PCOA_format)

time_finished = time.time()

print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")