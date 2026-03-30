import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from inspect import currentframe, getframeinfo

#my imports
from config import path_main_melange, stacked_order_color_file_name, stacked_graph_format, stacked_dpi, stacked_width
from config import stacked_graph_title, statistics_prefix, stacked_fontsize, stacked_fontstyle, stacked_ha
from config import stacked_rotation

time_start = time.time()

print(f"starting {os.path.basename(__file__)}...")

#paths
counts_files_path = os.path.normpath(os.path.join(path_main_melange, "Outputs/Statistics/presence_counter.csv"))
index_tree_order = os.path.normpath(os.path.join(os.path.join(path_main_melange,"Melange_extra/MF_statistics"), stacked_order_color_file_name))
output = os.path.normpath(os.path.join(path_main_melange, "Outputs/Statistics"))



#Make output if it doesn't exist
if os.path.exists(output):
    #Check if there are files with PCOA.svg termination
    bar_files = [file for file in os.listdir(output) if f"stacked_bar.{stacked_graph_format}" in file]
    if bar_files:
        for item in bar_files:
            os.remove(os.path.normpath(os.path.join(output,item))) #removes old files
else:
    os.makedirs(output) #if folder dont exist, make it

#df with presence_counter.csv
df = pd.read_csv(counts_files_path, index_col=0)

#Reindexes the barplot to match an order chosen by the user in tree_order.csv
try:
    df_index = pd.read_csv(index_tree_order, index_col="Tree_order")
    print(f"using the index in {index_tree_order}")
    df_reindexed = df.reindex(df_index.index)
except:
    print("index_tree_order not set or column name not 'Tree_order', proceeding unordered")
    df_reindexed = df_index

height_dict = {}
for column in df_reindexed.columns:
    clist = df_reindexed.loc[:,column]
    cdict = {
        column:clist
    }
    height_dict.update(cdict)

#making color maps
cat_len = len(height_dict)
if cat_len <= 10:
    cmap = plt.cm.tab10
elif cat_len > 10 and cat_len <= 10:
    cmap = plt.cm.tab20
else:
    frameinfo = getframeinfo(currentframe())
    print("Look I don't know what you're trying to do, but really think if you WANT to represent "
          "more than 20 categories in one stacked bar plot... \n"
          "... And then kindly figure out how to work with palettes with more than 20 colors, cause I really"
          " can't be bothered. Cherio!")
    print(f"script path: {frameinfo.filename}, code line: {frameinfo.lineno}")
    print("quitting....")
    quit()

#initialise the graph
fig, ax = plt.subplots(figsize=(30,6), constrained_layout=True)

#Make each bar
colors = cmap(np.linspace(0, 1, cat_len))
bottom = np.zeros(len(df_reindexed.index))
for (category, count), color in zip(height_dict.items(), colors):
    ax.bar(df_reindexed.index,
           count,
           stacked_width,
           bottom = bottom,
           label = category,
           color = color,
           edgecolor = "0.3",
           linewidth = 0.4)
    bottom += count

#Set the tick labels
labels_x = np.arange(len(df_reindexed.index))  # the label locations
ax.set_xticks(ticks = labels_x,
              labels = df_reindexed.index,
              minor = False,
              fontsize = stacked_fontsize,
              fontstyle = stacked_fontstyle,
              ha = stacked_ha,
              rotation = stacked_rotation)

#set tick colors based on "Color" columns of tree_order.csv
df_index["Color"] = df_index["Color"].fillna(value="#000000") #fill nan with black color hexadecimal
if "Color" in df_index.columns:
    for tick_label, color in zip(ax.get_xticklabels(), df_index.Color):
        tick_label.set_color(color)
else:
    print("Color column not found in tree_order.csv, proceeding with uncoloured labels")

ax.set_title(stacked_graph_title)
ax.legend(loc="upper right")
fig.savefig(fname=os.path.normpath(os.path.join(output, f"{statistics_prefix}_stacked_bar.{stacked_graph_format}")),
            dpi=stacked_dpi,
            format=stacked_graph_format)

time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")