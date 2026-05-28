###########################################
# Annotation Filters
###########################################

"""
Define annotation-based categories used to select ORFs/features from Melange

Rules:
- Each category (cat_n) can contain as many annotations as you want from only ONE annotation type.
  Examples:
    cat_1 -> PFAM entries
    cat_2 -> PFAM entries
    cat_3 -> CAZYmes entries

Supported annotation types :
    - COG
    - KO
    - MEROPS
    - PFAM
    - CAZymes

(These depending on their availability from your Melange results, this script does not support the spontaneaous generation hypothesis)

Variables:
    cat_n       : List of annotation identifiers.
    name_cat_n  : Human-readable category name. Used for FASTA and statistics outputs.
    type_n      : Annotation type corresponding to cat_n.

Example:
    cat_1 = ["PF00182", "PF06483"]
    name_cat_1 = "Endochitinases"
    type_1 = "PFAM"
    
    cat_2 = ["COG3325", "COG3469"]
    name_cat_2 = "GH18_endochitinases"
    type_2 = "COG" #case sensitive! double check!
    
    cat_3 = ["GH18"]
    name_cat_3 = "GH18_endo"
    type_3 = "CAZymes" #case sensitive! double check!
    
"""

"""
Additional categories can be freely added or removed as long as they follow
the cat_n / name_cat_n / type_n structure.
"""

cat_1 = ["PF00182", "PF06483","PF08329"]
name_cat_1 = "Endochitinases"
type_1 = "PFAM" ### Case sensitive! double check!


###########################################
# Secondary Annotation Filter
###########################################

"""
Optional second-pass filtering step.

After category matching, results can be filtered again by selecting
specific annotation values from one or more annotation columns.

Parameters:
    last_filter          : List of annotation identifiers to retain.
                           Leave as [] to disable this filter.

    last_filter_columns  : Annotation columns to search.
                           Same syntax as annotation_columns_target.


"""

last_filter = []  #Leave blank if you don't want this to run
last_filter_columns = ["PFAM"] #same syntax as annotation_columns_target above


###########################################
# Annotation Finder Options
###########################################

"""
Annotation columns from the MELANGE annotation table that will be searched.

Supported values:
    "COG"
    "KO"
    "MEROPS"
    "PFAM"
    "CAZy"

Multiple columns may be provided in any order.

Example:
    ["PFAM", "CAZy"]
"""

annotation_columns_target = ["PFAM"]

"""
Controls how genome names are extracted for the 'genome_prokka_feats' output.

The matched text is combined with the corresponding PROKKA feature ID.

Options:
    "no_change"
        Use filenames as they are encountered.

    Custom regex
        Extract genome names using a user-defined pattern.

Useful regexes:

Removes GCA/GCF numbers from the genome names 
-> r"GC[AF]_\\d+\\.\\d+_(.*?)_all_features\\.csv"

Captures only Melange-named files
-> r"(.*)_all_features\\.csv"
"""

name_change_condition = "no_change"

###########################################
# dbCAN Parser Options
###########################################

"""
If enabled, the following steps are executed:

    1. Process_dbcan_results
    2. parse_dbcan_to_melange_ORF

Results in incorporating user-generated dbCAN results 
into the melange_Orfs_per_genome table, and into the 
final overview and assembly summary tables.
"""

run_dbcan_parser = False

"""
Controls genome name extraction for dbCAN-derived outputs.

Used by:
    - process_dbcan_results.py

Affects genome names stored in the CAZymes_counts table.

Options:
    "no_change"
    Custom regex

Default:
    r"^(?:[^_]*_){2}(.*)$"
"""


###### Unsure if this does anything, please check!!!
dbcan_name_change = "^(?:[^_]*_){2}(.*)$"


###########################################
# FASTA Generation Options
###########################################

"""
If True:
    Creates a fasta file at the end of Melange_filter

If False:
    FASTA generation is skipped.
"""

run_make_fasta = True

"""
Category groups for FASTA generation.

Each list element produces one FASTA file.

Examples:
    ["cat_1"]
        One FASTA file containing only cat_1 matches.

    ["cat_1!cat_2!cat_3"]
        One FASTA file containing all three categories.

    ["cat_1!cat_2", "cat_3"]
        Two FASTA files, the first with hits from cat_1 and cat_2 and 
        the second with hits from cat_3

Naming:
    The FASTA filename is derived from the name_cat_n value of
    the first category in each group.
"""

make_fasta_groups = ["cat_1"]

"""
Sequence type used for FASTA generation.

Options:
    True  -> DNA sequences
    False -> Amino acid sequences
"""

select_by_dna = False

###########################################
# Statistics Options
###########################################

"""
Prefix added to statistics output files
(e.g. presence counters and stacked bar plots).

Not necessary to run
"""

statistics_prefix = "Meunier"

"""
If True, statistics files are stored between runs.
The Statistics/ folder gets wiped every run, but the 
directory below does not

Options:
    True
    False
"""

statistics_multiple_runs = True




#### Unsure if this actually does anything, please check!!
"""
Directory used to store statistics outputs when
statistics_multiple_runs is enabled. Placed inside Outputs


"""

save_directory = "Statistics_general"


###########################################
# Presence Counter Options
###########################################

"""
Run the presence_counter module. Counts how many
sequences are annotated in each group defined below
in each genome

Options:
    True  -> Execute presence counter
    False -> Skip presence counter
"""

run_presence_counter = True

"""
If True, outputs are additionally saved to:

    Outputs/Stats_perm

This directory is not cleaned between runs.

options:
    True
    False

"""

save_to_perm = False

"""
Category groups used by presence_counter.

Uses the same syntax as make_fasta_groups and is typically identical.

Example:
presence_counter_groups = ["cat_1!cat_2", "cat_3"]
        Two groups, the first with the count of hits from cat_1
        and cat_2 and the second with hits from cat_3
"""

presence_counter_groups = ["cat_1"]


###########################################
# Stacked Bar Plot Options
###########################################

"""
Generate a stacked bar plot from presence_counter results.

options:
    True
    False
"""

run_stacked_bar_plot = True

"""
Optional metadata file located in:
    Melange_filter/MF_statistics/

Allows setting the order of the stacked bar plot,
and setting the color for each label

Requires two columns columns:
    Tree_order:
        Controls sample order in the figure.

    Color:
        Controls label colors.
"""

stacked_order_color_file_name = ""

"""
Output figure format.

Examples:
    png
    pdf
    svg
    eps

See matplotlib.pyplot.savefig documentation for additional formats.
"""

stacked_graph_format = "svg"

"""
Figure resolution (in DPI).
"""
stacked_dpi = 300

"""
Width of individual bars. 0.4 is good for tight graphs
"""
stacked_width = 0.4

"""
Figure title.
"""
stacked_graph_title = "Endochitinase gene counts"

"""
Prefix used in output filenames.

Example:
    Meunier_stacked_bar.svg
"""

stacked_graph_name_prefix = "Meunier"


###########################################
# Stacked Bar Label Options
###########################################

"""
Font size for sample labels.
"""
stacked_fontsize = 10

"""
Font style.

Options:
    normal
    italic
    bold
"""

stacked_fontstyle = "italic"

"""
Alignment of the labels slabels.

Matplotlib values:
    left
    center
    right
"""

stacked_ha = "right"

"""
Label rotation angle in degrees.
"""

stacked_rotation = 45


###########################################
# PFAM Summary Table Options
###########################################

"""
Generate a summary table of PFAMs with count higher than the threshold.

Options:
    Integer value -> minimum number of hits required.
    Empty value   -> disable table generation.
"""

pfam_topn_list = 10


###########################################
# End of User-Editable Configuration
###########################################

"""
Do not modify.

These variables are required by downstream scripts to
retrieve configuration values from this file.
"""

import os
config_locals = locals().items()

# Common paths
path_main_melange = os.path.dirname(os.getcwd())

