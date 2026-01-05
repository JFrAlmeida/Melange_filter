import re
import os

#Objective is to hold variables to be share bewteen the melange processing scripts here
#no need to add file extensions (e.g. .fasta) to the names

#edit each category



########################## Edit this part of the script!!! ################################################

#--------------------------
#This specifies the target column of MELANGE annotation to select for
# Valid options are ["COG"], ["CAZymes"], ["KO"], ["MEROPS"], ["PFAM"], in any order, and any number and combination of
# them, as long as inside "" and separated by a comma; only COGS has been tested, but all exept PFAM should work
annotation_columns_target = ["COG"]
#--------------------------
#Variabes to tell program if you want to also include rows that contain specific auxiliary activities (AA), or
# Polyssacharide lyase (PL) enzymes; Carbohydrate esterases (CE) are not supported by MELANGE.
# Glicosyl Transferases (GT) are included in the CAZYmes column with gycosyl Hydrosyl (GH) enzymes;
#Valid options: ["none"], ["AA"],["PL"] and ["AA,"PL]
additional_cazymes = ["none"]
cat_cazy_AA = ["AA10"] #If you want specific AA cazy families, specify them here
name_cat_cazy_PL = ["PL1"] #If you want specific PL cazy families, specify them here

#---
# This variable specifies how the genome name should be edited for the result of this script; The default assumes your
# genome has the format GCA_XXXXXXXXX.X_whatever_strain_names_exist_all_features.csv;
#If nothing is speficied, or something other than the valid options it falls back to default
#Valid options are "default", "no_change", "Custom". For custom, add your own regex below to custom_regex
name_change_condition = "default"
# custom_regex = "REGEX CONDITION" #place your own regex condition here and uncomment this line
#---

#---

#Important, each cat can only take one type of annotation, only COGs, or only Cazy, or only MEROPS etc.

cat_1 = ["COG3325", "COG3469"] # GH18 endochitinases
name_cat_1 = "GH18_endochitinases" #will be the name of the file outputed
type_1 = "COG"

cat_2 = ["COG3179"] # GH19 endochitinases
name_cat_2 = "GH19_endochitinases" #will be the name of the file outputed
type_2 = "COG"

cat_3 = ["COG4724", "COG3525", "COG3979"] # exochitiniases
name_cat_3 = "Exochitiniases" #will be the name of the file outputed
type_3 = "COG"

cat_4 = ["COG3397"] # LPMOs
name_cat_4 = "LPMOs"
type_4 = "COG"

cat_5 = ["COG3394", "COG0726", "COG2861"] # Deacitylases
name_cat_5 = "Deacetylases" #will be the name of the file outputed
type_5 = "COG"

cat_6 = [] # empty
name_cat_6 = "empty" #will be the name of the file outputed
type_6 = "COG"

cat_7 = [] # empty
name_cat_7 = "empty" #will be the name of the file outputed
type_7 = "COG"

cat_8 = [] # empty
name_cat_8 = "empty" #will be the name of the file outputed
type_8 = "COG"

cat_9 = [] # empty
name_cat_9 = "empty" #will be the name of the file outputed
type_9 = "COG"

cat_10 = [] # empty
name_cat_10 = "empty" #will be the name of the file outputed
type_10 = "COG"

########################## Dont edit anything below!!!! ################################################

config_locals = locals().items()

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



#Code below! Dont touch!





























