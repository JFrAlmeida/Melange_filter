import re
import os
import time

#Config imports
from ..config import config_locals, annotation_columns_target, run_dbcan_parser

#Timer
time_start = time.time()

#Run dbcan_parser
if run_dbcan_parser:
    from dbcan_processing import process_dbcan_results, parse_dbcan_to_melange_ORF

#Makes a dictionary with the format cat_n : ["COGXX1","COGXX2" etc.], excludes empty lists
cats_with_stuff = {
    name: value
    for name, value in config_locals
    if name.startswith("cat_") and name[3:] and not value == []
}

#Makes a dictionary with the format name_cat_n : "GH18_endochitinase
cats_with_names = {
    name: value
    for name, value in config_locals
    if name.startswith("name_cat") and name[4:] and (not value == "" and not value == "empty")
}

#pattern and dict for finding the type of annotation (COG, CAZy, MEROPS, etc.) each category points to
pattern_type = re.compile(r"^type_(\d+)$")
types_temp = {}

#create a dict with the category annotation types
for name, value in config_locals:
    match = pattern_type.match(name)
    if value == "":
        continue
    if match:
        n = match.group(0)
        types_temp[n] = value

#Make a dict relating each cat_n to its type
cats_with_types = {
    name : value
    for name, value in zip(list(cats_with_stuff.keys()), list(types_temp.values()))
}

#List with the category labels (e.g. cat_1, cat_2 etc.)
stuff = list(cats_with_stuff.keys())

#List with the category name keys (e.g. name_cat_1, name_cat_2)
names = list(cats_with_names.keys())

#Make a list of what is inside all categories, maybe useless later, will see
all_cats =  [item
       for elements in cats_with_stuff.values()
       for item in elements
       ]


#Runtime calculation
time_finished = time.time()

print(f"{os.path.basename(__file__)},  took , {time_finished - time_start},  seconds to run \n proceeding to config_checks")


# print(cats_with_types)
# print(cats_with_stuff)