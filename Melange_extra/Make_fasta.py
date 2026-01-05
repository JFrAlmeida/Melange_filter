import os
import time
import pandas as pd

# Import variables from modules:
from config import name_change_condition, annotation_columns_target
from category_checks import cats_with_stuff, cats_with_types


time_start = time.time()


#Use the concatenated genomes table. Then search for each category and type.

#One by one search all feats in that cat, then do a copy of that table and export it with the name from name_cat_n

# use the function to generate the actual fasta file contained in the old fasta file.

#Implement a choice for user, string vs list and use .upper() to distinguish them within a try statement





#     for item in annotation_columns_target: #item is the target column ("COG" "KO" etc.)
#         for type in list(cats_with_types.keys()): #type here is the key to the type of annotation in config ("COG" "KO"
#             # etc.), its value needs to match item to proceed
#             print(f"catswith types: {cats_with_types[type]}, and item: {item}")
#             if cats_with_types[type] == item:
#                 for category in list(cats_with_stuff.keys()): #category is cat_1, cat_2 etc.
#                     for anno in cats_with_stuff[category]: #anno is the annotation targets, "COG132" "COG312" etc.
#                         # print(f"looking for {anno} in {item}")
#                         finder = genome[item].str.contains(anno, case=False, na=False) # True where match is found
#                         if sum(finder) != 0: # to avoid operations when no information is added
#                             if not filter_bool.size == 0:
#                                 filter_bool = filter_bool + finder # add finder to the boolean mask of other loops
# # inside the same genome. Builds a large boolean mask aggregating all True hits over every cat_
#                             else:
#                                 filter_bool = finder # if its not made yet, make it


time_finished = time.time()

print(f"{os.path.basename(__file__)}, took , {time_finished - time_start}, seconds to run")