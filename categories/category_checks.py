import re
import os
import time

#my imports
from categories.category_processing import cats_with_stuff, cats_with_types, cats_with_names
from config import annotation_columns_target

time_start = time.time()


#This script checks if the categories in config are well defined or not!

#This block checks if the number of categories to category names is the same or not as an initial quality check for the
# user; If the number is the same, and each category has its counterpart name, proceeds, otherwise prints a warning and
# quits

#List with the category labels (e.g. cat_1, cat_2 etc.)
stuff = list(cats_with_stuff.keys())

#List with the category name keys (e.g. name_cat_1, name_cat_2)
names = list(cats_with_names.keys())

if len(cats_with_stuff) == len(cats_with_names): #if number of categories of cats and name_cats is ==
    if  (one_stuff in one_name for
        all_stuff, all_names in
        zip(stuff, names)):  #if for every number of catN you have name_catN
        print("you have the same number of categories and category names!"
              " ദ്ദി(• ˕ •マ.ᐟ"
              "\n proceeding with checks...")

    else:
        print("your cats and categories are the same in number but have different names! please correct and run again! "
              "/ᐠﹷ ‸ ﹷ ᐟ\ﾉ"
              "\n Quitting...")
        quit()

elif len(cats_with_stuff) > len(cats_with_names):
    print("You seem to have more cats than name_cats, please name all your cats! "
          "/ᐠﹷ ‸ ﹷ ᐟ\ﾉ"
          "\n Quitting...")
    quit()

elif len(cats_with_names) > len(cats_with_stuff):
    print("You seem to have more named_cats than cats exist, please remove names or add more cats! "
          "/ᐠﹷ ‸ ﹷ ᐟ\ﾉ"
          "\n Quitting...")
    quit()




#this large block check if the categories were inserted correctly or not. It iterates through the total of annotations
# requested, and additional cazymes columns

#Need to make one for CAZy and the others
for item in annotation_columns_target:
    if annotation_columns_target == "COG":
        for key in cats_with_stuff.keys():
            for item in cats_with_stuff[key]:
                if not re.match(r"^COG\d+$", item): #works
                    if item == "":
                        print(f'An item in category {key} has the characters "" in it instead of a COG number, please eliminate'
                              ' them/correct it and run again! (๑ᵔ⤙ᵔ๑)'
                              '\n Quitting...')
                        quit()
                    if item == "COG": #works
                        print(f"An item in category {key} is missing the numbers after the 'COG' letters, please correct it and"
                              " run again! (๑ᵔ⤙ᵔ๑)"
                              "\n Quitting...")
                        quit()
                    if re.match(r"^[A-Za-z]+$", item): #works
                        print(f"An item in category {key} is composed only of letters other than COG, please correct and"
                              " run again! (๑ᵔ⤙ᵔ๑)"
                              "\n Quitting...")
                        quit()
                    if re.match(r"^\d+$", item): #works
                        print(f"OCategory {key} is composed only of numbers, please correct it and run again! (๑ᵔ⤙ᵔ๑)"
                              f"\n Quitting...")
                        quit()
                    else:
                        print("Something strange is going on with your categories, please verify them and run again! (๑ᵔ⤙ᵔ๑)")
                        quit()
    elif annotation_columns_target == "KO":
        for key in cats_with_stuff.keys():
            for item in cats_with_stuff[key]:
                if not re.match(r"^K\d+$", item):
                    if item == "":
                        print(
                            f'An item in category {key} has the characters "" in it instead of a K number, please eliminate'
                            ' them/correct it and run again! (๑ᵔ⤙ᵔ๑)'
                            '\n Quitting...')
                        quit()
                    if item == "K":
                        print(f"An item in category {key} is missing the numbers after the 'K' letter, please "
                              "correct it and run again! (๑ᵔ⤙ᵔ๑)"
                              "\n Quitting...")
                        quit()
                    if re.match(r"^[A-Za-z]+$", item):
                        print(
                            f"An item in category {key} is composed only of letters other than K, please correct and"
                            " run again! (๑ᵔ⤙ᵔ๑)"
                            "\n Quitting...")
                        quit()
                    if re.match(r"^\d+$", item):
                        print(f"OCategory {key} is composed only of numbers, please correct it and run again! (๑ᵔ⤙ᵔ๑)"
                              f"\n Quitting...")
                        quit()
                    else:
                        print(
                            "Something strange is going on with your categories, please verify them and run again! (๑ᵔ⤙ᵔ๑)")
                        quit()
    elif annotation_columns_target == "PFAM":
        for key in cats_with_stuff.keys():
            for item in cats_with_stuff[key]:
                if not re.search(r"^PF\d+$", item):
                    if item == "":
                        print(
                            f'An item in category {key} has the characters "" in it instead of a PFAM number, please eliminate'
                            ' them/correct it and run again! (๑ᵔ⤙ᵔ๑)'
                            '\n Quitting...')
                        quit()
                    if item == "PF":
                        print(f"An item in category {key} is missing the numbers after the 'PF' letters, please "
                              "correct it and run again! (๑ᵔ⤙ᵔ๑)"
                              "\n Quitting...")
                        quit()
                    if re.match(r"^[A-Za-z]+$", item):
                        print(
                            f"An item in category {key} is composed only of letters, please correct and"
                            " run again! (๑ᵔ⤙ᵔ๑)"
                            "\n Quitting...")
                        quit()
                    if re.match(r"^\d+$", item):
                        print(f"OCategory {key} is composed only of numbers, please correct it and run again! (๑ᵔ⤙ᵔ๑)"
                              f"\n Quitting...")
                        quit()
                    else:
                        print(
                            "Something strange is going on with your categories, please verify them and run again! (๑ᵔ⤙ᵔ๑)")
                        quit()



#REMEMBER TO MAKE THE annotation_columns_target into  A FOR LOOP

#Runtime calculation
time_finished = time.time()

print(f"{os.path.basename(__file__)} took {time_finished - time_start}  seconds to run")

print(f"Category contents are ok! (ฅ^•ﻌ•^ฅ)"
      f"\n proceeding to batch_COG_matcher ")


