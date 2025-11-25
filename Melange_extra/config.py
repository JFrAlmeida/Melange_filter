import re
import os

#Objective is to hold variables to be share bewteen the melange processing scripts here
#no need to add file extensions (e.g. .fasta) to the names

#edit each category



########################## Edit this part of the script!!! ################################################

#This specifies the target column of MELANGE annotation to select for
# Valid options are "COG", "CAZymes", "KO", "MEROPS", "PFAM"; only COGS has been tested
annotation_column_target = "COG"

cat_1 = ["COG3325", "COG3469"] # GH18 endochitinases
name_cat_1 = "GH18_endochitinases" #will be the name of the file outputed

cat_2 = ["COG3179"] # GH19 endochitinases
name_cat_2 = "GH19_endochitinases" #will be the name of the file outputed

cat_3 = ["COG4724", "COG3525", "COG3979"] # exochitiniases
name_cat_3 = "Exochitiniases" #will be the name of the file outputed

cat_4 = ["COG3397"] # LPMOs
name_cat_4 = "LPMOs"

cat_5 = ["COG3394", "COG0726", "COG2861"] # Deacitylases
name_cat_5 = "Deacetylases" #will be the name of the file outputed

cat_6 = [] # empty
name_cat_6 = "empty" #will be the name of the file outputed

cat_7 = [] # empty
name_cat_7 = "empty" #will be the name of the file outputed

cat_8 = [] # empty
name_cat_8 = "empty" #will be the name of the file outputed

cat_9 = [] # empty
name_cat_9 = "empty" #will be the name of the file outputed

cat_10 = [] # empty
name_cat_10 = "empty" #will be the name of the file outputed

########################## Dont edit anything below!!!! ################################################

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

#Makes a dictionary with the format cat_n : ["COGXX1","COGXX2" etc.], excludes empty lists
cats_with_stuff = {
    name: value
    for name, value in locals().items()
    if name.startswith("cat_") and name[3:] and not value == []
}

#Makes a dictionary with the format name_cat_n : "GH18_endochitinase
cats_with_names = {
    name: value
    for name, value in locals().items()
    if name.startswith("name_") and name[4:] and (not value == "" and not value == "empty")
}
# print(cats_with_names)

#List with the category labels (e.g. cat_1, cat_2 etc.)
stuff = list(cats_with_stuff.keys())

#List with the category name keys (e.g. name_cat_1, name_cat_2)
names = list(cats_with_names.keys())

#Make a list of what is inside all categories, maybe useless later, will see
all_cogs =  [item
       for elements in cats_with_stuff.values()
       for item in elements
       ]

#This block checks if the number of categories to category names is the same or not as an initial quality check for the
# user; If the number is the same, and each category has its counterpart name, proceeds, otherwise prints a warning and
# quits
if len(cats_with_stuff) == len(cats_with_names):
    if  (one_stuff in one_name for
        all_stuff, all_names in
        zip(stuff, names)):
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


print(f"Category contents are ok! Proceeding with {os.path.basename(__file__)}"
                  f"\n (ฅ^•ﻌ•^ฅ)")




























