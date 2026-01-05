import os
import re
import time
import shutil
import numpy as np
import pandas as pd





#-------------------  Functions below  -------------------#

#Takes a dict, key, separator and a string object to store results; Used in script to parse cazymes annotations
def Annotation_builder(dict, key, string, sep):
    for n in range(0, int(dict[key])):
        if string == "":
            string = string + key
        else:
            string = string + sep + key
    return string


#-------------------  Functions above  -------------------#


#Timer
time_start = time.time()

#-------------------  hard-set variables  -------------------#

#separator used in each cell to separate annotations (GH13+CBM6+PL4)
sep="+"

#Relative path to input folder
input_folder = os.path.normpath("dbcan_output/")

#Relative path to the desired Output folder
output_folder = "processed_dbcan/"

overview_notfound_counter = 0

############################## Code starts ##############################

#Functions

#pandas_df is the target df, column is your target column, conditions is a list of any numebrs to be passed to pandas .isin()
def pandas_rows_by_listofnumbers(pandas_df,column, conditions):
    finder = pandas_df[column].isin(conditions)
    new_df = pandas_df.loc[finder].copy()
    return new_df

#Check if the input folder exists and has stuff
if not os.path.exists(input_folder):
    print("your input folder was not detected, check the name! \n quitting...")
    quit()
elif len(os.listdir(input_folder)) == 0:
    print("your input folder is empty! \n quitting...")
    quit()

#Make necessary folders
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
else: #Clean out the output folder before running
    shutil.rmtree(output_folder)
    os.makedirs(output_folder)

for folder in os.listdir(input_folder):
    genome_name = folder
    genome_overview_path = os.path.join(genome_name, "overview.txt")
    overview_path = os.path.join(input_folder, genome_overview_path)


    try:
        #Catch events where overview cannot be found
        df = pd.read_csv(overview_path,
                         sep= "\t").drop(
            labels= "EC#",
            axis= "columns")
    except:
        print(f"overview.txt could not be found inside {genome_name}, skipping...")
        overview_notfound_counter = overview_notfound_counter + 1



    df_good_hits = pandas_rows_by_listofnumbers(df,"#ofTools",[3,2])


    #The following block cleans out the annotation columns of their weird annotations, standardizing them

    #HMMER column
    #This cleans out the things between parenthesis in HMMMER annotation
    df_good_hits["HMMER_clean"] = df_good_hits["HMMER"].str.replace(pat= "\([^()]*\)", repl= "", regex= True)

    #This cleans out the sub_family annontation of HMMER "PL7_5" to "PL7"
    df_good_hits["HMMER_clean"]= df_good_hits["HMMER_clean"].str.replace(pat= "_\d+", repl= "", regex= True)

    #dbCAN_sub column
    #cleans out the _enumbers
    df_good_hits["dbCAN_sub_clean"] = df_good_hits["dbCAN_sub"].str.replace(pat= "\_e\d+", repl= "", regex= True)

    #DIAMOND column
    df_good_hits["DIAMOND_clean"] = df_good_hits["DIAMOND"].str.replace(pat= "_\d+", repl= "", regex= True)
    df_good_hits["DIAMOND_clean"] = df_good_hits["DIAMOND_clean"].str.replace(pat= "\d+(?:\.\d+)+", repl= "", regex= True)



    #requires something in the position of "None" (last parameter); breaks if "" or Nan is given. couldn't be bothered
    # to figure out why, later in the code I just replace it with ""
    df_good_hits["inter_anno"] = np.where(
        (df_good_hits["HMMER_clean"] == df_good_hits["dbCAN_sub_clean"]) & (df_good_hits["HMMER_clean"] == df_good_hits["DIAMOND_clean"]),
        df_good_hits["HMMER_clean"],
        "None"
    )

    df_good_hits["inter_anno2"] = df_good_hits["inter_anno"].replace("None","")

    #dataframe slice where only lines where the 3 annotations are not the same
    df_not_equal_anno = df_good_hits[df_good_hits["inter_anno"] == "None"]

    # Make an object to hold the final series to be added back to the pandas df
    index_list_final_series = []
    data_list_final_series = []

    #This block will count how many instances of each thing exists in each tool, and compare it to the others.
    #Objective is to standerdize how many counts of each hit are to be considered, rule is at least two hits of the same things
    for item1,item2,item3, my_index in zip(df_not_equal_anno.HMMER_clean,
                                           df_not_equal_anno.dbCAN_sub_clean,
                                           df_not_equal_anno.DIAMOND_clean,
                                           df_not_equal_anno.index):

        # List with the row-wise entries, including index to merge back later
        row_list = [my_index, item1, item2, item3]

        # print(f"row_list: {row_list}")

        #This small block will make three lists with the annotations of each tool as a single item
        #HMMER_clean
        if "+" in item1:
            list_item1 = item1.split("+")
            # print(f"list_item1 '+': {list_item1}")
        else:
            list_item1 = [item1]
            # print(f"list_item1: {list_item1}")

        #dbCAN_sub_clean
        if "+" in item2:
            list_item2 = item2.split("+")
            # print(f"list_item2 ''+': {list_item2}")

        else:
            list_item2 = [item2]
            # print(f"list_item2: {list_item2}")


        #DIAMOND_clean
        if "+" in item3:
            list_item3 = item3.split("+")
            # print(f"list_item3 ''+': {list_item3}")
        else:
            list_item3 = [item3]
            # print(f"list_item3: {list_item3}")

        #Initial dicts to be used below
        dict1 = dict(); dict2 = dict(); dict3 = dict()

        # now lets create dictionaries with counts of each separate item
        for key in list_item1:
            if not key == "-":
                dict1.update({key : list_item1.count(key)})
            else:
                dict1 = {}

        for key in list_item2:
            if not key == "-":
                dict2.update({key : list_item2.count(key)})
            else:
                dict2 = {}

        for key in list_item3:
            if not key == "-":
                dict3.update({key : list_item3.count(key)})
            else:
                dict3 = {}

        # print(f"dict1: {dict1}, dict2: {dict2}, dict3: {dict3}")

        #append index to the index_list
        index_list_final_series.append(my_index)

        #This stores the outcome of string_end each pass through the key loop below, in case multiple keys exist
        interm_key_outcome = ""

        #Now to use those dictionaries to compare how many of which annotations exist in each tool; For example, if "CBM50"
        # is present in 2 copies in one annotation, and in 5 in another annotation, it counts as 2; This block will always
        # take the smaller number of annotations of each hit, and require at least 2 presences
        for key in sorted(list(set([*dict1, *dict2, *dict3])),reverse=True):
            # print(f"Index: {my_index}; Set selection {set([*dict1, *dict2, *dict3])}")
            #Check which dicts to use for the following checks
            try:
                use1 = dict1[key] != "-"
            except:
                use1 = False
            try:
                use2 = dict2[key] != "-"
            except:
                use2 = False
            try:
                use3 = dict3[key] != "-"
            except:
                use3 = False

            #This controls whether if statements down below are run or not
            solution_found = False

            #to build the actual final annotation column
            string_end = ""

            #Require presence in at least two annotation tools; keyindictn is a checker for if the key exists in that tool
            keyindict1 = key in dict1; keyindict2 = key in dict2; keyindict3 = key in dict3
            keycheck = keyindict1 + keyindict2 + keyindict3

            # print(f"keycheck is {keycheck}, keyindict1 is: {keyindict1}, keyindict2 is: {keyindict2}, keyindict3 is: {keyindict3}")

            #Only include that annotation if it is annotated by at elast two tools, otherwise skip it
            if keycheck >= 2:

                # If key in all 3 tools
                if keycheck == 3:
                    # print("I am on all 3 keys")

                    #In case where 2 annotation tools agree on the count of one annotation, but the third doesn't, take the
                    # case where 2 agree
                    if not solution_found: #These means basically, if a solution was found before dont run this code!
                        if dict1[key] == dict2[key] and (dict3[key] < dict1[key] or dict3[key] > dict1[key]):
                            string_end = Annotation_builder(dict1, key, string_end, sep=sep)
                            solution_found = True
                            # print("check 0.3")

                    if not solution_found: #These means basically, if a solution was found before dont run this code!
                        if dict2[key] == dict3[key] and (dict1[key] < dict2[key] or dict1[key] > dict2[key]):
                            string_end = Annotation_builder(dict2, key, string_end, sep=sep)
                            solution_found = True
                            # print("check 0.6")

                    if not solution_found:  # These means basically, if a solution was found before dont run this code!
                        if dict3[key] == dict1[key] and (dict2[key] < dict3[key] or dict2[key] > dict3[key]):
                            string_end = Annotation_builder(dict3, key, string_end, sep=sep)
                            solution_found = True
                            # print("check 0.9")

                    #In case no 2 tools have the same number of a given annotation
                    if not solution_found: #These means basically, if a solution was found before dont run this code!
                        #dict 2 v 3
                        if dict2[key] <= dict3[key]:
                            #then 2 v 1
                            if dict2[key] <= dict1[key]:
                                string_end = Annotation_builder(dict2, key, string_end, sep=sep)
                                solution_found = True
                                # print("check1")

                            else:
                                string_end = Annotation_builder(dict1, key, string_end, sep=sep)
                                solution_found = True
                                # print("check2")

                        else:
                            #else dict3 v 1
                            if dict3[key] <= dict1[key]:
                                string_end = Annotation_builder(dict3, key, string_end, sep=sep)
                                solution_found = True
                                # print("check3")

                            else:
                                string_end = Annotation_builder(dict1, key, string_end, sep=sep)
                                solution_found = True
                                # print("check4")

                #If keycheck is 2, then key must not be in one of the tools, and loop jumps here, do the comparison below

                if not solution_found:  # These means basically, if a solution was found before dont run this code!
                    # if its not in dict1
                    if not keyindict1:
                        if dict2[key] <= dict3[key]:
                            string_end = Annotation_builder(dict2, key, string_end, sep=sep)
                            solution_found = True
                            # print("check5")

                        else:
                            string_end = Annotation_builder(dict3, key, string_end, sep=sep)
                            solution_found = True
                            # print("check6")

                    #Now if its not on dict2
                    if not keyindict2:
                        if dict1[key] <= dict3[key]:
                            string_end = Annotation_builder(dict1,key,string_end,sep=sep)
                            solution_found = True
                            # print("check7")

                        else:
                            string_end = Annotation_builder(dict3,key,string_end,sep=sep)
                            solution_found = True
                            # print("check8")

                    # Now if its not on dict3
                    if not keyindict3:
                        if dict1[key] <= dict2[key]:
                            string_end = Annotation_builder(dict1,key,string_end,sep=sep)
                            solution_found = True
                            # print("check9")

                        else:
                            string_end = Annotation_builder(dict2,key,string_end,sep=sep)
                            solution_found = True
                            # print("check10")

            # print(f"this is the string: {string_end}, my index is: {my_index}")

            if solution_found: # only add to string if a solution has been found
                if interm_key_outcome != "":
                    interm_key_outcome = interm_key_outcome + sep + string_end
                else:
                    interm_key_outcome = string_end

                # print(f"this interm_key_outcome: {interm_key_outcome}")


        data_list_final_series.append(interm_key_outcome)

    #make a pandas series from the index list and the list built from interm_key_outcome
    anno_series = pd.Series(data=data_list_final_series, index=index_list_final_series, name="inter_ann3")

    #make it a column now
    df_good_hits["inter_anno3"] = anno_series

    #Concatenate both annotation columns
    df_good_hits["final_annotation"] = df_good_hits["inter_anno2"].astype(str) + df_good_hits["inter_anno3"].astype(str)

    #Clean columns, as Nan gets converted to str as well as an artifact
    df_good_hits["final_annotation"] = df_good_hits["final_annotation"].str.replace("nan", repl= "")

    #Clean out the intermediate columns for the final annotation
    df_good_hits = df_good_hits.drop(labels= ["inter_anno3",
                                              "inter_anno2",
                                              "inter_anno"],
                                     axis=1
                                     )

    #Clean out additionally the now outdated annotations
    df_good_hits_clean_absolut = df_good_hits.drop(labels= ["HMMER",
                                                            "dbCAN_sub",
                                                            "DIAMOND",
                                                            "#ofTools"],
                                                   axis=1
                                                   )
    #Export table with the old annotations as included
    # df_good_hits.to_csv("/home/jfa/Aquimarina_review/process_dbcan/almost_clean.tsv", sep= "\t", index=False)

    #Outgoing path
    export_path = os.path.normpath(output_folder + genome_name + "_dbcan_clean.csv")

    #Full cleaned table with only cleaned annotations and final annotation (and PROKKA features
    df_good_hits_clean_absolut.to_csv(export_path, sep= ",", index=False)


#print how many genomes did not have an overview file
if overview_notfound_counter > 0:
    print(f"The Overview.txt file, expected inside the folder of each genome as the output of dbcan was not found for "
          f"{overview_notfound_counter} genomes")

time_finished = time.time()

print(f"{os.path.basename(__file__)}  took , {time_finished - time_start},  seconds to run \n Thank you for using this script!! JFA")
