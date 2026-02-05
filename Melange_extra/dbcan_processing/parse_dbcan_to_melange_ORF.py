import os
import pandas as pd
import shutil
import time

#my imports
from config import path_main_melange

time_start = time.time()
print(f"starting {os.path.basename(__file__)}...")

#paths:
dbcan_folder_path = os.path.normpath(os.path.join(path_main_melange, "processed_dbcan"))
melange_ORF_path = os.path.normpath(os.path.join(path_main_melange, "Annotation_results/Orfs_per_genome"))
cazy_updated_path = os.path.normpath(os.path.join(path_main_melange, "Annotation_results/Orfs_per_genome_CAZy"))

#File terminations, in case it changes in the future:
dbcan_term = "_dbcan_clean.csv"
mel_term = "_all_features.csv"

#Make folders
if os.path.exists(cazy_updated_path):
    shutil.rmtree(cazy_updated_path)
    os.makedirs(cazy_updated_path)
else:
    os.makedirs(cazy_updated_path)

#QA check
if len(os.listdir(dbcan_folder_path)) != len(os.listdir(melange_ORF_path)):
    print("dbcan_folder",os.listdir(dbcan_folder_path))
    print("dbcan_folder", os.listdir(melange_ORF_path))
    print(f"it appears you have a different number of files in your Melange ORF folder -> {melange_ORF_path}"
          f"than in your dbcan folder -> {melange_ORF_path}; "
          f"\n sort them out then run again please!")
    print(f"your dbcan folder path is: {dbcan_folder_path}")
    print(f"your Melange ORF path is: {melange_ORF_path}")
    print("quitting...")
    quit()

#Loop goes through files in processed_dbcan
for file in os.listdir(dbcan_folder_path):

    #make the path for dbcan table
    dbcan_file_path = os.path.join(dbcan_folder_path,file)

    #create path to the corresponding file in Orfs_per_genome in melange results
    melange_file_name = file.replace(dbcan_term, mel_term)
    file_melange = os.path.normpath(os.path.join(melange_ORF_path,melange_file_name))

    #checks if the melange file exists, if not quits
    if not os.path.isfile(file_melange):
        print(f"the file {file_melange} exists in processed_dbcan, but not Orfs_per_genome, please ensure your dataset "
              f"is correct! \n quitting...")
        quit()

    #Load both files as
    df_dbcan = pd.read_csv(dbcan_file_path)
    df_mel = pd.read_csv(file_melange)

    #Move final_annotation and signalP to the end of the melange annotation
    df_mel["CAZymes"] = df_dbcan["final_annotation"]
    df_mel["Signalp"] = df_dbcan["Signalp"]
    df_mel = df_mel[['row_0', 'COG', 'KO', 'MEROPS','CAZymes','Signalp', 'PFAM']] #reorder

    #write it to cazy_updated_path
    df_mel.to_csv(os.path.join(cazy_updated_path, melange_file_name),sep= ",", index=False)

time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")