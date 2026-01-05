import os
import pandas as pd
import shutil


#purpose is to add the cazymes info created from process_dbcan_results into Melange ORFs table

#Get the set of files from one and another, compare if the same, then good, other spit out the ones whose name is
# different

#Import the dbcan files

#Import paths:
dbcan_folder_path = "./processed_dbcan"
melange_ORF_path = "./Annotation_results/Orfs_per_genome"
cazy_updated_path = "./Annotation_results/Orfs_per_genome_CAZy"

#File terminations, in case it changes in the future:
dbcan_term = "_dbcan_clean.csv"
mel_term = "_all_features.csv"

#Make folders
if os.path.exists("./Annotation_results/Orfs_per_genome_CAZy"):
    shutil.rmtree("./Annotation_results/Orfs_per_genome_CAZy")
    os.makedirs("./Annotation_results/Orfs_per_genome_CAZy")
else:
    os.makedirs("./Annotation_results/Orfs_per_genome_CAZy")



#hardset variables


#print(os.listdir(dbcan_folder_path))
#print(os.listdir(melange_ORF_path))

if len(os.listdir(dbcan_folder_path)) != len(os.listdir(melange_ORF_path)):
    print(f"it appears you have a different number of files in your Melange ORF folder -> {melange_ORF_path}"
          f"than in your dbcan folder -> {melange_ORF_path}; "
          f"\n sort them out then run again please!")
    quit()

#make a for loop from dbcan, first check if all items have the same number of stuff

#Loop goes through files in processed_dbcan
for file in os.listdir(dbcan_folder_path):

    #make the path for dbcan table
    dbcan_file_path = os.path.join(dbcan_folder_path,file)

    #create path to the corresponding file in Orfs_per_genome in melange results
    melange_file_name = file.replace(dbcan_term, mel_term)
    file_melange = os.path.join(melange_ORF_path,melange_file_name)

    #checks if the melange file exists, if not quits
    if not os.path.isfile(file_melange):
        print(f"the file {file_melange} exists in processed_dbcan, but not Orfs_per_genome, please ensure your dataset "
              f"is correct! \n quitting...")
        quit()

    #Load both files as
    df_dbcan = pd.read_csv(dbcan_file_path)
    df_mel = pd.read_csv(file_melange)

    #Move final_annotation to the end of the melange annotation
    df_mel["CAZymes"] = df_dbcan["final_annotation"]
    df_mel = df_mel[['row_0', 'COG', 'KO', 'MEROPS','CAZymes', 'PFAM']]


    #write it to cazy_updated_path
    df_mel.to_csv(os.path.join(cazy_updated_path, melange_file_name),sep= ",", index=False)

