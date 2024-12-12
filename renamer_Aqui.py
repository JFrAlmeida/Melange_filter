import os
import pandas as pd

# Load Excel file
excel_path = "/home/jfa/Aquimarina_review/Aquimarina_genomes/all_genomes/namer_2.xlsx"  # Update this path, point to where the file is
data = pd.read_excel(excel_path)

# Define the folder containing the subfolders
base_folder = "/home/jfa/Aquimarina_review/Aquimarina_genomes/all_genomes/data_sans_useless"  # Update this path

# Iterate through the Excel rows
for index in data.index:
    row = data.loc[index]  # Access the row using the index
    folder_path = row['Assembly Accession']  # Row with current folder names
    new_file_name = row['Name final']  # Row with new file names
    GCA_folder = os.path.join(base_folder, folder_path)

    if os.path.isdir(GCA_folder):
        #print("found folder" + folder_path)
        # Find the existing file in the folder
        files = os.listdir(GCA_folder)

        if len(files) == 1:  # Assumes only one file per folder
            old_file_path = os.path.join(GCA_folder, files[0])
            new_file_path = os.path.join(GCA_folder, new_file_name)

            #Check if ".fasta" is not in the path name, if it isnt, add ".fasta", and then rename, otherwise just rename
            if ".fasta" not in new_file_path:
                new_file_path_fasta = new_file_path + ".fasta"

                # Rename the file
                os.rename(old_file_path, new_file_path_fasta)
                print(f"Renamed {old_file_path} to {new_file_path_fasta}")
            else:
                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f"Renamed {old_file_path} to {new_file_path}")
        else:
            print(f"Multiple or no files in folder {GCA_folder}. Skipping...")
    else:
        print(f"Folder {GCA_folder} not found. Skipping...")
