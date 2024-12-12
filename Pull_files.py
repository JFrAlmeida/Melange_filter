import os
import shutil

# Path to the main folder containing subfolders
main_folder = '/home/jfa/Aquimarina_review/Aquimarina_genomes/all_genomes/data_sans_useless'

# Path to the target folder where all files should be moved
target_folder = '/home/jfa/Aquimarina_review/Aquimarina_genomes/all_genomes/all_together_now'

# Loop through the main folder to find subfolders
for subfolder in os.listdir(main_folder):
    subfolder_path = os.path.join(main_folder, subfolder)

    # List all files inside the subfolder
    files = os.listdir(subfolder_path)
        
        # Assuming there is only one file per subfolder
    for file in files:
        file_path = os.path.join(subfolder_path, file)

            # Move the file to the target folder
        shutil.copy(file_path, target_folder)
        print(f"copied: {file_path} to {target_folder}")
