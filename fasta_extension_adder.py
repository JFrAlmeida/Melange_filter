import os

# Specify the directory path
directory_path = "/home/jfa/Aquimarina_review/Aquimarina_genomes/all_genomes/all_together_now"  # Replace with your directory path

# Get a list of all files in the directory
for filename in os.listdir(directory_path):
    # Get the full path of the file
    full_file_path = os.path.join(directory_path, filename)

    # Create the new filename with .fasta extension
    new_filepath = full_file_path + ".fasta"
        
    # Rename the file
    os.rename(full_file_path, new_filepath)
    print(f"Renamed: {filename} -> {new_filepath}")

