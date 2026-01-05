import os
import re
import shutil
import time



#This script is meant to rename genomes (fasta) batch-downloaded from NCBI. Since the naming is a bit shite, this
# looks into the first line of the genome, and so the name of the first contig and renames the genome into the format:
# GCFnr_Name_Of_the_isolate_sp._strain.termination, copying the termination from the original file.

#This script does NOT delete files, only copies them





time_start = time.time()

#Add fucntions here
def get_first_line (gpath):
    with open(gpath) as f:
        first_line = f.readline().strip('\n')
        return first_line




#Place this script in the same folder as the "ncbi_dataset" folder obtained after extraction of the zip file downloaded from ncbi

#Default folder where each genome-containing-folder is located at for NCBI
path_to_data = "ncbi_dataset/data/"
secondary_path_to_data = "ncbi_dataset/ncbi_dataset/data/"

#folder where you want to output your genomes to
destination_path = "Genome_fasta_dump/"


#make the folder to output the genomes to
if os.path.exists(destination_path): #cleans up folder if it has a bunch of trash in it
    shutil.rmtree(destination_path)
    os.makedirs(destination_path)
else:
    os.makedirs(destination_path)



#This block finds the path to every genome and makes lists with it
if os.path.isdir(path_to_data):
    list_temp = os.listdir(path_to_data) #lists everything in data
    right_path = path_to_data
elif os.path.isdir(secondary_path_to_data): # some extractions put the folder in this structure
    list_temp = os.listdir(secondary_path_to_data) #lists everything in data
    right_path = secondary_path_to_data
else:
    print("path to your data could not be found, make sure it is present!n/ quitting...")
    quit()

list_GCA = [file
            for file in list_temp
            if re.match(".*GCA*", file)] #stores only the GCA numbers of the genomes

path_to_gcanr = [right_path + file
                    for file in list_temp
                    if re.match(".*GCA*", file)] #stores the path up to /data/GCAnumber

original_full_path = [thing + "/" + os.listdir(thing)[0]
                   for thing in path_to_gcanr] #now has the full path to each genome fasta file





#This block makes the new genome name, moves and renames each fasta file to the new name into Genome_fasta_dump
x = 0
for item in original_full_path:
    first_line = get_first_line(item)
    match = re.match("^[^ ]+\s+(.+?)(?=\s+\S+,)", first_line).group(1)
    match = match.replace(" ","_")
    match = match.replace(":", "")

    #Determine file extension
    file_ext = re.search("\.[^.]+$", item).group(0)

    #builds the path for the correctly named genome
    final_path = destination_path + list_GCA[x] + "_" + match + file_ext


    #moves the genome to the new path (renaming it in the process)
    shutil.copy2(item, final_path)
    x = x + 1

#nice goodbye message
print("Thank you for using my scripts, hope it helped (｡◕‿◕｡) -- JFA")

#Runtime calculation
time_finished = time.time()
print(os.path.basename(__file__), " took ", (time_finished - time_start), " seconds to run")



# the naming convention appears to be ENA | Calgstuff | more calgstuff.version bacteria_taxonomy assembly_name status_of_assembly, conting_nr


## Examples
#>ENA|CALGLV010000001|CALGLV010000001.1 uncultured Aquimarina sp. isolate L2_Bin_MAXBIN__015_sub_1 genome assembly, contig: c_000003857977
#>ENA|OMKF01000001|OMKF01000001.1 Aquimarina sp. Aq78 isolate Aquimarina sp. strain Aq78 assembly 01 genome assembly, contig: contig_1




