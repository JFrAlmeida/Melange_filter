from Bio import SeqIO

#looks into a pandas df to find a column with the second parameter as name; prints "ITS HERE" if found
def find_my_p_column(pandas_df, colname):
    if colname in list(pandas_df.columns.values):
        print(f"I found a col named {colname}!!!")
        return True
    else:
        return False

#counter incrementation. eats a count and a total. count is the progress, how many things has a loop gone through, total
# is the total of things that exists. returns (updated counter, percentage of progress) as a tuple.
def increment_counter (count, total):
    counter = count + 1
    counter_percent = round((counter / total * 100), 2)
    return (counter, counter_percent)


#Function to extract sequences from FASTA files, returns a list of sequences found with the name seqs_id
def extract_sequences(fasta_file, seqs_id):
    seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return [seq_dict.get(seq_id) for seq_id in seqs_id]  # Return sequence if ID exists

#this uses pandas .str.contains() to check for a regex condition and returns a copy of the passed df selected by condition
def pandas_rows_contains(pandas_df,column, condition):
    finder = pandas_df[column].str.contains(condition, case=False, regex=True, na=False)
    new_df = pandas_df.loc[finder].copy()
    return new_df

#pandas_df is the target df, column is your target column, conditions is a list of any numebrs to be passed to pandas .isin()
def pandas_rows_by_listofnumbers(pandas_df,column, conditions):
    finder = pandas_df[column].isin(conditions)
    new_df = pandas_df.loc[finder].copy()
    return new_df

#Takes a dict, key, separator and a string object to store results; Used in script to parse cazymes annotations
def Annotation_builder(dict, key, string, sep):
    for n in range(0, int(dict[key])):
        if string == "":
            string = string + key
        else:
            string = string + sep + key
    return string

##### Usage of extrac_sequences example:
# genome_selected.loc[:, "dna_sequence"] = extract_sequences(ffn_path, genome_selected.loc[:, "prokka_features"])
#       Where genome_selected is a pandas dataframe, and "prokka_feature" is a column containing the prokka number of
#       features to pull the DNA sequences for; Works for amino acids and for anything that is text

