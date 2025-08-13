# Global imports
import os
import argparse
from Bio import SeqIO

# Script information - Written in Python 3.9.12 - May 2023
__author__ = "Guillermo Carrillo Martin & Ricardo Fong Zazueta"
__maintainer__ = "Guillermo Carrillo Martin"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script searches for duplicate and fragmentary sequences between records 
in a multi-fasta and removes that redundancy. It writes two outputs, one with 
the filtered multi-fasta (without redundancy) and another multi-fasta with the 
records removed. The way it detects redundancy is based on the sequence (DNA or 
protein), and it can occur in two different situations:

1. Both records have the same sequence. In this case, the removed record
    is the first to be found in the multi-fasta. For example, among the following 
    records, 'A' will be removed from the multi-fasta:

    A. MPIVCLGLLVFGLT
    B. MPIVCLGLLVFGLT

2. One record has a sequence that is a substring from another. In this case, the
    removed record is the one containing the substring, in other words, the one with 
    the shortest sequence length. For example, among the following records, 'B'
    will be removed from the multi-fasta:

    A. MPIVCLGLLVFGLT
    B.   IVCLGL
"""

def main():

    print("# Removing redundant records")
    input_file_path, output_path, output_folder_name = parser()

    # Create an output folder to store the results
    if not os.path.exists(f"{output_path}/{output_folder_name}"):
        os.mkdir(f"{output_path}/{output_folder_name}")

    # Import the multi-fasta file as a list of records
    fasta_record_list = list(SeqIO.parse(input_file_path, "fasta"))

    # Remove duplicate and substring records
    no_duplicate_list, duplicate_list, dup_count = remove_duplicate_records(fasta_record_list)
    no_redundant_list, redundant_list, substring_count = remove_substring_records(no_duplicate_list, duplicate_list)

    # Write the non-redundand database and the (removed) redundant records
    with open(f"{output_path}/{output_folder_name}/filtered_database.fasta", "w") as output_fasta:
        SeqIO.write(no_redundant_list, output_fasta, "fasta")
    with open(f"{output_path}/{output_folder_name}/redundant_records.fasta", "w") as output_fasta:
        SeqIO.write(redundant_list, output_fasta, "fasta")

    # Print the number of duplicate and substring records removed
    print(f"   {dup_count} records with the same sequence removed")
    print(f"   {substring_count} fragment records removed")
    
def parser():
    """
    This function parses the required arguments from the terminal to the python script.
    
    #OUTPUT
    - input_file_path (string); The path to the input multi-fasta database
    - output_path (string); The directory path to write the results folder.
    - output_folder_name (string); The name of the folder to store the results.
    """
    # Set the arguments to run the program from the command line (parser)
    parser = argparse.ArgumentParser(description="This script removes duplicate and fragmentary sequences between records in a multi-fasta")
    parser.add_argument("--input-path", dest="input_path", type=str, help="The path to the input multi-fasta file", required=True, nargs=1)
    parser.add_argument("--output-path", dest="output_path", type=str, help="The directory path to write the results folder (default: working directory)", required=False, default=["."],  nargs=1)
    parser.add_argument("--output-folder-name", dest="output_folder_name", type=str, help="The name of the folder to store the results (default: fasta_remove_redundancy)", required=False, default=["fasta_remove_redundancy"], nargs=1)

    args = parser.parse_args()

    input_file_path = os.path.realpath(args.input_path[0])
    output_path = os.path.realpath(args.output_path[0])
    output_folder_name = args.output_folder_name[0]

    return input_file_path, output_path, output_folder_name

def remove_duplicate_records(fasta_record_list):
    """ 
    This function keeps one copy of each duplicated record in a multi-fasta file. Two records are
    considered as duplicated if both of them have the exact same sequence.
    
    #INPUT
    - fasta_record_list (list); A nested list with all the multi-fasta records, parsed 
      by the SeqIO module.
    #OUTPUT
    - no_duplicate_list (list); A nested list, in SeqIO format, without exact duplicate records.
    - duplicate_list (list); A nested list, in SeqIO format, containing the removed records.
    - dup_count (integer); Number of exact duplicates removed.
    """
    no_duplicate_dic = {}
    duplicate_list = []
    dup_count = 0

    # Store each sequence into dictionary. Only once per each different sequence
    for record in fasta_record_list:

        if str(record.seq) not in no_duplicate_dic:
            no_duplicate_dic[str(record.seq)] = record
        
        else:
            dup_count += 1
            duplicate_list.append(record)

    # Turn the dictionary values into a list
    no_duplicate_list = list(no_duplicate_dic.values())

    return no_duplicate_list, duplicate_list, dup_count

def remove_substring_records(no_duplicate_list, duplicate_list):
    """
    This function removes all records whose sequence is a substring of another record's sequence.
    To do so, sequences are sorted by length, and compared to the ones in the list having more
    length. This way, we avoid extra comparisons, as a protein with higher length can not be a 
    substring of a sorter protein.
    
    #INPUT
    - no_duplicate_list (list); A nested list, in SeqIO format, without exact duplicate records.
    - duplicate_list (list); A nested list, in SeqIO format, containing the removed duplicate 
      records.
    #OUTPUT
    - no_redundant_list (list); A nested list, in SeqIO format, without exact duplicate or 
      substring records.
    - redundant_list (list); A nested list, in SeqIO format, containing the all the 
      duplicate/substring removed records.
    - substring_count (integer); Number of substring records removed.
    """

    redundant_list = duplicate_list.copy()
    no_redundant_list = []
    substring_count = 0

    # Sort the sequences by length
    no_duplicate_list_sorted = sorted(no_duplicate_list, key=lambda record: len(record.seq))
    sequences_list_sorted  = [str(record.seq) for record in no_duplicate_list_sorted]

    # Subset the substring/non-substring sequences in different lists
    for index, sequence_record in enumerate(sequences_list_sorted):

        # Concatenating all next values and checking for existence
        if str(sequence_record) not in ', '.join(sequences_list_sorted[index + 1:]):
            no_redundant_list.append(no_duplicate_list_sorted[index])
        
        else:
            redundant_list.append(no_duplicate_list_sorted[index])
            substring_count += 1

    # Reorder the records based on the previous database order
    original_record_order = {str(record): index for index, record in enumerate(no_duplicate_list)}
    no_redundant_list = sorted(no_redundant_list, key=lambda x: original_record_order[str(x)])

    return no_redundant_list, redundant_list, substring_count

main()