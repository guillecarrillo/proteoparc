# Global imports
import re
import os
import argparse

# Script information - Written in Python 3.9.12 - May 2023
__author__ = "Guillermo Carrillo Martín & Ricardo Fong Zazueta"
__maintainer__ = "Guillermo Carrillo Martín"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script searches for duplicate sequences between records in a multi-fasta 
and removes that redundancy. It creates two outputs, one with the filtered 
multi-fasta (without redundancy) and another multi-fasta with the records removed. 
The way it detects redundancy is based on the sequence (DNA or protein), and 
it can occur in two different situations:

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
    B. IVCLGL

------------------------------------------------------------------------------------    

The script has three main sections:

1. Append the multi-fasta content into a Python list. To do so, firstly the 
    multi-fasta file is converted into another multi-fasta but with the sequences
    of each record collapsed in one row. Then, that new multi-fasta is iterated to 
    append the content into a Pyhton list. 
2. Remove the redundant records from the list and write the outputs. 
3. Print the number of removed records. Specifically, it prints in the terminal
    the number of exact duplicated records removed and the number of substrings
    records removed.
"""

#0. Set the arguments to run the program from the command line (parser)
parser = argparse.ArgumentParser(description="A script that removes two records with redundant sequences in a FASTA file.")
parser.add_argument("--input-path", dest="input_path", type=str, help="the path to the input multi-fasta file", required=True, nargs=1)
parser.add_argument("--output-folder", dest="folder_name", type=str, help="the name of the output folder", required=True, nargs=1)
parser.add_argument("--output-name", dest="multi_fasta", type=str, help="the name of the filtered multi-fasta (default='filtered_mf.fasta')", required=False, nargs=1, default="/filtered_mf.fasta")

args = parser.parse_args()

input_file_path = args.input_path[0]
output_folder_path = args.folder_name[0]
output_mf_name = args.multi_fasta[0]

#1. Append the multi-fasta content into a Python list
#1.1 Remove breaks ('\n') at the end of the sequence lines
def fasta_oneliner(input_path, output_path): 
    fasta_file = open(input_path, "rt") 
    output = open(output_path, "wt")

    file_line = fasta_file.readline() #Read and wrie the first header in the output multi-fasta
    output.write(file_line)

    for line in fasta_file: #For each remaining row...
        
        match = re.search(r"^>",line) #If the row is a header...
        if match:
            output.write("\n" + line) #Write a page break and then the header
        
        elif not match: #If the line is part of a sequence
            output.write(line.replace("\n","")) #Write the sequence without the page break

    fasta_file.close()
    output.close()

fasta_oneliner(input_file_path, "./oneliner.fasta")

#1.2 Append the records to a nested list: [[header1, sequence1], [header2, sequence2], ... ]
record_list = []

with open("./oneliner.fasta", "rt") as multifasta_file: 
   
    while True: #Iterate each row in the oneliner file
        tag = multifasta_file.readline().strip() #Header
        
        if not tag: 
            break #If there are no more headers, finish the loop
        
        sequence = multifasta_file.readline().strip() #Sequence
        
        record_list.append([tag, sequence]) #Append the header and the sequence

os.remove("./oneliner.fasta") #Remove the temporal sequence-oneliner file

#2. Remove the redundant records from the list and write the outputs
#Create a new list to remove the sequences
record_list_removed = record_list.copy()

#Create variables to count the number of removed records
records_equal_count = 0 
records_contained_count = 0

"""
The comparison is made between a record and the following ones to avoid 
repetitive comparations.
"""

#2.1 Remove the duplicates in the new list
print("# Removing duplicate records")

for record1 in record_list_removed: #For each record...
    for record2 in record_list_removed[record_list_removed.index(record1) + 1:]: #And for each following record...
        
        if record1[1] == record2[1]: #Compares if they are equal
            record_list_removed.remove(record1) #If so, remove the first record and move to the next record
            records_equal_count += 1
            break

        elif record1[1] in record2[1]: #Compares if the first record is a substring of the second one
            record_list_removed.remove(record1) #If so, remove the substring record and move to the next record
            records_contained_count += 1
            break

        elif record2[1] in record1[1]: #Compares if the second record is a substring of the first one
            record_list_removed.remove(record2) #If so, remove the substring record and continue with the comparison
            records_contained_count += 1

#2.2 Write the records in a new file
if not os.path.exists(output_folder_path): #Create the output folder if it does not exist
    os.mkdir(output_folder_path)

#2.2.1 Write in a new file the filtered multi-fasta
with open(output_folder_path + "/" + output_mf_name, "wt") as output:
    for record in record_list_removed: #For each record...

        output.write(record[0] + "\n") #Write the header

        sequence_split = re.findall(".{1,60}", record[1]) #Split the sequence into rows of 60 amino acids
        for row in sequence_split:
            output.write(row + "\n")

#2.2.2. Write in a new file the removed records
with open(output_folder_path + "/duplicate_records.fasta", "wt") as duplicate_records_file:
    for record in record_list:
        if record not in record_list_removed: #For each record that is not in the filtered list...
            
            duplicate_records_file.write(record[0] + "\n") #Write the header

            sequence_split = re.findall(".{1,60}", record[1]) #Split the sequence into rows of 60 amino acids
            for row in sequence_split:
                duplicate_records_file.write(row + "\n")

#3. Print the number of removed records
print("   %d records with the same sequence removed" % records_equal_count)
print("   %d fragment records removed" % records_contained_count)