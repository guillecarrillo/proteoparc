# Global imports
import os
import re
import shutil
import pandas as pd
import argparse
from datetime import datetime

# Script information - Written in Python 3.9.12 - May 2023
__author__ = "Guillermo Carrillo Martín"
__maintainer__ = "Guillermo Carrillo Martín"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script reads a protein multi-fasta database and creates a folder with metadata 
information. The seven files created by the script are:

1. summary.txt; A text file that contains a summary of all the metadata information
    retrieved from the database. It also shows the paths to all the other files.

2. records_info.csv; A CSV file that contains all the information present in each 
    record header. See the README.md file for a detailed description of the header 
    information.

3. databases_employed.csv; A CSV file that contains the number of the different
    sources (databases) used to build the multi-fasta. 

4. genes_retrieved.csv; A CSV file that contains the number of genes retrieved.

5. genes_NOT_retrieved.txt; A text file that contains the genes not retrieved. This
    file is only created if the list of genes used to build the multi-fasta is 
    specified.

6. species_retrieved.csv; A CSV file that contains the number of species retrieved, 
    indicating the scientific name and the TaxID.

7. species_genes.csv; A CSV file that shows the number of genes retrieved for each 
    species.

----------------------------------------------------------------------------------

The script has two main sections:

1. Calculate all the metadata information and stats from the multi-fasta. Also, 
    creates all the files previously mentioned except the metadata_readme.txt file.

2. Create the metadata_readme.txt file and write there all the information previously
    calculated.
"""

#0. Set the arguments to run the program from the command line (parser)
parser = argparse.ArgumentParser(description="A script that creates metadata information for a database created by uniparc_db.py")
parser.add_argument("--input-path", dest="input_path", type=str, help="the path to the input database", required=True, nargs=1)
parser.add_argument("--output-folder", dest="folder_name", type=str, help="the name of the output folder with the metadata information", required=True, nargs=1)
parser.add_argument("--genes", dest="gene_list", type=str, help="the path to the list of genes (not required)", required=False, nargs=1)

args = parser.parse_args()

database_path = args.input_path[0]
metadata_folder_path = args.folder_name[0]
gene_list_path = str(args.gene_list[0])

#1. Creating stats
#Create the metadata folder
if os.path.exists(metadata_folder_path):
    shutil.rmtree(metadata_folder_path)

os.mkdir(metadata_folder_path)

#Current time and date
now = datetime.now()
dat_time = now.strftime("%d/%m/%Y - %H:%M:%S")

#Create the gene list if the path is specified
if gene_list_path != "None":

    #Append in a list the genes in the gene list
    gene_list = []
    with open(gene_list_path, "rt") as gene_list_file: #Open the gene list file
        for gene in gene_list_file:
            gene_list.append(gene.strip().replace("\n", "")) #Append the genes in the file to a list

#1.1 Stats about the protein records
#1.1.1 Number of records; records with gene name, scientific name, and TaxID
record_num = os.popen(f"grep -E -c '^>' {database_path}").read() 
record_with_gene = os.popen(f"cat {database_path} | grep -E '^>' | grep -c 'GN='").read()
record_with_taxname = os.popen(f"cat {database_path} | grep -E '^>' | grep -c 'OS='").read()
record_with_taxid = os.popen(f"cat {database_path} | grep -E '^>' | grep -c 'OX='").read()

#1.1.2 Create records_info.csv
#a. Create a function that retrieves a header tag value using a regular expression
def retrieve_tag(regular_expression, header):
    tag_match = re.search(r"" + regular_expression, header)
    
    if tag_match:
        tag = tag_match.group(1)
    
    elif not tag_match:
        tag = ""

    return tag

#b. Iterate through the multi-fasta 
with open(database_path, "rt") as database_file:
    records = [] 

    for line in database_file:  
        if line.startswith(">"): #For each header in the multi-fasta
            
            #Retrieve the values of the following parameters using the retrieve tag function
            upi = retrieve_tag("\|([^|]+)\|", line)  #Unic identifier
            database = retrieve_tag("^>([^|]+)", line) #Database
            gene = retrieve_tag("GN=(\w+)", line).upper() #Gene
            species = retrieve_tag("OS=(.*?)\sOX=", line) #Species
            taxid = retrieve_tag("OX=(\w+)", line) #TaxID
            update = retrieve_tag("\|([\d-]+)", line) #Update date

            header_info = "{},{},{},{},{},{}".format(upi,database,gene,species,taxid,update)
            records.append(header_info)

#c. Print all the information in a CSV file
with open(metadata_folder_path + "/records_info.csv", "wt") as records_info_file:
    records_info_file.write("Unic Identifier,Database,Gene,Species,TaxID,Last Update\n")
    for header in records:
        records_info_file.write(header + "\n")
        
#d. Read the record CSV file as a Panda data frame for future calculations
record_df = pd.read_csv(metadata_folder_path + "/records_info.csv")

#1.2 Stats about the employed databases
databases_count_sr = record_df["Database"].value_counts() #A Panda series counting the values of "Databases"
databases_count_sr.to_csv(metadata_folder_path + "/databases_employed.csv", index_label="Database", header=["Count"])

databases_num = databases_count_sr.count() #Count the number of databases

#1.3 Stats about the genes found
#1.3.1 genes retrieved
record_df["Gene"] = record_df["Gene"].str.upper() #Capitalize all genes
gene_count_sr = record_df["Gene"].value_counts() #A Panda series counting the values of "Databases"
gene_count_sr.to_csv(metadata_folder_path + "/genes_retrieved.csv", index_label="Gene", header=["Count"])

genes_num = gene_count_sr.count() #Count the number of genes

# 1.3.2 Genes not retrieved, if we are using a gene list
if gene_list_path != "None":

    #a. Append the genes found in a list
    gene_found_ls = list(record_df["Gene"].unique())

    #b. Compare the genes_found with the gene_list and create genes_not_found_list
    gene_not_found_ls = []
    for gene in gene_list:
        if gene not in gene_found_ls:
            gene_not_found_ls.append(gene)

    #c. Write the genes not retrieved in a text file
    with open(metadata_folder_path + "/genes_NOT_retrieved.txt", "wt") as genes_not_retrieved:
        for gene in gene_not_found_ls:
            genes_not_retrieved.write(gene + "\n")

#1.4 Stats about the species in the database
#1.4.1 Species retrieved
species_count_sr = record_df[["Species", "TaxID"]].value_counts() #A Panda series counting the values of "Databases"
species_count_sr.to_csv(metadata_folder_path + "/species_retrieved.csv", index_label=["Species", "TaxID"], header=["Count"])

species_num = species_count_sr.count() #Count the number of species

#1.4.2 Genes per species

"""
In the gene per species table, the genes that are not found in a species are written
in the table with a 0 in the count column. The set of genes used to fill the table
in this way depends if the database was built using a gene list or not. In the first
case, the genes used to create the combinations of species, TaxID, and genes are the 
ones in the gene list. In the second case, the genes are the ones found in the 
database.
"""

with open(database_path, "rt") as database:
    species_genes_dic = {} #Create a dictionary to count the genes per specie

    for line in database:
        if line.startswith(">"): #For each header in the database
            
            #a. Count the different combinations of specie, TaxID, and gene
            specie_match = re.search(r"OS=(.*?)\sOX=", line)
            taxid_match = re.search(r"OX=(\w+)", line)
            gene_match = re.search(r"GN=(\w+)", line)
            
            if specie_match and taxid_match and gene_match:
                species = specie_match.group(1)
                taxid = taxid_match.group(1)
                gene = gene_match.group(1).upper()

                trio = "{},{},{},".format(species,taxid,gene)
                
                if trio not in species_genes_dic:
                    species_genes_dic[trio] = 1
                elif trio in species_genes_dic:
                    species_genes_dic[trio] += 1
    
    #b. Create a list with all the species + TaxID found in the database
    names_found_ls = list(record_df["Species"].unique())
    tax_id_found_ls = list(record_df["TaxID"].unique())
    
    species_found_ls = []

    for i in range(0,len(names_found_ls)):
        species_found_ls.append("{},{},".format(names_found_ls[i], tax_id_found_ls[i]))

    #c. Fill the table with the combinations not found
    if gene_list_path != "None": #If the database is built using a gene list...

        for species in species_found_ls: 
            for gene in gene_list: #Use the gene list to generate the combinations

                #If a combination is not in the database, set it as 0 counts
                if species + gene + "," not in species_genes_dic:
                    species_genes_dic[species + gene + ","] = "0"  

    if gene_list_path == "None": #If the database is NOT built using a gene list...

        gene_found_ls = list(record_df["Gene"].unique())

        for species in species_found_ls:
            for gene in gene_found_ls: #Use the gene founds to generate the combination
                
                #If a combination is not in the database, set it as 0 counts
                if species + str(gene) + "," not in species_genes_dic:
                    species_genes_dic[species + str(gene)] = ",0"  

# d. Create the CSV file
with open(metadata_folder_path + "/species_genes.csv", "wt") as species_genes:
    species_genes.write("Species,TaxID,Gene,Count\n")
    for header, count in species_genes_dic.items():
        species_genes.write(header + str(count) + "\n")

#e. Sort the CSV file first by the specie name and then numerically by the number of counts
species_genes_table = pd.read_csv(metadata_folder_path + "/species_genes.csv")
species_genes_table.sort_values(["Species", "Count"], axis=0,ascending=[True, False], inplace=True)

species_genes_table.dropna(subset="Gene", inplace=True) #
species_genes_table.to_csv(metadata_folder_path + "/species_genes.csv", index=False)

#2. Write the metadata information in a file
with open(metadata_folder_path + "/summary.txt", "wt") as metadata:

    #Metadata file header
    metadata.write("## METADATA for the database '{}'".format(database_path) + "\n" + \
                   "{}".format(dat_time) + "\n" + "\n")

    #2.1 Warnings
    warnings = 0 
                   
    metadata.write("# WARNINGS\n")
    if record_num != record_with_gene:
        warnings += 1
        metadata.write("Check gene names (GN=)\n")
    
    if record_num != record_with_taxname:
        warnings += 1
        metadata.write("Check taxa names (OS=)\n")

    if record_num != record_with_taxid:
        warnings += 1
        metadata.write("Check taxIDs (OX=)\n")

    if gene_list_path != "None":
        if len(gene_found_ls) != len(gene_list):
            warnings += 1
            metadata.write("{} genes not found\n".format((len(gene_list)) - (len(gene_found_ls))))

    if warnings == 0:
        metadata.write("No warnings\n")

    #2.2 Stats about the protein records
    metadata.write("\n# Stats about the protein records\n" + \
                    "records: {}".format(record_num) + \
                    "records_with_gene: {}/{}".format(record_with_gene[:-1],record_num) + \
                    "records_with_taxname: {}/{}".format(record_with_taxname[:-1],record_num) + \
                    "records_with_taxid: {}/{}".format(record_with_taxid[:-1],record_num) + 
                    "all_records_path: {}\n".format(metadata_folder_path + "/records_info.csv"))

    #2.3 Stats about the employed databases
    metadata.write("\n# Stats about the employed databases\n" + \
                    "databases_employed: {}\n".format(databases_num) + \
                    "databases_employed_path: {}\n\n".format(metadata_folder_path + "/databases_employed.csv"))

    #2.4 Stats about the genes found
    if gene_list_path != "None":            
        metadata.write("# Stats about the genes found\n" + \
                        "genes_retrieved: {}/{}\n".format(genes_num,str(len(gene_list))) + \
                        "genes_retrieved_path: {}\n".format(metadata_folder_path + "/genes_retrieved.csv") + \
                        "genes_not_retrieved_path: {}\n".format(metadata_folder_path + "/genes_NOT_retrieved.txt"))
    
    elif gene_list_path == "None":
        metadata.write("# Stats about the genes found\n" + \
                        "genes_retrieved: {}\n".format(genes_num) + \
                        "genes_retrieved_path: {}\n".format(metadata_folder_path + "/genes_retrieved.csv"))

    #2.5 Stats about the species in the database
    metadata.write("\n# Stats about the species in the database\n" + \
                    "species_retrieved: {}\n".format(species_num) + \
                    "species_retrieved_path: {}\n".format(metadata_folder_path + "/species_retrieved.csv") + \
                    "species_and_genes_path: {}".format(metadata_folder_path + "/species_genes.csv"))                    