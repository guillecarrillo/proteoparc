# Global imports
import os
import argparse
import requests

# Script information - Written in Python 3.9.12 - June 2023
__author__ = "Guillermo Carrillo Martín"
__maintainer__ = "Guillermo Carrillo Martín"
__email__ = "guillermo.carrillo@upf.edu"

"""
This pipeline combines some scripts to create a protein multi-fasta database
meant to be used in mass spectrometry software for protein identification. 

The pipeline is based on three modules:

    1. Download module. Creates a multi-fasta database file with the proteins 
    that fulfil the search requirements.

    2. Post processing module. Does optional modifications to the multi-fasta 
    database, like removing redundant records with exact or substring sequences 
    or removing the word 'isoform' from the record headers.

    3. Metadata module. Generates some tables and files with metadata information 
    about the database, like the number of species retrieved or the genes not found 
    during the search.
"""

def main():
    RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENES, do_remove_duplicates, do_remove_isoform_word = parser()

    if not internet_on():
        print("ERROR: No internet connection")
        exit(0)
    
    #1. Download module
    """
    if there is no database file after the first module (no proteins found), 
    delete the resuls folder and stop the pipeline
    """
    download_proteins(RESULTS_FOLDER, TAX_ID, GENES)
    
    if not os.path.exists(RESULTS_FOLDER + "/unfiltered_database.fasta"):
        os.system(f"rm -r {RESULTS_FOLDER}")
        print("ERROR: NO PROTEINS FOUND")
        exit(0)

    #2. Post processing module
    if do_remove_duplicates:
        remove_duplicates(RESULTS_FOLDER, DATABASE_NAME)   
    elif not remove_duplicates:
        os.system(f"mv {RESULTS_FOLDER}/unfiltered_database.fasta {RESULTS_FOLDER}/{DATABASE_NAME}")

    if do_remove_isoform_word:
        remove_isoform_word(RESULTS_FOLDER, DATABASE_NAME)

    #3. Metadata module
    produce_metadata(RESULTS_FOLDER, DATABASE_NAME, GENES)

#0. Set the arguments to run the program from the command line (parser)
def parser():
    parser = argparse.ArgumentParser(description="A pipeline that creates a protein multi-fasta database using a TaxID")
    parser.add_argument("--project", "-p", dest="project", type=str, help="The name of the project", required=True, nargs=1)
    parser.add_argument("--tax-id", "-t", dest="TaxID", type=int, help="the TaxID used to search for proteins", required=True, nargs=1)
    parser.add_argument("--genes", "-g", dest="gene_list", type=str, help="the path to the list of genes (not required)", required=False, nargs=1)
    parser.add_argument("--path", dest="path", type=str, help="The path to locate the folder with the result (default: ./)", required=False, default="./", nargs=1)
    parser.add_argument("--remove-duplicates", dest="remove_duplicates", action=argparse.BooleanOptionalAction, help="Specify if the duplicate records are removed", default=True, required=False)
    parser.add_argument("--remove-isoform-word", dest="remove_isoform_word", action=argparse.BooleanOptionalAction, help="Specify if the word 'isoform' is removed from the headers", default=True, required=False)

    # Recovering the arguments 
    args = parser.parse_args()

    project = args.project[0]
    path = args.path[0]
    if path[-1] != "/":
        path += "/" 

    RESULTS_FOLDER = path + project
    DATABASE_NAME = project + "_database.fasta"

    TAX_ID = args.TaxID[0]
    if args.gene_list:
        GENES = str(args.gene_list[0])
    else:
        GENES = None

    do_remove_duplicates = args.remove_duplicates
    do_remove_isoform_word = args.remove_isoform_word

    return RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENES, do_remove_duplicates, do_remove_isoform_word

# Check if there is internet connection
def internet_on():
    try:
        requests.get('http://www.google.com', timeout=5)
        return True
    except requests.ConnectionError as err: 
        return False

#1. Download module 

def download_proteins(RESULTS_FOLDER, TAX_ID, GENES):
    # Create the folder to store the database
    if os.path.exists(RESULTS_FOLDER):
        os.system(f"rm -r {RESULTS_FOLDER}")
    os.mkdir(RESULTS_FOLDER)

    # Download the proteins
    os.system(f"python3 scripts/uniparc_db.py \
            --output-name {RESULTS_FOLDER}/unfiltered_database.fasta \
            --tax-id {TAX_ID} \
            --genes {GENES}")

#2. Post processing module 
# Remove duplicate records
def remove_duplicates(RESULTS_FOLDER, DATABASE_NAME):
    os.system(f"python3 scripts/remove_dup_records.py \
            --input-path {RESULTS_FOLDER}/unfiltered_database.fasta \
            --output-folder {RESULTS_FOLDER}/fasta \
            --output-name filtered_database.fasta")
    os.system(f"mv {RESULTS_FOLDER}/unfiltered_database.fasta {RESULTS_FOLDER}/fasta")
    os.system(f"mv {RESULTS_FOLDER}/fasta/filtered_database.fasta {RESULTS_FOLDER}/{DATABASE_NAME}")

# Remove "isoform" word from the header
def remove_isoform_word(RESULTS_FOLDER, DATABASE_NAME):
    os.system(f"sed 's/ isoform//g' {RESULTS_FOLDER}/{DATABASE_NAME} > {RESULTS_FOLDER}/sed.fasta")
    os.system(f"rm {RESULTS_FOLDER}/{DATABASE_NAME}")
    os.system(f"mv {RESULTS_FOLDER}/sed.fasta {RESULTS_FOLDER}/{DATABASE_NAME}")

#3. Metadata module
def produce_metadata(RESULTS_FOLDER, DATABASE_NAME, GENES):
    print("# Creating metadata information")
    os.system(f"python3 scripts/metadata_db.py \
            --input-path {RESULTS_FOLDER}/{DATABASE_NAME} \
            --output-folder {RESULTS_FOLDER}/metadata \
            --genes {GENES}")

main()