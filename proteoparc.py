# Global imports
import os
import argparse
import requests
import subprocess

# Script information - Written in Python 3.9.12 - June 2023
__author__ = "Guillermo Carrillo Martin"
__maintainer__ = "Guillermo Carrillo Martin"
__email__ = "guillermo.carrillo@upf.edu"

"""
This software concatenates some scripts to create a reference protein multi-fasta 
database designed to be used in LC-MS/MS protein identification. The performance 
of this pipeline can be split into three different steps:

    1. Download. Creates a multi-fasta database file with the proteins 
    that fulfil the search requirements.

    2. Post-processing. Does optional modifications to the multi-fasta 
    database: 
     - Remove redundant records with exact or substring sequences.
     - Align each protein record by the gene name.

    3. Metadata. Generates some CSV tables and files with metadata information 
    about the database, like the number of species retrieved or the genes not found 
    during the search.
"""

def main():

    # Parse the input variables from the terminal
    RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, do_remove_redundancy, do_align_database = parser()
    script_directory_path = f"{os.path.dirname(__file__)}/scripts"

    # Interrupt the execution if the user is not connected to internet
    if not internet_on():
        print("ERROR: No internet connection detected")
        exit(0)
    
    # DOWNLOAD STEP
    if os.path.exists(RESULTS_FOLDER):
        os.system(f"rm -r {RESULTS_FOLDER}")
    os.mkdir(RESULTS_FOLDER)

    download_proteins(RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, script_directory_path)
    
    # Delete the results folder if no proteins were downloaded
    if not os.path.exists(f"{RESULTS_FOLDER}/{DATABASE_NAME}"):
        os.system(f"rm -r {RESULTS_FOLDER}")
        print("ERROR: NO PROTEINS FOUND")
        exit(0)

    # POST-PROCESSING STEP
    if do_remove_redundancy:
        remove_redundancy(RESULTS_FOLDER, DATABASE_NAME, script_directory_path)
    if do_align_database:
        align_database_per_gene(RESULTS_FOLDER, DATABASE_NAME, script_directory_path)

    # METADATA STEP
    produce_metadata(RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, script_directory_path)
    plot_metadata(RESULTS_FOLDER, script_directory_path)

def parser():
    """
    This function parses the required arguments from the terminal to the python script.
    
    #OUTPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    - DATABASE_NAME (String); The name of the database file and folder.
    - TAX_ID (Integer); The TaxID number employed to construct the multi-fasta database.
    - GENE_LIST (String); The path to the gene list employed to construct the multi-fasta.
      database. If no gene list is specified, the variable is assigned as None.
    - do_remove_redundancy (Boolean); A boolean indicator to indicate if the 
      'remove redundancy' process happens.
    - do_align_database (Boolean); A boolean indicator to indicate if the 
      'align database' process happens.
    """

    # Setting up the parser
    parser = argparse.ArgumentParser(description="A pipeline to generate protein multi-fasta databases using a TaxID")
    parser.add_argument("--project", "-p", dest="project", type=str, help="The name of the project", required=True, nargs=1)
    parser.add_argument("--output-path", dest="output_path", type=str, help="The path to write the result folder (default: working directory)", required=False, default=["."], nargs=1)
    parser.add_argument("--tax-id", "-t", dest="taxid", type=int, help="The TaxID number employed to construct the multi-fasta database", required=True, nargs=1)
    parser.add_argument("--genes", "-g", dest="gene_list", type=str, help="The path to the list of genes (not mandatory)", required=False, nargs=1)
    parser.add_argument("--remove-redundancy", dest="remove_redundancy", action=argparse.BooleanOptionalAction, help="Specify if the redundant records are removed", default=True, required=False)
    parser.add_argument("--align-database", dest="align_database", action=argparse.BooleanOptionalAction, help="Specify if the database is also aligned per protein", default=True, required=False)

    # Recovering the arguments 
    args = parser.parse_args()

    project_name = args.project[0]
    output_path = os.path.realpath(args.output_path[0])

    RESULTS_FOLDER = f"{os.path.realpath(output_path)}/{project_name}"
    DATABASE_NAME = project_name + "_database.fasta"

    TAX_ID = args.taxid[0]
    if args.gene_list:
        GENE_LIST = str(args.gene_list[0])
    else:
        GENE_LIST = None

    do_remove_redundancy = args.remove_redundancy
    do_align_database = args.align_database

    return RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, do_remove_redundancy, do_align_database

def internet_on():
    """
    This function detects if the user's computer is connected to internet
    or not.

    #OUPUT
    - Boolean (True or False).
    """

    try:
        requests.get('http://www.google.com', timeout=5)
        return True
    except requests.ConnectionError as err: 
        return False

def download_proteins(RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, script_directory_path):
    """
    This function creates a multi-fasta protein database using the UniParc archive, a 
    non-redundant repository that contains all the proteins sequenced or predicted in 
    UniProt, NCBI, and other repositories. The search is focused on a specific taxonomic 
    group by a NCBI TaxID and can be restricted to a certain group of genes, indicated by
    a text file. Other specificities, such as the description of the protein header, can 
    be seen in the README.md file
    
    #INPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    - DATABASE_NAME (String); The name of the database file and folder.
    - TAX_ID (Integer); The TaxID number employed to construct the multi-fasta database.
    - GENE_LIST (String); The path to the gene list employed to construct the multi-fasta.
      database. If no gene list is specified, the variable is assigned as None.
    #WRITE OUTPUT
    - {DATABASE_NAME}.fasta; A multi-fasta protein database constructed following the
      search criteria for certain species and genes.
    """

    # Download the proteins
    if GENE_LIST:
        command_line_uniparc = f"python3 -u {script_directory_path}/uniparc_download.py \
                --output-path {RESULTS_FOLDER} \
                --output-name {DATABASE_NAME}\
                --tax-id {TAX_ID} \
                --genes {GENE_LIST}"
        subprocess.run(command_line_uniparc, shell=True)
    
    elif not GENE_LIST:
        command_line_uniparc = f"python3 -u {script_directory_path}/uniparc_download.py \
                --output-path {RESULTS_FOLDER} \
                --output-name {DATABASE_NAME}\
                --tax-id {TAX_ID}"
        subprocess.run(command_line_uniparc, shell=True)

def remove_redundancy(RESULTS_FOLDER, DATABASE_NAME, script_directory_path):
    """
    This function removes duplicate and substring records from the  
    database multi-fasta file, based on the amino acid sequence.

    #INPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    - DATABASE_NAME (String); The name of the database file and folder.
    #WRITE OUTPUT
    - fasta_remove_redundancy/unfiltered_database.fasta; The unfiltered version 
      of the multi-fasta protein database. 
    - fasta_remove_redundancy/redundant_records.fasta; A multi-fasta file containig
      all the removed records (due to redundancy).
    - {DATABASE_NAME}.fasta; A non-redundant version of the multi-fasta protein database
    """
    
    # Remove redundant records
    command_line_remove_redundancy = f"python3 -u {script_directory_path}/remove_redundant_records.py \
            --input-path {RESULTS_FOLDER}/{DATABASE_NAME} \
            --output-path {RESULTS_FOLDER} \
            --output-folder-name fasta_remove_redundancy"
    subprocess.run(command_line_remove_redundancy, shell=True)
    
    # Move the unifiltered database and the redundant records to the 'fasta_remove_redundancy' folder
    os.system(f"mv {RESULTS_FOLDER}/{DATABASE_NAME} {RESULTS_FOLDER}/fasta_remove_redundancy/unfiltered_database.fasta")
    os.system(f"mv {RESULTS_FOLDER}/fasta_remove_redundancy/filtered_database.fasta {RESULTS_FOLDER}/{DATABASE_NAME}")

def align_database_per_gene(RESULTS_FOLDER, DATABASE_NAME, script_directory_path):
    """
    This function generates an aligned multi-fasta file per each different 
    gene present in a multi-fasta. To do so, the header format should indicate 
    the gene name before the string "GN=". 'mafft' v7.525 has been to choosen 
    tool to build each alignment.

    #INPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    - DATABASE_NAME (String); The name of the database file and folder.
    #WRITE OUTPUT
    - aligned_database/{gene_name}_aligned.fasta; An aligned multi-fasta file in 
      mafft format per each gene present in the protein database.
    """

    # Align the database per each gene name
    align_database_command_line = f"python3 -u {script_directory_path}/align_database_per_gene.py \
               --input-path {RESULTS_FOLDER}/{DATABASE_NAME} \
               --output-path {RESULTS_FOLDER} \
               --output-folder-name alignment_per_gene"
    subprocess.run(align_database_command_line, shell=True)

def produce_metadata(RESULTS_FOLDER, DATABASE_NAME, TAX_ID, GENE_LIST, script_directory_path):
    """
    This function generates metadata files with information about a multi-fasta protein 
    database outputed from the uniparc_download.py script. The software only generates 
    the "genes_NOT_retrieved.csv" file if a gene list has been specified. 

    #INPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    - DATABASE_NAME (String); The name of the database file and folder.
    - TAX_ID (Integer); The TaxID number employed to construct the multi-fasta database.
    - GENE_LIST (String); The path to the gene list employed to construct the multi-fasta.
      database. If no gene list is specified, the variable is assigned as None.
    #WRITE OUTPUT
    - metadata/summary.txt; A text file with a summary of all the metadata information 
      retrieved from the multi-fasta database. It also shows the paths to all the other files.
    - metadata/genes_NOT_retrieved.txt; A text file that contains the genes not retrieved. This
      file is only created if the list of genes used to build the multi-fasta is 
      specified.
    - metadata/records_info.csv; A CSV file that contains all the information present in each 
      record header. See the README.md file for a detailed description of the header 
      information.
    - metadata/genes_retrieved.csv; A CSV file that contains the number of genes retrieved.
    - metadata/species_retrieved.csv; A CSV file that contains the number of species retrieved, 
      indicating the scientific name and the TaxID.
    - metadata/repositories_employed.csv; A CSV file that contains the number of the different
      sources (repositories) used to build the multi-fasta. 
    - metadata/species_genes.csv; A CSV file that shows the number of genes retrieved for each 
      species.
    """

    # Generate metadata files, whether there is or not a gene list
    if GENE_LIST:
        metadata_command_line = f"python3 -u {script_directory_path}/metadata_proteoparc.py \
                --input-path {RESULTS_FOLDER}/{DATABASE_NAME} \
                --output-path {RESULTS_FOLDER} \
                --output-folder-name metadata \
                --genes {GENE_LIST} \
                --tax-id {TAX_ID}"
        subprocess.run(metadata_command_line, shell=True)

    elif not GENE_LIST:
        metadata_command_line = f"python3 -u {script_directory_path}/metadata_proteoparc.py \
                --input-path {RESULTS_FOLDER}/{DATABASE_NAME} \
                --output-path {RESULTS_FOLDER} \
                --output-folder-name metadata \
                --tax-id {TAX_ID}"
        subprocess.run(metadata_command_line, shell=True)

def plot_metadata(RESULTS_FOLDER, script_directory_path):
    """
    This function plots the frequency and presence of the species and genes
    retrieved in the multi-fasta database.

    #INPUT
    - RESULTS_FOLDER (String); The folder's absolute path to write the output.
    #WRITE OUTPUT
    - metadata/plots/species_per_gene_barplot.jpg; A barplot representing the
      number of diferent protein records per each specie and gene name.
    - metadata/plots/species_per_gene_grid.jpg; A grid plot representing the
      presence or absence of each gene name per species.
    """

    # Build barplot
    barplot_command_line = f"Rscript {script_directory_path}/proteoparc_barplot.R \
              {RESULTS_FOLDER}/metadata/species_genes.csv"
    subprocess.run(barplot_command_line, shell=True)
    
    # Build grid
    grid_command_line = f"Rscript {script_directory_path}/proteoparc_grid.R \
              {RESULTS_FOLDER}/metadata/species_genes.csv"
    subprocess.run(grid_command_line, shell=True)
    
    # Generate plots folder
    if not os.path.exists(f"{RESULTS_FOLDER}/metadata/plots"):
        os.mkdir(f"{RESULTS_FOLDER}/metadata/plots")

    os.system(f"mv {RESULTS_FOLDER}/metadata/*.jpg {RESULTS_FOLDER}/metadata/plots")

main()