# Global imports
import os
import argparse
import subprocess
from Bio import SeqIO

# Script information - Written in Python 3.9.12 - May 2024
__author__ = "Guillermo Carrillo Martin"
__maintainer__ = "Guillermo Carrillo Martin"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script generates an aligned multi-fasta file per each different gene name present
in a multi-fasta. To do so, the header format should indicate the gene name before 
the string "GN=". mafft v7.525 has been to choosen to build each alignment, under the
--auto parameter.
"""

def main():

    print("# Aligning database per gene name")
    fasta_real_path, output_folder_realpath = parser()

    if not os.path.exists(output_folder_realpath):
        os.mkdir(output_folder_realpath)

    # Split the multi-fasta into different files based on gene names
    per_gene_fasta_path_list = split_fasta_per_gene(fasta_real_path, output_folder_realpath)

    # Align each different gene-multi-fasta
    for gene_fasta_path in per_gene_fasta_path_list:
        multi_fasta_aligner(gene_fasta_path)
        os.remove(gene_fasta_path)

def parser():
    """
    This function parses the required arguments from the terminal to the python script.
    
    #OUTPUT
    - fasta_real_path (string); The path to the input multi-fasta file.
    - output_folder_path (string); The path to an existing folder where the alignments will be stored.
    """
    parser = argparse.ArgumentParser(description="A script to align a multi-fasta file per each annotated gene")
    parser.add_argument("--input-path", dest="input_path", type=str, help="The path to the input multi-fasta", required=True, nargs=1)
    parser.add_argument("--output-path", dest="output_path", type=str, help="The path to write the folder storing the alignments (default: working directory)", required=False, default=["."], nargs=1)
    parser.add_argument("--output-folder-name", dest="output_folder_name", type=str, help="The name of the folder storing the alignments (default: alignment_per_gene)", required=False, default=["alignment_per_gene"], nargs=1)

    args = parser.parse_args()

    fasta_real_path = os.path.realpath(args.input_path[0])
    output_path = os.path.realpath(args.output_path[0])
    output_folder_name = args.output_folder_name[0]

    output_folder_realpath = f"{output_path}/{output_folder_name}"

    return fasta_real_path, output_folder_realpath

def split_fasta_per_gene(fasta_real_path, output_folder_realpath):
    """
    This function splits a multi-fasta into different gene multi-fasta (temporary) files.
    To do so, the header format should indicate the gene name before the string "GN=".
    
    #INPUT
    - fasta_real_path (string); The path to the input multi-fasta file.
    - output_folder_path (string); The path to an existing folder where the alignments will be stored.
    #OUTPUT
    - per_gene_fasta_path_list (list); A list containing the path of each gene multi-fasta.
    #WRITE OUTPUT
    - {gene}.temp; A multi-fasta file with all the records belonging to a certain gene.
    """

    gene_records_dic = {}
    per_gene_fasta_path_list = []

    # Read each record from the input multi-fasta file
    for record in SeqIO.parse(fasta_real_path, "fasta"):
        
        # Skip the record if it does not have a gene name
        if not "GN=" in record.description:
            continue

        # Extract gene name
        gene_name = record.description.split("GN=")[1].split()[0]
        
        # Sort the record in a dictionary. The gene name of the header acts as the key
        if gene_name not in gene_records_dic:
            gene_records_dic[gene_name] = []
        gene_records_dic[gene_name].append(record)

    # Write in a multi-fasta the sorted records per each gene name and store the file path into a list
    for gene_name, record_list in gene_records_dic.items():
        per_gene_fasta_path_list.append(f"{output_folder_realpath}/{gene_name}.temp")
        SeqIO.write(record_list, f"{output_folder_realpath}/{gene_name}.temp", "fasta")

    return per_gene_fasta_path_list

def multi_fasta_aligner(multi_fasta_path):
    """
    This function aligns a multi-fasta file by the mafft software. It has been tested using
    mafft v7.525.
    
    #INPUT
    - multi_fasta_path (string); The path to the multi-fasta file.
    #WRITE OUTPUT
    - {gene}_aligned.fasta; An aligned multi-fasta file in mafft alignment format.
    """

    alignment_path = multi_fasta_path.replace(".temp", "_aligned.fasta")
    subprocess.run(f"mafft --auto --quiet {multi_fasta_path} > {alignment_path}", shell=True)

main()