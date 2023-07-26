# Global imports
import json
import re
import os
import requests
import argparse
from requests.adapters import HTTPAdapter, Retry

# Script information - Written in Python 3.9.12 - May 2023
__author__ = "Guillermo Carrillo Martín"
__maintainer__ = "Guillermo Carrillo Martín"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script creates a protein database (multi-fasta) using the UniParc 
archive, a non-redundant database that contains all the proteins sequenced 
or predicted in UniProt, NCBI, and other resources. The search is focused 
on a specific taxonomic group (indicated using the NCBI tax ID) and can be 
reduced to a certain group of genes, indicated in text file. Other specificities,
such as the description of the protein header, can be seen in the README.md file.

The script has four main sections:

1. Importing the genes from the gene_list file to a Python list (if a path to
    the gene list is specified).
2. Recovering information about the TaxID using the Uniprot API. Specifically, 
    the clade name that corresponds to the TaxID number and a list of the 
    descendent clades (e.g. Homo sapiens [TaxID=9606] is a descendent clade of 
    Hominidae [TaxID=9604]). 
3. Defining a function that takes a JSON object and turns it into FASTA format.
    This function is designed specifically to retrieve information from the UniParc 
    JSON objects and also creating a specific FASTA header (see README.md file).
4. Using the UniParc API to retrieve the proteins from that database. The protein 
    records are downloaded in JSON format and then converted into FASTA using the
    function described in section 3. If a gene list is specified, the search is 
    iterated for each gene.
"""

#0. Set the arguments to run the program from the command line (parser)
parser = argparse.ArgumentParser(description="A script that creates some metadata information from a database build using uniparc_db.py")
parser.add_argument("--output-name", dest="database", type=str, help="the path to the database output file", required=True, nargs=1)
parser.add_argument("--tax-id", dest="TaxID", type=int, help="the TaxID used to search for proteins", required=True, nargs=1)
parser.add_argument("--genes", dest="gene_list", type=str, help="the path to the list of genes (not required)", required=False, nargs=1)

args = parser.parse_args()

db_name = args.database[0]
tax_id = args.TaxID[0]
gene_list_path = str(args.gene_list[0])

#1. Create the gene list (only if the path is indicated)
if gene_list_path != "None":
    gene_list = []

    with open(gene_list_path, "rt") as gene_list_file:
        for gene in gene_list_file:
            gene_list.append(gene.strip().replace("\n", "")) #Append the genes in the file to a list

#2. Recover TaxID information
"""
In UniParc, all the proteins with the same sequence have the same unique
identifier even though they belong to different species. So, in the end, if a 
protein is conserved through many species it appears in the database as
one record (and the proper unique identifier) associated with many species.
Because of this, sometimes the query retrieves a random specie instead of one
corresponding to our lineage. To fix this error, we create a list with all the 
descendent's clade's TaxID to evaluate if the association belongs or not to our 
TaxID. If it founds a mismatch, there will be the following modifications in the 
header:

OS= TaxID corresponding name (for example, Proboscidea)
OX= TaxID number (for example, 9779)
"""

#2.1 Build API link for downloading in UniProt the list of taxonomic descendent species contained in the TaxID
url_tax_id_descendent = "https://rest.uniprot.org/taxonomy/stream?format=list&query=%28%28ancestor%3A{}%29%29".format(str(tax_id))

tax_id_download = str(requests.get(url_tax_id_descendent).text) #Downloads a string with all the TaxID separated by "\n"
tax_id_list = tax_id_download.split("\n") 
tax_id_list.pop() #Remove the last element of the list (is empty)
tax_id_list = [int(taxid) for taxid in tax_id_list] #Convert strings into integers

#2.2 TaxID organism name
url_org_name = "https://rest.uniprot.org/taxonomy/stream?format=json&query=%28%28tax_id%3A{}%29%29".format(str(tax_id))

org_name_str = str(requests.get(url_org_name).text) #Downloads a string with the TaxID information in JSON format
org_name_dic = json.loads(org_name_str) #Turn the string into a dictionary (JSON format)

org_name = org_name_dic["results"][0]["scientificName"] #Retrieve the TaxID corresponding name

#3. Convert JSON format to FASTA
def json_to_fasta(record, gene_name=None):

    #Retrieves all the header information from the JSON object
    header = ">" + record["uniParcCrossReferences"][0]["database"] + "|" + record["uniParcId"] + "|" + \
    record["uniParcCrossReferences"][0]["lastUpdated"]+ " "
    
    if "proteinName" in record["uniParcCrossReferences"][0]:
        header += record["uniParcCrossReferences"][0]["proteinName"] + " "

    for database in record["uniParcCrossReferences"]:

        if "organism" in database:          
            if database["organism"]["taxonId"] not in tax_id_list:
                header += "OS=" + org_name + " " + \
                "OX=" + str(tax_id) + " "
            else:
                header += "OS=" + database["organism"]["scientificName"] + " " + \
                "OX=" + str(database["organism"]["taxonId"]) + " "
            
            break
                    
    if gene_name != None:
        header += "GN=" + gene_name
    elif gene_name == None:
        if "geneName" in record["uniParcCrossReferences"][0]:
            header += "GN=" + record["uniParcCrossReferences"][0]["geneName"]

    #Retrieves the sequence from the JSON object
    sequence = record["sequence"]["value"]

    return (header, sequence) #Return the header in FASTA format and the protein sequence

#4. Programmatic acces to UniParc through the API
"""
The download process through the API is done by pagination. This is a
special type of download that retrieves records in batches of 500, allowing
for faster and more efficient downloading.  
"""

if os.path.exists(db_name): #Remove the previous database with the same name to avoid problems
    os.remove(db_name)

#4.1 Programmatic access setup via pagination (see http://uniprot.org/help/programmatic_access for the code details)
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

#4.2.1 If there is a list of genes, iterate the download through it and create the database FASTA file to append the protein records
progress = 0

if gene_list_path != "None":
    print("# Downloading the proteins indicated in '%s'" % gene_list_path)

    for gene in gene_list:

        #Build the URL to download the proper proteins via API
        url = "https://rest.uniprot.org/uniparc/search?compressed=false&format=json&query=%28%28gene%3A{}%29%20AND%20%28taxonomy_id%3A{}%29%29&size=500".format(gene, str(tax_id))

        with open(db_name, "at") as file:
            for batch, total in get_batch(url): #For each search...
                if int(total) > 0: #If there is more than one record found
                    
                    json_retrieve = batch.json() #Interprets the download as a JSON object

                    for protein in json_retrieve["results"]: #For each protein record...
                        
                        retrieve = json_to_fasta(protein, gene) #Turn the record into a FASTA
                        
                        header = retrieve[0].replace(',', '').replace(';', '').replace(':', '')
                        file.write(header + "\n")
                        
                        sequence_split = re.findall(".{1,60}", retrieve[1]) #Split the sequence into rows of 60 amino acids
                        for row in sequence_split:
                            file.write(row + "\n")
                    
                    progress += int(total)
    
    if progress == 0:
        os.remove(db_name)
        print("0 proteins found")
        exit(0)

    print("   %d proteins downloaded" % progress)

#4.2.2 If there is no list of genes, download the full proteome and create the database FASTA file to append the protein records
if gene_list_path == "None":
    print("# Downloading the full proteome")

    #Build the URL to download the proper proteins via API
    url = "https://rest.uniprot.org/uniparc/search?compressed=false&format=json&query=%28%28taxonomy_id%3A{}%29%29&size=500".format(str(tax_id))

    with open(db_name, "wt") as file:
        for batch, total in get_batch(url): #For each search...
            
            json_retrieve = batch.json() #Interprets the download as a JSON object
            
            for protein in json_retrieve["results"]: #For each protein record...
                
                retrieve = json_to_fasta(protein) #Turn the record into a FASTA
                
                header = retrieve[0].replace(',', '').replace(';', '').replace(':', '')
                file.write(header + "\n")
                        
                sequence_split = re.findall(".{1,60}", retrieve[1]) #Split the sequence into rows of 60 amino acids
                for row in sequence_split:
                    file.write(row + "\n")
            
            progress += len(json_retrieve["results"])
            print(f'   {progress} / {total}')
    
    print("   %d proteins downloaded" % progress)