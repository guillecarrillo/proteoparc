# Global imports
import re
import os
import argparse
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Script information - Written in Python 3.9.12 - May 2023
__author__ = "Guillermo Carrillo Martin"
__maintainer__ = "Guillermo Carrillo Martin"
__email__ = "guillermo.carrillo@upf.edu"

"""
This script creates a multi-fasta protein database using the UniParc archive, a 
non-redundant archive containing all the proteins sequenced or predicted in 
UniProt, NCBI, and other repositories. The search is focused on a specific taxonomic 
group by a TaxID and can be restricted to a certain group of genes, indicated by
a text file. Proteins from the FusionGDB repository are not being downloaded, as those 
peptides have a synthetic origin.

Other specificities, such as the description of the protein header, can 
be seen in the README.md file.
"""

def main():
    output_path, output_name, tax_id, gene_list_path, gene_list = parser()

    # Generate a list with all the descendents taxid clades of the input taxid
    tax_id_descendent_list = api_taxid_descendent_list(tax_id)

    # Download the protein records and store the results in a SeqRecord list
    if gene_list:
        print(f"# Downloading the proteins indicated in '{gene_list_path}'")
        records_fasta_list, records_total_count = per_gene_protein_downloader(gene_list, tax_id, tax_id_descendent_list)
    elif not gene_list:
        print("# Downloading the whole proteome")
        records_fasta_list, records_total_count = whole_proteome_downloader(tax_id, tax_id_descendent_list)

    # Print an error if there no records were downloaded
    if records_total_count == 0:
        print("EXIT: No proteins found with the set conditions")
        exit(0)

    # Write the multi-fasta file and print the number of proteins downloaded
    with open(f"{output_path}/{output_name}", "w") as output_fasta_file:
        SeqIO.write(records_fasta_list, output_fasta_file, "fasta")
    print(f"   {records_total_count} proteins downloaded")

def parser():
    """
    This function parses the required arguments from the terminal to the python script.
    
    #OUTPUT
    - db_output_path (string); The file's absolute path to output the multi-fasta database.
    - tax_id (integer); The TaxID number employed to construct the multi-fasta database.
    - gene_list_path (string); The path to the gene list employed to construct the multi-fasta.
      database. If no gene list is specified, the variable is assigned as None.
    - gene_list (list); A list with all the gene names present in the input gene list file.
      If no gene list is specified, the variable is assigned as None.
    """
    parser = argparse.ArgumentParser(description="This script creates some metadata information from a database build using uniparc_download.py")
    parser.add_argument("--output-path", dest="output_path", type=str, help="The path to write the multi-fasta database (default: working directory)", required=False, default=["."], nargs=1)
    parser.add_argument("--output-name", dest="output_name", type=str, help="The name of the multi-fasta database (default: database.fasta)", required=False, default=["database.fasta"], nargs=1)
    parser.add_argument("--tax-id", dest="TaxID", type=int, help="The TaxID number employed to construct the multi-fasta database", required=True, nargs=1)
    parser.add_argument("--genes", dest="gene_list", type=str, help="The path to the list of genes (not mandatory)", required=False, default=[None], nargs=1)

    args = parser.parse_args()

    output_path = os.path.realpath(args.output_path[0])
    output_name = args.output_name[0]
    tax_id = args.TaxID[0]
 
    if args.gene_list[0]:
        gene_list_path = str(args.gene_list[0])
    elif not args.gene_list[0]:
        gene_list_path = None

    # Parse the gene list file as a python list
    if gene_list_path:
        gene_list = []
        with open(gene_list_path, "rt") as gene_list_file:
            for gene in gene_list_file:
                if gene == "\n":
                    continue
                gene_list.append(gene.strip().replace("\n", ""))
    
    elif not gene_list_path:
        gene_list = None

    return output_path, output_name, tax_id, gene_list_path, gene_list

def api_taxid_descendent_list(tax_id):
    """
    This function employs UniProt API to generate a list with all the descendent
    TaxID clades of the input TaxID. As UniParc is a non-redundant repositorie, 
    all the proteins with the same sequence have the same 'unique identifier (UPI)'
    even though they belong to different species. So, in the end, if a protein is fully
    conserved through a clade it appears in the database as one record tagged with 
    many species. Because of this, sometimes the query retrieves metadata from a random 
    specie instead of the one corresponding to our desired linage. To fix this error, we 
    create a list with all the descendent's clade's TaxID to evaluate and filter if the 
    retrieved metadata belongs or not to our TaxID linage. (e.g. Homo sapiens [TaxID=9606]
    is a descendent clade of Hominidae [TaxID=9604])
    
    #INPUT
    - tax_id (integer); The TaxID number employed to construct the multi-fasta database.
    #OUTPUT
    - tax_id_descendent_list (list); A list with all the descendent TaxID clades of the 
      input TaxID
    """
    url_tax_id_descendent = f"https://rest.uniprot.org/taxonomy/stream?format=list&query=%28%28ancestor%3A{str(tax_id)}%29%29"

    # Download a string with all the TaxID separated by "\n"
    tax_id_download = str(requests.get(url_tax_id_descendent).text)
    tax_id_descendent_list = tax_id_download.split("\n") 

    # Remove the last element of the list (empty value)
    tax_id_descendent_list.pop()

    # Convert the list values from strings to integers
    tax_id_descendent_list = [int(taxid) for taxid in tax_id_descendent_list]
    tax_id_descendent_list.append(int(tax_id))

    return tax_id_descendent_list

def api_get_json_batch(batch_url):
    """
    This function systematically downloads, as a JSON object, all the protein 
    records found in UniParc under the criteria set by an URL. The download process
    through UniParc API is done by pagination, a feature to retrieve information in 
    batches of 500 records, allowing for faster and more efficient downloading.
    see http://uniprot.org/help/programmatic_access for more code details.
    
    #INPUT
    - batch_url (string); A link indicating the search criteria to focus the protein
      search. It is generated by the criteria set in the UniProt API documentation:
      https://www.uniprot.org/help/api_queries. 
    #OUTPUT
    - records_batch_json (dictionary); A batch of protein records (sequence + metadata) with
      500 entries maximum.
    - record_count (integer); The number of total records downloaded.
    """
    # Set up the API downloader
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    while batch_url:
        records_batch = session.get(batch_url)
        records_batch.raise_for_status()
        records_batch_json = records_batch.json()
        record_count = int(records_batch.headers["x-total-results"])
        
        yield records_batch_json, record_count

        # If there are more than 500 records, retrieve the link for the next batch of records 
        if "Link" in records_batch.headers:
            match = re_next_link.match(records_batch.headers["Link"])
            batch_url = match.group(1)
        elif not "Link" in records_batch.headers:
            batch_url = None

def json_to_fasta(protein_record, tax_id_descendent_list, gene_name=None):
    """
    This function turns a UniParc JSON protein record into fasta format, stored as
    a header and the corresponding amino acid sequence. Check the README.md file 
    to know more about the metadata stored in the header.

    #INPUT
    - protein_record (dictionary) The protein record sequence and metadata stored as
      a JSON object.
    - tax_id_descendent_list (list); A list with all the descendent TaxID clades of the 
      input TaxID.
    - gene_name (string); The gene name assigned to the protein record. If it has not been
      set in advance, this variable has a None boolean value.
    #OUTPUT
    - header (string); The formatted fasta header of the protein record.
    - sequence (string); The amino acid sequence of the protein record.
    """
    # Iterate for each UniParc record in the JSON file to find the correct metadata
    for database in protein_record["uniParcCrossReferences"]:
        
        # Discard incorrect UniParc metadata records
        if not "organism" in database:
            continue
        if not database["organism"]["taxonId"] in tax_id_descendent_list:
            continue

        # Retrieve the sequence metadata from the header
        if (not gene_name) and ("geneName" in database):
            gene_name = database["geneName"].upper()

        repository = database["database"]
        upi_tag = protein_record["uniParcId"] + (f"_{gene_name}" if gene_name else "")
        last_update = database["lastUpdated"]
        protein_name = f"{database['proteinName']}" if "proteinName" in database else None            
        specie = database["organism"]["scientificName"]
        taxid = str(database["organism"]["taxonId"])
        sequence_version = str(database["versionI"])

        # Construct the header
        header = f"{repository}|{upi_tag}|{last_update}"
        header += f" {protein_name}" if protein_name else ""
        header += f" OS={specie} OX={taxid}"
        header += f" GN={gene_name}" if gene_name else ""
        header += f" SV={sequence_version}"

        break

    # Retrieve the sequence
    sequence = protein_record["sequence"]["value"]

    return header, sequence

def per_gene_protein_downloader(gene_list, tax_id, tax_id_descendent_list):
    """
    This function downloads all the protein records in UniParc associated with 
    a clade (tax_id) and a list of genes (gene_list). To do so, an API link is
    generated per each different gene, and then protein records are downloadaded 
    in JSON format to be converted into fasta format. All fasta records are stored 
    in a python list to be easily outputed.

    #INPUT
    - gene_list (list); A list with all the gene names present in the input gene list file.
      If no gene list is specified, the variable is assigned as None.
    - tax_id (integer); The TaxID number employed to construct the multi-fasta 
      database.
    - tax_id_descendent_list (list); A list with all the descendent TaxID clades 
      of the input TaxID.
    #OUTPUT
    - records_fasta_list (list); A list including all the downloaded fasta records 
      in biopython SeqRecord format.
    - records_total_count (integer); Total count of downloaded records.
    """

    records_total_count = 0
    records_fasta_list = []

    for gene_name in gene_list:
        uniparc_api_url = f"https://rest.uniprot.org/uniparc/search?compressed=false&format=json&query=%28%28gene%3A{gene_name}%29%20AND%20%28taxonomy_id%3A{str(tax_id)}%29%20NOT%20%28database%3AFusionGDB%29%29&size=500"
        
        # Download each protein sequence found by the link specifications, in batches of 500 records
        for records_batch_json, record_count in api_get_json_batch(uniparc_api_url):
            
            # Skip the batch if there are no records fulfilling the query requirements
            if record_count == 0:
                continue

            # Append each protein record within the batch to a python list in fasta format
            for protein_record in records_batch_json["results"]:
                
                # Turn the protein JSON record into fasta format
                header, sequence = json_to_fasta(protein_record, tax_id_descendent_list, gene_name)
                header = header.replace(',', '').replace(';', '').replace(':', '')
                
                # Create Biopython fasta record object
                record_fasta = SeqRecord(Seq(sequence), id=header, description="")
                records_fasta_list.append(record_fasta)
        
            # Count the number of records downloaded per batch
            records_total_count += record_count

    return records_fasta_list, records_total_count

def whole_proteome_downloader(tax_id, tax_id_descendent_list):
    """
    This function downloads all the protein records in UniParc associated with 
    a clade (tax_id). To do so, an API link is generated following the API 
    requitements. Then, protein records are downloadaded in JSON format to be 
    converted into fasta format. All fasta records are stored in a python list 
    to be easily outputed.

    #INPUT
    - tax_id (integer); The TaxID number employed to construct the multi-fasta 
      database.
    - tax_id_descendent_list (list); A list with all the descendent TaxID clades 
      of the input TaxID.
    #OUTPUT
    - records_fasta_list (list); A list including all the downloaded fasta records 
      in biopython SeqRecord format.
    - records_total_count (integer); Total count of downloaded records.
    """
    records_total_count = 0
    records_fasta_list = []

    uniparc_api_url = f"https://rest.uniprot.org/uniparc/search?compressed=false&format=json&query=%28%28taxonomy_id%3A{str(tax_id)}%29%20NOT%20%28database%3AFusionGDB%29%29&size=500"
    
    # Download each protein sequence found by the link specifications, in batches of 500 records
    for records_batch_json, record_count in api_get_json_batch(uniparc_api_url):

        # Skip the batch if there are no records fulfilling the query requirements
        if record_count == 0:
            continue
        
        # Append each protein record within the batch to a python list in fasta format
        for protein_record in records_batch_json["results"]:
            
            # Turn the record into a fasta
            header, sequence = json_to_fasta(protein_record, tax_id_descendent_list, gene_name=None)
            header = header.replace(',', '').replace(';', '').replace(':', '')
            
            # Create Biopython fasta record object
            record_fasta = SeqRecord(Seq(sequence), id=header, description="")
            records_fasta_list.append(record_fasta)
        
        # Count the number of records downloaded per batch and print the download progress
        records_total_count += len(records_batch_json["results"])
        print(f"   {records_total_count}/{str(record_count)}")

    return records_fasta_list, records_total_count

main()