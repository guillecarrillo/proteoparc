# Global imports
import re
import os
import argparse
import requests
import json
from requests.adapters import HTTPAdapter, Retry
from concurrent.futures import as_completed
from requests_futures.sessions import FuturesSession
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
peptides have a synthetic origin. Additionally, 3 JSON files are ouput to store the extra 
metadata of the repositories, species, and TaxID of each record in the database. The 
objective of these files is to store the cases of records associated with more than one
repository, species, or TaxID.

Other specificities, such as the description of the protein header, can 
be seen in the README.md file.
"""

def main():
	
	output_path, output_name, tax_id, gene_list_path, gene_list = parser()

	# Download the protein records IDs (UPI) and store them in a list
	upi_id_list = []

	if gene_list:
		upi_gene_dic = {} # Create a dictionary to store the gene name for each UPI ID

		print(f"# Downloading the proteins indicated in '{gene_list_path}'")
		for gene_name in gene_list:
			# Download the UPI IDs for each gene name and assign the gene name in a dictionary
			upi_id_batch_list = api_get_uniparc_record_id_list(tax_id, gene_name)
			upi_gene_batch_dic = {key: gene_name for key in upi_id_batch_list}

			upi_id_list.extend(upi_id_batch_list)
			upi_gene_dic.update(upi_gene_batch_dic)

	elif not gene_list:
		print("# Downloading the whole proteome")
		upi_id_list = api_get_uniparc_record_id_list(tax_id)
		upi_gene_dic = None

	# Download each JSON record based on the UPI ID list
	json_record_list, records_total_count = api_get_json_record_list(upi_id_list)

	if records_total_count == 0:
		print("EXIT: No proteins found with the set conditions")
		exit(0)

	# Generate a list with all the descendents taxid clades of the input taxid to filter the json records
	tax_id_descendent_list = api_get_taxid_descendent_list(tax_id)

	# Turn the JSON records into a multi-fasta file
	records_fasta_list, extra_metadata_repos_dic, extra_metadata_species_dic, extra_metadata_taxid_dic = json_to_fasta(json_record_list, tax_id_descendent_list, upi_gene_dic)

	# Write the multi-fasta file and print the number of proteins downloaded
	with open(f"{output_path}/{output_name}", "w") as output_fasta_file:
		SeqIO.write(records_fasta_list, output_fasta_file, "fasta")

	# Write the extra metadata files
	with open(f"{output_path}/.repos_metadata.json", "w") as output_repos_metadata_file:
		json.dump(extra_metadata_repos_dic, output_repos_metadata_file, indent=4)

	with open(f"{output_path}/.species_metadata.json", "w") as output_species_metadata_file:
		json.dump(extra_metadata_species_dic, output_species_metadata_file, indent=4)

	with open(f"{output_path}/.taxid_metadata.json", "w") as output_taxid_metadata_file:
		json.dump(extra_metadata_taxid_dic, output_taxid_metadata_file, indent=4)

def parser():
	"""
	This function parses the required arguments from the terminal to the python script.
	
	#OUTPUT
	- output_path (string); The absolute folder path to write the multi-fasta database.
	- output_name (string); The name of the multi-fasta database.
	- tax_id (integer); The TaxID number employed to construct the multi-fasta database.
	- gene_list_path (string); The absolute path to the gene list employed to build the multi-fasta.
		database. If no gene list is specified, the variable is assigned as None.
	- gene_list (list); A list with all the gene names present in the input gene list file.
		If no gene list is specified, the variable is assigned as None.
	"""
	parser = argparse.ArgumentParser(description="This script generates a multi-fasta database from the UniParc archive. The search is focused on a specific taxonomic group by a TaxID and can be restricted to a certain group of genes, indicated by a text file")
	parser.add_argument("--output-path", dest="output_path", type=str, help="The folder path to write the multi-fasta database (default: working directory)", required=False, default=["."], nargs=1)
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
			gene_list = [gene.strip() for gene in gene_list_file if gene.strip()]
	elif not gene_list_path:
		gene_list = None

	return output_path, output_name, tax_id, gene_list_path, gene_list

def api_get_uniparc_record_id_list(tax_id, gene_name=None):
	"""
	This function employs the UniProt API to retrieve a list of UniParc IDs (UPI)
	for a specific taxonomic group, defined by a TaxID. If a gene name is provided,
	the function will filter the results to include only proteins associated with 
	that gene.
	#INPUT
	- tax_id (integer); The TaxID number employed to construct the multi-fasta database.
	- gene_name (string); The gene name employed to filter the UniParc records. If no gene name
		is specified, the variable is assigned as None.
	#OUTPUT
	- upi_id_batch_list (list); A list with all the UniParc IDs (UPI) of the proteins
		associated with the input TaxID and gene name.
	"""
	
	upi_id_batch_list = []

	# Set the URL for the UniParc API request
	if gene_name:
		url_record_upi_id = f"https://rest.uniprot.org/uniparc/stream?format=list&query=%28%28gene%3A{gene_name}%29+AND+%28taxonomy_id%3A{str(tax_id)}%29%29&size=500"
	elif not gene_name:
		url_record_upi_id = f"https://rest.uniprot.org/uniparc/search?format=list&query=%28taxonomy_id%3A{str(tax_id)}%29&size=500"

	# Set up the batch API downloader
	re_next_link = re.compile(r'<(.+)>; rel="next"')
	retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
	session = requests.Session()
	session.mount("https://", HTTPAdapter(max_retries=retries))

	# Download in batches the UniParc IDs that fulfill the query conditions
	while url_record_upi_id:
		id_batch_500 = session.get(url_record_upi_id)
		id_batch_500.raise_for_status()

		upi_id_batch_500_list = [upi_id for upi_id in id_batch_500.text.split("\n") if upi_id]
		upi_id_batch_list.extend(upi_id_batch_500_list)

		# If there are more than 500 records, retrieve the link for the next batch of records 
		if "Link" in id_batch_500.headers:
			match = re_next_link.match(id_batch_500.headers["Link"])
			url_record_upi_id = match.group(1)
		elif not "Link" in id_batch_500.headers:
			url_record_upi_id = None

	return upi_id_batch_list

def api_get_json_record_list(upi_id_list):
	"""
	This function employs the UniProt API to download the JSON records of the
	UniParc IDs (UPI) provided in the input list. The function uses parallel requests
	to speed up the download process, and it returns a list of JSON records along with
	the total number of records downloaded.
	#INPUT
	- upi_id_list (list); A list with all the UniParc IDs (UPI) of the proteins.
	#OUTPUT
	- json_records_list (list); A list with all the JSON records of the UniParc IDs.
	- records_total_count (integer); The total number of records downloaded.
	"""

	json_records_list = []
	total_records = len(upi_id_list)
	counter = 0

	print(f"   {total_records} records will be downloaded")

	# Paralelize JSON download
	with FuturesSession() as session:

		# Create the download futures
		download_list = [session.get(f"https://rest.uniprot.org/uniparc/{upi_id}.json") for upi_id in upi_id_list]
		
		# Process the results as they complete
		for download_json in as_completed(download_list):
			json_record = download_json.result().json()
			json_records_list.append(json_record)
			counter += 1

			# Print the download progress
			if counter % 500 == 0:
				print(f"{counter}/{total_records} proteins downloaded", end="\r")

	records_total_count = len(json_records_list)

	return json_records_list, records_total_count

def api_get_taxid_descendent_list(tax_id):
	"""
	This function employs the UniProt API to generate a list with all the descendent
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

def json_to_fasta(json_record_list, tax_id_descendent_list, upi_gene_dic=None):
	"""
	This function parses a list of JSON records from UniParc into a list of
	SeqRecord objects, which can be written to a multi-fasta file. It also collects
	extra metadata about the repositories, species, and TaxID associated with each 
	UniParc ID.
	#INPUT
	- json_record_list (list); A list with all the JSON records of the UniParc IDs.
	- tax_id_descendent_list (list); A list with all the descendent TaxID clades of the
		input TaxID.
	- upi_gene_dic (dictionary); A dictionary with the gene name for each UPI ID. If no gene
		is specified, the variable is assigned as None.
	#OUTPUT
	- records_fasta_list (list); A list with all the SeqRecord objects created from the JSON records.
	- extra_metadata_repos_dic (dictionary); A dictionary with the repositories associated with each UPI ID.
	- extra_metadata_species_dic (dictionary); A dictionary with the species associated with each UPI ID.
	- extra_metadata_taxid_dic (dictionary); A dictionary with the TaxID associated with each UPI ID.
	"""

	extra_metadata_repos_dic = {}
	extra_metadata_species_dic = {}
	extra_metadata_taxid_dic = {}

	records_fasta_list = []
	
	# Iterate for each UniParc record in the JSON file
	for json_protein_record in json_record_list:

		upi_tag = json_protein_record["uniParcId"]
		seen_repos = set()
		seen_species = set()

		upi_repos = []
		upi_species = []
		upi_taxid = []

		for repository_metadata in json_protein_record["uniParcCrossReferences"]:
			
			# Discard incorrect UniParc metadata records
			if not "organism" in repository_metadata:
				continue
			if not repository_metadata["organism"]["taxonId"] in tax_id_descendent_list:
				continue
			
			# Collect the repository name in a unique and ordered way
			if repository_metadata["database"] not in seen_repos:
				seen_repos.add(repository_metadata["database"])
				upi_repos.append(repository_metadata["database"])

			# Collect the species and taxid name in a unique and ordered way
			if repository_metadata["organism"]["scientificName"] not in seen_species:
				seen_species.add(repository_metadata["organism"]["scientificName"])
				upi_species.append(repository_metadata["organism"]["scientificName"])
				upi_taxid.append(repository_metadata["organism"]["taxonId"])
		
		extra_metadata_repos_dic[upi_tag] = upi_repos
		extra_metadata_species_dic[upi_tag] = upi_species
		extra_metadata_taxid_dic[upi_tag] = upi_taxid

		# Iterate for each repository metadata in the record to find the correct metadata
		for repository_metadata in json_protein_record["uniParcCrossReferences"]:
				
			# Discard incorrect UniParc metadata records
			if not "organism" in repository_metadata:
				continue
			if not repository_metadata["organism"]["taxonId"] in tax_id_descendent_list:
				continue
			
			# Retrieve gene name
			if upi_gene_dic:
				gene_name = upi_gene_dic[upi_tag]

			elif (not upi_gene_dic) and ("geneName" in repository_metadata):
				gene_name = repository_metadata["geneName"].upper()

			else:
				gene_name = None

			repository = repository_metadata["database"]
			last_update = repository_metadata["lastUpdated"]
			protein_name = f"{repository_metadata['proteinName']}" if "proteinName" in repository_metadata else None            
			specie = repository_metadata["organism"]["scientificName"]
			taxid = str(repository_metadata["organism"]["taxonId"])
			sequence_version = str(repository_metadata["versionI"])

			# Construct the header
			header = f"{repository}|{upi_tag}"
			header += f"_{gene_name}" if upi_gene_dic else ""
			header += f"|{last_update}"
			header += f" {protein_name}" if protein_name else ""
			header += f" OS={specie} OX={taxid}"
			header += f" GN={gene_name}" if gene_name else ""
			header += f" SV={sequence_version}"

			break

		# Retrieve the sequence
		sequence = json_protein_record["sequence"]["value"]
		
		# Ensemble the fasta record
		protein_record_fasta = SeqRecord(Seq(sequence), id=header, description="")
		records_fasta_list.append(protein_record_fasta)

	return records_fasta_list, extra_metadata_repos_dic, extra_metadata_species_dic, extra_metadata_taxid_dic

main()