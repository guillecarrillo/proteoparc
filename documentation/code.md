# Code overview ProteoParc 1.0

This document explains the code and outputs of all the scripts that compose ProteoParc, a pipeline for the creation of protein databases focused on paleoprotein mass spectrometry identification. This software was entirely written in Python v3.9.12 along with the following packages:

-   numpy v1.25.0
-   pandas v2.0.3
-   requests v2.31.0

## ProteoParc output

### Protein multi-fasta database

The main output of this pipeline is a **multi-fasta protein database**. The header of each record has no commas (,), semicolons (;) or colons(:) to avoid interferences when creating csv files. Also, the architecture of the header is based on the UniProt format but including some changes:

``` textinfo
>Database|Unique_identifier|Last_update [Name_of_the_protein] OS=Specie OX=TaxID [GN=Gene_name] SV=Sequence_version
```

| Field                   | Description                                                                                                                                                                                                   |
|--------------------------|----------------------------------------------|
| **Database**            | Source database of the record. The list of databases that conform UniParc can be consulted [here](https://www.uniprot.org/help/uniparc).                                                                      |
| **Unique Identifier**   | Primary accession number corresponding to the UniParc record.                                                                                                                                                 |
| **Last update**         | Date of the last time the UniParc record was modified. Changes can be due to typos correction or metadata alteration and not just to sequence modification.                                                   |
| **Name of the protein** | Short name to describe the protein. This field can be omitted due to lack of information and poor quality annotations.                                                                                        |
| **Specie**              | Scientific name of the organism associated with the protein.                                                                                                                                                  |
| **TaxID**               | TaxID of the specie                                                                                                                                                                                           |
| **Gene name**           | Abbreviation in capital letters of the gene that codifies the protein. If the download was done without using a gene list, this field can be omitted due to lack of information and poor quality annotations. |
| **Sequence version**    | Number that identifies the version number of the sequence. In contrast to the last update field, this number only changes if the sequence is modified.                                                        |

Below, you can find the first header of the [enamel Proboscidea database](../documentation/example/proboscidea_enamelome/proboscidea_enamelome_database.fasta) as an example:

``` textinfo
>RefSeq|UPI002116250B|2023-03-08 albumin X1 OS=Elephas maximus indicus OX=99487 GN=ALB SV=1
```

As it is explained in the [tutorial](../documentation/tutorial.md), the protein download can be focused on some specific proteins using a gene list to restrict the search. If this option is executed, all the records should have a gene name indicated in the header. Otherwise, if the full proteome of a taxon is downloaded, it can happen that some records would not have an associated gene name due to bad annotated proteins.

### Metadata files

The **metadata module creates seven files** in a metadata folder about database information. This output files will change depending on whether the search was restricted to a list of genes or not. In short, the seven metadata files are:

1.  **summary.txt**; A text file that contains a summary of all the metadata information retrieved from the multi-fasta. It also shows the paths to all the other metadata files and indicates a warning message if there are records with missing fields, like the gene name or the species.

2.  **records_info.csv**; A CSV file that contains all the information present in each record header. See the previous section for a detailed description of the header information.

3.  **databases_employed.csv**; A CSV file that contains the number of the different sources (databases) used to build the multi-fasta.

4.  **genes_retrieved.csv**; A CSV file that shows the genes retrieved and counts the number of records with each gene.

5.  **genes_NOT_retrieved.txt**; A text file that contains the genes not retrieved. This file is only created if the multi-fasta was built using a list of genes.

6.  **species_retrieved.csv**; A CSV file that contains the number of species retrieved, indicating the scientific name and the TaxID. It also counts the number of records with each species. Species and subspecies are count as different taxa in order to calculate the metadata info.

7.  **species_genes.csv**; A CSV file that shows the genes retrieved for each species. In this table, the genes that are not found in a species are set as 0 in the count column. The set of genes used to fill the table in this way depends if the database was built using a gene list or not. In the first case, the genes used to create the combinations of species, TaxID, and genes are the ones in the gene list. In the second case, the genes are the ones found in the database.In other words, when the protein search is done using a gene list, this table counts the abundance of each gene present in the gene list per species. Otherwise, if the search was not restricted, the table counts the abundance of each gene retrieved per species.

It is important to mention that the information present in the metadata files can be underestimated due to the UniParc non-redundant architecture. If two different databases have the exact same sequence, it is only stored once in UniParc, mixing the metadata information. This can also happen with different species, so if a protein is really conserve through evolution and two different species share the exact same sequence, the UniParc record that belongs to that protein will have more than one specie in the metadata field.

At the end, one UniParc record can be associated with more than one database, species, or gene due to metadata mixing. As this pipeline only recovers one result per field in the fasta header -the one that is set in first position in UniParc- some information can be lose in the process. As a recommendation, **take the metadata information as the minimum knowledge of your database**.

### Other optional outputs

Additionally, if the **remove duplicates option is set up**, the pipeline creates a folder named "fasta" that contains all the removed records (duplicate_records.fasta) and the unfiltered database (unfiltered_database.fasta). In this case, the printed output indicates the number of the exact duplicates (same sequence) and partial duplicates (fragment records) removed. For a more detailed explanation of the way the removal is done, check the remove_dup_records.py code description in this document.

Also, if the remove "isoform" word option is set up, only the final database file will have that word removed.

## ProteoParc code

The ProteoParc pipeline consist of four Python scripts. Three of them (uniparc_db.py, remove_dup_records.py, and metadata_db.py) perform specific actions in each module, while the other one (proteoparc.py) works as an script that merges and concatenates this steps to create the usage of the pipeline.

### 1. proteoparc.py

This script combines three previous scripts plus some new code lines to **execute the ProteoParc pipeline** and create a multi-fasta database meant to be used in mass spectrometry software for protein identification. The execution of this script can be seen in the [tutorial](../documentation/tutorial.md).

The performance of this pipeline can be split in **three different modules**:

1.1. **Download module**. Creates the result folder and a multi-fasta database file with the proteins that fulfil the search requirements. Here, the uniparc_db.py script is used, and if there is no database file after the download, because no proteins where found due to erroneous TaxID or gene list, the result folder is deleted and the execution of the pipeline is finished.

1.2. **Post processing module**. Does optional modifications to the multi-fasta database, like removing redundant records with exact or substring sequences (using remove_dup_records.py) or removing the word 'isoform' from the record headers.

1.3. **Metadata module**. Generates some tables and files with metadata information about the database, like the number of species retrieved or the genes not found during the search. To know more about metadata module output, go to the *Metadata files* section.

### 2. uniparc_db.py

This script creates a protein database (multi-fasta) using the UniParc archive, a non-redundant database that contains all the proteins sequenced or predicted in UniProt, NCBI, and other resources. The search is focused on a specific taxonomic group (indicated using the NCBI tax ID) and can be reduced to a certain group of genes, indicated in text file.

The script has four main sections:

2.1. Importing the genes from the gene_list file to a Python list (if a path to the gene list is specified).

2.2. Recovering information about the TaxID using the Uniprot API. Specifically, a list of the descendent clades that belongs to the selected TaxID (e.g. Homo sapiens [TaxID=9606] is a descendent clade of Hominidae [TaxID=9604]).

2.3. Defining a function that takes a JSON object and turns it into FASTA format. This function is designed specifically to retrieve information from the UniParc JSON objects and also creating a specific FASTA header (see README.md file).

2.4. Using the UniParc API to retrieve the proteins from that database. The protein records are downloaded in JSON format and then converted into FASTA using the function described in section 3. If a gene list is specified, the search is iterated for each gene.

The download process through the API is done by pagination. This is a special type of download that retrieves records in batches of 500, allowing for faster and more efficient performance.

### 3. remove_dup_records.py

This script **searches for duplicate sequences** between records in a multi-fasta and **removes that redundancy**. It creates two files, one with the filtered multi-fasta (without redundancy) and another multi-fasta with the records removed. The way it detects redundancy is based on the sequence (DNA or protein), and it can occur in two different situations:

I.  Both records have the same sequence. In this case, the removed record is the first to be found in the multi-fasta. For example, among the following records, 'A' will be removed from the multi-fasta:

``` texinfo
A. MPIVCLGLLVFGLT
B. MPIVCLGLLVFGLT
```

II. One record has a sequence that is a substring from another. In this case, the removed record is the one containing the substring, in other words, the one with the shortest sequence length. For example, among the following records, 'B' will be removed from the multi-fasta:

``` texinfo
A. MPIVCLGLLVFGLT
B. IVCLGL
```

The script has three main sections:

3.1. Append the multi-fasta content into a Python list. To do so, firstly the multi-fasta file is converted into another multi-fasta but with the sequences of each record collapsed in one row. Then, that new multi-fasta is iterated to append the content into a Pyhton list.

3.2. Remove the redundant records from the list and write the outputs.

3.3. Print the number of removed records. Specifically, it prints in the terminal the number of exact duplicated records removed and the number of substrings records removed.

### 4. metadata_db.py

This script reads a protein multi-fasta database and **creates a folder with metadata information**. It has two main sections:

4.1. Calculate all the metadata information and stats from the multi-fasta. Also, creates all the files previously mentioned except the metadata_readme.txt file.

4.2. Create the metadata_readme.txt file and write there all the information previously calculated.
