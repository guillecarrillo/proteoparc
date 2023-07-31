# Code overview ProteoParc 1.0

This document explains the code and outputs of all the scripts that compose ProteoParc, a pipeline for the creation of protein databases focused on paleoprotein mass spectrometry identification. This software was entirely written in Python v3.9.12 along with the following packages:

-   numpy v1.25.0
-   pandas v2.0.3
-   requests v2.31.0

## ProteoParc output

### Protein multi-fasta database

The main output of this pipeline is a multi-fasta protein database. The header of each record has no commas (,), semicolons (;) or colons(:) to avoid interferences when creating csv files, also, the architecture of the header itself is based on the UniProt format but including some changes:

``` textinfo
>Database|Unique_identifier|Last_update [Name_of_the_protein] OS=Specie OX=TaxID [GN=Gene_name] SV=Sequence_version
```

| Field                   | Description                                                                                                                                                                                                   |
|----------------------|-------------------------------------------------|
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

As can be seen in the [tutorial](../documentation/tutorial.md), the protein download can be focused in some specific proteins using a gene list to restrict the search. If this option is executed, all the records should have a gene name specified in the header. Otherwise, if the full proteome of a taxon is downloaded because there is no gene restriction, it can happen that a record has not an associated gene name due to bad annotated proteins.

### Metadata files

The metadata module creates seven files including metadata information about the generated database. This output files will change depending on whether the search was restricted to a list of genes or not. In short, the seven metadata files are:

1.  **summary.txt**; A text file that contains a summary of all the metadata information retrieved from the multi-fasta. It also shows the paths to all the other metadata files and indicates a warning message if there are records with missing fields, like the gene name or the species.

2.  **records_info.csv**; A CSV file that contains all the information present in each record header. See the previous section for a detailed description of the header information.

3.  **databases_employed.csv**; A CSV file that contains the number of the different sources (databases) used to build the multi-fasta.

4.  **genes_retrieved.csv**; A CSV file that shows the genes retrieved and counts the number of records with each gene.

5.  **genes_NOT_retrieved.txt**; A text file that contains the genes not retrieved. This file is only created if the multi-fasta was built using a list of genes.

6.  **species_retrieved.csv**; A CSV file that contains the number of species retrieved, indicating the scientific name and the TaxID. It also counts the number of records with each species.

7.  **species_genes.csv**; A CSV file that shows the genes retrieved for each species. Explain difference between gene list or full proteome.

metadata = explain outputs and infraestimation. Species and subspecies are count as different taxa to the metadata info.

### Other optional outputs

Additionally, if remove duplicates: fasta folder with unfiltered database and records removed. If remove isoform word: only in the final database, not in the unfiltered one.

## ProteoParc code

The ProteoParc pipeline consist of four Python scripts. Three of them perform specific actions in each module, being the other one an script that merges and concatenates this steps to create the usage of the pipeline.

### 1. proteoparc.py

script that merges the other three. Options mentioned in LINK TO TUTORIAL.
