# Code overview ProteoParc 1.0

This document explains the code and output of all the scripts that compose ProteoParc, a pipeline for the creation of protein databases focused on paleoprotein mass spectrometry identification. This software was entirely written in Python v3.9.12 along with the following packages:

-   numpy v1.25.0
-   pandas v2.0.3
-   requests v2.31.0

## ProteoParc output

### Protein multi-fasta database

The main output of this pipeline is the multi-fasta protein database. This file contains all the protein records in the standardized fasta format. The header format of each protein record is based on the UniProt one but including some changes. In addition, commas (,), semicolons (;) and colons(:) are removed from the header to avoid interferences when creating csv files. The representation of the header format is:

``` textinfo
>Database|Unique_identifier|Last_update Name_of_the_protein OS=Specie OX=TaxID GN=Gene_name SV=Sequence_version
```

| Field               | Col2 |
|---------------------|------|
| Database            |      |
| Unique Identifier   |      |
| Last update         |      |
| Name of the protein |      |
| Specie              |      |
| TaxID               |      |
| Sequence version    |      |

sequence version = more info <https://www.uniprot.org/help/uniparc>

Below, you can find the first header of the database stored in documentation/example/proboscidea_enamelome as an example

``` textinfo
>RefSeq|UPI002116250B|2023-03-08 albumin X1 OS=Elephas maximus indicus OX=99487 GN=ALB
```

### Metadata tables

metadata = explain outputs and infraestimation

### Other optional outputs

-   Additionally, if remove duplicates: fasta folder with unfiltered database and records removed. If remove isoform word: only in the final database, not in the unfiltered one.

## ProteoParc code

All the scripts are properly commented within each script.
