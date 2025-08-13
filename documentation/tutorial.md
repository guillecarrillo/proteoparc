# ProteoParc 1.2 tutorial

## Requirements
1.  **macOS or Linux operating system.** This software was thought to be employed in macOS and Linux operating systems. It won't probably work on Windows unless an alternative terminal is used.

2.  **Internet connection.** During the performance, it is important to have a stable internet connection. Otherwise, ProteoParc will stop and print an error.

3.  **Anaconda**. Anaconda/Conda is an environment management tool to quickly install, run, and update packages and their dependencies. It can be easily downloaded [here](https://www.anaconda.com/download). We extremely recommend rebooting your computer after the download. To check if Conda is already installed, type `conda` in the terminal.

    Disclaimer! This pipeline has been tested using conda v24.5.0

4.  **Git**. Git is a repository management software to work with file versions and download software from GitHub. It can be easily downloaded [here](https://git-scm.com). As in the previous step, you can check if Git is installed on your computer by typing `git` in the terminal.

## Set-up
To **download ProteoParc,** navigate through the terminal ([change directory](https://www.cyberciti.biz/faq/how-to-change-directory-in-linux-terminal/)) to the folder you want to store the pipeline. Then, clone the GitHub repository as it is shown in the next command.

``` bash
cd /Users/user_name/Documents/software # Example path
git clone https://github.com/guillecarrillo/proteoparc
```

To set-up the package requirements, change directory to `/proteoparc` and **create a conda virtual environment** by typing the following command.

``` bash
cd /Users/user_name/Documents/software/proteoparc # Example path
conda env create -f .set_up.yml # This action can take some minutes
```

This environment has all the package requirements to run the pipeline. If the previous command **brings an error**, try to install each package from scratch as it is shown below.

``` bash
conda create -n proteoparc-v1.2 python=3.12.3
conda activate proteoparc-v1.2

# Bioconda and python modules
conda install bioconda::mafft==7.525
pip install pandas==2.2.2
pip install requests==2.32.3
pip install bio==1.7.1
pip install requests_futures==1.0.2

# R packages
conda install conda-forge::r-ggplot2==3.4.4
conda install conda-forge::r-dplyr==1.1.4
conda install conda-forge::r-reshape2==1.4.4
conda install conda-forge::r-tidyr==1.3.1
```

The conda environment must be activated each time before running ProteoParc. If it is not, type the following command:

``` bash
conda activate proteoparc-v1.2
```

## Execution time! Make it simple
Before running ProtepParc, the user must **select a TaxID**. This number, which is shared between UniProt and NCBI, works as a unique identifier representing a taxonomic category. For instance, *Homo sapiens* (species) TaxID is 9606 and Mammalia (class) is 40674. To find a TaxID, you can search both in NCBI or [UniProt taxonomy browser](https://www.uniprot.org/taxonomy). We extremely recommend the UniProt search engine as it has a more intuitive interface. In general, the higher the scale of the TaxID the more proteins you will download. For instance, 9779 TaxID will download all Proboscidea (elephants, order) proteins, while 9783 will only download *Elephas maximus* (Indian elephant, species) proteins.

Now it's time to execute the pipeline, typing `python3 proteoparc.py` + the desired options. There are only two mandatory arguments, which are:

-   `-p` or `--project` to indicate the name of the project

-   `-t` or `--tax-id` to indicate the TaxID

A simple and common execution can be:

``` bash
python3 proteoparc.py -p mammalia_enamelome -t 40674 -g enamelome.txt
```

In this example, ProteoParc is generating a database with enamel proteins for all the mammalian species. The argument `-g` indicates the path to [enamelome.txt](../documentation/example/enamelome.txt), a list of 15 enamel gene names to restrict the search scope.

## Software arguments
**Mandatory arguments**

-   `-p` or `--project`. Indicates the name of the project.

-   `-t` or `--tax-id`. Indicates the TaxID.

**Optional arguments**

-   `-g` or `--genes`. Indicates a path to a text file containing a list of genes (one per row). The search will be focused on proteins annotated under those gene names. An example of a gene list file can be seen in [enamelome.txt](../documentation/example/enamelome.txt).

-   `--output-path`. Indicates a path to store the output folder. If this argument is not indicated, the result folder will be generated in the current working directory.

-   `--remove-redundancy` \| `--no-remove-redundancy`. Remove redundant records from the final database. By default, this action is switched on. Two types of redundant records are removed throught this process:

         I.  Both records have the same sequence. In this case, the removed record is the first to be found in the multi-fasta. For example, among the following records, 'A' will be removed from the multi-fasta:

         A. MPIVCLGLLVFGLT
         B. MPIVCLGLLVFGLT

         II. One record has a sequence that is a substring from another. In this case, the removed record is the substring, in other words, the one with the shortest sequence length. For example, among the following records, 'B' will be removed from the multi-fasta:

         A. MPIVCLGLLVFGLT
         B.   IVCLGL

-   `--align-database` \| `--no-align-database`. Generate an alignment multi-fasta file per each different gene name. By default, this action is switched on.

-   `--ignore-json` \| `--no-ignore-json`. Ignore the JSON files to generate the metadata values. This way, there will be only one repository, species and TaxID in each record metadata.

## Output example
Two example outputs can be found in the documentation/[example](../documentation/example) directory.

**1.- proboscidea_enamelome**. A protein database including 15 enamel genes from all the Proboscidea (elephants) species.

``` bash
python3 proteoparc.py -p proboscidea_enamelome -t 9779 -g enamelome.txt
```

Printed output (August 13th 2025)
``` texinfo
# Downloading the proteins indicated in 'enamelome.txt'
   66 records will be downloaded
# Removing redundant records
   0 records with the same sequence removed
   6 fragment records removed
# Aligning database per gene name
# Generating metadata information
```

**2.- mammuthus_proteome**. A protein database including the the whole available proteome of *Mammuthus* genus (mammoths)

``` bash
python3 proteoparc.py -p mammuthus_proteome -t 37348
```

Printed output (August 13th 2025)
``` texinfo
# Downloading the whole proteome
   1250 records will be downloaded
# Removing redundant recordsd
   0 records with the same sequence removed
   64 fragment records removed
# Aligning database per gene name
# Generating metadata information
[1] "WARNING: Barplot might be messy due to a high number of gene names or species"
[1] "WARNING: Grid plot might be messy due to a high number of gene names or species"
```

Notice that a warning message is printed due to the high number of gene names recovered. This can mess up the metadata plots, making them unintelligible.
