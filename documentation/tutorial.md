# ProteoParc 1.1 tutorial

## Requirements

1.  **macOS or Linux operating system.** This software was thought to be employed in macOS and Linux operating systems. It won't probably work on Windows unless an alternative terminal is used.

2.  **Internet connection.** During the performance, it is important to have a stable internet connection. Otherwise, ProteoParc will stop and print an error.

3.  **Anaconda**. Anaconda/Conda is an environment management tool to quickly install, run, and update packages and their dependencies. It can be easily downloaded [here](https://www.anaconda.com/download). We extremely recommend rebooting your computer after the download. To check if Conda is already installed, type `conda` in the terminal.

    Disclaimer! This pipeline has been tested using conda v24.5.0

4.  **Git**. Git is a repository management software to work with file versions and download software from GitHub. It can be easily downloaded [here](https://git-scm.com). As in the previous step, you can check if Git is installed on your computer by typing `git` in the terminal.

## Set-up

To **download ProteoParc,** navigate through the terminal ([change directory](https://www.cyberciti.biz/faq/how-to-change-directory-in-linux-terminal/)) to the folder you want to store the pipeline. Then, type `git clone` the GitHub repository as it is shown in the next command.

``` bash
cd /Users/user_name/Documents/software # This path is an example
git clone https://github.com/guillecarrillo/proteoparc
```

To set-up the package requirements, navigate to `/proteoparc` and **create a conda virtual environment** by typing the following command.

``` bash
cd /Users/user_name/Documents/software/proteoparc # This path is an example
conda env create -f .set_up.yml # This action can take some minutes
```

This environment has all the package requirements to run the pipeline. If the previous command **brings an error**, try to install each package from scratch as it is shown below.

``` bash
conda create -n proteoparc python=3.12.3
conda activate proteoparc

# Bioconda and python modules
conda install bioconda::mafft==7.525
pip install pandas==2.2.2
pip install requests==2.32.3
pip install bio==1.7.1

# R packages
conda install conda-forge::r-ggplot2==3.4.4
conda install conda-forge::r-dplyr==1.1.4
conda install conda-forge::r-reshape2==1.4.4
conda install conda-forge::r-tidyr==1.3.1
```

The conda environment must be activated each time before running ProteoParc. If it is not, type the following command in the terminal:

``` bash
conda activate proteoparc
```

## Execution time! Make it simple

Before running ProtepParc, the user must **select a TaxID**. This number, which is shared between UniProt and NCBI, works as a unique identifier representing a taxonomic category. For instance, *Homo sapiens* (species) TaxID is 9606 and Mammalia (class) is 40674. To find a TaxID, you can search both in NCBI or [UniProt taxonomy browser](https://www.uniprot.org/taxonomy). We extremely recommend the UniProt one as it has a more intuitive interface. In general, the higher the scale of the TaxID the more proteins you will download. For instance, 9779 TaxID will download all Proboscidea (elephants, order) proteins, while 9783 will only download *Elephas maximus* (Indian elephant, species) proteins.

Now it's time to execute the pipeline, typing `python3 proteoparc.py` plus the desired options. There are only two mandatory arguments, which are:

-   `-p` or `--project` to indicate the name of the project

-   `-t` or `--tax-id` to indicate the TaxID

A simple and common execution can be:

``` bash
python3 proteoparc.py -p mammalia_enamelome -t 40674 -g enamelome.txt
```

In this example, ProteoParc is generating a database with enamel proteins for all the mammalian species. The argument `-g` indicates the path to [enamelome.txt](../documentation/example/enamelome.txt), a list of 15 enamel gene names to focus the search.

## Software arguments

**Mandatory arguments**

-   `-p` or `--project`. Indicate the name of the project.

-   `-t` or `--tax-id`. Indicate the TaxID.

**Optional arguments**

-   `-g` or `--genes`. Indicate a path to a text file containing a list of genes (one per row). The search will be focused on proteins annotated under those gene names. Look at [enamelome.txt](../documentation/example/enamelome.txt) for the correct way of formatting a gene list text file.

-   `--output-path`. Indicate a path to store the output folder. If this argument is not indicated, the result folder will be generated in the current working directory.

-   `--remove-redundancy` \| `--no-remove-redundancy`. Remove redundant records from the final database. By default, this action is switched on. Two types of redundant records are taken into account under this process:

    -   More than two records with the same sequence (one copy is maintained).

    -   One record with a sequence being a sub-string of another record's sequence (The sub-string sequence record is removed).

-   `--align-database` \| `--no-align-database`. Generate an alignment multi-fasta file per each different gene name. By default, this action is switched on.

## Output example

Two example outputs can be found in the documentation/[example](../documentation/example) folder

**1.- proboscidea_enamelome**. A protein download focused on 15 enamel proteins for all the Proboscidea (elephants) species.

``` bash
python3 proteoparc.py -p proboscidea_enamelome -t 9779 -g enamelome.txt
```

The printed output was:

``` texinfo
# Downloading the proteins indicated in 'enamelome.txt'
   59 proteins downloaded
# Removing redundant records
   0 records with the same sequence removed
   6 fragment records removed
# Aligning database per gene name
# Generating metadata information
```

**2.- mammuthus_proteome**. A protein download focused on the whole available proteome of *Mammuthus* genus (mammoths)

``` bash
python3 proteoparc.py -p mammuthus_proteome -t 37348
```

The printed output was:

``` texinfo
# Downloading the whole proteome
   500/1249
   1000/1249
   1249/1249
   1249 proteins downloaded
# Removing redundant records
   0 records with the same sequence removed
   63 fragment records removed
# Aligning database per gene name
# Generating metadata information
[1] "WARNING: Barplot might be messy due to a high number of gene names or species"
[1] "WARNING: Grid plot might be messy due to a high number of gene names or species"
```

Notice a warning message is printed due to a high number of gene names recovered. This can mess up the metadata plots, making them unintelligible.
