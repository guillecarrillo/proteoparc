# ProteoParc 1.0 tutorial

This document describes the usage of ProteoParc, a pipeline for the creation of protein databases focused on paleoprotein mass spectrometry identification.

## Requirements

1.  **macOS or Linux operating system.** This software was meant to be used in macOS and Linux operating systems, so it won't probably work on Windows operating system unless an alternative terminal or computer cluster is used.

2.  **Internet connection.** During the usage of the pipeline it's important to have a stable internet connection. The download time will be affected by this and the performance of the pipeline will be affected by this.

3.  **Anaconda software**. Anaconda (or just Conda) is an environment management system that allows you to quickly install, run, and update packages and their dependencies. The user can easily download this software [here](https://www.anaconda.com/download), we extremely recommend to reboot the system after the download. You can check if Conda is installed in your computer just typing `conda` in your terminal.

4.  **Git software**. Git is a software that allows the user to work with file repositories and download software from GitHub. It can be easily downloaded from [here](https://git-scm.com). As in the previous step, you can check if Git is installed in your computer just typing `git` in your terminal.

## Set up

First of all, the user needs to **download the pipeline**. To do so, open the terminal and navigate to the folder you want to store the pipeline using the [`cd` command](https://www.cyberciti.biz/faq/how-to-change-directory-in-linux-terminal/). As an example, imagine you want to store the pipeline in a folder named *pipelines* inside the *Documents* folder, in this case you should type:

``` bash
# This path is an example, the user should writte it's own persolanized path
cd /Users/user_name/Documents/pipeline
```

Then, type the following in the terminal:

``` bash
git clone https://github.com/guillecarrillo/proteoparc
```

Once the pipeline is downloaded in your computer, change again your working directory to the proteoparc folder using the `cd` command. Following the previous example, you can type:

``` bash
# This path is an example, the user should writte it's own persolanized path
cd /Users/user_name/Documents/pipeline/proteoparc
```

Now, **create the conda virtual environment** typing the following command. This action can take some minutes.

``` bash
conda env create -f .set_up.yml
```

This environment will have all the software needed to run the pipeline, which mainly consists of Python 3.9 plus some packages such as *pandas* or *requests*. If the previous command **brings an error**, try this instead:

``` bash
conda create -n proteoparc python=3.9
conda activate proteoparc
pip install pandas
pip install requests
conda deactivate
```

These steps only have to be done once, so now you are ready to use the pipeline!

## Execution time! Make it simple

Firstly, you need to **select a TaxID** to focus your protein search. This number, which is shared between UniProt and NCBI, works as a unique identifier that represents a taxonomic category, for instance, *Homo sapiens* TaxID is 9606 and Mammalia (class) is 40674. To find the TaxID that corresponds to your desired taxa, you can search both in NCBI or [UniProt taxonomy browser](https://www.uniprot.org/taxonomy), but we extremely recommend the UniProt one as it has a more intuitive interface. The higher the scale of the TaxID the more proteins you will download. For instance, using 9779 TaxID will download all Proboscidea (elephants, order) proteins, and using 9783 will only download *Elephas maximus* (Indian elephant) proteins.

The next step is to **activate the conda environment**, to do so, type the following command in the terminal:

``` bash
conda activate proteoparc
```

Remember to activate the environment every time before using the pipeline. You can deactivate the environment by typing `conda deactivate` if you want to stop working with ProteoParc.

Now it's time to execute the pipeline, to do so, **your working directory must be in the proteoparc folder**. Now, type `python3 proteoparc.py` plus the required and optional arguments. The only mandatory arguments are:

-   `-p` or `--project` to indicate the name of the project

-   `-t` or `--tax-id` to indicate the TaxID

A simple execution can be:

``` bash
python3 proteoparc.py -p mammalia_enamelome -t 40674 -g enamelome.txt
```

In this case, we are downloading enamel proteins for all the mammalian species. The argument `-g` indicates the path to a list of 15 enamel genes to focus the search. Remember that you can see all the output explanation in the [**Code Overview**](../documentation/code.md) manual.

## All the pipeline arguments

More options can be added to the pipeline to perform a more personalized search.

**Mandatory arguments**

| Command             | Description                      |
|---------------------|----------------------------------|
| `-p` or `--project` | Indicate the name of the project |
| `-t` or `--tax-id`  | Indicate the TaxID               |

**Optional arguments**

-   `-g` or `--genes`. Indicate a path to a text file containing a list of genes (one per row) and focus the search to only the proteins that come from those genes. In the next link you can find an example of a valid [gene list](../documentation/example/enamelome.txt) text file.

-   `--path`. Indicate a path to store the output folder. If this argument is not indicated, the result folder will be created in the current working directory.

-   `--remove-duplicates` \| `--no-remove-duplicates`. Specify if the duplicate records are removed during the post-processing module. By default, the duplicates will be removed even if `--remove-duplicates` is not specified. This option is designed to simplify the database and reduce computational cost during protein identification.

-   `--remove-isoform-word` \| `--no-remove-isoform-word`. Specify if the 'isoform' word is removed from the multi-fasta headers during the post-processing module. By default, the word 'isoform' will be removed even if `--remove-duplicates` is not specified. This option was designed because of an observed output interference in some protein identification software.

Remember, you can check easily all the pipeline options by typing `python3 proteoparc.py -h`

## Output example

In the [example](../documentation/example) folder, you can find two database examples created using different parameters. These examples were done in July 2023, so future performances with the same framework can change the output due to an UniParc update.

1.  **Quick download**. A protein search focused on enamel Proboscidea (elephants) proteins. The gene list used contained 15 enamel gene names.

``` bash
python3 proteoparc.py -p proboscidea_enamelome -t 9779 -g enamelome.txt
```

The printed output was:

``` texinfo
# Downloading the proteins indicated in 'enamelome.txt'
   59 proteins downloaded
# Removing duplicate records
   0 records with the same sequence removed
   5 fragment records removed
# Creating metadata information
```

2.  **Full proteome**. Downloading the full available proteome of the *Mammuthus* genus (mammoths)

``` bash
python3 proteoparc.py -p mammuthus_proteome -t 37348
```

The printed output was:

``` texinfo
# Downloading the full proteome
   500 / 1249
   1000 / 1249
   1249 / 1249
   1249 proteins downloaded
# Removing duplicate records
   0 records with the same sequence removed
   61 fragment records removed
# Creating metadata information
```
