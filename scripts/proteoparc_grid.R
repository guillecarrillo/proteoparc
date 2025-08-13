# This script takes the species_genes.csv file (ouputed from proteoparc.py)
# and generates a grid representing the presence or absence of a protein 
# per each species in the database.

# Programmed by Guillermo Carrillo Martín - May 2024

library(ggplot2)
library(reshape2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)

# Parse the input database path
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript proteoparc_grid.R <species_genes.csv> [<gene_list.txt>]")
}

species_genes_path <- normalizePath(args[1]) # Input the species_genes.csv file
if (length(args) == 2) {
  gene_list_path <- args[2]  # user provided a path
} else {
  gene_list_path <- NULL  # no gene list provided
}

# Read files
species_gene_df <- read.csv(species_genes_path)

if (is.null(gene_list_path)) { # If no gene list is provided, keep all genes from data frame
  gene_list <- unique(species_gene_df$Gene)
} else { # If gene list, read it as a vector of characters
  gene_list <- readLines(gene_list_path)
}

# Remove the "no gene" records from the gene list and the dataframe
if ("no gene" %in% gene_list) {
  gene_list <- gene_list[gene_list != "no gene"]
}

species_gene_df <- species_gene_df[species_gene_df$Gene != "no gene", ]

# Complete data frame with all possible combinations, setting count as 0
species_gene_df <- species_gene_df %>% complete(Gene = gene_list,
                                                Species = unique(species_gene_df$Species), 
                                                fill = list(Count = 0))

# Create the presence or absence column 
species_gene_df <- species_gene_df %>%
  mutate(presence_or_absence = ifelse(Count > 0, "Presence", "Absence"))

# Print a warning if more than 35 genes or 15 species are indicated
if (n_distinct(species_gene_df$Gene) > 35 | n_distinct(species_gene_df$Species) > 15){
  print("WARNING: Grid plot might be messy due to a high number of gene names or species")
}

# Create the heat map
species_gene_grid <- ggplot(species_gene_df, aes(x = Gene, y = Species, fill = presence_or_absence)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  scale_alpha_identity() +
  scale_fill_manual(values = c("Absence" = "white", "Presence" = "royalblue")) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Protein presence in database") +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    text = element_text(size=16)
  )

# Output the plot in as a jpg file
output_path <- paste(dirname(species_genes_path), "species_per_gene_grid.jpg", sep = "/")
ggsave(output_path, plot = species_gene_grid, dpi = 300, width = 11, height = 6)