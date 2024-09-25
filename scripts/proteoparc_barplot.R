# This script takes the species_genes.csv file (ouputed from proteoparc.py)
# and generates a regular barplot representing the number of different protein
# isoforms per gene name and species.

# Programmed by Guillermo Carrillo Martín - May 2024

library(ggplot2)
library(dplyr, warn.conflicts = FALSE)

# Parse the input database path
species_genes_path <- normalizePath(commandArgs(trailingOnly = TRUE)[1]) # Input the species_genes.csv file
species_gene_df <- read.csv(species_genes_path)

# Summarize the 'Count' by 'Gene' and 'Species'
species_gene_reordered_df <- species_gene_df %>%
  group_by(Gene, Species) %>%
  summarise(Total_Count = sum(Count), .groups = 'drop') %>%
  ungroup() %>%
  group_by(Gene) %>%
  mutate(Total_Gene_Count = sum(Total_Count)) %>%
  ungroup() %>%
  arrange(desc(Total_Gene_Count))

# Create a new column to know if the value in "Gene" is a gene name or an "na" value.
species_gene_reordered_df <- species_gene_reordered_df %>%
  mutate(gene_or_na = ifelse(grepl("no gene", Gene), "na_value", "gene_name"))

# Print a warning if more than 35 genes or 15 species are indicated
if (n_distinct(species_gene_reordered_df$Gene) > 35 | n_distinct(species_gene_reordered_df$Species) > 15){
  print("WARNING: Barplot might be messy due to a high number of gene names or species")
}

# Create the bar plot
species_gene_barplot <- ggplot(species_gene_reordered_df, aes(x = reorder(Gene, -Total_Gene_Count), y = Total_Count, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Protein variations count in database",
    fill = "Species"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size=16)
  ) +
  facet_grid(. ~ gene_or_na, scales = "free", space = "free")

# Output the plot in as a jpg file
output_path <- paste(dirname(species_genes_path), "species_per_gene_barplot.jpg", sep = "/")
ggsave(output_path, plot = species_gene_barplot, dpi = 300, width = 11, height = 6)