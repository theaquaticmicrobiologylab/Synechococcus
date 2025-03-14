# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(gtools)
library(ggtree)
library(ggtext)
setwd("/Users/niu_mwh/Library/CloudStorage/GoogleDrive-hensonmw@gmail.com/My\ Drive/Syn_salinity_EM")
# Load data
#counts_df <- read.csv("merged_counts_derep97.txt", sep="\t", check.names=FALSE)  # Read merged count table
counts_df <- read.csv("merged_counts_UpdatedAllGenomes.txt", sep="\t", check.names=FALSE)
#metadata_df <- read.csv("Syn_salinity_2025_v3.csv")  # Read genome metadata
metadata_df <- read.csv("Syn_salinity_2025.csv")  # Read genome metadata

# Rename "Reference" to "genome_name" to match metadata
colnames(counts_df)[1] <- "genome_name"

# Step 1a: Aggregate contigs to genome level
counts_df$genome_name <- gsub("_00000[0-9]+", "", counts_df$genome_name)  # Remove contig numbers
agg_counts <- counts_df %>%
  group_by(genome_name) %>%
  summarise(across(everything(), sum, na.rm = TRUE))  # Sum counts per genome
agg_counts <- agg_counts[agg_counts$genome_name != "MV0610", ]
agg_counts_long <- agg_counts %>%
  pivot_longer(cols = -c(genome_name, Length), names_to = "sample", values_to = "count")
#write.csv(agg_counts, "agg_rawcounts.csv")

agg_counts_long <- agg_counts_long %>%
  group_by(sample) %>%
  mutate(total_reads_per_sample=sum(count), RPKM= (count / (Length / 1000)) / (total_reads_per_sample / 1e6))
# Step 3: Extract salinity from sample names
agg_counts_long_wSalinity <- agg_counts_long %>%
  filter(total_reads_per_sample > 50) %>%
  group_by(sample) %>%
  mutate(Salinity = as.numeric(sapply(strsplit(sample, "_"), function(x) x[2])))
#agg_counts_long_wSalinity$RPKM[is.na(agg_counts_long_wSalinity$RPKM)] <- 0 
agg_counts_long_wSalinity_filtered <- agg_counts_long_wSalinity %>%
  group_by(genome_name) %>%
  mutate(count_RPKM=sum(RPKM > 2), samples=n(), cutoff_value=(samples*0.05)) %>%
  filter(count_RPKM >= cutoff_value)

# Step 4: Determine peak salinity & salinity range for each genome

# Find peak salinity per genome
peak_salinity <- agg_counts_long_wSalinity_filtered %>%
  group_by(genome_name) %>%
  filter(RPKM == max(RPKM)) %>%
  select(genome_name, peak_salinity = Salinity)

# Identify salinity range before significant decline
salinity_range <- agg_counts_long_wSalinity %>%
  group_by(genome_name) %>%
  filter(RPKM >= (max(RPKM) * 0.25)) %>%  # Threshold: 50% of max RPKM
  summarise(min_salinity = min(Salinity), max_salinity = max(Salinity), min_RPKM=min(RPKM),
            max_RPKM=max(RPKM), median_RPKM=median(RPKM), average_RPKM=mean(RPKM), mean_sal=mean(Salinity), median_Sal=median(Salinity),
            ) 

# Merge results to classify genomes
final_classification <- peak_salinity %>%
  left_join(salinity_range, by = "genome_name") %>%
  mutate(category = case_when(
    min_salinity >= 18 & max_salinity >= 20 ~ "Marine",  # Both min & max above 20
    min_salinity >= 0 & min_salinity <= 10 & max_salinity >= 30 & max_salinity <= 40 ~ "Brackish",  # Range between 0-40
    max_salinity < 6 ~ "Freshwater",  # Max salinity below 5
    TRUE ~ "Brackish"  # Default category
  ))
merged_df <- final_classification %>%
  inner_join(metadata_df, by="genome_name")  %>%
  filter(median_RPKM >= 5) %>%
  select(genome_name, Taxonomy, peak_salinity, max_salinity, max_RPKM, median_RPKM, median_Sal,category, Source_salinity)

agg_counts_long_wSalinity_filtered <- agg_counts_long_wSalinity %>%
  group_by(genome_name) %>%
  mutate(count_RPKM=sum(RPKM > 5), samples=n(), cutoff_value=(samples*0.05)) %>%
  filter(count_RPKM >= cutoff_value)
agg_count_merged<- agg_counts_long_wSalinity_filtered %>%
  inner_join(metadata_df, by="genome_name") %>%
  inner_join(final_classification, by="genome_name") %>%
  mutate(match=factor((category==Source_salinity)))
agg_count_merged <- agg_count_merged[order(as.numeric(as.character(agg_count_merged$Subcluster))), ]

# Ensure facet labels are formatted, bold for MV isolates
agg_count_merged$facet_labels <- with(agg_count_merged, 
                                      ifelse(grepl("^MV", Short.Name), 
                                             sprintf("<b>%s</b>", Short.Name),  # Bold for MV isolates
                                             Short.Name)
)
agg_count_merged$facet_labels<-factor(agg_count_merged$facet_labels, levels = unique(agg_count_merged$facet_labels))



# Plot with modified facet labels
data <- ggplot(agg_count_merged, aes(Salinity, RPKM, color = category, shape = Source_salinity), size = 4) +
  geom_point() +
  scale_color_manual(values = c("#E69F00", "#009E73","#56B4E9")) +
  geom_smooth(se = FALSE, color = "black") +  # Ensure smooth line is black
  facet_wrap(~facet_labels, scales = "free_y") +  # Use the custom formatted labels
  geom_vline(xintercept = c(5, 18, 35)) +
  theme_bw() +
  theme(strip.text = element_markdown(size = 8))  # Enable markdown formatting for facet labels

data

# Save output
write.csv(agg_count_merged, "RPKM_classifications.csv", row.names=FALSE)



# Load the metadata file
metadata_redo <- read.csv("RPKM_reclaiffication.csv")
gene_clusters_redo <- read.delim("reclassified_SaluniqKO.tsv", sep='\t', stringsAsFactors=FALSE)

# Select relevant metadata columns
metadata_columns <- c("genome_name", "Taxonomy","RPKM", "OriginationSource", "Source_salinity", "peak_salinity", "min_salinity", "max_salinity", 
                      "min_RPKM", "max_RPKM", "median_RPKM", "average_RPKM", "mean_sal", "median_Sal", "category", "match")
metadata_V2 <- metadata_redo %>% select(all_of(metadata_columns))
metadata_V2_unique <- metadata_V2 %>% distinct(genome_name, Taxonomy, Source_salinity,category, match, .keep_all = FALSE)
metadata_3 <- metadata_V2 %>% distinct(genome_name, peak_salinity, Taxonomy, Source_salinity,category, .keep_all = FALSE)
merged_V2 <- left_join(metadata_df, metadata_3, by="genome_name")
merged_Cluster_peak<- merged_V2[!is.na(merged_V2$peak_salinity), ] %>%
  group_by(Subcluster) %>%
  summarise(median_SC_Sal=median(peak_salinity), mean=mean(peak_salinity))

cluster_peakavg<-metadata_3%>%
  group_by(Subcluster) %<%
  summarise(median=mean(peak_salinity)
            
            
# Count genes per KOfam per genome
kofam_counts <- gene_clusters_redo %>% 
  group_by(genome_name, KOfam_ACC, KOfam) %>% 
  summarise(gene_count = n(), .groups = 'drop')
# Merge with metadata
merged <- left_join(kofam_counts, metadata_V2_unique, by="genome_name")

#write.csv(merged, "Reclassified_KOcounts.csv")

heatmap_data <- merged %>% 
  select(genome_name, KOfam_ACC, gene_count) %>% 
  pivot_wider(names_from = KOfam_ACC, values_from = gene_count, values_fill = list(gene_count = 0)) %>% 
  column_to_rownames("genome_name")
#Run if you want only taxonomy as heatmap column 
# Check if Source_salinity is available for annotation
if (!"Source_salinity" %in% colnames(merged)) {
  stop("Error: Source_salinity column is missing in merged dataset")
}

# Create row annotation for Source_salinity
row_annotations <- merged %>% 
  select(genome_name, Source_salinity, category) %>% 
  distinct() %>% 
  column_to_rownames("genome_name")


# Generate color palettes for annotations
heatmap_colors <- colorRampPalette(c("white", "orange", "red"))(100)
salinity_colors <- c("#56B4E9","#E69F00","#009E73")  # Source_salinity
category_colors <- c("#56B4E9","#E69F00","#009E73")  # Category
# Assign colors to annotation values
names(salinity_colors) <- unique(row_annotations$Source_salinity)
names(category_colors) <- unique(row_annotations$category)

annotation_colors <- list(
  Source_salinity = salinity_colors,
  category = category_colors
)


# Generate heatmap
pheatmap(
  as.matrix(heatmap_data), 
  annotation_row = row_annotations, #%>% select(-PercComplete),  # Only keep relevant annotations
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  main = "KO Gene Counts Heatmap",
  color = heatmap_colors
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tibble)



# Read TSV files
metadata <- read.csv("RPKM_reclaiffication.csv")
gene_clusters_redo <- read.delim("reclassified_SaluniqKO.tsv", sep='\t', stringsAsFactors=FALSE)
all_metadata<-read.csv("Syn_salinity_2025.csv")
all_gene_clusters<-read.delim("SynSal_2025_pan_gene_clusters_summary.tsv", sep='\t', stringsAsFactors=FALSE)
# Count genes per KOfam per genome
kofam_counts <- all_gene_clusters %>% 
  group_by(genome_name, gene_cluster_id, KOfam_ACC, KOfam) %>% 
  summarise(gene_count = n(), .groups = 'drop')
merged <- left_join(kofam_counts, all_metadata, by=c("genome_name"))
#metadata_redo <- read.delim(metadata_redo_file, sep='\t', stringsAsFactors=FALSE)

# Select relevant metadata columns
metadata_columns <- c("genome_name", "RPKM", "Source_salinity", "peak_salinity", "min_salinity", "max_salinity", 
                      "min_RPKM", "max_RPKM", "median_RPKM", "average_RPKM", "mean_sal", "median_Sal", "category", "match","PercComplete", "Taxonomy")
metadata <- metadata %>% select(all_of(metadata_columns))

# Extract unique values from metadata_redo based on genome_id, category, and match
metadata_redo_unique <- metadata %>% distinct(genome_name, category, match,.keep_all = TRUE)

# Count genes per KOfam per genome
kofam_counts <- gene_clusters_redo %>% 
  group_by(genome_name, KOfam, COG20_PATHWAY_ACC) %>% 
  summarise(gene_count = n(), .groups = 'drop')

# Merge with unique metadata_redo
merged <- left_join(kofam_counts, metadata_redo_unique, by=c("genome_name"))

#
KOdata<-read.delim("Salinity_gene_KOlist.txt", sep='\t') %>%
  rename(KOfam=KO)
merged<-left_join(merged, KOdata, by="KOfam")
# Prepare data for heatmap
heatmap_data <- merged %>% 
  select(genome_name, Taxonomy, gene_count) %>% 
  pivot_wider(names_from = genome_name, values_from = gene_count, values_fill = list(gene_count = 0)) %>% 
  column_to_rownames("genome_name")

# Check if Source_salinity is available for annotation
if (!"Source_salinity" %in% colnames(merged)) {
  stop("Error: Source_salinity column is missing in merged dataset")
}

# Create row annotation for Source_salinity
row_annotations <- merged %>% 
  select(genome_name, Source_salinity, category, PercComplete) %>% 
  distinct() %>% 
  column_to_rownames("genome_name")

# Generate color palettes for annotations
heatmap_colors <- colorRampPalette(c("white", "orange", "red"))(100)
salinity_colors <- c("#56B4E9","#E69F00","#009E73")  # Source_salinity
category_colors <- c("#56B4E9","#E69F00","#009E73")  # Category
# Assign colors to annotation values
names(salinity_colors) <- unique(row_annotations$Source_salinity)
names(category_colors) <- unique(row_annotations$category)

annotation_colors <- list(
  Source_salinity = salinity_colors,
  category = category_colors
)


# Generate heatmap
pheatmap(
  as.matrix(heatmap_data), 
  annotation_row = row_annotations %>% select(-PercComplete),  # Only keep relevant annotations
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  main = "KO Gene Counts Heatmap",
  color = heatmap_colors
)


library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(ape)
tree<-read.tree("iqtree-out_bootstrapped/iqtree_SynSal_derep2025_out.treefile")
tree <- root(tree, outgroup = "GCF_000011385", resolve.root = TRUE)
metadata<-read.csv("Syn_salinity_2025.csv")

metadata$Taxonomy <- as.character(metadata$Taxonomy)
tree$tip.label <- as.character(tree$tip.label)
metadata <- metadata %>% filter(genome_name %in% tree$tip.label)

p <- ggtree(tree) +
  geom_tiplab(size = 3) +  # Label branches with taxonomy
  geom_tippoint(aes(color = metadata$Source_Salinity), size = 3) +  # Color by salinity
  new_scale_fill() +  # Add new layer for additional data
  geom_facet(geom = geom_bar, aes(x = GC), stat = "identity", data = metadata, panel = "% GC Content") +
  geom_facet(geom = geom_bar, aes(x = EstGenomeSize), stat = "identity", data = metadata, panel = "Genome Size") +
  theme_tree2()  # Use a nice theme

p <- ggtree(tree) %<+% metadata + 
  geom_tiplab(aes(label = Taxonomy), size = 3)  # Label the tips with Taxon names

# Add GC Content as a heatmap
# Plot the tree
p <- ggtree(tree) +
  geom_tiplab(size = 3) +  # Label branches with taxonomy
  geom_tippoint(aes(color = metadata$Source_Salinity), size = 3) +    # Color tips by salinity
  
  new_scale_fill() +  # Add new fill scale for %GC heatmap
  geom_fruit(
    data = metadata,
    geom = geom_tile,  # Heatmap for %GC
    mapping = aes(y = genome_name, fill = GC)
  ) + scale_fill_gradient(low = "blue", high = "red", name = "% GC") +
  
  new_scale_fill() +  # Add another fill scale for Genome Size heatmap
  geom_fruit(
    data = metadata,
    geom = geom_tile,  # Heatmap for Genome Size
    mapping = aes(y = genome_name, fill = EstGenomeSize),
    pwidth = 0.15
  ) + scale_fill_gradient(low = "yellow", high = "purple", name = "Genome Size (Mb)") +
  
  theme_tree2()  # Apply nice theme

p <- ggtree(tree) +  
  geom_fruit(
    data = metadata,
    geom = geom_tile,  # Color by Source Salinity categories
    mapping = aes(y = genome_name, fill = Source_salinity),
    pwidth = 0.15
  ) + scale_fill_manual(values=c("Blue", "Green", "Red", "Gold", "Pink", "Navy", "Yellow", "Gray"))+
  # Add the heatmaps first
  geom_fruit(
    data = metadata,
    geom = geom_tile,  # Heatmap for %GC
    mapping = aes(y = genome_name, fill = GC)
  ) + 
  scale_fill_gradient(low = "blue", high = "red", name = "% GC") +
  
  new_scale_fill() +  # Add new scale for Genome Size heatmap
  geom_fruit(
    data = metadata,
    geom = geom_tile,  # Heatmap for Genome Size
    mapping = aes(y = genome_name, fill = EstGenomeSize),
    pwidth = 0.15
  ) +
  scale_fill_gradient(low = "yellow", high = "purple", name = "Genome Size (Mb)") +

  
  # Add the tip labels and tip points last
  geom_tiplab(aes(color = metadata$Source_Salinity), size = 3) +  # Label branches with taxonomy and color by salinity
  geom_tippoint(aes(color = metadata$Source_Salinity), size = 3) +  # Color tips by salinity
  
  theme_tree2()  # Apply nice theme

# Load required libraries
library(readxl)
library(ggplot2)
library(ggrepel)

# Load the data from the Excel file
df <- read_excel("Manuscript/Fig/MV0715_proteomics_plotting.xlsx")


# Convert relevant columns to numeric
df$Diff <- as.numeric(df$Diff)  # Log fold-change
df$abs_Z <- as.numeric(df$abs_Z)  # Significance measure (Z-score)


# Create a threshold for significance (adjust if needed)
significance_threshold <- 1.96  # Commonly used Z-score cutoff

# Add a new column for significance classification
df$Significance <- ifelse(df$abs_Z > significance_threshold, "Significant", "Not Significant")

# Create the volcano plot
ggplot(df, aes(x = Diff, y = abs_Z)) +
  geom_point(alpha = 0.7) +  # Scatter plot with transparency
  labs(title = "Volcano Plot of Proteomics Data",
       x = "Log Fold-Change (Diff)",
       y = "Z-score (abs_Z)") +
  theme_minimal() +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "blue") +# Threshold line
  geom_text_repel(aes(label = ifelse(abs_Z > 0.9, `KO/COG/pFam`, "")), 
                  size = 3, max.overlaps = 15) + # Label only those with Z > 0.9
  geom_text_repel(aes(label = ifelse(abs(Diff) > 0.75, `KO/COG/pFam`, "")), 
                  size = 3, max.overlaps = 15)  # Label only those with Z > 0.75
