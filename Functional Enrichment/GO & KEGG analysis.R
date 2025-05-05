BiocManager::install("enrichplot")


# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)   # Database for human gene annotations
library(enrichplot)      # For enrichment visualization
library(ggplot2)         # For general plotting
library(ReactomePA)
library(reactome.db)
library(ggpubr)

# Define the list of key genes
key_genes <- c("ADM", "CCL2", "AXL", "TLR7", "AIM2", "BATF2", "GBP1", "IL15RA", "IRF1", "TRIM21")

# Convert gene symbols to ENTREZ IDs
gene_ids <- bitr(key_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

#######----Perform GO enrichment analysis (Biological Process)
go_enrichment <- enrichGO(gene          = gene_ids$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "BP",  # Biological Process (Can be "MF" or "CC" as well)
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
head(go_enrichment)

# Plot GO enrichment results
dotplot(go_enrichment, showCategory=10) + theme_minimal() +
  theme(text = element_text(family = "Arial", size = 14, face = "bold"),  # Change font and size
        axis.text = element_text(family = "Times", size = 12, face = "bold.italic"), 
        axis.title = element_text(family = "Arial", size = 14, face = "bold"),
        plot.title = element_text(family = "Helvetica", size = 16, face = "bold", hjust = 0.5)) +
  labs(y = "Gene Ontology Terms")  # Add Y-axis title


##save in high dimension PNG image
ggsave("GO_enrichment_plot.png", 
       plot = last_plot(), 
       width = 8, height = 6, 
       dpi = 300)
##save in high dimension PNG image
ggsave("GO_enrichment_plot.tiff", 
       plot = last_plot(), 
       width = 8, height = 6, 
       dpi = 600, compression = "lzw")


#####################---Perform KEGG pathway enrichment analysis---############
kegg_enrichment <- enrichKEGG(gene         = gene_ids$ENTREZID,
                              organism     = 'hsa',  # 'hsa' for human
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.05)
head(kegg_enrichment)

# Plot KEGG pathway enrichment results
dotplot(kegg_enrichment, showCategory=10, title="KEGG Pathway Enrichment") + theme_minimal()

#####################--REACTOME enrichment analysis---#########################

# Install missing dependency (reactome.db)
BiocManager::install("reactome.db") 
#the reactome.db pathway was installed using source package. Downloading source package directly from website and adding in 'install packages' section.

BiocManager::install("ReactomePA")
library(ReactomePA)
library(dplyr)

reactome_results <- enrichPathway(gene = gene_ids$ENTREZID, 
                                  organism = "human", 
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.05, 
                                  readable = TRUE)  # Convert IDs to gene symbols
head(reactome_results)
dotplot(reactome_results, showCategory=10, title="Reactome Pathway Enrichment") + theme_minimal()


################## KEGG and REACTOME combined enrichment analysis ###################

# Add source column to each result
kegg_df <- as.data.frame(kegg_enrichment)
kegg_df$Source <- "KEGG"
reactome_df <- as.data.frame(reactome_results)
reactome_df$Source <- "Reactome"

# Combine both enrichment results
combined_df <- bind_rows(kegg_df, reactome_df)

# Optional: Filter top pathways by adjusted p-value
top_combined <- combined_df %>%
  group_by(Source) %>%
  slice_min(order_by = p.adjust, n = 10)

# Plot
ggplot(top_combined, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Pathway Enrichment of hub genes",
       x = "Pathways",
       y = "-log10(adj. p-value)") +
  theme_minimal() +
  scale_fill_manual(values = c("KEGG" = "#1f77b4", "Reactome" = "#ff7f0e")) +
  theme(
    text = element_text(family = "Arial", size = 14, face = "bold"),  # Overall text
    axis.text.y = element_text(family = "Times", size = 12, face = "bold.italic"),  # Y-axis labels
    axis.text.x = element_text(family = "Times", size = 12, face = "bold"),         # X-axis labels
    axis.title = element_text(family = "Arial", size = 14, face = "bold"),          # Axis titles
    plot.title = element_text(family = "Helvetica", size = 16, face = "bold", hjust = 0.5)  # Centered title
  )

ggsave("Pathway_Enrichment_Hub_Genes.png", width = 12, height = 4, dpi = 600, units = "in")
ggsave("Pathway_Enrichment_Hub_Genes.tiff", width = 12, height = 4, dpi = 600, units = "in")



