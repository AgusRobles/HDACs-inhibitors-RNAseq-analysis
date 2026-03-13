library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(dplyr)
library(enrichplot)
library(pheatmap)
library(scales)

dds <- readRDS("data/DEseq_analysis.rds")

res966 <- results(dds, contrast = c("treatment", "RGFP966", "vehicle"))



# ORA analysis -----------------------------------------------------------------

universe_genes <- rownames(res966_df[!is.na(res966_df$padj), ])     # define the universe as the genes obtained from the DESeq analysis with adjusted p-value

DE_up <- subset(res966_df, padj < 0.05 & log2FoldChange > 0) %>%    # create a data frame with up-regulated genes only
  rownames()

DE_down <- subset(res966_df, padj < 0.05 & log2FoldChange < 0) %>%  # create a data frame with down-regulated genes only
  rownames()

# Up-regulated genes
resGO_up <- enrichGO(
  gene = DE_up, 
  universe = universe_genes, 
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP", 
  pvalueCutoff = 1, 
  qvalueCutoff = 1)

resGO_upTable <- as.data.frame(resGO_up)

n_11 = resGO_upTable$Count                                  # number of DE genes in the particular gene sets (GO terms)
n_10 = length(intersect(resGO_up@gene, resGO_up@universe))  # total number of DE genes
n_01 = as.numeric(gsub("/.*$", "", resGO_upTable$BgRatio))  # total number of genes in the gene sets
n = length(resGO_up@universe)                               # total number of genes (intersection between universe_genes and collection from org.Mm.eg.db)

hyper_mean = n_01 * n_10/n                              # mean of hypergeometric distribution   

n_02 = n - n_01                                         # total number of genes that are not in the specific gene sets
n_20 = n - n_10                                         # total number of genes that are not DE
hyper_var = n_01 * n_10 / n * n_20 * n_02 / n / (n - 1) # standard deviation

resGO_upTable <- resGO_upTable %>% 
  mutate(
    DE_Ratio = n_11/n_01,                            # proportion of DE genes IN the geneset
    GS_size = n_01,                                  # size of gene sets
    log2_Enrichment = log2( (n_11/n_10)/(n_01/n) ),  # log2 enrichment
    zScore = (n_11 - hyper_mean)/sqrt(hyper_var)     # z-score
  )


# Down-regulated genes
resGO_down <- enrichGO(
  gene = DE_down, 
  universe = universe_genes, 
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP", 
  pvalueCutoff = 1, 
  qvalueCutoff = 1)

resGO_downTable <- as.data.frame(resGO_down)

n_11 = resGO_downTable$Count
n_10 = length(intersect(resGO_down@gene, resGO_down@universe))
n_01 = as.numeric(gsub("/.*$", "", resGO_downTable$BgRatio))
n = length(resGO_down@universe)

hyper_mean = n_01 * n_10/n

n_02 = n - n_01
n_20 = n - n_10
hyper_var = n_01 * n_10 / n * n_20 * n_02 / n / (n - 1)

resGO_downTable <- resGO_downTable %>% 
  mutate(
    DE_Ratio = n_11/n_01,                             # Proportion of DE genes IN the geneset
    GS_size = n_01,                                   # size of gene sets
    log2_Enrichment = log2( (n_11/n_10)/(n_01/n) ), 
    zScore = -(n_11 - hyper_mean)/sqrt(hyper_var)
  )

resGO_combined <- rbind(resGO_downTable, resGO_upTable)

# Visualizations ---------------------------------------------------------------

resGO_top20 <- rbind(
  resGO_downTable %>% 
  arrange(p.adjust) %>%                                         # orders list of genes based on p.adjust 
  slice(1:10), 
  resGO_upTable %>% 
    arrange(p.adjust) %>%                                         # orders list of genes based on p.adjust 
    slice(1:10)) %>%                                               # pick top 20 terms by p.adjust
  mutate(GO = factor(Description, levels = rev(Description)),   # make GO term a factor to control y-axis order
         log10_padj = -log10(p.adjust))   

# Bar plot
total_bar <- ggplot(resGO_top20, aes(x = log10_padj, y = GO, fill = zScore)) +  
  geom_col() +  
  scale_fill_gradient2(
    low = "#1a80bb", mid = "grey70", high = "#ea801c", midpoint = 0
  ) +
  labs(
    x = expression(-log[10]("adjusted p-value")),
    y = "",
    fill = "z-score"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_blank()
  ) +
  geom_text(aes(label = label_scientific(digits = 2)(p.adjust)), 
            x = max(resGO_top20$log10_padj) * 1.05,
            hjust = 1)
total_bar
ggsave("output/RGFP966 20 most over-represented GO terms.pdf", total_bar, width = 8, height = 7)


# create a list with GO terms significantly enriched related to neuronal function 

selected_terms <- c(
  "GO:0061564",  # axon development
  "GO:0007409",  # axonogenesis
  "GO:0010975",  # regulation of neuron projection development
  "GO:0099173",  # postsynapse organization
  "GO:0035418"   # protein localization to synapse
)

# Filter the results based on the selected GO terms 
neuronal_terms <- resGO_upTable %>%                  
  filter(ID %in% selected_terms & p.adjust < .05)

neuronal_genes <- neuronal_terms$geneID %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

neuronal_DE <- res966[intersect(neuronal_genes, rownames(res966)), ] %>%
  as.data.frame() %>%
  rownames_to_column("gene")

neuronal_DE <- neuronal_DE %>%
  filter(padj < 0.05) %>% 
  arrange(padj)

write.csv(neuronal_DE, file = "output/neuronal GO terms DE genes RGFP966.csv")

go_gene_table <- neuronal_terms %>%    # create a table with each gene and its associated GO term in rows
  select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/")

gene_counts <- go_gene_table %>%       # count the number of terms each gene appears in
  count(geneID, sort = TRUE)

gene_counts %>% filter(n > 1)          # filter the genes that appear in more than one GO term


ggplot(go_gene_table, aes(x = Description, y = geneID)) +
  geom_tile(fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.background = element_blank()
  ) +
  labs(x = "GO term", y = "Gene")

# Volcano plot with DESeq results ----------------------------------------------

res966_df <- as.data.frame(res966) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%                           # eliminates genes without adjusted p-value 
  mutate(
    significance = case_when(                        # add a column "significance"
      padj < 0.05 & log2FoldChange > 0 ~ "Up",       # genes with a significant positive log2 fold change are labeled "Up" (up-regulated)
      padj < 0.05 & log2FoldChange < 0 ~ "Down",     # genes with a significant negative log2 fold change are labeled "Down" (down-regulated)
      TRUE ~ "Not significant"
    )
  )

genes_to_label <- neuronal_DE %>%      # creates a list of the top 10 DE genes in order of significance
  arrange(padj) %>%                                         
  slice(1:10) %>% 
  pull(gene)

volcano <- ggplot(res966_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.8, size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(                                                               
    data = subset(res963_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 4,
    max.overlaps = 50,
    box.padding = 0.5
  ) +
  scale_color_manual(values = c(
    "Up" = "#ea801c",
    "Down" = "#1a80bb",
    "Not significant" = "grey70"
  ), labels = c("Down-regulated", "Unchanged", "Up-regulated")) +
  
  xlim(-9, 9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 1), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value", 
    color = ""
  )
volcano
ggsave("output/RGFP966 volcano plot.pdf", volcano, width = 6, height = 4.5)
