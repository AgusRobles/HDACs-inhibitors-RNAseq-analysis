suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
  library(ExploreModelMatrix)
  library(cowplot)
  library(ComplexHeatmap)
  library(apeglm)
})

se <- readRDS("data/norm_data_se2.rds")

dds <- DESeqDataSet(se, design = ~ treatment)

# Re-order the levels to put 'vehicle' as the reference level
levels(dds$treatment)
dds$treatment <- relevel(dds$treatment, ref = "vehicle")

# Differential expression analysis
dds <- DESeq(dds) 

# Results of DE analysis -------------------------------------------------------
results(dds)
resultsNames(dds)

res963 <- results(dds, contrast = c("treatment", "RGFP963", "vehicle"))
summary(res963)

res966 <- results(dds, contrast = c("treatment", "RGFP966", "vehicle"))
summary(res966)

# Plots ------------------------------------------------------------------------

plotMA(res963)

# This removes low expression and highly variable genes
res963Lfc <- lfcShrink(dds, coef = "treatment_RGFP963_vs_vehicle", res = res963)
plotMA(res963Lfc)

# This is another way of filtering low expression genes
res963NotFiltered <- results(dds,
                              contrast = c("treatment", "RGFP963", "vehicle"), 
                              independentFiltering = FALSE)
summary(res963)
summary(res963NotFiltered)
plotMA(res963NotFiltered)


plotMA(res966)

res966Lfc <- lfcShrink(dds, coef = "treatment_RGFP966_vs_vehicle", res = res966)
plotMA(res966Lfc)

res966NotFiltered <- results(dds,
                                  contrast = c("treatment", "RGFP966", "vehicle"), 
                                  independentFiltering = FALSE)
summary(res966)
summary(res966NotFiltered)
plotMA(res966NotFiltered)


# Visualizations

# Transform counts
vsd <- vst(dds, blind = TRUE)

##### RGFP966 #####
# Get top 10 DE genes
genes_966 <- res966[order(res966$pvalue), ] %>% 
  head(10) %>% 
  rownames()

# Create the dataframe with the counts of the top 10 DE genes
heatmapData_966 <- assay(vsd)[genes_966, ]

# Scale counts for visualization
heatmapData_966 <- t(scale(t(heatmapData_966)))

# Define treatment colors
treatment_colors <- c(
  "vehicle"  = "#F5D000",   
  "RGFP963"  = "#96C65C",   
  "RGFP966"  = "#1BB6AF"    
)

# Create a dataframe with the list of treatments
col_annot_df <- data.frame(
  Treatment = colData(vsd)$treatment
)
# Chenges row names from numbers to the sample names 
rownames(col_annot_df) <- colnames(vsd)


ha_col <- HeatmapAnnotation(
  Treatment = col_annot_df$Treatment,
  col       = list(Treatment = treatment_colors),
  annotation_name_side = "right",         
  show_annotation_name = F
)


ht_966 <- Heatmap(
  heatmapData_966,
  name = "Z-score",                       
  
  # Clustering
  cluster_rows    = T,
  cluster_columns = F,  
  show_row_dend   = F,
  
  # Labels and side
  row_names_side     = "left",
  column_names_side  = "bottom",          
  show_column_names  = TRUE,              
  column_names_gp    = gpar(fontsize = 10),
  row_names_gp       = gpar(fontsize = 10),
  #width = unit(10, "cm"),
  
  # Annotations
  bottom_annotation = ha_col,            
  
  # Layout tweaks
  heatmap_legend_param = list(
    title = "Z-score",
    title_position = "topcenter"
    )
  )

pdf("output/Heatmap DE genes RGFP966.pdf", ht_966, width = 4, height = 6)

draw(
  ht_966,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = T
)


dev.off()

##### RGFP963 #####
# Get top 10 DE genes
genes_963 <- res963[order(res963$pvalue), ] |>
  head(10) |>
  rownames()

heatmapData_963 <- assay(vsd)[genes_963, ]

# Scale counts for visualization
heatmapData_963 <- t(scale(t(heatmapData_963)))


ha_col <- HeatmapAnnotation(
  Treatment = col_annot_df$Treatment,
  col       = list(Treatment = treatment_colors),
  annotation_name_side = "right",         
  show_annotation_name = F
)


ht_963 <- Heatmap(
  heatmapData_963,
  name = "Z-score",                       
  
  # Clustering
  cluster_rows    = T,
  cluster_columns = F,  
  show_row_dend   = F,
  
  # Labels and side
  row_names_side     = "left",
  column_names_side  = "bottom",          
  show_column_names  = TRUE,              
  column_names_gp    = gpar(fontsize = 10),
  row_names_gp       = gpar(fontsize = 10),
  #width = unit(10, "cm"),
  
  # Annotations
  bottom_annotation = ha_col,            
  
  # Layout tweaks
  heatmap_legend_param = list(
    title = "Z-score",
    title_position = "topcenter"
  )
)

pdf("output/Heatmap DE genes RGFP963.pdf", ht_963, width = 4, height = 6)

draw(
  ht_963,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = T
)

dev.off()


res_RGFP963 <- as.data.frame(res963)

write.csv(res_RGFP963, file = "output/RGFP963_vs_veh.csv")

res_RGFP966 <- as.data.frame(res966)

write.csv(res_RGFP966, file = "output/RGFP966_vs_veh.csv")

saveRDS(dds, "data/DEseq_analysis.rds")
saveRDS(res963, "data/RGFP963 results")
saveRDS(res966, "data/RGFP966 results")


