suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(tidySummarizedExperiment)
  library(ComplexHeatmap)
  library(RColorBrewer)

})


counts <- read.csv("data/norm_count.csv", 
                   row.names = 1)
dim(counts)

sample_data <- read.csv("data/sample_data.csv", 
                        row.names = 1)

counts <- round(counts)  # rounds the values to be used as counts
range(counts)            # min should be 0, no negatives

all.equal(colnames(counts), rownames(sample_data)) # check the sample names are the same for both datasets

# Create the Summarized Experiment 
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  colData = sample_data
)

head(assay(se))
dim(assay(se))
colData(se)

# Filter genes with low number of counts
nrow(se)
se <- se[rowSums(assay(se, "counts")) > 5, ] # Remove genes/rows that do not have > 5 total counts
nrow(se)

# Exploratory analysis ---------------------------------------------------------

# library size
se$libSize <-  colSums(assay(se))

treatment_colors <- c(
  "vehicle"  = "#F5D000",   
  "RGFP963"  = "#96C65C",   
  "RGFP966"  = "#1BB6AF"    
)

# Bar graph of library sizes
libSize <- colData(se) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = row.names(sample_data), y = libSize / 1e6, fill = treatment)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  scale_fill_manual(values = treatment_colors) +
  labs(x = "Sample", 
       y = "Total count in millions", 
       fill = "Treatment")
libSize

ggsave("Library size.pdf", libSize, width = 6, height = 6)

# Create the DESeq data set
dds <- DESeqDataSet(se, design = ~ treatment)
dds <- estimateSizeFactors(dds)

# Plot the size factors against library size and look for any patterns by group
sizeFact <- ggplot(data.frame(libSize = colSums(assay(dds)),
                  sizeFactor = sizeFactors(dds),
                  Group = colnames(dds)),
       aes(x = libSize / 1e6, y = sizeFactor, col = colnames(dds))) + 
  geom_point(size = 5) + 
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(x = "Library size", 
       y = "Size factor", 
       color = "Sample")

sizeFact

ggsave("Size factor vs library size.pdf", sizeFact, width = 6, height = 6)

assay(dds)

# Relationship mean-sd 
meanSdPlot(assay(dds), ranks = FALSE) 

# Data transformation
vsd <- vst(dds, blind = TRUE)

# Relationship mean-sd 
meanSdPlot(assay(vsd), ranks = FALSE)

# Heatmap
dst <- dist(t(assay(vsd)))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

Heatmap(
  as.matrix(dst), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(dst),
  cluster_columns = hclust(dst),
  bottom_annotation = columnAnnotation(
    treatment = dds$treatment,
    col = list(treatment = treatment_colors))
)


# PCA
pcaData <- plotPCA(vsd, intgroup = c("treatment"),
                           returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_plot <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = treatment), size = 5) +
  geom_text(
    aes(label = rownames(pcaData)),
    size = 4,
    vjust = -1.2,          
    hjust = 0.5             
  ) +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  scale_color_manual(values = treatment_colors) +
  labs(color = "Treatment")

PCA_plot
ggsave("PCA.pdf", PCA_plot, width = 8, height = 8)

# descriptive statistics -------------------------------------------------------

desc <-  se %>%
  summarise(mean_size = mean(libSize), 
            sd_size = sd(libSize))
  
limits <- desc %>% 
  mutate(upper_limit = mean_size + 2 * sd_size, 
   lower_limit = mean_size - 2 * sd_size)

se$libSize - limits$upper_limit # m6 library size is higher than 2 sd from the mean library size
limits$lower_limit - se$libSize # There aren't any outliers on the lower side 

# ------------------------------------------------------------------------------

# Removes m6 sample
se2 <-  se[ , colnames(se) != "m6"]
assay(se2)

# Exploratory analysis without m6 ----------------------------------------------
se2$libSize <-  colSums(assay(se2))

# Library size
libSize2 <- colData(se2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>%
  ggplot(aes(x = Sample, y = libSize / 1e6, fill = treatment)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(x = "Sample", 
       y = "Total count in millions", 
       fill = "Treatment") + 
  scale_fill_manual(values = treatment_colors)

libSize2
ggsave("Library size without outlier.pdf", libSize2, width = 6, height = 6)

# Create the DESeq data set
dds2 <- DESeqDataSet(se2, design = ~ treatment)
dds2 <- estimateSizeFactors(dds2)


# Size factor vs library size
sizeFact2 <- ggplot(data.frame(libSize = colSums(assay(dds2)),
                  sizeFactor = sizeFactors(dds2),
                  Group = colnames(dds2)),
       aes(x = libSize / 1e6, y = sizeFactor, col = colnames(dds2))) + 
  geom_point(size = 5) + 
  theme_minimal() + 
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(x = "Library size", 
       y = "Size factor", 
       color = "Sample")

sizeFact2

ggsave("Size factor vs library size without outlier.pdf", sizeFact2, width = 6, height = 6)


assay(dds2)

# Relationship mean-sd 
meanSdPlot(assay(dds2), ranks = FALSE) 

# Pruebo transformar los datos 
vsd2 <- vst(dds2, blind = TRUE)

# Relationship mean-sd 
meanSdPlot(assay(vsd2), ranks = FALSE)

# Heatmap
dst2 <- dist(t(assay(vsd2)))
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
Heatmap(
  as.matrix(dst2), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(dst2),
  cluster_columns = hclust(dst2),
  bottom_annotation = columnAnnotation(
    treatment = dds2$treatment,
    col = list(treatment = treatment_colors))
)

# PCA
pcaData2 <- plotPCA(vsd2, intgroup = c("treatment"),
                           returnData = TRUE)

percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
PCA_plot2 <- ggplot(pcaData2, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = treatment), size = 5) +
  geom_text(
    aes(label = rownames(pcaData2)),
    size = 3.5,
    vjust = -1.2,           
    hjust = 0.5             
  ) +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  scale_color_manual(values = treatment_colors) +
  labs(color = "Treatment")

PCA_plot2
ggsave("PCA without outlier.pdf", PCA_plot2, width = 8, height = 8)


saveRDS(se2, "data/norm_data_se2.rds")
