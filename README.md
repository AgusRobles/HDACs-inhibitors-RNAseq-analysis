# HDACs-inhibitors-RNAseq-analysis
RNAseq analysis workflow in R using Bioconductor. The project re-analyzes mouse brain transcriptomic data from Rumbaugh, Gavin et al. “Pharmacological Selectivity Within Class I Histone Deacetylases Predicts Effects on Synaptic Function and Memory Rescue.” as a bioinformatics learning exercise.

## Disclaimer

This repository was created as a bioinformatics learning and portfolio exercise.

The RNA-seq analysis is based on data reported in the paper Rumbaugh, Gavin et al. “Pharmacological Selectivity Within Class I Histone Deacetylases Predicts Effects on Synaptic Function and Memory Rescue.” Neuropsychopharmacology : official publication of the American College of Neuropsychopharmacology vol. 40,10 (2015): 2307-16. doi:10.1038/npp.2015.93

The expression table provided in the publication contains normalized values rather than raw sequencing counts. 
For the purpose of practicing a typical RNA-seq workflow with DESeq2, these values were rounded and treated as counts.

Because of this preprocessing step, the differential expression and enrichment results presented here should not be considered a faithful reproduction of the original study and are intended only for methodological demonstration.
