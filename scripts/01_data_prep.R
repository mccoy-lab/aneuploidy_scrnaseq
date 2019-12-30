list_of_packages <- c("BiocStyle", "biomaRt", "devtools", "dplyr", "gplots", "gridExtra", 
                      "here", "lme4", "MultiAssayExperiment", "readxl", "Rtsne", "scater", 
                      "scploid", "scran", "tidyr", "umap")

# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE)

# Source my collection of functions
# Easily install and load packages
install_and_load_packages <- function(pkg){
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
  if(length(new_packages))
    install.packages(new_packages, dependencies = TRUE)
  sapply(list_of_packages, require, character.only = TRUE)
}

install_and_load_packages(list_of_packages)

# Data were acquired for the human embryo (EMTAB3929) scRNA-seq aneuploidy project as follows:  
# (1) EMTAB3929 data (PMID 27062923) were downloaded on 11/6/2018 from http://imlspenticton.uzh.ch:3838/conquer/.  
# (a) MultiAssay Experiment (EMTAB3929.rds)  
# (b) MultiQC report  
# (c) Scater report  
# (d) Salmon archive (EMTAB3929_salmo.tar and EMTAB3929 folder)  
# (2) Supplementary files from Griffiths et al., 2017 were forked on 11/06/2018 from MarioniLab/Aneuploidy2017 on Github: https://github.com/MarioniLab/Aneuploidy2017.  
# (3) Data from Griffiths et al., 2017 were downloaded on 11/14/2018 using the shell script supplied in their supplementary files (sh get_data.sh)  

# Load EMTAB3929 data
emtab3929_meta <- readRDS(here("RawData/EMTAB3929.rds"))

# (1) Remove version numbers from Ensembl gene IDs
rownames(emtab3929_meta@ExperimentList@listData$gene@assays$data$count) <- sapply(strsplit(rownames(emtab3929_meta@ExperimentList@listData$gene@assays$data$count), "\\."), `[`, 1)
dim(emtab3929_meta@ExperimentList@listData$gene@assays$data$count) # 65218  1529

# (2) Obtain annotation for reference genome GRCh38.84 (which is GRCh38.p5 Ensembl 84:Mar 2016); keep autosomal genes only
human_ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         host = "mar2016.archive.ensembl.org", 
                         path = "/biomart/martservice", 
                         dataset = "hsapiens_gene_ensembl")

annotation <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "gene_biotype"),
  mart = human_ensembl,
  values = as.character(rownames(emtab3929_meta@ExperimentList@listData$gene@assays$data$count)),
  filters = "ensembl_gene_id"
)
annotation <- annotation [annotation $chromosome_name %in% 1:22, ] # 56400 genes

# (3) Cell lineage information. Compile sample metasheet.
cell_lineage_data <- read_xlsx(here("external_metadata/stirparo2018_tableS4.xlsx"), sheet = 1)
cell_lineage_data <- cell_lineage_data[cell_lineage_data$Study == "Petropoulos et al., 2016 (ERP012552)", ] # 1,481 cells
cell_lineage_data$Cell <- gsub("_", ".", cell_lineage_data$Cell)
cell_lineage_data$EStage <- sapply(strsplit(cell_lineage_data$Embryo, "_"), "[", 1)

metasheet <- cell_lineage_data[, c(2:6, 8:10)]
colnames(metasheet)[2] <- "Sample"
metasheet$Cell <- sapply(strsplit(metasheet$Sample, "\\."), tail, n = 1)
rownames(metasheet) <- metasheet$Sample
metasheet <- metasheet[, c(2, 8, 3, 1, 9, 4:7)]

salmon_summary <- emtab3929_meta@metadata$salmon_summary[, c(1, 6:8)]
salmon_summary$sample <- gsub("_", ".", salmon_summary$sample)

metasheet <- inner_join(metasheet, salmon_summary, by = c("Sample" = "sample"))
colnames(metasheet)[colnames(metasheet) == "num_processed"] <- "Processed Reads" 
colnames(metasheet)[colnames(metasheet) == "num_mapped"] <- "Mapped Reads"
colnames(metasheet)[colnames(metasheet) == "percent_mapped"] <- "Percent Mapped"
metasheet$`Revised lineage (this study)` <- gsub("epiblast", "Epiblast", metasheet$`Revised lineage (this study)`)
metasheet$`Revised lineage (this study)` <- gsub("Inner cell mass", "ICM", metasheet$`Revised lineage (this study)`)
metasheet$`Revised lineage (this study)` <- gsub("intermediate", "Intermediate", metasheet$`Revised lineage (this study)`)
metasheet$`Revised lineage (this study)` <- gsub("primitive_endoderm", "Primitive Endoderm", metasheet$`Revised lineage (this study)`)
metasheet$`Revised lineage (this study)` <- gsub("trophectoderm", "Trophectoderm", metasheet$`Revised lineage (this study)`)
metasheet$`Revised lineage (this study)` <- gsub("undefined", "Undefined", metasheet$`Revised lineage (this study)`)

# (4) Modify EMTAB gene expression matrix to contain only information for genes in `annotation` and samples with cell lineage information, and then convert counts into CPM and apply gene expression filter
emtab3929_counts <- emtab3929_meta@ExperimentList@listData$gene@assays$data$count[annotation$ensembl_gene_id, ]
colnames(emtab3929_counts) <- gsub("_", ".", colnames(emtab3929_counts))
emtab3929_counts <- emtab3929_counts[, colnames(emtab3929_counts) %in% metasheet$Sample]

emtab3929_cpm <- edgeR::cpm(emtab3929_counts, normalized.lib.sizes = TRUE, log = FALSE)
emtab3929_cpm <- emtab3929_cpm[, colnames(emtab3929_cpm) %in% metasheet$Sample]
dim(emtab3929_cpm) # 56,400 genes and 1,481 cells
emtab3929_log2cpm <- log2(emtab3929_cpm + 1)

# Apply gene filter used in Griffiths analysis.
gene_filter <- apply(emtab3929_cpm, 1, median) > 50
filtered_cpm <- emtab3929_cpm[gene_filter, ] 
dim(filtered_cpm) # 2,991 genes and 1,481 cells
filtered_log2cpm <- log2(filtered_cpm + 1)

filtered_counts <- emtab3929_counts[rownames(emtab3929_counts) %in% rownames(filtered_cpm), 
                                    colnames(emtab3929_counts) %in% colnames(filtered_cpm)]
dim(filtered_counts) # 2,991 genes and 1,481 cells

# Save objects for easy uploading in the future.

save(metasheet, emtab3929_counts, filtered_counts, emtab3929_cpm, filtered_cpm, filtered_log2cpm, annotation, 
     file = here("ProcessedData/EMTAB3929_DataPrep.RData"))

# To keep all analyses consistent, the following modifications to the original EMTAB3929 data were made:  
# (1) Removed version numbers from EMTAB3929 Ensembl IDs  
# (2) Included GRCh38.p5 Ensembl 84:Mar2016 (same as GRCh38.84) reference gene annotation. Kept only information for autosomal        genes  
# (3) Included cell lineage information from Stirparo et al., 2018 (Table S4)  
# (4) Applied gene filter of median CPM > 50 used in Griffiths analysis  
# **The final data used in downstream analyses contain 2,991 genes and 1,481 cells from a total of 88 embryos.**

devtools::session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.6.1 (2019-07-05)
# os       macOS Mojave 10.14.6        
# system   x86_64, darwin15.6.0        
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       America/New_York            
# date     2019-12-30                  
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date       lib source                                    
# AnnotationDbi          1.48.0    2019-10-29 [1] Bioconductor                              
# askpass                1.1       2019-01-13 [1] CRAN (R 3.6.0)                            
# assertthat             0.2.1     2019-03-21 [1] CRAN (R 3.6.0)                            
# backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.0)                            
# beeswarm               0.2.3     2016-04-25 [1] CRAN (R 3.6.0)                            
# Biobase              * 2.46.0    2019-10-29 [1] Bioconductor                              
# BiocFileCache          1.10.2    2019-11-08 [1] Bioconductor                              
# BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor                              
# BiocManager            1.30.10   2019-11-16 [1] CRAN (R 3.6.0)                            
# BiocNeighbors          1.4.1     2019-11-01 [1] Bioconductor                              
# BiocParallel         * 1.20.0    2019-10-29 [1] Bioconductor                              
# BiocSingular           1.2.0     2019-10-29 [1] Bioconductor                              
# BiocStyle            * 2.14.0    2019-10-29 [1] Bioconductor                              
# biomaRt              * 2.42.0    2019-10-29 [1] Bioconductor                              
# bit                    1.1-14    2018-05-29 [1] CRAN (R 3.6.0)                            
# bit64                  0.9-7     2017-05-08 [1] CRAN (R 3.6.0)                            
# bitops                 1.0-6     2013-08-17 [1] CRAN (R 3.6.0)                            
# blob                   1.2.0     2019-07-09 [1] CRAN (R 3.6.0)                            
# boot                   1.3-23    2019-07-05 [1] CRAN (R 3.6.0)                            
# callr                  3.3.2     2019-09-22 [1] CRAN (R 3.6.0)                            
# caTools                1.17.1.2  2019-03-06 [1] CRAN (R 3.6.0)                            
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.0)                            
# cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.0)                            
# colorspace             1.4-1     2019-03-18 [1] CRAN (R 3.6.0)                            
# crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.0)                            
# curl                   4.2       2019-09-24 [1] CRAN (R 3.6.0)                            
# DBI                    1.0.0     2018-05-02 [1] CRAN (R 3.6.0)                            
# dbplyr                 1.4.2     2019-06-17 [1] CRAN (R 3.6.0)                            
# DelayedArray         * 0.12.0    2019-10-29 [1] Bioconductor                              
# DelayedMatrixStats     1.8.0     2019-10-29 [1] Bioconductor                              
# desc                   1.2.0     2018-05-01 [1] CRAN (R 3.6.0)                            
# devtools             * 2.2.1     2019-09-24 [1] CRAN (R 3.6.0)                            
# digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.0)                            
# dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.0)                            
# dqrng                  0.2.1     2019-05-17 [1] CRAN (R 3.6.0)                            
# edgeR                  3.28.0    2019-10-29 [1] Bioconductor                              
# ellipsis               0.3.0     2019-09-20 [1] CRAN (R 3.6.0)                            
# evaluate               0.14      2019-05-28 [1] CRAN (R 3.6.0)                            
# fs                     1.3.1     2019-05-06 [1] CRAN (R 3.6.0)                            
# gdata                  2.18.0    2017-06-06 [1] CRAN (R 3.6.0)                            
# GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor                              
# GenomeInfoDbData       1.2.2     2019-11-26 [1] Bioconductor                              
# GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor                              
# ggbeeswarm             0.6.0     2017-08-07 [1] CRAN (R 3.6.0)                            
# ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.0)                            
# glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.0)                            
# gplots               * 3.0.1.1   2019-01-27 [1] CRAN (R 3.6.0)                            
# gridExtra            * 2.3       2017-09-09 [1] CRAN (R 3.6.0)                            
# gtable                 0.3.0     2019-03-25 [1] CRAN (R 3.6.0)                            
# gtools                 3.8.1     2018-06-26 [1] CRAN (R 3.6.0)                            
# here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.0)                            
# hms                    0.5.2     2019-10-30 [1] CRAN (R 3.6.0)                            
# htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.0)                            
# httr                   1.4.1     2019-08-05 [1] CRAN (R 3.6.0)                            
# igraph                 1.2.4.1   2019-04-22 [1] CRAN (R 3.6.0)                            
# IRanges              * 2.20.1    2019-11-20 [1] Bioconductor                              
# irlba                  2.3.3     2019-02-05 [1] CRAN (R 3.6.0)                            
# jsonlite               1.6       2018-12-07 [1] CRAN (R 3.6.0)                            
# KernSmooth             2.23-16   2019-10-15 [1] CRAN (R 3.6.0)                            
# knitr                  1.26      2019-11-12 [1] CRAN (R 3.6.0)                            
# lattice                0.20-38   2018-11-04 [1] CRAN (R 3.6.1)                            
# lazyeval               0.2.2     2019-03-15 [1] CRAN (R 3.6.0)                            
# lifecycle              0.1.0     2019-08-01 [1] CRAN (R 3.6.0)                            
# limma                  3.42.0    2019-10-29 [1] Bioconductor                              
# lme4                 * 1.1-21    2019-03-05 [1] CRAN (R 3.6.0)                            
# locfit                 1.5-9.1   2013-04-20 [1] CRAN (R 3.6.0)                            
# magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.0)                            
# MASS                   7.3-51.4  2019-03-31 [1] CRAN (R 3.6.1)                            
# Matrix               * 1.2-17    2019-03-22 [1] CRAN (R 3.6.1)                            
# matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.0)                            
# memoise                1.1.0     2017-04-21 [1] CRAN (R 3.6.0)                            
# minqa                  1.2.4     2014-10-09 [1] CRAN (R 3.6.0)                            
# MultiAssayExperiment * 1.12.0    2019-10-29 [1] Bioconductor                              
# munsell                0.5.0     2018-06-12 [1] CRAN (R 3.6.0)                            
# nlme                   3.1-142   2019-11-07 [1] CRAN (R 3.6.0)                            
# nloptr                 1.2.1     2018-10-03 [1] CRAN (R 3.6.0)                            
# openssl                1.4.1     2019-07-18 [1] CRAN (R 3.6.0)                            
# pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.0)                            
# pkgbuild               1.0.6     2019-10-09 [1] CRAN (R 3.6.0)                            
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.0)                            
# pkgload                1.0.2     2018-10-29 [1] CRAN (R 3.6.0)                            
# prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.6.0)                            
# processx               3.4.1     2019-07-18 [1] CRAN (R 3.6.0)                            
# progress               1.2.2     2019-05-16 [1] CRAN (R 3.6.0)                            
# ps                     1.3.0     2018-12-21 [1] CRAN (R 3.6.0)                            
# purrr                  0.3.3     2019-10-18 [1] CRAN (R 3.6.0)                            
# R6                     2.4.1     2019-11-12 [1] CRAN (R 3.6.0)                            
# rappdirs               0.3.1     2016-03-28 [1] CRAN (R 3.6.0)                            
# Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.0)                            
# RCurl                  1.95-4.12 2019-03-04 [1] CRAN (R 3.6.0)                            
# readxl               * 1.3.1     2019-03-13 [1] CRAN (R 3.6.0)                            
# remotes                2.1.0     2019-06-24 [1] CRAN (R 3.6.0)                            
# reticulate             1.13      2019-07-24 [1] CRAN (R 3.6.0)                            
# rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.0)                            
# rmarkdown              1.17      2019-11-13 [1] CRAN (R 3.6.0)                            
# rprojroot              1.3-2     2018-01-03 [1] CRAN (R 3.6.0)                            
# RSQLite                2.1.2     2019-07-24 [1] CRAN (R 3.6.0)                            
# rstudioapi             0.10      2019-03-19 [1] CRAN (R 3.6.0)                            
# rsvd                   1.0.2     2019-07-29 [1] CRAN (R 3.6.0)                            
# Rtsne                * 0.15      2018-11-10 [1] CRAN (R 3.6.0)                            
# S4Vectors            * 0.24.0    2019-10-29 [1] Bioconductor                              
# scales                 1.1.0     2019-11-18 [1] CRAN (R 3.6.0)                            
# scater               * 1.14.4    2019-11-18 [1] Bioconductor                              
# scploid              * 0.9       2019-11-26 [1] Github (MarioniLab/Aneuploidy2017@286c064)
# scran                * 1.14.5    2019-11-19 [1] Bioconductor                              
# sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.6.0)                            
# SingleCellExperiment * 1.8.0     2019-10-29 [1] Bioconductor                              
# statmod                1.4.32    2019-05-29 [1] CRAN (R 3.6.0)                            
# stringi                1.4.3     2019-03-12 [1] CRAN (R 3.6.0)                            
# stringr                1.4.0     2019-02-10 [1] CRAN (R 3.6.0)                            
# SummarizedExperiment * 1.16.0    2019-10-29 [1] Bioconductor                              
# testthat               2.3.0     2019-11-05 [1] CRAN (R 3.6.0)                            
# tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.0)                            
# tidyr                * 1.0.0     2019-09-11 [1] CRAN (R 3.6.0)                            
# tidyselect             0.2.5     2018-10-11 [1] CRAN (R 3.6.0)                            
# umap                 * 0.2.3.1   2019-08-21 [1] CRAN (R 3.6.0)                            
# usethis              * 1.5.1     2019-07-04 [1] CRAN (R 3.6.0)                            
# vctrs                  0.2.0     2019-07-05 [1] CRAN (R 3.6.0)                            
# vipor                  0.4.5     2017-03-22 [1] CRAN (R 3.6.0)                            
# viridis                0.5.1     2018-03-29 [1] CRAN (R 3.6.0)                            
# viridisLite            0.3.0     2018-02-01 [1] CRAN (R 3.6.0)                            
# withr                  2.1.2     2018-03-15 [1] CRAN (R 3.6.0)                            
# xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.0)                            
# XML                    3.98-1.20 2019-06-06 [1] CRAN (R 3.6.0)                            
# XVector                0.26.0    2019-10-29 [1] Bioconductor                              
# yaml                   2.2.0     2018-07-25 [1] CRAN (R 3.6.0)                            
# zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.0)                            
# zlibbioc               1.32.0    2019-10-29 [1] Bioconductor                              
# 
# [1] /Library/Frameworks/R.framework/Versions/3.6/Resources/library