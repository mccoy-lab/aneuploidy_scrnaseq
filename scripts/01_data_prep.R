list_of_packages <- c("BiocStyle", "biomaRt", "devtools", "dplyr", "here", "MultiAssayExperiment", 
                      "readxl", "scran", "tidyr")

# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE)

# install and load packages
# https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
install.packages.auto <- function(x) { 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    eval(parse(text = sprintf("BiocManager::install(\"%s\")", x, update = FALSE)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

lapply(list_of_packages, function(x) {message(x); install.packages.auto(x)})

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

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.6.1 (2019-07-05)
# os       CentOS Linux 7 (Core)       
# system   x86_64, linux-gnu           
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       America/New_York            
# date     2020-04-22                  
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source        
# AnnotationDbi          1.48.0   2019-10-29 [1] Bioconductor  
# askpass                1.1      2019-01-13 [2] CRAN (R 3.6.1)
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 3.6.1)
# backports              1.1.5    2019-10-02 [1] CRAN (R 3.6.1)
# beeswarm               0.2.3    2016-04-25 [1] CRAN (R 3.6.1)
# Biobase              * 2.46.0   2019-10-29 [1] Bioconductor  
# BiocFileCache          1.10.2   2019-11-08 [1] Bioconductor  
# BiocGenerics         * 0.32.0   2019-10-29 [1] Bioconductor  
# BiocManager            1.30.10  2019-11-16 [1] CRAN (R 3.6.1)
# BiocNeighbors          1.4.1    2019-11-01 [1] Bioconductor  
# BiocParallel         * 1.20.1   2019-12-21 [1] Bioconductor  
# BiocSingular           1.2.1    2019-12-23 [1] Bioconductor  
# BiocStyle            * 2.14.4   2020-01-09 [1] Bioconductor  
# biomaRt              * 2.42.0   2019-10-29 [1] Bioconductor  
# bit                    1.1-14   2018-05-29 [2] CRAN (R 3.6.1)
# bit64                  0.9-7    2017-05-08 [2] CRAN (R 3.6.1)
# bitops                 1.0-6    2013-08-17 [1] CRAN (R 3.6.1)
# blob                   1.2.0    2019-07-09 [2] CRAN (R 3.6.1)
# callr                  3.4.2    2020-02-12 [1] CRAN (R 3.6.1)
# cellranger             1.1.0    2016-07-27 [1] CRAN (R 3.6.1)
# cli                    2.0.1    2020-01-08 [1] CRAN (R 3.6.1)
# colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 3.6.1)
# curl                   4.3      2019-12-02 [1] CRAN (R 3.6.1)
# DBI                    1.1.0    2019-12-15 [1] CRAN (R 3.6.1)
# dbplyr                 1.4.2    2019-06-17 [1] CRAN (R 3.6.1)
# DelayedArray         * 0.12.2   2020-01-06 [1] Bioconductor  
# DelayedMatrixStats     1.8.0    2019-10-29 [1] Bioconductor  
# desc                   1.2.0    2018-05-01 [1] CRAN (R 3.6.1)
# devtools             * 2.2.1    2019-09-24 [1] CRAN (R 3.6.1)
# digest                 0.6.24   2020-02-12 [1] CRAN (R 3.6.1)
# dplyr                * 0.8.4    2020-01-31 [1] CRAN (R 3.6.1)
# dqrng                  0.2.1    2019-05-17 [1] CRAN (R 3.6.1)
# edgeR                  3.28.0   2019-10-29 [1] Bioconductor  
# ellipsis               0.3.0    2019-09-20 [1] CRAN (R 3.6.1)
# evaluate               0.14     2019-05-28 [1] CRAN (R 3.6.1)
# fansi                  0.4.0    2018-10-05 [2] CRAN (R 3.6.1)
# fs                     1.3.1    2019-05-06 [1] CRAN (R 3.6.1)
# GenomeInfoDb         * 1.22.0   2019-10-29 [1] Bioconductor  
# GenomeInfoDbData       1.2.2    2019-12-07 [1] Bioconductor  
# GenomicRanges        * 1.38.0   2019-10-29 [1] Bioconductor  
# ggbeeswarm             0.6.0    2017-08-07 [1] CRAN (R 3.6.1)
# ggplot2                3.2.1    2019-08-10 [1] CRAN (R 3.6.1)
# glue                   1.3.1    2019-03-12 [2] CRAN (R 3.6.1)
# gridExtra              2.3      2017-09-09 [1] CRAN (R 3.6.1)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
# here                 * 0.1      2017-05-28 [1] CRAN (R 3.6.1)
# hms                    0.5.3    2020-01-08 [1] CRAN (R 3.6.1)
# htmltools              0.4.0    2019-10-04 [1] CRAN (R 3.6.1)
# httr                   1.4.1    2019-08-05 [1] CRAN (R 3.6.1)
# igraph                 1.2.4.2  2019-11-27 [1] CRAN (R 3.6.1)
# IRanges              * 2.20.2   2020-01-13 [1] Bioconductor  
# irlba                  2.3.3    2019-02-05 [1] CRAN (R 3.6.1)
# knitr                  1.28     2020-02-06 [1] CRAN (R 3.6.1)
# lattice                0.20-38  2018-11-04 [2] CRAN (R 3.6.1)
# lazyeval               0.2.2    2019-03-15 [2] CRAN (R 3.6.1)
# lifecycle              0.1.0    2019-08-01 [2] CRAN (R 3.6.1)
# limma                  3.42.2   2020-02-03 [1] Bioconductor  
# locfit                 1.5-9.1  2013-04-20 [2] CRAN (R 3.6.1)
# magrittr               1.5      2014-11-22 [2] CRAN (R 3.6.1)
# Matrix                 1.2-17   2019-03-22 [2] CRAN (R 3.6.1)
# matrixStats          * 0.55.0   2019-09-07 [1] CRAN (R 3.6.1)
# memoise                1.1.0    2017-04-21 [2] CRAN (R 3.6.1)
# MultiAssayExperiment * 1.12.2   2020-01-21 [1] Bioconductor  
# munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
# openssl                1.4.1    2019-07-18 [2] CRAN (R 3.6.1)
# pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.1)
# pkgbuild               1.0.6    2019-10-09 [1] CRAN (R 3.6.1)
# pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.1)
# pkgload                1.0.2    2018-10-29 [1] CRAN (R 3.6.1)
# plyr                   1.8.5    2019-12-10 [1] CRAN (R 3.6.1)
# prettyunits            1.0.2    2015-07-13 [2] CRAN (R 3.6.1)
# processx               3.4.2    2020-02-09 [1] CRAN (R 3.6.1)
# progress               1.2.2    2019-05-16 [1] CRAN (R 3.6.1)
# ps                     1.3.1    2020-02-12 [1] CRAN (R 3.6.1)
# purrr                  0.3.3    2019-10-18 [1] CRAN (R 3.6.1)
# R6                     2.4.1    2019-11-12 [1] CRAN (R 3.6.1)
# rappdirs               0.3.1    2016-03-28 [1] CRAN (R 3.6.1)
# Rcpp                   1.0.3    2019-11-08 [1] CRAN (R 3.6.1)
# RCurl                  1.98-1.1 2020-01-19 [1] CRAN (R 3.6.1)
# readxl               * 1.3.1    2019-03-13 [1] CRAN (R 3.6.1)
# remotes                2.1.0    2019-06-24 [1] CRAN (R 3.6.1)
# reshape2               1.4.3    2017-12-11 [2] CRAN (R 3.6.1)
# rlang                  0.4.4    2020-01-28 [1] CRAN (R 3.6.1)
# rmarkdown              2.1      2020-01-20 [1] CRAN (R 3.6.1)
# rprojroot              1.3-2    2018-01-03 [1] CRAN (R 3.6.1)
# RSQLite                2.2.0    2020-01-07 [1] CRAN (R 3.6.1)
# rstudioapi             0.11     2020-02-07 [1] CRAN (R 3.6.1)
# rsvd                   1.0.2    2019-07-29 [1] CRAN (R 3.6.1)
# S4Vectors            * 0.24.3   2020-01-18 [1] Bioconductor  
# scales                 1.1.0    2019-11-18 [1] CRAN (R 3.6.1)
# scater                 1.14.6   2019-12-16 [1] Bioconductor  
# scran                * 1.14.6   2020-02-03 [1] Bioconductor  
# sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 3.6.1)
# SingleCellExperiment * 1.8.0    2019-10-29 [1] Bioconductor  
# statmod                1.4.33   2020-01-10 [1] CRAN (R 3.6.1)
# stringi                1.4.3    2019-03-12 [2] CRAN (R 3.6.1)
# stringr                1.4.0    2019-02-10 [2] CRAN (R 3.6.1)
# SummarizedExperiment * 1.16.1   2019-12-19 [1] Bioconductor  
# testthat               2.3.1    2019-12-01 [1] CRAN (R 3.6.1)
# tibble                 2.1.3    2019-06-06 [1] CRAN (R 3.6.1)
# tidyr                * 1.0.2    2020-01-24 [1] CRAN (R 3.6.1)
# tidyselect             0.2.5    2018-10-11 [2] CRAN (R 3.6.1)
# usethis              * 1.5.1    2019-07-04 [1] CRAN (R 3.6.1)
# vctrs                  0.2.2    2020-01-24 [1] CRAN (R 3.6.1)
# vipor                  0.4.5    2017-03-22 [1] CRAN (R 3.6.1)
# viridis                0.5.1    2018-03-29 [1] CRAN (R 3.6.1)
# viridisLite            0.3.0    2018-02-01 [2] CRAN (R 3.6.1)
# withr                  2.1.2    2018-03-15 [2] CRAN (R 3.6.1)
# xfun                   0.12     2020-01-13 [1] CRAN (R 3.6.1)
# XML                    3.99-0.3 2020-01-20 [1] CRAN (R 3.6.1)
# XVector                0.26.0   2019-10-29 [1] Bioconductor  
# yaml                   2.2.1    2020-02-01 [1] CRAN (R 3.6.1)
# zlibbioc               1.32.0   2019-10-29 [1] Bioconductor  
# 
# [1] /home-net/home-4/rmccoy22@jhu.edu/R/x86_64-pc-linux-gnu-library/3.6/gcc/5.5
# [2] /software/apps/R/3.6.1/gcc/5.5.0/lib64/R/library
