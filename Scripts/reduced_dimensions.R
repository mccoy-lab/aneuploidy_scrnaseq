list_of_packages <- c("biomaRt", "BiocStyle", "broom", "cowplot", "data.table", "devtools", 
                      "ggdendro", "ggrepel", "gmodels", "gplots", "gridExtra", "lme4", "MAST", "metap", 
                      "MultiAssayExperiment", "mppa", "plyr", "readxl", "Rtsne", "scales", "scater", 
                      "scran", "stringr", "tidyr", "TreeBH", "tools", "umap", "zoo", 
                      "scploid", "dplyr", "margins")

# Source my collection of functions
source("~/Desktop/AneuploidyProject/UsefulFunctions.R")

install_and_load_packages(list_of_packages)

# Resolve masked functions
# getCounts <- scploid::getCounts

# Ensure clean R environment
rm(list=ls())
# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE) 

### dimension reduction visualization

library(monocle3)
library(SummarizedExperiment)
library(MultiAssayExperiment)

emtab3929 <- readRDS("~/Desktop/AneuploidyProject/RawData/emtab3929/EMTAB3929.rds")
emtab3929_gene <- experiments(emtab3929)[["gene"]]
assays(emtab3929_gene)[["count"]][1:10, 1:10]

expr_matrix <- assays(emtab3929_gene)[["count"]]
colnames(expr_matrix) <- gsub("_", ".", colnames(expr_matrix))
expr_matrix <- expr_matrix[, unique(results$cell)]

sample_sheet <- results[!duplicated(cell)] %>%
  setnames(., "cell", "sampleNames") %>%
  as.data.frame()
rownames(sample_sheet) <- sample_sheet$sampleNames
gene_annotation <- data.table(gene_short_name = rowRanges(emtab3929@ExperimentList[[1]])$symbol,
                              ensembl_id = rowRanges(emtab3929@ExperimentList[[1]])$gene) %>% 
  as.data.frame()
rownames(gene_annotation) <- rownames(expr_matrix)

cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = sample_sheet, 
                         gene_metadata = gene_annotation)

rm(emtab3929)
rm(emtab3929_gene)
rm(expr_matrix)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "embryo", residual_model_formula_str = "~ mapped_reads")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "tSNE", max_components = 3)
cds <- reduce_dimension(cds, reduction_method = "UMAP", max_components = 3)

## Step 4: Plot the data
cds$lineage <- factor(cds$lineage, levels = c("Undefined", "ICM", "Trophectoderm", "Intermediate", "Epiblast", "Primitive Endoderm"))

umap_12 <- plot_cells(cds, x = 1, y = 2, color_cells_by = "lineage", cell_size = 1, 
                      label_groups_by_cluster = FALSE, show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Dark2")

umap_23 <- plot_cells(cds, x = 2, y = 3, color_cells_by = "lineage", cell_size = 1, 
                      label_groups_by_cluster = FALSE, show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2")

plot_grid(umap_12, umap_13, rel_widths = c(0.7, 1))

cds$sig_cell <- factor(cds$sig_cell, levels = c("0", "1"))

cds$is_aneuploid <- as.character(cds$sig_cell)
cds$is_aneuploid <- revalue(cds$is_aneuploid, c("0" = "Euploid", "1" = "Aneuploid"))

umap_aneuploid_12 <- plot_cells(cds, x = 1, y = 2, color_cells_by = "is_aneuploid", cell_size = 1, 
                             label_groups_by_cluster = FALSE, show_trajectory_graph = FALSE,
                             label_cell_groups = FALSE) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#4e79a7", "#f28e2b"))

umap_aneuploid_23 <- plot_cells(cds, x = 2, y = 3, color_cells_by = "is_aneuploid", cell_size = 1, 
                                label_groups_by_cluster = FALSE, show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#4e79a7", "#f28e2b"))

grid_left <- plot_grid(umap_12, umap_aneuploid_12, by_celltype_plot, align = "v", axis = "lr", nrow = 3, 
                       labels = c("A", "C", "E"), rel_heights = c(0.8, 0.8, 1))
grid_right <- plot_grid(umap_23, umap_aneuploid_23, enrichment_plot, align = "v", axis = "lr", nrow = 3,
                       labels = c("B", "D", "F"), rel_heights = c(0.8, 0.8, 1))
plot_grid(grid_left, grid_right, ncol = 2, rel_widths = c(0.6, 1))

plot_cells_3d(cds, color_cells_by = "lineage", cell_size = 50, 
              color_palette = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'),
              show_trajectory_graph = FALSE)
