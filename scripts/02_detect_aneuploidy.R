list_of_packages <- c("BiocStyle", "biomaRt", "broom", "cowplot", "data.table", "devtools", 
                      "dplyr", "ggdendro", "ggrepel", "gmodels", "gplots", "gridExtra", 
                      "here", "lme4", "margins", "MAST", "metap", "mixtools", "monocle3", 
                      "mppa", "MultiAssayExperiment", "MultiAssayExperiment", "plyr", 
                      "readxl", "Rtsne", "scales", "scater", "scploid", "scran", "SMITE", 
                      "stringr", "SummarizedExperiment", "survcomp", "tidyr", "tools", 
                      "TreeBH", "umap", "zoo")

# Easily install and load packages
install_and_load_packages <- function(pkg){
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
  if(length(new_packages))
    install.packages(new_packages, dependencies = TRUE)
  sapply(list_of_packages, require, character.only = TRUE)
}

install_and_load_packages(list_of_packages)

# over-ride masked functions
here <- here::here
summarize <- dplyr::summarize

# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE) 

# Load processed EMTAB3929 data and check data
load(here("ProcessedData/EMTAB3929_DataPrep.RData"))
dim(annotation) # 56,400 genes
dim(filtered_counts) # 2,991 genes and 1,481 cells
dim(filtered_cpm) # 2,991 genes and 1,481 cells
dim(filtered_log2cpm) # 2,991 genes and 1,481 cells
dim(metasheet) # 1,481 samples
length(unique(metasheet$Embryo)) # 88 embryos

emtab3929_counts <- emtab3929_counts[, colnames(emtab3929_counts) %in% metasheet$Sample]
emtab3929_counts <- emtab3929_counts[order(rownames(emtab3929_counts)), order(colnames(emtab3929_counts))]
annotation <- annotation[order(annotation$ensembl_gene_id),] %>%
  as.data.table()
metasheet <- metasheet[order(metasheet$Sample),]

# to split data, use annotations from Stirparo et al. 2018
metasheet$Group <- paste0(metasheet$EStage, "_", metasheet$`Revised lineage (this study)`)
metasheet <- metasheet %>%
  as.data.table()
metasheet[, qc_pass := (`Percent Mapped` > quantile(`Percent Mapped`, 0.1) & `Mapped Reads` > quantile(`Mapped Reads`, 0.1))]
metasheet[, Embryo := gsub("_", ".", Embryo)]

## create scploid object
ploidytest_dt <- makeAneu(counts = emtab3929_counts[, metasheet[qc_pass == TRUE]$Sample],
                          genes = annotation$ensembl_gene_id,
                          chrs = annotation$chromosome_name,
                          cellNames = metasheet[qc_pass == TRUE]$Sample,
                          cellGroups = metasheet[qc_pass == TRUE]$Group) # split data by EStage and cell type

# data split into 13 subsets, one for each EStage_celltype combination
spt <- splitCellsByGroup(ploidytest_dt) 

# drop groups that fail scploid QC
# excluded stages include E3, where MZT is not yet complete
ploidytest_qc <- do.call(rbind, lapply(spt, assessMetrics)) %>%
  as.data.table(keep.rownames = TRUE)
groups_qc_pass <- ploidytest_qc[!grepl("Poor", residuals) & !grepl("Poor", ngenes) & !grepl("Poor", zeros)]$rn

qc_melt <- melt(ploidytest_qc, id.vars = c("rn")) 
qc_melt[, c("stage", "cell_type") := tstrsplit(rn, "_", fixed = TRUE)]
qc_melt[, group :=  paste(stage, cell_type, sep = " ")]
setorder(qc_melt, stage, cell_type)
qc_melt$group <- factor(qc_melt$group, levels = rev(unique(qc_melt$group)))
qc_melt$value <- factor(qc_melt$value, levels = c("Poor quality", "Exercise Caution", "Good quality"))
qc_melt$cell_type <- factor(qc_melt$cell_type, levels = rev(c("Undefined", "ICM", "Trophectoderm", "Epiblast", "Primitive Endoderm", "Intermediate")))
qc_melt[variable == "residuals", variable := "resid."]
qc_melt[variable == "ngenes", variable := "genes"]

qc_fig <- ggplot(data = qc_melt, aes(x = variable, y = cell_type, fill = value)) +
  theme_minimal() +
  geom_tile(size = 2, color = "white") +
  scale_fill_manual(values = c("#e41a1c", "#e6ab02", "#1b9e77"), name = "") +
  xlab("QC metric") +
  ylab("") +
  facet_grid(. ~ stage)

# visualize with PCA
pc <- prcomp(t(filtered_cpm), center = TRUE, scale = TRUE)
eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
pca_dims <- data.table(colnames(filtered_cpm), pc$x[, 1:4]) %>%
  setnames(., c("Sample", "PC1", "PC2", "PC3", "PC4"))
pca_dims <- merge(pca_dims, metasheet, "Sample")
pca_dims$`Revised lineage (this study)` <- factor(pca_dims$`Revised lineage (this study)`, levels = unique(pca_dims$`Revised lineage (this study)`))
pca_dims$EStage <- factor(pca_dims$EStage, levels = unique(pca_dims$EStage))

pca_fig <- ggplot(data = pca_dims[qc_pass == TRUE], aes(x = PC1, y = PC3, color = `Revised lineage (this study)`, shape = EStage)) +
  geom_point(size = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'), name = "Cell type", drop = FALSE) +
  scale_shape_discrete(name = "Stage", drop = FALSE) +
  theme_classic()

plot_grid(pca_fig, qc_fig, labels = c("A", "B"), ncol = 1, rel_heights = c(1, 0.75))

# subset data to stage/cell-type groups that pass QC
spt <- spt[names(spt) %in% groups_qc_pass]

# run scploid
expression_results <- do.call(rbind, lapply(spt, calcAneu)) %>%
  as.data.table()
expression_results[, chr := paste0("chr", chr)]

# add cell type information from Stirparo et al.
metasheet <- metasheet %>% 
  as.data.table()
metadata <- metasheet[, c("Sample", "Embryo", "EStage", "Stage", "Revised lineage (this study)")] %>%
  setnames(., c("cell", "embryo", "EStage", "stage", "lineage"))
expression_results <- merge(expression_results, metadata, "cell")
setnames(expression_results, "z", "scploid_z")
setnames(expression_results, "score", "scploid_score")
setnames(expression_results, "p", "scploid_p")

## add ASE data

if (file.exists(here("results/ase_by_chr.txt"))) {
  ase_by_chr <- fread(here("results/ase_by_chr.txt"))
} else {
  file_list <- list.files(here("results/ase_tables"), pattern = "*.table", full.names = TRUE)
  mappings <- fread(here("results/ase_tables/E-MTAB-3929.sdrf.txt"))[, c(1, 31)] %>%
    setnames(., c("cell", "accession"))
  read_ase <- function(file_name, metadata) {
    id <- sub('\\..*$', '', basename(file_name))
    dt <- fread(file_name)
    dt[, cell := mappings[accession == id]$cell]
    dt[, accession := id]
    return(dt)
  }
  ase <- do.call(rbind, lapply(file_list, function(x) read_ase(x, mappings)))
  ase[, snp_id := paste(contig, position, sep = "_")]
  ase[, embryo := sub("^(.*)[.].*", '\\1', cell)]
  ase[, embryo := gsub("_", ".", embryo)]
  ase[, cell := gsub("_", ".", cell)]
  ase[, embryo_snp_id := paste(embryo, snp_id, sep = "_")]
  ase[, minCount := pmin(refCount, altCount)]
  ase[, maxCount := pmax(refCount, altCount)]
  
  # summarize ASE per cell-chromosome
  ase_by_chr <- group_by(ase, cell, contig) %>%
    summarize(., allelic_ratio = sum(minCount) / sum(totalCount), 
              min_count_sum = sum(minCount), 
              total_reads = sum(totalCount)) %>%
    as.data.table()
  fwrite(ase_by_chr, file = here("results/ase_by_chr.txt"), 
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# merge with expression data
ase_by_chr[, chrom := paste(cell, contig, sep = ".")]
expression_results[, chrom := paste(cell, chr, sep = ".")]
results <- merge(expression_results, 
                 ase_by_chr[, c("chrom", "allelic_ratio", "total_reads")], 
                 "chrom", all.x = TRUE)

# get mapped reads metadata
coverage_metadata <- metasheet[, c("Sample", "Mapped Reads")] %>%
  setnames(., c("cell", "mapped_reads"))
results <- merge(results, coverage_metadata, "cell")

# if no SNPs discovered for a chromosome (despite sufficient coverage), set allelic ratio to 0
results[is.na(allelic_ratio), allelic_ratio := 0]
results[is.na(allelic_ratio), total_reads := 0]

# correct allelic imbalance based on converage
results[, resid_allelic_ratio := resid(lm(data = results, formula = allelic_ratio ~ mapped_reads))]

# estimate variance and compute allelic imbalance z-scores
ase_iqr <- quantile(results$allelic_ratio, c(0.25, 0.75))
iqr_indices <- which(results$allelic_ratio > ase_iqr[1] & results$allelic_ratio < ase_iqr[2])
m_iqr <- results$resid_allelic_ratio[iqr_indices]
null_allelic_ratio <- mean(m_iqr)
null_var_allelic_ratio <- {IQR(results$resid_allelic_ratio)/(2 * qnorm(.75))} ^ 2
results[, ase_z := (resid_allelic_ratio - null_allelic_ratio) / sqrt(null_var_allelic_ratio)]
results[, ase_p := pnorm(ase_z)]

# examine overlap between two signatures of aneuploidy
cor.test(results$scploid_p, results$ase_p, method = "kendall")

scploid_ase_overlap <- function(pval_threshold) {
  fisher_exact_results <- group_by(results, ase_p < pval_threshold, scploid_p < pval_threshold) %>% 
    summarize(., n()) %>% 
    as.data.table() %>%
    setnames(., c("ase", "expression", "n")) %>%
    dcast.data.table(formula = ase ~ expression, value.var = "n") %>%
    dplyr::select("FALSE", "TRUE") %>%
    fisher.test() %>%
    tidy() %>%
    mutate(., threshold = pval_threshold) %>%
    as.data.table()
  return(fisher_exact_results)
}

enrichment_results <- do.call(rbind, lapply(c(3.2e-1, 1e-1, 3.2e-2, 1e-2, 3.2e-3, 1e-3, 
                                              3.2e-4, 1e-4), function(x) try(scploid_ase_overlap(x))))

ggplot(data = enrichment_results, aes(x = threshold, y = estimate, ymin = conf.low, ymax = conf.high)) + 
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_x_log10() +
  geom_hline(yintercept = 1, lty = "dotted") +
  xlab("P-value threshold") +
  ylab("Odds Ratio")

fisher_wrapper <- function(pval_1, pval_2, wt_1 = 1, wt_2 = 1) {
  return(tryCatch(combine.test(p = c(pval_1, pval_2), weight = c(wt_1, wt_2), method = "fisher"), error = function(e) NA))
}

# impose effect size threshold, consistent with Griffiths et al.; set p-values to 1
results[, scploid_effect_p := scploid_p]
results[(scploid_score > 0.8 & scploid_score < 1.2), scploid_effect_p := 1]
results[, fisher_p := mapply(fisher_wrapper, results$ase_p, results$scploid_effect_p)]

# control the FDR with TreeBH
vary_fdr <- function(fdr, results_dt) {
  message(fdr)
  sc_groups <- as.matrix(results[, c("embryo", "cell", "chrom")])
  calls <- suppressWarnings(get_TreeBH_selections(results$fisher_p,
                                                  sc_groups,
                                                  q = c(fdr, fdr, fdr)))
  
  results_dt[, sig_embryo := calls[, 1]]
  results_dt[, sig_cell := calls[, 2]]
  results_dt[, sig_chrom := calls[, 3]]
  
  results_dt[cell %in% unique(results_dt[sig_cell == 1]$cell), sig_cell := 1]
  results_dt[embryo %in% unique(results_dt[sig_embryo == 1]$embryo), sig_embryo := 1]
  
  fdr_table <- data.table(
    table(results_dt[!duplicated(chrom)]$sig_chrom),
    table(results_dt[!duplicated(cell)]$sig_cell),
    table(results_dt[!duplicated(embryo)]$sig_embryo),
    FDR = fdr) %>%
    setnames(c("ploidy", "n_chromosomes", "ploidy_2", "n_cells", "ploidy_3", "n_embryos", "FDR"))
  return(fdr_table)
}

sig_by_fdr <- do.call(rbind, 
                      lapply(sort(unique(c(rev(1*10^-(1:5)), 5e-3, seq(1, 10, 2) * 10^-2, seq(0.2, 0.5, 0.1)))), 
                             function(x) vary_fdr(x, results)))

fdr_plot <- ggplot(data = sig_by_fdr[ploidy == 1 & FDR <= 0.5]) +
  geom_line(aes(x = FDR, y = n_chromosomes / sum(sig_by_fdr[1:2,]$n_chromosomes), color = "Chromosomes"), lwd = 1) +
  geom_line(aes(x = FDR, y = n_cells / sum(sig_by_fdr[1:2,]$n_cells), color = "Cells"), lwd = 1) +
  geom_line(aes(x = FDR, y = n_embryos / sum(sig_by_fdr[1:2,]$n_embryos), color = "Embryos"), lwd = 1) +
  ylab("Prop. aneuploid") +
  xlab("FDR") +
  theme_classic() +
  theme(legend.position = "none") +
  ylim(0, 1) +
  xlim(0, 0.75) + 
  annotate(geom = "text", x = 0.52, y = 0.91, label = "Embryos", hjust = "left") +
  annotate(geom = "text", x = 0.52, y = 0.51, label = "Cells", hjust = "left") +
  annotate(geom = "text", x = 0.52, y = 0.11, label = "Chromosomes", hjust = "left") +
  scale_color_brewer(palette = "Dark2")

# set FDR = 1%
fdr <- 0.01
sc_groups <- as.matrix(results[, c("embryo", "cell", "chrom")])
calls <- suppressWarnings(get_TreeBH_selections(results$fisher_p,
                                                sc_groups,
                                                q = c(fdr, fdr, fdr)))

results[, sig_embryo := calls[, 1]]
results[, sig_cell := calls[, 2]]
results[, sig_chrom := calls[, 3]]

results[cell %in% unique(results[sig_cell == 1]$cell), sig_cell := 1]
results[embryo %in% unique(results[sig_embryo == 1]$embryo), sig_embryo := 1]

length(unique(results[sig_chrom == 1]$chrom))
length(unique(results[sig_chrom == 1]$cell))
length(unique(results[sig_chrom == 1]$embryo))

cell_fraction <- group_by(results[!duplicated(cell)], EStage, embryo) %>%
  summarize(., prop_aneuploid = mean(sig_cell), n = n()) %>%
  as.data.table()

cell_fraction_plot <- ggplot(data = cell_fraction) +
  geom_histogram(aes(x = prop_aneuploid), bins = 30) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Prop. aneuploid cells") +
  ylab("Number of embryos")

chr_fraction <- group_by(results, EStage, embryo, chr) %>%
  summarize(., prop_aneuploid = mean(sig_chrom), n = n()) %>%
  as.data.table()

length(unique(chr_fraction[prop_aneuploid >= 0.75]$embryo))
length(unique(chr_fraction[prop_aneuploid > 0 & prop_aneuploid < 0.75]$embryo))
sum(unique(chr_fraction[prop_aneuploid >= 0.75]$embryo) %in% 
      unique(chr_fraction[prop_aneuploid > 0 & prop_aneuploid < 0.75]$embryo))

plot_grid(fdr_plot, cell_fraction_plot, labels = c('A', 'B'))

# calculate aneuploidies per chromosome
aneuploid_by_chr <- group_by(results, chr) %>%
  summarize(aneuploid_cells = sum(sig_chrom == 1), euploid_cells = sum(sig_chrom == 0)) %>%
  as.data.table()
aneuploid_by_chr$chr <- factor(aneuploid_by_chr$chr, paste0("chr", c(1:22, "X")))
aneuploid_by_chr[, chr_numeric := gsub("chr", "", chr)]
aneuploid_by_chr$chr_numeric <- factor(aneuploid_by_chr$chr_numeric, c(1:22, "X"))

gencode <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz")
genes_per_chr <- data.table(table(gencode[V3 == "gene" & grepl("protein_coding", V9)][, 1])) %>%
  setnames(., c("chr", "n_genes"))
chrom_lengths <- fread("https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes") %>%
  setnames(., c("chr", "len"))

aneuploid_by_chr <- merge(chrom_lengths, merge(aneuploid_by_chr, genes_per_chr, "chr"), "chr")

cor.test(aneuploid_by_chr$aneuploid_cells, aneuploid_by_chr$len)
cor.test(aneuploid_by_chr$aneuploid_cells, aneuploid_by_chr$n_genes)

by_chrom_plot <- ggplot(data = aneuploid_by_chr, aes(x = n_genes , y = aneuploid_cells, label = chr)) +
  theme_classic() +
  xlab("Number of protein-coding genes") +
  ylab("Number of aneuploid cells") +
  geom_point() +
  geom_label_repel(size = 4) +
  ylim(0, 125) +
  xlim(0, 2250)

m1 <- glmer(data = results, formula = (sig_chrom == 1) ~ (1 | embryo / cell) + (1 | lineage) + chr, family = binomial, nAGQ = 0)
m0 <- glmer(data = results, formula = (sig_chrom == 1) ~ (1 | embryo / cell) + (1 | lineage), family = binomial, nAGQ = 0)
anova(m1, m0, test = "Chisq") # embryo-specific models; to average over levels of random effect, need to extract average marginal effects as below
mx <- margins(m1, type = "response", variables = "chr")
b <- summary(mx)
cov_mat <- attr(mx, "vcov")
k <- diag(nrow = ncol(cov_mat))
kvar <- t(k) %*% cov_mat %*% k
kb <- k %*% b$AME
my_chi <- t(kb) %*% solve(kvar) %*% kb
pchisq(my_chi[1, 1], df = ncol(cov_mat), lower.tail = F)

# cluster with k-means, assign clusters to monosomy and trisomy
km <- results[sig_chrom == 1, c("scploid_z", "ase_z")] %>%
  kmeans(centers = 2)
km_clusters <- results[sig_chrom == 1, c("chrom", "scploid_z", "ase_z")]
km_clusters[, cluster := km$cluster]

results <- merge(results, km_clusters[, c("chrom", "cluster")], "chrom", all.x = TRUE)

results[, ploidy := 2]
if (mean(results[cluster == 1]$ase_z) < mean(results[cluster == 2]$ase_z)) {
  results[cluster == 1, ploidy := 1]
  results[cluster == 2, ploidy := 3] 
} else {
  results[cluster == 2, ploidy := 1]
  results[cluster == 1, ploidy := 3] 
}

ggplot() +
  geom_point(data = results[ploidy == 2], aes(x = scploid_z, y = ase_z, col = "disomy"), size = 0.3) +
  geom_point(data = results[ploidy == 1], aes(x = scploid_z, y = ase_z, col = "monosomy"), size = 0.3) +
  geom_point(data = results[ploidy == 3], aes(x = scploid_z, y = ase_z, col = "trisomy"), size = 0.3) +
  theme_classic() +
  scale_color_manual(name = "", values = c("gray", "blue", "red")) +
  xlab("scploid z-score") +
  ylab("Allelic imbalance z-score")

# plot aneuploidy heatmaps for individual embryos

results[, chr_num := gsub("chr", "", chr)]
results$chr_num <- factor(results$chr_num, levels = 1:22)
results[, log_fisher_p := -log10(fisher_p)]
results[, monosomy := scploid_z < 0 | ase_z < -3]

plot_mca <- function(embryo_id, combine = FALSE, legend = TRUE, cluster_method = "ward.D2", dist_method = "euclidean") {
  heatmap_data <- pivot_wider(results[(embryo == embryo_id), c("chr_num", "cell", "ploidy")], names_from = chr_num, values_from = ploidy)
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$cell
  
  distance.row <- dist(heatmap_matrix, method = dist_method)
  cluster.row <- hclust(distance.row, method = cluster_method)
  
  dendrogram <- ggplot(segment(dendro_data(cluster.row))) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0.2, 0)) + 
    theme_dendro() + 
    scale_x_reverse()
  
  sample_order <- rev(cluster.row$labels[cluster.row$order])
  
  dt_to_plot <- results[embryo == embryo_id]
  dt_to_plot$cell <- factor(dt_to_plot$cell, levels = sample_order)
  
  exp_heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = scploid_z)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradientn(name = "Z-score", colors = c("blue", rep("white", 3), "red"), limits = c(-5, 5), na.value = 1, oob = squish) +
    xlab("Chromosome") +
    ylab("Cell") +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  
  ase_heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = ase_z)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradientn(name = "Z-score", colors = c("blue", rep("white", 3), "red"), limits = c(-3, 3), na.value = 1, oob = squish) +
    xlab("Chromosome") +
    ylab("Cell") +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  
  heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = monosomy, alpha = log_fisher_p + 1)) +
    geom_tile() +
    theme_bw() +
    scale_fill_manual(name = "", values = c("red", "blue"), labels = c("Trisomy", "Monosomy")) +
    scale_alpha(name = "-log10(p)", range = c(0, 1), limits = c(2, 6), na.value = 1, oob = squish) +
    xlab("Chromosome") +
    ylab("") +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(5.5, 5.5, 5.5, -3), "pt"), panel.grid = element_blank())
  
  if (combine == TRUE) {
    plot_grid(dendrogram, heatmap, align = "h", axis = "b", rel_widths = c(0.3, 1), scale = c(1, 0.95))
  } else {
    plot_grid(exp_heatmap, ase_heatmap, align = "h", axis = "b", rel_widths = c(1, 1))
  }
}

plot_mca_sig <- function(embryo_id) {
  heatmap_data <- pivot_wider(results[(embryo == embryo_id), c("chr_num", "cell", "ploidy")], names_from = chr_num, values_from = ploidy)
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$cell
  
  distance.row <- dist(heatmap_matrix, method = "euclidean") # same parameters as honeyBADGER
  cluster.row <- hclust(distance.row, method = "ward.D") # same parameters as honeyBADGER
  
  dendrogram <- ggplot(segment(dendro_data(cluster.row))) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + 
    scale_y_reverse(expand=c(0.2, 0)) + 
    theme_dendro() + 
    scale_x_reverse()
  
  sample_order <- rev(cluster.row$labels[cluster.row$order])
  
  dt_to_plot <- results[embryo == embryo_id]
  dt_to_plot$cell <- factor(dt_to_plot$cell, levels = sample_order)
  dt_to_plot[, polarized_log_fisher_p := log_fisher_p]
  dt_to_plot[monosomy == TRUE, polarized_log_fisher_p := -1 * log_fisher_p]
  dt_to_plot[allelic_ratio < 0.1, polarized_log_fisher_p := -1 * log_fisher_p]
  
  dt_to_plot$ploidy_factor <- factor(dt_to_plot$ploidy, levels = c(1, 2, 3))
  
  heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = ploidy_factor)) +
    geom_tile() +
    theme_bw() +
    scale_fill_manual(name = "Ploidy", values = c("blue", "white", "red"), drop = FALSE) +
    xlab("Chromosome") +
    ylab("") +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(5.5, 5.5, 5.5, 0), "pt"))
  
  plot_grid(dendrogram, heatmap, align = "h", rel_widths = c(0.3, 1), rel_heights = c(1.1, 1))
}

for (embryo_id in unique(results$embryo)) {
  pdf(file = paste0("~/Downloads/mosaic_aneuploidy_plots/", embryo_id, ".pdf"), height = 4, width = 6)
  try(print(plot_mca(embryo_id, combine = TRUE, cluster_method = "ward.D2", dist_method = "euclidean")))
  dev.off()
}

plot_a <- plot_mca("E7.3", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_b <- plot_mca("E5.13", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_c <- plot_mca("E7.17", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_d <- plot_mca("E7.5", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_grid(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'))

for (embryo_id in unique(results$embryo)) {
  pdf(file = paste0("~/Downloads/mosaic_aneuploidy_plots/", embryo_id, ".pdf"), height = 4, width = 6)
  try(print(plot_mca_sig(embryo_id)))
  dev.off()
}

### statistical models of cell-type-specificity

results[, num_estage := scale(as.numeric(gsub("E", "", EStage)) - 3)]

m1 <- glmer(data = results[!duplicated(cell)], 
            formula = sig_cell ~ num_estage + (1 | embryo) + lineage, 
            family = binomial)

mx <- margins(m1, type = "response", variables = "lineage")
b <- summary(mx)
cov_mat <- attr(mx, "vcov")
k <- diag(nrow = ncol(cov_mat))
kvar <- t(k) %*% cov_mat %*% k
kb <- k %*% b$AME
my_chi <- t(kb) %*% solve(kvar) %*% kb
pchisq(my_chi[1, 1], df = ncol(cov_mat), lower.tail = F)

results[, is_trophectoderm := lineage == "Trophectoderm"]
te_enrich <- glmer(data = results[!duplicated(cell) & lineage != "Undefined"], 
                   formula = sig_cell ~ num_estage + (1 | embryo)  + is_trophectoderm, 
                   family = binomial)
summary(margins(te_enrich))

std <- function(x) sd(x) / sqrt(length(x))

aneuploid_by_lineage <- group_by(results[!duplicated(cell)], lineage) %>%
  summarize(., n_euploid = sum(sig_cell == 0), n_aneuploid = sum(sig_cell == 1), total = n())

aneuploid_by_lineage_ci <- mapply(function(x, y) tidy(prop.test(x, y)), 
                                  aneuploid_by_lineage$n_aneuploid, aneuploid_by_lineage$total) %>%
  t() %>%
  as.data.table()

aneuploid_by_lineage <- cbind(aneuploid_by_lineage, aneuploid_by_lineage_ci)
aneuploid_by_lineage[, estimate := as.numeric(estimate)]
aneuploid_by_lineage[, conf.low := as.numeric(conf.low)]
aneuploid_by_lineage[, conf.high := as.numeric(conf.high)]

aneuploid_by_lineage$lineage <- str_wrap(aneuploid_by_lineage$lineage, width = 10)

aneuploid_by_lineage$lineage <- factor(aneuploid_by_lineage$lineage, 
                                       levels = c("Undefined", "ICM", "Trophectoderm",
                                                  "Intermediate",
                                                  "Epiblast", "Primitive\nEndoderm"))

by_celltype_plot <- ggplot(data = aneuploid_by_lineage, 
                           aes(x = lineage, y = estimate, ymin = conf.low, ymax = conf.high, fill = lineage)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.25) +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Prop. aneuploid cells (95% CI)") +
  xlab("Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylim(0, 1)

results$lineage <- factor(results$lineage, ordered = FALSE)
results$lineage <- relevel(results$lineage, ref = "Trophectoderm")

enrich_model <- glmer(data = results[!duplicated(cell)], 
                      formula = sig_cell ~ (1 | embryo) + lineage, 
                      family = binomial)

enrich_coef <- summary(margins(enrich_model, type = "response", data = results[!duplicated(cell)])) %>%
  as.data.table() %>%
  setnames(., "factor", "term")

enrich_coef[, lineage := gsub("lineage", "", term)]
enrich_coef$lineage <- str_wrap(enrich_coef$lineage, width = 10)

enrich_coef$lineage <- factor(enrich_coef$lineage, 
                              levels = c("Undefined", "ICM",
                                         "Intermediate",
                                         "Epiblast", "Primitive\nEndoderm"))

enrichment_plot <- ggplot(data = enrich_coef, aes(x = lineage, y = AME, 
                                                  ymin = lower, ymax = upper,
                                                  color = lineage)) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  ylab("AME vs. trophectoderm (95% CI)") +
  xlab("Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  geom_hline(yintercept = 0, lty = "dashed", color = "gray") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a", "#66a61e", "#e6ab02"))


### dimension reduction visualization

emtab3929 <- readRDS(here("/RawData/emtab3929/EMTAB3929.rds"))
emtab3929_gene <- experiments(emtab3929)[["gene"]]
# assays(emtab3929_gene)[["count"]][1:10, 1:10]

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
