list_of_packages <- c("BiocStyle", "biomaRt", "broom", "cowplot", "data.table", "devtools", 
                      "dplyr", "ggdendro", "ggrepel", "gmodels", "gplots", "gridExtra", 
                      "here", "lme4", "margins", "mixtools", "monocle3", "mppa", "MultiAssayExperiment", 
                      "MultiAssayExperiment", "plyr", "readxl", "Rtsne", "scales", 
                      "scater", "scploid", "scran", "stringr", "survcomp", "tidyr", 
                      "tools", "TreeBH", "umap", "zoo")

# https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
install.packages.auto <- function(x) { 
  if (isTRUE(x %in% .packages(all.available = TRUE))) { 
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

# install scploid
if (!("scploid" %in% .packages(all.available = TRUE)))
  devtools::install_github("MarioniLab/Aneuploidy2017", subdir = "package")

# install TreeBH
if (!("TreeBH" %in% .packages(all.available = TRUE)))
  install.packages("https://odin.mdacc.tmc.edu/~cbpeterson/TreeBH_1.0.tar.gz", repos = NULL, type = "source")

# install monocle3
if (!("monocle3" %in% .packages(all.available = TRUE))) {
  monocle_prereqs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor')
  lapply(monocle_prereqs, function(x) {message(x); install.packages.auto(x)})
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
}

# install and load other packages
lapply(list_of_packages, function(x) {message(x); install.packages.auto(x)})

# over-ride masked functions
here <- here::here
summarize <- dplyr::summarize

# Change global default setting so every data frame created will not auto-convert to factors unless explicitly instructed
options(stringsAsFactors = FALSE) 

# Load data from mouse 8-cell stage G&T-seq and human Trisomy 21 G&T-seq
load(here("RawData/proc_data/emb8_data.RData"))
load(here("RawData/proc_data/tris_data.RData"))

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

# apply methods from Griffiths et al. (https://github.com/MarioniLab/Aneuploidy2017) to assess gene expression variance 
test_median <- 50
sample_size <- 1000
set.seed(12983)

## human trisomy 21 G&T-seq data
t21_trisomy <- getCPM(splitCellsByGroup(t21)[["T21"]])
t21_norm <- getCPM(splitCellsByGroup(t21)[["Diploid"]])

t21_trisomy_HEXP <- t21_trisomy[apply(t21_trisomy, 1, median) > test_median, ]
t21_norm_HEXP <- t21_norm[apply(t21_norm, 1, median) > test_median, ]

t21_trisomy_mean <- apply(t21_trisomy_HEXP, 1, mean)
t21_trisomy_sd <- apply(t21_trisomy_HEXP, 1, sd)
t21_norm_mean <- apply(t21_norm_HEXP, 1, mean)
t21_norm_sd <- apply(t21_norm_HEXP, 1, sd)

t21_trisomy_lm <- lm(log10(t21_trisomy_sd) ~ log10(t21_trisomy_mean))
t21_norm_lm <- lm(log10(t21_norm_sd) ~ log10(t21_norm_mean))

t21_trisomy_sample <- sample(length(t21_trisomy_mean), sample_size)
t21_norm_sample <- sample(length(t21_norm_mean), sample_size)

## mouse embryo G&T-seq data
emb8_reversine <- getCPM(splitCellsByGroup(emb8)[["Reversine"]])
emb8_control <- getCPM(splitCellsByGroup(emb8)[["Control"]])

emb8_reversine <- emb8_reversine[rowMedians(emb8_reversine) > test_median, ]
emb8_control <- emb8_control[rowMedians(emb8_control) > test_median,]

emb8_reversine_mean <- rowMeans(emb8_reversine)
emb8_reversine_sd <- apply(emb8_reversine, 1, sd)
emb8_control_mean <- rowMeans(emb8_control)
emb8_control_sd <- apply(emb8_control, 1, sd)

emb8_reversine_lm <- lm(log10(emb8_reversine_sd) ~ log10(emb8_reversine_mean))
emb8_control_lm <- lm(log10(emb8_control_sd) ~ log10(emb8_control_mean))

emb8_reversine_sample <- sample(length(emb8_reversine_mean), sample_size)
emb8_control_sample <- sample(length(emb8_control_mean), sample_size)

## human EMTAB3929 scRNA-seq data
group <- list()
group_mean <- list()
group_sd <- list()
group_lm <- list()
group_sample <- list()

sample_mean <- list()
sample_sd <- list()

for (k in 1:length(spt)) {
  
  group[[k]] <- spt[[k]]@cpm
  group[[k]] <- group[[k]][rowMedians(group[[k]]) > test_median, ]
  
  group_mean[[k]] <- rowMeans(group[[k]])
  group_sd[[k]] <- apply(group[[k]], 1, sd)
  group_lm[[k]] <- lm(log10(group_sd[[k]]) ~ log10(group_mean[[k]]))
  group_sample[[k]] <- sample(length(group_mean[[k]]), sample_size)
  
  sample_mean[[k]] <- group_mean[[k]][group_sample[[k]]]
  sample_sd[[k]] <- group_sd[[k]][group_sample[[k]]]
  
}

means <- do.call("cbind", sample_mean) %>% as.data.frame()
colnames(means) <- names(spt)
means$t21_trisomy <- t21_trisomy_mean[t21_trisomy_sample]
means$t21_normal <- t21_norm_mean[t21_norm_sample]
means$emb8_reversine <- emb8_reversine_mean[emb8_reversine_sample]
means$emb8_control <- emb8_control_mean[emb8_control_sample]

variances <- do.call("cbind", sample_sd) %>% as.data.frame()
colnames(variances) <- names(spt)
variances$t21_trisomy <- t21_trisomy_sd[t21_trisomy_sample]
variances$t21_normal <- t21_norm_sd[t21_norm_sample]
variances$emb8_reversine <- emb8_reversine_sd[emb8_reversine_sample]
variances$emb8_control <- emb8_control_sd[emb8_control_sample]

means_long <- reshape2::melt(means)
variances_long <- reshape2::melt(variances)

mean_variance <- data.table(
  sample = means_long$variable,
  mean = means_long$value,
  variance = variances_long$value
)
mean_variance[, c("stage", "lineage") := tstrsplit(sample, "_", fixed = TRUE)]

intercept_slope <- t(sapply(list(group_lm[[1]], 
                                 group_lm[[2]],
                                 group_lm[[3]],
                                 group_lm[[4]],
                                 group_lm[[5]],
                                 group_lm[[6]],
                                 group_lm[[7]],
                                 group_lm[[8]],
                                 group_lm[[9]],
                                 group_lm[[10]],
                                 group_lm[[11]],
                                 t21_trisomy_lm, 
                                 t21_norm_lm, 
                                 emb8_reversine_lm, 
                                 emb8_control_lm), coef)) %>%
  as.data.table()

intercept_slope[, sample := colnames(variances)]
intercept_slope[, c("stage", "lineage") := tstrsplit(sample, "_", fixed = TRUE)]

plot_mean_sd <- function(day, mean_variance_dt, intercept_slope_dt) {
  mean_variance_dt$facet_title <- day
  intercept_slope_dt$facet_title <- day
  mean_sd_plot <- ggplot(mean_variance_dt[stage == day | (sample %in% c("t21_trisomy", "t21_normal", "emb8_reversine", "emb8_control"))], 
                         aes(x = mean, y = variance, col = sample)) + 
    geom_point(size = 0.75, alpha = 0.5) +  
    scale_x_log10() + scale_y_log10() +
    theme_classic() +
    scale_color_manual(values = c("t21_trisomy" = "grey34", "t21_normal" = "grey34",
                                  "emb8_reversine" = "grey77", "emb8_control" = "grey77",
                                  "E4_Undefined" = "#1b9e77",
                                  "E5_ICM" = "#d95f02",
                                  "E5_Trophectoderm" = "#7570b3",
                                  "E5_Undefined" = "#1b9e77",
                                  "E6_Trophectoderm" = "#7570b3",
                                  "E6_Primitive Endoderm" = "#e6ab02",
                                  "E6_Epiblast" = "#66a61e",
                                  "E7_Trophectoderm" = "#7570b3",
                                  "E7_Intermediate" = "#e7298a",
                                  "E7_Epiblast" = "#66a61e",
                                  "E7_Primitive Endoderm" = "#e6ab02"), name = "") +
    geom_abline(data = intercept_slope_dt[stage == day | (sample %in% c("t21_trisomy", "t21_normal", "emb8_reversine", "emb8_control"))], 
                aes(slope = `log10(group_mean[[k]])`, intercept = `(Intercept)`, col = sample), alpha = 0.8, lwd = 1) +
    labs(x = expression("log"[10]*"(mean)"), y = expression("log"[10]* "(standard deviation)")) +
    facet_wrap(~ facet_title)
  
  return(mean_sd_plot)
}

mean_sd_plots <- lapply(paste0("E", 4:7), function(x) plot_mean_sd(x, mean_variance, intercept_slope))

plot_grid(plotlist = mean_sd_plots, ncol = 2, align = "hv", axis = "tblr")

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
  theme_classic() +
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
  pdf(file = paste0(here("/results/mosaic_aneuploidy_plots/"), embryo_id, ".pdf"), height = 4, width = 6)
  try(print(plot_mca(embryo_id, combine = TRUE, cluster_method = "ward.D2", dist_method = "euclidean")))
  dev.off()
}

plot_a <- plot_mca("E7.3", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_b <- plot_mca("E5.13", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_c <- plot_mca("E7.17", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_d <- plot_mca("E7.5", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_grid(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'))

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

emtab3929 <- readRDS(here("/RawData/EMTAB3929.rds"))
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

fwrite(results, here("results/aneuploidy_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

###

devtools::session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date       lib source                                     
# abind                  1.4-5      2016-07-21 [1] CRAN (R 3.6.0)                             
# acepack                1.4.1      2016-10-29 [1] CRAN (R 3.6.0)                             
# AnnotationDbi          1.48.0     2019-10-29 [1] Bioconductor                               
# askpass                1.1        2019-01-13 [1] CRAN (R 3.6.0)                             
# assertthat             0.2.1      2019-03-21 [1] CRAN (R 3.6.0)                             
# backports              1.1.5      2019-10-02 [1] CRAN (R 3.6.0)                             
# base64enc              0.1-3      2015-07-28 [1] CRAN (R 3.6.0)                             
# batchelor              1.2.2      2019-11-07 [1] Bioconductor                               
# beeswarm               0.2.3      2016-04-25 [1] CRAN (R 3.6.0)                             
# BiasedUrn              1.07       2015-12-28 [1] CRAN (R 3.6.0)                             
# bibtex                 0.4.2      2017-06-30 [1] CRAN (R 3.6.0)                             
# Biobase              * 2.46.0     2019-10-29 [1] Bioconductor                               
# BiocFileCache          1.10.2     2019-11-08 [1] Bioconductor                               
# BiocGenerics         * 0.32.0     2019-10-29 [1] Bioconductor                               
# BiocManager            1.30.10    2019-11-16 [1] CRAN (R 3.6.0)                             
# BiocNeighbors          1.4.1      2019-11-01 [1] Bioconductor                               
# BiocParallel         * 1.20.0     2019-10-29 [1] Bioconductor                               
# BiocSingular           1.2.0      2019-10-29 [1] Bioconductor                               
# BiocStyle            * 2.14.0     2019-10-29 [1] Bioconductor                               
# biomaRt              * 2.42.0     2019-10-29 [1] Bioconductor                               
# BioNet                 1.46.0     2019-10-29 [1] Bioconductor                               
# Biostrings             2.54.0     2019-10-29 [1] Bioconductor                               
# bit                    1.1-14     2018-05-29 [1] CRAN (R 3.6.0)                             
# bit64                  0.9-7      2017-05-08 [1] CRAN (R 3.6.0)                             
# bitops                 1.0-6      2013-08-17 [1] CRAN (R 3.6.0)                             
# blob                   1.2.0      2019-07-09 [1] CRAN (R 3.6.0)                             
# boot                   1.3-23     2019-07-05 [1] CRAN (R 3.6.0)                             
# bootstrap              2019.6     2019-06-17 [1] CRAN (R 3.6.0)                             
# broom                * 0.5.2      2019-04-07 [1] CRAN (R 3.6.0)                             
# callr                  3.3.2      2019-09-22 [1] CRAN (R 3.6.0)                             
# caTools                1.17.1.2   2019-03-06 [1] CRAN (R 3.6.0)                             
# cellranger             1.1.0      2016-07-27 [1] CRAN (R 3.6.0)                             
# checkmate              1.9.4      2019-07-04 [1] CRAN (R 3.6.0)                             
# cli                    1.1.0      2019-03-19 [1] CRAN (R 3.6.0)                             
# cluster                2.1.0      2019-06-19 [1] CRAN (R 3.6.1)                             
# codetools              0.2-16     2018-12-24 [1] CRAN (R 3.6.1)                             
# colorspace             1.4-1      2019-03-18 [1] CRAN (R 3.6.0)                             
# cowplot              * 1.0.0      2019-07-11 [1] CRAN (R 3.6.0)                             
# crayon                 1.3.4      2017-09-16 [1] CRAN (R 3.6.0)                             
# curl                   4.2        2019-09-24 [1] CRAN (R 3.6.0)                             
# data.table           * 1.12.6     2019-10-18 [1] CRAN (R 3.6.0)                             
# DBI                    1.0.0      2018-05-02 [1] CRAN (R 3.6.0)                             
# dbplyr                 1.4.2      2019-06-17 [1] CRAN (R 3.6.0)                             
# DelayedArray         * 0.12.0     2019-10-29 [1] Bioconductor                               
# DelayedMatrixStats     1.8.0      2019-10-29 [1] Bioconductor                               
# desc                   1.2.0      2018-05-01 [1] CRAN (R 3.6.0)                             
# devtools             * 2.2.1      2019-09-24 [1] CRAN (R 3.6.0)                             
# digest                 0.6.23     2019-11-23 [1] CRAN (R 3.6.0)                             
# dplyr                * 0.8.3      2019-07-04 [1] CRAN (R 3.6.0)                             
# dqrng                  0.2.1      2019-05-17 [1] CRAN (R 3.6.0)                             
# edgeR                  3.28.0     2019-10-29 [1] Bioconductor                               
# ellipsis               0.3.0      2019-09-20 [1] CRAN (R 3.6.0)                             
# evaluate               0.14       2019-05-28 [1] CRAN (R 3.6.0)                             
# farver                 2.0.1      2019-11-13 [1] CRAN (R 3.6.0)                             
# foreign                0.8-72     2019-08-02 [1] CRAN (R 3.6.0)                             
# Formula                1.2-3      2018-05-03 [1] CRAN (R 3.6.0)                             
# fs                     1.3.1      2019-05-06 [1] CRAN (R 3.6.0)                             
# gbRd                   0.4-11     2012-10-01 [1] CRAN (R 3.6.0)                             
# gdata                  2.18.0     2017-06-06 [1] CRAN (R 3.6.0)                             
# geneLenDataBase        1.22.0     2019-11-05 [1] Bioconductor                               
# generics               0.0.2      2018-11-29 [1] CRAN (R 3.6.0)                             
# GenomeInfoDb         * 1.22.0     2019-10-29 [1] Bioconductor                               
# GenomeInfoDbData       1.2.2      2019-11-26 [1] Bioconductor                               
# GenomicAlignments      1.22.1     2019-11-12 [1] Bioconductor                               
# GenomicFeatures        1.38.0     2019-10-29 [1] Bioconductor                               
# GenomicRanges        * 1.38.0     2019-10-29 [1] Bioconductor                               
# ggbeeswarm             0.6.0      2017-08-07 [1] CRAN (R 3.6.0)                             
# ggdendro             * 0.1-20     2016-04-27 [1] CRAN (R 3.6.0)                             
# ggplot2              * 3.2.1      2019-08-10 [1] CRAN (R 3.6.0)                             
# ggrepel              * 0.8.1      2019-05-07 [1] CRAN (R 3.6.0)                             
# glue                   1.3.1      2019-03-12 [1] CRAN (R 3.6.0)                             
# gmodels              * 2.18.1     2018-06-25 [1] CRAN (R 3.6.0)                             
# GO.db                  3.10.0     2019-12-20 [1] Bioconductor                               
# goseq                  1.38.0     2019-10-29 [1] Bioconductor                               
# gplots               * 3.0.1.1    2019-01-27 [1] CRAN (R 3.6.0)                             
# gridExtra            * 2.3        2017-09-09 [1] CRAN (R 3.6.0)                             
# gtable                 0.3.0      2019-03-25 [1] CRAN (R 3.6.0)                             
# gtools                 3.8.1      2018-06-26 [1] CRAN (R 3.6.0)                             
# here                 * 0.1        2017-05-28 [1] CRAN (R 3.6.0)                             
# Hmisc                  4.3-0      2019-11-07 [1] CRAN (R 3.6.0)                             
# hms                    0.5.2      2019-10-30 [1] CRAN (R 3.6.0)                             
# htmlTable              1.13.2     2019-09-22 [1] CRAN (R 3.6.0)                             
# htmltools              0.4.0      2019-10-04 [1] CRAN (R 3.6.0)                             
# htmlwidgets            1.5.1      2019-10-08 [1] CRAN (R 3.6.0)                             
# httr                   1.4.1      2019-08-05 [1] CRAN (R 3.6.0)                             
# igraph                 1.2.4.1    2019-04-22 [1] CRAN (R 3.6.0)                             
# IRanges              * 2.20.1     2019-11-20 [1] Bioconductor                               
# irlba                  2.3.3      2019-02-05 [1] CRAN (R 3.6.0)                             
# jsonlite               1.6        2018-12-07 [1] CRAN (R 3.6.0)                             
# KEGG.db                3.2.3      2019-12-20 [1] Bioconductor                               
# KernSmooth             2.23-16    2019-10-15 [1] CRAN (R 3.6.0)                             
# knitr                  1.26       2019-11-12 [1] CRAN (R 3.6.0)                             
# labeling               0.3        2014-08-23 [1] CRAN (R 3.6.0)                             
# lattice                0.20-38    2018-11-04 [1] CRAN (R 3.6.1)                             
# latticeExtra           0.6-28     2016-02-09 [1] CRAN (R 3.6.0)                             
# lava                   1.6.6      2019-08-01 [1] CRAN (R 3.6.0)                             
# lazyeval               0.2.2      2019-03-15 [1] CRAN (R 3.6.0)                             
# lifecycle              0.1.0      2019-08-01 [1] CRAN (R 3.6.0)                             
# limma                  3.42.0     2019-10-29 [1] Bioconductor                               
# lme4                 * 1.1-21     2019-03-05 [1] CRAN (R 3.6.0)                             
# locfit                 1.5-9.1    2013-04-20 [1] CRAN (R 3.6.0)                             
# magrittr               1.5        2014-11-22 [1] CRAN (R 3.6.0)                             
# margins              * 0.3.23     2018-05-22 [1] CRAN (R 3.6.0)                             
# MASS                   7.3-51.4   2019-03-31 [1] CRAN (R 3.6.1)                             
# MAST                 * 1.12.0     2019-10-29 [1] Bioconductor                               
# Matrix               * 1.2-17     2019-03-22 [1] CRAN (R 3.6.1)                             
# matrixStats          * 0.55.0     2019-09-07 [1] CRAN (R 3.6.0)                             
# memoise                1.1.0      2017-04-21 [1] CRAN (R 3.6.0)                             
# metap                * 1.1        2019-02-06 [1] CRAN (R 3.6.0)                             
# mgcv                   1.8-31     2019-11-09 [1] CRAN (R 3.6.0)                             
# minqa                  1.2.4      2014-10-09 [1] CRAN (R 3.6.0)                             
# mixtools             * 1.1.0      2017-03-10 [1] CRAN (R 3.6.0)                             
# monocle3             * 0.2.0      2019-11-27 [1] Github (cole-trapnell-lab/monocle3@9becd94)
# mppa                 * 1.0        2014-08-23 [1] CRAN (R 3.6.0)                             
# MultiAssayExperiment * 1.12.0     2019-10-29 [1] Bioconductor                               
# munsell                0.5.0      2018-06-12 [1] CRAN (R 3.6.0)                             
# nlme                   3.1-142    2019-11-07 [1] CRAN (R 3.6.0)                             
# nloptr                 1.2.1      2018-10-03 [1] CRAN (R 3.6.0)                             
# nnet                   7.3-12     2016-02-02 [1] CRAN (R 3.6.1)                             
# openssl                1.4.1      2019-07-18 [1] CRAN (R 3.6.0)                             
# org.Hs.eg.db           3.10.0     2019-12-20 [1] Bioconductor                               
# pillar                 1.4.2      2019-06-29 [1] CRAN (R 3.6.0)                             
# pkgbuild               1.0.6      2019-10-09 [1] CRAN (R 3.6.0)                             
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 3.6.0)                             
# pkgload                1.0.2      2018-10-29 [1] CRAN (R 3.6.0)                             
# plyr                 * 1.8.4      2016-06-08 [1] CRAN (R 3.6.0)                             
# prediction             0.3.14     2019-06-17 [1] CRAN (R 3.6.0)                             
# prettyunits            1.0.2      2015-07-13 [1] CRAN (R 3.6.0)                             
# processx               3.4.1      2019-07-18 [1] CRAN (R 3.6.0)                             
# prodlim              * 2019.11.13 2019-11-17 [1] CRAN (R 3.6.0)                             
# progress               1.2.2      2019-05-16 [1] CRAN (R 3.6.0)                             
# ps                     1.3.0      2018-12-21 [1] CRAN (R 3.6.0)                             
# purrr                  0.3.3      2019-10-18 [1] CRAN (R 3.6.0)                             
# qvalue               * 2.18.0     2019-10-29 [1] Bioconductor                               
# R.methodsS3            1.7.1      2016-02-16 [1] CRAN (R 3.6.0)                             
# R.oo                   1.23.0     2019-11-03 [1] CRAN (R 3.6.0)                             
# R.utils                2.9.0      2019-06-13 [1] CRAN (R 3.6.0)                             
# R6                     2.4.1      2019-11-12 [1] CRAN (R 3.6.0)                             
# rappdirs               0.3.1      2016-03-28 [1] CRAN (R 3.6.0)                             
# RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 3.6.0)                             
# Rcpp                   1.0.3      2019-11-08 [1] CRAN (R 3.6.0)                             
# RcppAnnoy              0.0.14     2019-11-12 [1] CRAN (R 3.6.0)                             
# RcppParallel           4.4.4      2019-09-27 [1] CRAN (R 3.6.0)                             
# RCurl                  1.95-4.12  2019-03-04 [1] CRAN (R 3.6.0)                             
# Rdpack                 0.11-0     2019-04-14 [1] CRAN (R 3.6.0)                             
# reactome.db            1.70.0     2019-12-20 [1] Bioconductor                               
# readxl               * 1.3.1      2019-03-13 [1] CRAN (R 3.6.0)                             
# remotes                2.1.0      2019-06-24 [1] CRAN (R 3.6.0)                             
# reshape2               1.4.3      2017-12-11 [1] CRAN (R 3.6.0)                             
# reticulate             1.13       2019-07-24 [1] CRAN (R 3.6.0)                             
# rlang                  0.4.2      2019-11-23 [1] CRAN (R 3.6.0)                             
# rmarkdown              1.17       2019-11-13 [1] CRAN (R 3.6.0)                             
# rmeta                  3.0        2018-03-20 [1] CRAN (R 3.6.0)                             
# rpart                  4.1-15     2019-04-12 [1] CRAN (R 3.6.1)                             
# rprojroot              1.3-2      2018-01-03 [1] CRAN (R 3.6.0)                             
# Rsamtools              2.2.1      2019-11-06 [1] Bioconductor                               
# RSpectra               0.15-0     2019-06-11 [1] CRAN (R 3.6.0)                             
# RSQLite                2.1.2      2019-07-24 [1] CRAN (R 3.6.0)                             
# rstudioapi             0.10       2019-03-19 [1] CRAN (R 3.6.0)                             
# rsvd                   1.0.2      2019-07-29 [1] CRAN (R 3.6.0)                             
# rtracklayer            1.46.0     2019-10-29 [1] Bioconductor                               
# Rtsne                * 0.15       2018-11-10 [1] CRAN (R 3.6.0)                             
# S4Vectors            * 0.24.0     2019-10-29 [1] Bioconductor                               
# scales               * 1.1.0      2019-11-18 [1] CRAN (R 3.6.0)                             
# scater               * 1.14.4     2019-11-18 [1] Bioconductor                               
# scploid              * 0.9        2019-11-26 [1] Github (MarioniLab/Aneuploidy2017@286c064) 
# scran                * 1.14.5     2019-11-19 [1] Bioconductor                               
# segmented              1.1-0      2019-12-10 [1] CRAN (R 3.6.0)                             
# sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 3.6.0)                             
# SingleCellExperiment * 1.8.0      2019-10-29 [1] Bioconductor                               
# SMITE                * 1.14.0     2019-10-29 [1] Bioconductor                               
# statmod                1.4.32     2019-05-29 [1] CRAN (R 3.6.0)                             
# stringi                1.4.3      2019-03-12 [1] CRAN (R 3.6.0)                             
# stringr              * 1.4.0      2019-02-10 [1] CRAN (R 3.6.0)                             
# SummarizedExperiment * 1.16.0     2019-10-29 [1] Bioconductor                               
# SuppDists              1.1-9.4    2016-09-23 [1] CRAN (R 3.6.0)                             
# survcomp             * 1.36.0     2019-10-29 [1] Bioconductor                               
# survival             * 3.1-7      2019-11-09 [1] CRAN (R 3.6.0)                             
# survivalROC            1.0.3      2013-01-13 [1] CRAN (R 3.6.0)                             
# testthat               2.3.0      2019-11-05 [1] CRAN (R 3.6.0)                             
# tibble                 2.1.3      2019-06-06 [1] CRAN (R 3.6.0)                             
# tidyr                * 1.0.0      2019-09-11 [1] CRAN (R 3.6.0)                             
# tidyselect             0.2.5      2018-10-11 [1] CRAN (R 3.6.0)                             
# TreeBH               * 1.0        2019-11-26 [1] local                                      
# umap                 * 0.2.3.1    2019-08-21 [1] CRAN (R 3.6.0)                             
# usethis              * 1.5.1      2019-07-04 [1] CRAN (R 3.6.0)                             
# uwot                   0.1.4      2019-09-23 [1] CRAN (R 3.6.0)                             
# vctrs                  0.2.0      2019-07-05 [1] CRAN (R 3.6.0)                             
# vipor                  0.4.5      2017-03-22 [1] CRAN (R 3.6.0)                             
# viridis                0.5.1      2018-03-29 [1] CRAN (R 3.6.0)                             
# viridisLite            0.3.0      2018-02-01 [1] CRAN (R 3.6.0)                             
# withr                  2.1.2      2018-03-15 [1] CRAN (R 3.6.0)                             
# xfun                   0.11       2019-11-12 [1] CRAN (R 3.6.0)                             
# XML                    3.98-1.20  2019-06-06 [1] CRAN (R 3.6.0)                             
# XVector                0.26.0     2019-10-29 [1] Bioconductor                               
# yaml                   2.2.0      2018-07-25 [1] CRAN (R 3.6.0)                             
# zeallot                0.1.0      2018-01-28 [1] CRAN (R 3.6.0)                             
# zlibbioc               1.32.0     2019-10-29 [1] Bioconductor                               
# zoo                  * 1.8-6      2019-05-28 [1] CRAN (R 3.6.0)                             
# 
# [1] /Library/Frameworks/R.framework/Versions/3.6/Resources/library


