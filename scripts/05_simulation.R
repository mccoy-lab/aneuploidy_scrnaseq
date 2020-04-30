list_of_packages <- c("BiocStyle", "broom", "cowplot", "data.table", "devtools", 
                      "dplyr", "ggrepel", "gmodels", "gplots", "gridExtra", 
                      "here", "mixtools", "mppa", "MultiAssayExperiment", 
                      "MultiAssayExperiment", "plyr", "readxl", "scales", 
                      "scater", "scploid", "scran", "stringr", "survcomp",
                      "tidyr", "tools", "zoo")

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

# install and load other packages
lapply(list_of_packages, function(x) {message(x); install.packages.auto(x)})

# over-ride masked functions
here <- here::here
summarize <- dplyr::summarize

# load gene Petropoulos et al. (2016) data

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

# load inferred aneuploidies
results <- fread(here("results/aneuploidy_results.txt")) # load results data

# examine allic ratios in real data
ggplot(data = results, aes(x = allelic_ratio, fill = factor(ploidy))) +
  geom_histogram(bins = 100) +
  facet_grid(factor(ploidy) ~ ., scales = "free_y")

# function to get parameters of a beta distribution given
# mean and variance
get_alpha_beta <- function(mu, var){
  b <- ( mu*(1-mu)*(1-mu)/var ) - (1-mu)
  a <- b*mu/(1-mu)
  c(a, b)
}

# function to simulate new allelic ratios of given ploidy
simulate_allelic_ratio <- function(results_dt, ploidy_sim = 2, n_sim = 1000) {
  
  mu_input <- mean(results_dt[ploidy == ploidy_sim]$allelic_ratio)
  var_input <- var(results_dt[ploidy == ploidy_sim]$allelic_ratio)

  beta_params <- get_alpha_beta(mu_input, var_input)

  sim_results <- data.table(ploidy = ploidy_sim,
                            sim_allelic_ratio = rbeta(n = n_sim, beta_params[1], beta_params[2]))
  return(sim_results)
  
}

# function to simulate dosage effects of aneuploidy, adapted from scploid
# code is identical, but can now work with gene_map as a data.table
subsample_simulate_aneuploidy <- function(counts, downsample_frac = 1, overdisperse_factor = 1, gene_map = NULL, known_aneuploidy = NULL) {
  
  if (is.null(gene_map) | is.null(known_aneuploidy)) {
    stop("Specify gene_map and known_aneuploidy")
  }
  
  # get the fit, allocate memory, get library sizes
  fit <- get_fit(counts)
  new_matrix = matrix(NA, ncol = ncol(counts), nrow = nrow(counts),
                      dimnames = list(rownames(counts), colnames(counts)))
  libs <- apply(counts, 2, sum)
  
  # make a ploidy state matrix
  ploidy_mat <- matrix(1, ncol = ncol(counts), nrow = nrow(counts),
                       dimnames = list(rownames(counts), colnames(counts)))
  for (i in 1:nrow(known_aneuploidy)) {
    chr <- as.numeric(as.character(known_aneuploidy$chr[i]))
    cell <- as.character(known_aneuploidy$cell[i])
    factor_aneu <- ifelse(known_aneuploidy$monosomy[i], 0.5, 1.5)
    
    target_genes <- gene_map[chromosome_name == chr]$ensembl_gene_id
    
    ploidy_mat[which(rownames(ploidy_mat) %in% target_genes), which(colnames(ploidy_mat) == cell)] <- factor_aneu
  }
  
  #estimate mu from the non-aneuploid cells
  normal_cells <- counts[, !colnames(counts)%in%known_aneuploidy$cell ]
  mu <- apply(normal_cells, 1, mean)
  #make a matrix of mus
  mu_matrix <- matrix(mu, byrow = FALSE, ncol = ncol(counts), nrow = nrow(counts),
                      dimnames = list(rownames(counts), colnames(counts)))
  
  #correct mus to maintain library size differences
  lib_ratio <- libs/mean(libs)
  mu_matrix <- sweep(mu_matrix, 2, lib_ratio, "*")
  
  #adjust the mean expression by the ploidy matrix
  mu_matrix <- mu_matrix * ploidy_mat

  #make the phi matrix, which will be the same dimensions and 
  #dimnames as the mu matrix
  phi_matrix <- mu_matrix
  
  #calculate phis from the fit of the adjusted mus
  for(col in 1:ncol(phi_matrix)) {
    #reassign columns of this matrix as we progress through it
    phi_matrix[, col] = 10 ^ predict(object = fit, newdata = log10(data.frame(x = phi_matrix[, col])))
  }
  
  genes_with_na <- apply(phi_matrix, 1, function(x) any(is.na(x)))
  #where phi is NA, we want instead to give it the value of the minimum considered mean expression in the fit
  phi_matrix[is.na(phi_matrix)] <- 10 ^ predict(object = fit, newdata = data.frame(x=log10(min(rowMeans(counts)[!genes_with_na]))))
  
  #sample from the matrix of mus, inferring phi from our fit for each mu separately
  for(cell in 1:ncol(counts)) {
    new_matrix[, cell] = rnbinom(n = nrow(mu_matrix), mu = mu_matrix[, cell] * downsample_frac, 
                                 size = 1 / (phi_matrix[, cell] * overdisperse_factor) )
  }
  return(new_matrix)
}

# overdispersion estimate, taken from scploid
od_est <- function(vec) {
  cv = sd(vec) / mean(vec)
  return(max((cv ^ 2 - (1 / mean(vec))), 0))
}

# function to fit the mean variance relationship, taken from scploid
get_fit <- function(counts_matrix){
  od <- apply(counts_matrix, 1, od_est)
  mean <- apply(counts_matrix, 1, mean)
  
  y <- log10(od)
  x <- log10(mean)
  
  #these genes really mess up the fit - remove them as they are not meaningful for our CPM>50 considerations
  remove <- which(is.na(y) | is.na(x) | is.infinite(y) | is.infinite(x) | y < -5)
  y <- y[-remove]
  x <- x[-remove]
  
  fit <- loess(y ~ x)
  return(fit)
}


### simulate "known" aneuploidies in the same quantities as our inferred aneuploidies

## simulate gene expression

# define cell groups over which to stratify simulation
results[, cell_group := paste(EStage, lineage, sep = "_")]

# wrapper function to peform simulation for a given cell group
sim_by_celltype <- function(results_dt, ploidytest_input, group_stratum) {
  
  message(paste("Simulating gene expression counts for group:", group_stratum))
  
  cell_list <- results_dt[!duplicated(cell) & cell_group == group_stratum]$cell
  
  known_aneuploidy_sim <- data.table(cell = results_dt[ploidy != 2 & cell %in% cell_list]$cell, 
                                     chr = results_dt[ploidy != 2 & cell %in% cell_list]$chr_num,
                                     monosomy = (results_dt[ploidy != 2 & cell %in% cell_list]$ploidy == 1))
  
  gene_map_sim <- data.table(ensembl_gene_id = rownames(ploidytest_input@counts),
                             chromosome_name = ploidytest_input@chrs)
  
  sim_counts <- subsample_simulate_aneuploidy(counts = ploidytest_input@counts[, cell_list], 
                                              overdisperse_factor = 1,
                                              gene_map = gene_map_sim,
                                              known_aneuploidy = known_aneuploidy_sim)
  
  return(sim_counts)
  
}

set.seed(12983)

# run gene expression simulation
sim_counts_stratified <- lapply(unique(results$cell_group), function(x) sim_by_celltype(results, ploidytest_dt, x))
sim_counts <- do.call(cbind, sim_counts_stratified)

# set NA counts to 0
sim_counts[is.nan(sim_counts)] <- 0

# counts to CPM
sim_counts_cpm <- sweep(sim_counts, 2, colSums(sim_counts), "/") * 1e6

# get the cell groups in the new order of the matrix columns
sim_cell_groups <- data.table(cell = colnames(sim_counts_cpm)) %>%
  .[, index := .I] %>%
  merge(., results[!duplicated(cell), c("cell", "cell_group")], by = "cell") %>%
  setorder(., index)

## infer aneuploidy in simulated data

gene_map <- data.table(ensembl_gene_id = rownames(ploidytest_dt@counts),
                       chromosome_name = ploidytest_dt@chrs)

# create scploid object
ploidytest_sim <- makeAneu(counts = sim_counts_cpm,
                           genes = rownames(sim_counts_cpm),
                           chrs = gene_map$chromosome_name,
                           cellNames = colnames(sim_counts_cpm),
                           cellGroups = sim_cell_groups$cell_group)

# data split into EStage_celltype subsets
spt_sim <- splitCellsByGroup(ploidytest_sim)

# run scploid
scploid_sim_results <- do.call(rbind, lapply(spt_sim, calcAneu)) %>%
  as.data.table() %>%
  setnames(., "z", "scploid_z") %>%
  setnames(., "score", "scploid_score") %>%
  setnames(., "p", "scploid_p")

known_aneuploidy_sim <- data.table(cell = results[ploidy != 2]$cell, 
                                   chr = results[ploidy != 2]$chr_num,
                                   monosomy = (results[ploidy != 2]$ploidy == 1))

scploid_sim_results[, chrom := paste(cell, chr, sep = "_")]
known_aneuploidy_sim[, chrom := paste(cell, chr, sep = "_")]

scploid_sim_results <- merge(scploid_sim_results, 
                             known_aneuploidy_sim[, -c("chr", "cell")], 
                             by = "chrom", all.x = TRUE, 
                             allow.cartesian = TRUE)

ggplot(data = scploid_sim_results, aes(x = scploid_p, fill = monosomy)) +
  geom_histogram() +
  facet_grid(monosomy ~ .)

# simulate allelic imbalance

simulated_ase <- rbind(simulate_allelic_ratio(results, 1, n = nrow(results[ploidy == 1])),
                       simulate_allelic_ratio(results, 2, n = nrow(results[ploidy == 2])),
                       simulate_allelic_ratio(results, 3, n = nrow(results[ploidy == 3])))

ggplot(data = simulated_ase, aes(x = sim_allelic_ratio, fill = factor(ploidy))) +
  geom_histogram(bins = 100) +
  facet_grid(factor(ploidy) ~ ., scales = "free_y")

monosomy_ar <- simulate_allelic_ratio(results, 1, n = nrow(scploid_sim_results[monosomy == TRUE]))
disomy_ar <- simulate_allelic_ratio(results, 2, n = nrow(scploid_sim_results[is.na(monosomy)]))
trisomy_ar <- simulate_allelic_ratio(results, 3, n = nrow(scploid_sim_results[monosomy == FALSE]))

scploid_sim_results[, allelic_ratio := as.numeric(NA)]
scploid_sim_results[monosomy == TRUE, allelic_ratio := monosomy_ar$sim_allelic_ratio]
scploid_sim_results[is.na(monosomy), allelic_ratio := disomy_ar$sim_allelic_ratio]
scploid_sim_results[monosomy == FALSE, allelic_ratio := trisomy_ar$sim_allelic_ratio]

# estimate variance and compute allelic imbalance z-scores
ase_iqr <- quantile(scploid_sim_results$allelic_ratio, c(0.25, 0.75))
iqr_indices <- which(scploid_sim_results$allelic_ratio > ase_iqr[1] & scploid_sim_results$allelic_ratio < ase_iqr[2])
m_iqr <- scploid_sim_results$allelic_ratio[iqr_indices]
null_allelic_ratio <- mean(m_iqr)
null_var_allelic_ratio <- {IQR(scploid_sim_results$allelic_ratio)/(2 * qnorm(.75))} ^ 2
scploid_sim_results[, ase_z := (allelic_ratio - null_allelic_ratio) / sqrt(null_var_allelic_ratio)]
scploid_sim_results[, ase_p := pnorm(ase_z)]

cor.test(scploid_sim_results$scploid_p, scploid_sim_results$ase_p, method = "kendall")

ggplot(data = scploid_sim_results, aes(x = scploid_z, y = ase_z, color = monosomy)) +
  geom_point()

fisher_wrapper <- function(pval_1, pval_2, wt_1 = 1, wt_2 = 1) {
  return(tryCatch(combine.test(p = c(pval_1, pval_2), weight = c(wt_1, wt_2), method = "fisher"), error = function(e) NA))
}

# impose effect size threshold, consistent with Griffiths et al.; set p-values to 1
scploid_sim_results[, scploid_effect_p := scploid_p]
scploid_sim_results[(scploid_score > 0.8 & scploid_score < 1.2), scploid_effect_p := 1]
scploid_sim_results[, scploid_effect_p_adj := p.adjust(scploid_effect_p, method = "BH")]
scploid_sim_results[, fisher_p := mapply(fisher_wrapper, scploid_sim_results$ase_p, scploid_sim_results$scploid_effect_p)]
scploid_sim_results[, fisher_p_adj := p.adjust(fisher_p, method = "BH")]

# cluster with k-means, assign clusters to monosomy and trisomy
km <- scploid_sim_results[fisher_p_adj < 0.01, c("scploid_z", "ase_z")] %>%
  kmeans(centers = 2)
km_clusters <- scploid_sim_results[fisher_p_adj < 0.01, c("chrom", "scploid_z", "ase_z")]
km_clusters[, cluster := km$cluster]

scploid_sim_results <- merge(scploid_sim_results, km_clusters[, c("chrom", "cluster")], "chrom", all.x = TRUE)

scploid_sim_results[, ploidy := 2]
if (mean(scploid_sim_results[cluster == 1]$ase_z) < mean(scploid_sim_results[cluster == 2]$ase_z)) {
  scploid_sim_results[cluster == 1, ploidy := 1]
  scploid_sim_results[cluster == 2, ploidy := 3] 
} else {
  scploid_sim_results[cluster == 2, ploidy := 1]
  scploid_sim_results[cluster == 1, ploidy := 3] 
}

## test performance of expression signature alone

known_aneuploidy_sim[, ploidy := as.numeric(NA)]
known_aneuploidy_sim[monosomy == TRUE, ploidy := 1]
known_aneuploidy_sim[monosomy == FALSE, ploidy := 3]

test_performance <- function(nominal_fdr, known_aneuploidy_dt, sim_results) {

  hit_chr <- paste(sim_results[scploid_effect_p_adj < nominal_fdr]$cell, 
                   sim_results[scploid_effect_p_adj < nominal_fdr]$chr,
                   sim_results[scploid_effect_p_adj < nominal_fdr]$monosomy)
  
  true_chr <- paste(known_aneuploidy_dt$cell, known_aneuploidy_dt$chr, known_aneuploidy_dt$monosomy)
  TP = sum(hit_chr %in% true_chr)
  FP = sum(!hit_chr %in% true_chr)
  FN = sum(!true_chr %in% hit_chr)
  TN = nrow(sim_results) - (TP + FP + FN)
  
  r1 <- data.table(signature = "expression",
                   sensitivity = TP/(TP + FN), 
                   precision = TP/(TP + FP), 
                   fdr = 1 - (TP/(TP + FP)), 
                   specificity = TN/(TN + FP), 
                   accuracy = (TP + TN)/(TP + TN + FP + FN), 
                   f1 = 2 * TP/(2 * TP + FP + FN), 
                   fpr = FP/(FP + TN))
  
  # test performance when allelic imbalance signature is added
  hit_chr <- paste(sim_results[fisher_p_adj < nominal_fdr]$cell, 
                   sim_results[fisher_p_adj < nominal_fdr]$chr,
                   sim_results[fisher_p_adj < nominal_fdr]$ploidy)
  
  true_chr <- paste(known_aneuploidy_dt$cell, known_aneuploidy_dt$chr, known_aneuploidy_dt$ploidy)
  TP = sum(hit_chr %in% true_chr)
  FP = sum(!hit_chr %in% true_chr)
  FN = sum(!true_chr %in% hit_chr)
  TN = nrow(sim_results) - (TP + FP + FN)
  
  r2 <- data.table(signature = "expression + allelic imbalance",
                   sensitivity = TP/(TP + FN), 
                   precision = TP/(TP + FP), 
                   fdr = 1 - (TP/(TP + FP)), 
                   specificity = TN/(TN + FP), 
                   accuracy = (TP + TN)/(TP + TN + FP + FN), 
                   f1 = 2 * TP/(2 * TP + FP + FN), 
                   fpr = FP/(FP + TN))
  
  performance_comparison <- rbind(r1, r2) %>%
    .[, nom_fdr := nominal_fdr]
  
  return(performance_comparison)
}

fdr_0.001 <- test_performance(0.001, known_aneuploidy_sim, scploid_sim_results)
fdr_0.005 <- test_performance(0.005, known_aneuploidy_sim, scploid_sim_results)
fdr_0.01 <- test_performance(0.01, known_aneuploidy_sim, scploid_sim_results)
fdr_0.05 <- test_performance(0.05, known_aneuploidy_sim, scploid_sim_results)

rbind(fdr_0.001, fdr_0.005, fdr_0.01, fdr_0.05)

###

devtools::session_info()
  
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.6.1 (2019-07-05)
# os       macOS Catalina 10.15.3      
# system   x86_64, darwin15.6.0        
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       America/New_York            
# date     2020-04-30                  
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date       lib source                                    
# assertthat             0.2.1      2019-03-21 [1] CRAN (R 3.6.0)                            
# backports              1.1.5      2019-10-02 [1] CRAN (R 3.6.0)                            
# beeswarm               0.2.3      2016-04-25 [1] CRAN (R 3.6.0)                            
# Biobase              * 2.46.0     2019-10-29 [1] Bioconductor                              
# BiocGenerics         * 0.32.0     2019-10-29 [1] Bioconductor                              
# BiocManager            1.30.10    2019-11-16 [1] CRAN (R 3.6.0)                            
# BiocNeighbors          1.4.1      2019-11-01 [1] Bioconductor                              
# BiocParallel         * 1.20.0     2019-10-29 [1] Bioconductor                              
# BiocSingular           1.2.0      2019-10-29 [1] Bioconductor                              
# BiocStyle            * 2.14.0     2019-10-29 [1] Bioconductor                              
# bitops                 1.0-6      2013-08-17 [1] CRAN (R 3.6.0)                            
# bootstrap              2019.6     2019-06-17 [1] CRAN (R 3.6.0)                            
# broom                * 0.5.2      2019-04-07 [1] CRAN (R 3.6.0)                            
# callr                  3.3.2      2019-09-22 [1] CRAN (R 3.6.0)                            
# caTools                1.17.1.2   2019-03-06 [1] CRAN (R 3.6.0)                            
# cellranger             1.1.0      2016-07-27 [1] CRAN (R 3.6.0)                            
# cli                    2.0.1      2020-01-08 [1] CRAN (R 3.6.0)                            
# colorspace             1.4-1      2019-03-18 [1] CRAN (R 3.6.0)                            
# cowplot              * 1.0.0      2019-07-11 [1] CRAN (R 3.6.0)                            
# crayon                 1.3.4      2017-09-16 [1] CRAN (R 3.6.0)                            
# data.table           * 1.12.6     2019-10-18 [1] CRAN (R 3.6.0)                            
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
# fansi                  0.4.1      2020-01-08 [1] CRAN (R 3.6.0)                            
# farver                 2.0.3      2020-01-16 [1] CRAN (R 3.6.0)                            
# fs                     1.3.1      2019-05-06 [1] CRAN (R 3.6.0)                            
# gdata                  2.18.0     2017-06-06 [1] CRAN (R 3.6.0)                            
# generics               0.0.2      2018-11-29 [1] CRAN (R 3.6.0)                            
# GenomeInfoDb         * 1.22.0     2019-10-29 [1] Bioconductor                              
# GenomeInfoDbData       1.2.2      2019-11-26 [1] Bioconductor                              
# GenomicRanges        * 1.38.0     2019-10-29 [1] Bioconductor                              
# ggbeeswarm             0.6.0      2017-08-07 [1] CRAN (R 3.6.0)                            
# ggplot2              * 3.2.1      2019-08-10 [1] CRAN (R 3.6.0)                            
# ggrepel              * 0.8.1      2019-05-07 [1] CRAN (R 3.6.0)                            
# glue                   1.3.1      2019-03-12 [1] CRAN (R 3.6.0)                            
# gmodels              * 2.18.1     2018-06-25 [1] CRAN (R 3.6.0)                            
# gplots               * 3.0.1.1    2019-01-27 [1] CRAN (R 3.6.0)                            
# gridExtra            * 2.3        2017-09-09 [1] CRAN (R 3.6.0)                            
# gtable                 0.3.0      2019-03-25 [1] CRAN (R 3.6.0)                            
# gtools                 3.8.1      2018-06-26 [1] CRAN (R 3.6.0)                            
# here                 * 0.1        2017-05-28 [1] CRAN (R 3.6.0)                            
# htmltools              0.4.0      2019-10-04 [1] CRAN (R 3.6.0)                            
# igraph                 1.2.4.1    2019-04-22 [1] CRAN (R 3.6.0)                            
# IRanges              * 2.20.1     2019-11-20 [1] Bioconductor                              
# irlba                  2.3.3      2019-02-05 [1] CRAN (R 3.6.0)                            
# KernSmooth             2.23-16    2019-10-15 [1] CRAN (R 3.6.0)                            
# knitr                  1.26       2019-11-12 [1] CRAN (R 3.6.0)                            
# labeling               0.3        2014-08-23 [1] CRAN (R 3.6.0)                            
# lattice                0.20-38    2018-11-04 [1] CRAN (R 3.6.1)                            
# lava                   1.6.6      2019-08-01 [1] CRAN (R 3.6.0)                            
# lazyeval               0.2.2      2019-03-15 [1] CRAN (R 3.6.0)                            
# lifecycle              0.1.0      2019-08-01 [1] CRAN (R 3.6.0)                            
# limma                  3.42.0     2019-10-29 [1] Bioconductor                              
# locfit                 1.5-9.1    2013-04-20 [1] CRAN (R 3.6.0)                            
# magrittr               1.5        2014-11-22 [1] CRAN (R 3.6.0)                            
# MASS                   7.3-51.4   2019-03-31 [1] CRAN (R 3.6.1)                            
# Matrix                 1.2-17     2019-03-22 [1] CRAN (R 3.6.1)                            
# matrixStats          * 0.55.0     2019-09-07 [1] CRAN (R 3.6.0)                            
# memoise                1.1.0      2017-04-21 [1] CRAN (R 3.6.0)                            
# mixtools             * 1.1.0      2017-03-10 [1] CRAN (R 3.6.0)                            
# mppa                 * 1.0        2014-08-23 [1] CRAN (R 3.6.0)                            
# MultiAssayExperiment * 1.12.0     2019-10-29 [1] Bioconductor                              
# munsell                0.5.0      2018-06-12 [1] CRAN (R 3.6.0)                            
# nlme                   3.1-142    2019-11-07 [1] CRAN (R 3.6.0)                            
# pillar                 1.4.3      2019-12-20 [1] CRAN (R 3.6.0)                            
# pkgbuild               1.0.6      2019-10-09 [1] CRAN (R 3.6.0)                            
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 3.6.0)                            
# pkgload                1.0.2      2018-10-29 [1] CRAN (R 3.6.0)                            
# plyr                 * 1.8.5      2019-12-10 [1] CRAN (R 3.6.0)                            
# prettyunits            1.0.2      2015-07-13 [1] CRAN (R 3.6.0)                            
# processx               3.4.1      2019-07-18 [1] CRAN (R 3.6.0)                            
# prodlim              * 2019.11.13 2019-11-17 [1] CRAN (R 3.6.0)                            
# ps                     1.3.0      2018-12-21 [1] CRAN (R 3.6.0)                            
# purrr                  0.3.3      2019-10-18 [1] CRAN (R 3.6.0)                            
# R6                     2.4.1      2019-11-12 [1] CRAN (R 3.6.0)                            
# Rcpp                   1.0.3      2019-11-08 [1] CRAN (R 3.6.0)                            
# RCurl                  1.95-4.12  2019-03-04 [1] CRAN (R 3.6.0)                            
# readxl               * 1.3.1      2019-03-13 [1] CRAN (R 3.6.0)                            
# remotes                2.1.0      2019-06-24 [1] CRAN (R 3.6.0)                            
# reshape2               1.4.3      2017-12-11 [1] CRAN (R 3.6.0)                            
# rlang                  0.4.4      2020-01-28 [1] CRAN (R 3.6.0)                            
# rmarkdown              1.17       2019-11-13 [1] CRAN (R 3.6.0)                            
# rmeta                  3.0        2018-03-20 [1] CRAN (R 3.6.0)                            
# rprojroot              1.3-2      2018-01-03 [1] CRAN (R 3.6.0)                            
# rstudioapi             0.10       2019-03-19 [1] CRAN (R 3.6.0)                            
# rsvd                   1.0.2      2019-07-29 [1] CRAN (R 3.6.0)                            
# S4Vectors            * 0.24.0     2019-10-29 [1] Bioconductor                              
# scales               * 1.1.0      2019-11-18 [1] CRAN (R 3.6.0)                            
# scater               * 1.14.4     2019-11-18 [1] Bioconductor                              
# scploid              * 0.9        2019-11-26 [1] Github (MarioniLab/Aneuploidy2017@286c064)
# scran                * 1.14.5     2019-11-19 [1] Bioconductor                              
# segmented              1.1-0      2019-12-10 [1] CRAN (R 3.6.0)                            
# sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 3.6.0)                            
# SingleCellExperiment * 1.8.0      2019-10-29 [1] Bioconductor                              
# statmod                1.4.32     2019-05-29 [1] CRAN (R 3.6.0)                            
# stringi                1.4.5      2020-01-11 [1] CRAN (R 3.6.0)                            
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
# usethis              * 1.5.1      2019-07-04 [1] CRAN (R 3.6.0)                            
# vctrs                  0.2.2      2020-01-24 [1] CRAN (R 3.6.0)                            
# vipor                  0.4.5      2017-03-22 [1] CRAN (R 3.6.0)                            
# viridis                0.5.1      2018-03-29 [1] CRAN (R 3.6.0)                            
# viridisLite            0.3.0      2018-02-01 [1] CRAN (R 3.6.0)                            
# withr                  2.1.2      2018-03-15 [1] CRAN (R 3.6.0)                            
# xfun                   0.11       2019-11-12 [1] CRAN (R 3.6.0)                            
# XVector                0.26.0     2019-10-29 [1] Bioconductor                              
# yaml                   2.2.0      2018-07-25 [1] CRAN (R 3.6.0)                            
# zlibbioc               1.32.0     2019-10-29 [1] Bioconductor                              
# zoo                  * 1.8-6      2019-05-28 [1] CRAN (R 3.6.0)                            
# 
# [1] /Library/Frameworks/R.framework/Versions/3.6/Resources/library
  
  
