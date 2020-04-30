list_of_packages <- c("broom", "car", "cowplot", "data.table", "dplyr", "emmeans", "fgsea",
                      "GenomicRanges", "ggbeeswarm", "ggplot2", "ggrepel", "here", "lme4", 
                      "margins", "msigdbr", "MultiAssayExperiment", "pbmcapply", "purrr",
                      "qvalue", "rtracklayer", "SCnorm", "SingleCellExperiment", "viridis")

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

here <- here::here

# http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/EMTAB3929.rds
emtab3929 <- readRDS(here("RawData/EMTAB3929.rds"))
results <- fread(here("results/aneuploidy_results.txt")) # load results data

emtab3929_gene <- experiments(emtab3929)[["gene"]]
emtab3929_count <- assays(emtab3929_gene)$count
colnames(emtab3929_count) <- gsub("_", ".", colnames(emtab3929_count))
# subset cells to those in which aneuploidies were called
emtab3929_count <- emtab3929_count[, unique(results$cell)]

if (file.exists(here("ProcessedData/emtab_sce.rds"))) {
  emtab_sce <- readRDS(here("ProcessedData/emtab_sce.rds"))
} else {
  emtab_sce <- SingleCellExperiment(assays = list(counts = emtab3929_count))
  
  scnorm_groups <- data.table(cell = colnames(assays(emtab_sce)$counts))
  scnorm_groups[, index := .I]
  scnorm_groups <- merge(scnorm_groups, results[!duplicated(cell)], "cell")
  
  # store cell annotations in colData
  colData(emtab_sce)$embryo <- scnorm_groups$embryo
  colData(emtab_sce)$EStage <- scnorm_groups$EStage
  colData(emtab_sce)$lineage <- scnorm_groups$lineage
  colData(emtab_sce)$num_estage <- scnorm_groups$num_estage
  colData(emtab_sce)$is_aneuploid <- scnorm_groups$sig_cell == 1
  
  # store gene information in RowData
  rowData(emtab_sce) <- rowData(emtab3929_gene)
  
  # set spike-ins; function is deprecated by SingleCellExperiment, but SCnorm still requires this
  isSpike(emtab_sce, type = "ERCC") <- grep("ERCC", rownames(assays(emtab_sce)$counts))
  
  par(mfrow = c(2, 2))
  # ultimately decided not to use spike-in data
  # see supplementary note S1 https://media.nature.com/original/nature-assets/nmeth/journal/v14/n6/extref/nmeth.4263-S1.pdf
  emtab_sce <- SCnorm(emtab_sce, Conditions = colData(emtab_sce)$lineage,
                      PrintProgressPlots = TRUE, NCores = 24, useSpikes = FALSE)
  saveRDS(emtab_sce, file = here("ProcessedData/emtab_sce.rds"))
}

# add chromosome location for each gene
# get reference file from conquer
conquer_ref <- readRDS(here("external_metadata/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds"))$gene_granges %>%
  as.data.table()
rowData(emtab_sce)$seqnames <- conquer_ref$seqnames

# function to get normalized counts table for a gene, merged with aneuploidy data
# filter cells with aneuploidies affecting the same chromosome
get_filtered_counts <- function(sce_object, results_data, gene_index) {
  
  cds_gene <- data.table(cell = names(assays(sce_object)$normcounts[gene_index,]), 
                         counts = assays(sce_object)$normcounts[gene_index,])
  
  if (nrow(cds_gene[counts >= 1]) < (nrow(cds_gene) / 2)) {
    stop("Too few cells with count >= 1.")
  }
  
  gene_id <- as.character(rowData(sce_object)$gene[gene_index])
  gene_symbol <- as.character(rowData(sce_object)$symbol[gene_index])
  gene_chr <- paste0("chr", as.character(rowData(sce_object)$seqnames[gene_index]))
  
  cds_gene <- suppressWarnings(merge(cds_gene, 
                                     colData(sce_object) %>% as.data.table(., keep.rownames = "cell"),
                                     "cell"))
  
  # drop cells with aneuploidies affecting the same chromosome as the location of the gene
  # being tested for differential expression
  aneuploid_cells <- unique(results_data[sig_chrom == 1 & chr == gene_chr]$cell)
  cds_gene <- cds_gene[!(cell %in% aneuploid_cells)]
  
  return(list(cds_gene, gene_id, gene_symbol))
  
}

nb_dge <- function(filtered_counts) {
  
  cds_gene <- filtered_counts[[1]]
  gene_id <- filtered_counts[[2]]
  gene_symbol <- filtered_counts[[3]]
  
  m1 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + (1 | lineage) + num_estage + is_aneuploid, 
                 nAGQ = 0)
  
  m2 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + ( 1 + is_aneuploid | lineage) + num_estage + is_aneuploid, 
                 nAGQ = 0)
  
  m2_m1_anova <- anova(m2, m1, test = "Chisq")
  m2_m1_anova <- suppressWarnings(tidy(m2_m1_anova))
  
  if (m2_m1_anova[2,]$p.value < 0.05) {
    
    results_summary <- suppressWarnings(tidy(m2)) %>%
      mutate(gene_id = gene_id) %>%
      mutate(gene_symbol = gene_symbol) %>%
      filter(group == "fixed") %>%
      as.data.table()
    
    results_summary[, model := "m2"]
    results_summary[, anova_pval := m2_m1_anova[2,]$p.value]
    
    results_margins <- suppressWarnings(summary(margins(m2))) %>%
      as.data.table()
    results_summary[, AME := as.numeric(NA)]
    results_summary[, AME_SE := as.numeric(NA)]
    results_summary[, AME_P := as.numeric(NA)]
    results_summary[term == "num_estage", AME := results_margins[factor == "num_estage"]$AME]
    results_summary[term == "num_estage", AME_SE := results_margins[factor == "num_estage"]$SE]
    results_summary[term == "num_estage", AME_P := results_margins[factor == "num_estage"]$p]
    results_summary[term == "is_aneuploidTRUE", AME := results_margins[factor == "is_aneuploid"]$AME]
    results_summary[term == "is_aneuploidTRUE", AME_SE := results_margins[factor == "is_aneuploid"]$SE]
    results_summary[term == "is_aneuploidTRUE", AME_P := results_margins[factor == "is_aneuploid"]$p]
    
  } else if (m2_m1_anova[2,]$p.value >= 0.05) {
    
    results_summary <- suppressWarnings(tidy(m1)) %>%
      mutate(gene_id = gene_id) %>%
      mutate(gene_symbol = gene_symbol) %>%
      filter(group == "fixed") %>%
      as.data.table()
    
    results_summary[, model := "m1"]
    results_summary[, anova_pval := m2_m1_anova[2,]$p.value]
    
    results_margins <- suppressWarnings(summary(margins(m1))) %>%
      as.data.table()
    results_summary[, AME := as.numeric(NA)]
    results_summary[, AME_SE := as.numeric(NA)]
    results_summary[, AME_P := as.numeric(NA)]
    results_summary[term == "num_estage", AME := results_margins[factor == "num_estage"]$AME]
    results_summary[term == "num_estage", AME_SE := results_margins[factor == "num_estage"]$SE]
    results_summary[term == "num_estage", AME_P := results_margins[factor == "num_estage"]$p]
    results_summary[term == "is_aneuploidTRUE", AME := results_margins[factor == "is_aneuploid"]$AME]
    results_summary[term == "is_aneuploidTRUE", AME_SE := results_margins[factor == "is_aneuploid"]$SE]
    results_summary[term == "is_aneuploidTRUE", AME_P := results_margins[factor == "is_aneuploid"]$p]
  }
  return(results_summary)
}

dge_wrapper <- function(sce_object, results_data, gene_index) {
  filtered_counts <- get_filtered_counts(sce_object, results_data, gene_index)
  results_summary <- nb_dge(filtered_counts)
  return(results_summary)
}

# cell-type-specific interaction model
nb_dge_interaction <- function(filtered_counts) {
  
  cds_gene <- filtered_counts[[1]]
  gene_id <- filtered_counts[[2]]
  gene_symbol <- filtered_counts[[3]]
  
  m1 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + num_estage + is_aneuploid + lineage, 
                 nAGQ = 0)
  
  m2 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + num_estage + is_aneuploid * lineage, 
                 nAGQ = 0)
  
  m2_anova <- Anova(m2, type = 3) %>%
    tidy() %>%
    as.data.table() %>%
    .[, gene_id := gene_id] %>%
    .[, gene_symbol := gene_symbol]
  
  m2_joint_test <- joint_tests(m2) %>%
    as.data.table()
  
  ame_main <- emmeans(m2, revpairwise ~ is_aneuploid | lineage, data = cds_gene)$contrasts %>%
    as.data.table() %>%
    .[, gene_id := gene_id] %>%
    .[, gene_symbol := gene_symbol] %>%
    .[, joint_test_pval := m2_joint_test[3,]$p.value]
  
  ame_interaction <- contrast(emmeans(m2, pairwise ~ is_aneuploid * lineage, data = cds_gene)[[1]], 
                              interaction = c("revpairwise")) %>%
    as.data.table() %>%
    .[, gene_id := gene_id] %>%
    .[, gene_symbol := gene_symbol] %>%
    .[, joint_test_pval := m2_joint_test[3,]$p.value]
  
  return(list(m2_anova, ame_main, ame_interaction))
}

dge_interaction_wrapper <- function(sce_object, results_data, gene_index) {
  filtered_counts <- get_filtered_counts(sce_object, results_data, gene_index)
  results_list <- nb_dge_interaction(filtered_counts)
  return(results_list)
}

# only test genes with expression in more than half of cells
hi_ex_indices <- unname(which(rowSums(assays(emtab_sce)$normcounts > 1) > (nrow(colData(emtab_sce)) / 2)))

if (file.exists(here("results/dge_dt.txt"))) {
  dge_dt <- fread(here("results/dge_dt.txt"))
} else {
  # this code captures warnings output by lme4
  dge_results <- do.call(list, pbmclapply(hi_ex_indices, 
                                          function(x) {
                                            r <- tryCatch(
                                              withCallingHandlers(
                                                {
                                                  error_text <- "No error."
                                                  list(value = dge_wrapper(emtab_sce, results, x), error_text = error_text)
                                                }, 
                                                warning = function(e) {
                                                  error_text <<- trimws(paste0("WARNING: ", e))
                                                  invokeRestart("muffleWarning")
                                                }), 
                                              error = function(e) {
                                                return(list(value = NA, error_text = trimws(paste0("ERROR: ", e))))
                                              }, 
                                              finally = {
                                              }
                                            )
                                            return(r)
                                          }, mc.cores = 24, ignore.interactive = getOption("ignore.interactive", T)))
  
  warning_text <- sapply(dge_results, function(x) x[2])
  no_warning_models <- which(grepl("No error.", warning_text))
  dge_output <- sapply(dge_results, function(x) x[1])
  no_warning_output <- dge_output[no_warning_models]
  
  dge_dt <- rbindlist(no_warning_output[unlist(lapply(no_warning_output, is.data.table))]) %>%
    setorder(., p.value)
  
  dge_dt[term == "is_aneuploidTRUE" & !grepl("ERCC", gene_id)]
  
  fwrite(dge_dt, file = here("results/dge_dt.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

dge_aneuploidy <- dge_dt[term == "is_aneuploidTRUE"] %>%
  filter(!(grepl("ERCC", gene_id))) %>%
  setorder(., statistic) %>%
  as.data.table()

dge_aneuploidy[, q.value := qvalue(dge_aneuploidy$p.value)$qvalues]
nrow(dge_aneuploidy[q.value < 0.05])

# interaction model
dge_results_interaction <- do.call(list, pbmclapply(hi_ex_indices, 
                                                    function(x) {
                                                      r <- tryCatch(
                                                        withCallingHandlers(
                                                          {
                                                            error_text <- "No error."
                                                            list(value = dge_interaction_wrapper(emtab_sce, results, x), error_text = error_text)
                                                          }, 
                                                          warning = function(e) {
                                                            error_text <<- trimws(paste0("WARNING: ", e))
                                                            invokeRestart("muffleWarning")
                                                          }), 
                                                        error = function(e) {
                                                          return(list(value = NA, error_text = trimws(paste0("ERROR: ", e))))
                                                        }, 
                                                        finally = {
                                                        }
                                                      )
                                                      return(r)
                                                    }, mc.cores = 48, ignore.interactive = getOption("ignore.interactive", T)))

no_warning_models <- which(grepl("No error.", map(dge_results_interaction, 2)))

cell_type_specific_anova <- sapply(map(dge_results_interaction, 1), function(x) x[1])[no_warning_models] %>% 
  rbindlist() %>%
  .[term == "is_aneuploid:lineage"] %>%
  setorder(p.value, gene_id) %>%
  .[, q.value := qvalue(p.value)$qvalues]

nrow(cell_type_specific_anova[q.value < 0.05])

cell_type_specific_main <- sapply(map(dge_results_interaction, 1), function(x) x[2])[no_warning_models] %>% 
  rbindlist() %>%
  setorder(joint_test_pval, gene_id)

cell_type_specific_interaction <- sapply(map(dge_results_interaction, 1), function(x) x[3])[no_warning_models] %>% 
  rbindlist() %>%
  setorder(joint_test_pval, gene_id)

# plot interaction
plot_interaction <- function(sce_object, results_data, gene_id) {
  
  gene_index <- which(rowData(sce_object)$gene == gene_id)
  
  filtered_counts <- get_filtered_counts(sce_object, results_data, gene_index)
  
  cds_gene <- filtered_counts[[1]]
  gene_id <- filtered_counts[[2]]
  gene_symbol <- filtered_counts[[3]]
  
  m1 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + num_estage + is_aneuploid + lineage, 
                 nAGQ = 0)
  
  m2 <- glmer.nb(data = cds_gene, 
                 formula = round(counts + 1) ~ (1 | embryo) + num_estage + is_aneuploid * lineage, 
                 nAGQ = 0)
  
  emmip_plot <- emmip(m2, lineage ~ is_aneuploid, CIs = TRUE, type = "response") +
    theme_classic() +
    scale_x_discrete(labels = c('Euploid','Aneuploid')) +
    scale_color_manual(values = c('#66a61e', '#d95f02','#e7298a','#e6ab02', '#7570b3', '#1b9e77'), name = "Cell type") +
    xlab("") +
    ylab("Predicted expression (normalized counts)") +
    ggtitle(gene_symbol) +
    theme(plot.title = element_text(face = "italic"))
  
  m2_joint_test <- joint_tests(m2) %>%
    as.data.table()
  
  ame_main <- emmeans(m2, revpairwise ~ is_aneuploid | lineage, data = cds_gene)$contrasts %>%
    as.data.table()
  
  ame_interaction <- contrast(emmeans(m2, pairwise ~ is_aneuploid * lineage, data = cds_gene)[[1]], 
                              interaction = c("revpairwise")) %>%
    as.data.table()
  
  return(list(m2_joint_test, ame_main, ame_interaction, emmip_plot))
}

plot_interaction(emtab_sce, results, cell_type_specific_anova[1,]$gene_id)[[4]]
plot_interaction(emtab_sce, results, "ENSG00000107485.15")[[4]]

## enrichment analysis

# get hallmark gene sets
m_df <- msigdbr(species = "Homo sapiens", category = "H")

m_list <- m_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

# rank on test statistic (i.e., signed p-value)
ranks <- dge_aneuploidy$statistic
names(ranks) <- dge_aneuploidy$gene_symbol

fgsea_results <- fgseaMultilevel(m_list, stats = ranks, minSize = 10, maxSize = 500, nproc = 1) %>%
  setorder(pval, NES)

gsea_all <- ggplot(fgsea_results[padj < 0.05], 
                   aes(x = reorder(gsub("_", " ", gsub("HALLMARK_", "", pathway)), NES), y = NES, fill = -log10(pval))) +
  geom_col() +
  coord_flip() +
  labs(x= "Pathway", y = "Normalized Enrichment Score") + 
  theme_minimal() +
  scale_fill_viridis(limits = c(0, 50), name = expression(-log[10](italic(p))), end = 0.9)

plotEnrichment <- function(pathway, stats, gseaParam = 1, ticksSize = 0.2, line_color = "green") {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = line_color, 
                                                      size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                               linetype = "dashed") + geom_hline(yintercept = min(bottoms), 
                                                                                                                 colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                                                                   colour = "black") + geom_line(color = line_color) + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
                                                               y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score")
  g
}

gsea_sub1 <- plotEnrichment(m_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranks, line_color = "#d95f02") + 
  labs(title = "TNFA SIGNALING VIA NFKB", size = 2) +
  theme(plot.title = element_text(size = 10)) +
  ylab("Enrichment Score") +
  xlab("Gene Rank")

gsea_sub2 <- plotEnrichment(m_list[["HALLMARK_MYC_TARGETS_V1"]], ranks, line_color = "#1b9e77") + 
  labs(title = "MYC TARGETS V1") +
  theme(plot.title = element_text(size = 10)) +
  ylab("Enrichment Score") +
  xlab("Gene Rank")

gsea_sub3 <- plotEnrichment(m_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks, line_color = "#7570b3") + 
  scale_color_manual(values = "red") +
  labs(title = "OXIDATIVE PHOSPHORYLATION") +
  theme(plot.title = element_text(size = 10)) +
  ylab("Enrichment Score") +
  xlab("Gene Rank")

le_s1 <- unlist(fgsea_results[pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]$leadingEdge)
le_s2 <- unlist(fgsea_results[pathway == "HALLMARK_MYC_TARGETS_V1"]$leadingEdge)
le_s3 <- unlist(fgsea_results[pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]$leadingEdge)

dge_aneuploidy[, pathway_label := as.character(NA)]
dge_aneuploidy[gene_symbol %in% le_s1, pathway_label := "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
dge_aneuploidy[gene_symbol %in% le_s2, pathway_label := "HALLMARK_MYC_TARGETS_V1"]
dge_aneuploidy[gene_symbol %in% le_s3, pathway_label := "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]

volcano <- ggplot() +
  geom_point(data = dge_aneuploidy[is.na(pathway_label)], 
             aes(x = estimate, y = -log10(p.value)), color = "gray") +
  geom_point(data = dge_aneuploidy[pathway_label == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"], 
             aes(x = estimate, y = -log10(p.value)), color = "#d95f02") +
  geom_point(data = dge_aneuploidy[pathway_label == "HALLMARK_MYC_TARGETS_V1"], 
             aes(x = estimate, y = -log10(p.value)), color = "#1b9e77") +
  geom_point(data = dge_aneuploidy[pathway_label == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"], 
             aes(x = estimate, y = -log10(p.value)), color = "#7570b3") +
  theme_classic() +
  scale_color_manual(name = "", 
                     values = c("gray", "#d95f02", "#7570b3", "#e7298a"),
                     labels = c("", "TNF signaling via  NF-kB", "Myc targets v1", "Oxidative phosphorylation")) +
  ylab(expression(-log[10](italic(p)))) +
  xlab(expression(hat(italic("\u03B2")))) +
  ylim(0, 18) +
  xlim(-1.1, 1.5) +
  theme(legend.position = c(0.27, 0.93),
        legend.background = element_blank()) +
  geom_text_repel(data = dge_aneuploidy[p.value < 1e-8], 
                  aes(label = gene_symbol, x = estimate, y = -log10(p.value)), size = 3, fontface = "italic")

plot_grid(plot_grid(gsea_all, volcano, labels = c("A", "B"), ncol = 1, rel_heights = c(0.7, 1)),
          plot_grid(gsea_sub1, gsea_sub2, gsea_sub3, ncol = 1, labels = c("C", "D", "E")), 
          ncol = 2, rel_widths = c(1, 0.5))

plot_dge <- function(sce_object, results_data, gene_symbol) {
  
  cds_gene <- data.table(cell = names(assays(sce_object)$normcounts[which(rowData(sce_object)$symbol == gene_symbol),]), 
                         counts = assays(sce_object)$normcounts[which(rowData(sce_object)$symbol == gene_symbol),])
  
  cds_gene <- suppressWarnings(merge(cds_gene, 
                                     colData(sce_object) %>% as.data.table(., keep.rownames = "cell"),
                                     "cell"))
  
  gene_chr <- paste0("chr", conquer_ref[symbol == gene_symbol]$seqnames)
  aneuploid_cells <- unique(results_data[sig_chrom == 1 & chr == gene_chr]$cell)
  
  summary(glmer.nb(data = cds_gene, formula = round(counts + 1) ~ (1 | embryo) + (1 | lineage) + num_estage + is_aneuploid, nAGQ = 0))
  
  cds_gene$lineage <- factor(cds_gene$lineage, 
                             levels = c("Undefined", "ICM", "Trophectoderm", "Intermediate", "Epiblast", "Primitive Endoderm"))
  
  p1 <- ggplot(data = cds_gene[counts != 0], aes(x = EStage, y = (counts + 1), color = is_aneuploid)) +
    theme_classic() +
    geom_beeswarm(dodge.width = 1, size = 0.2) +
    #geom_violin() +
    #scale_y_log10() +
    facet_wrap(~ lineage, nrow = 3) +
    scale_color_manual(values = c("#4e79a7", "#f28e2b"), labels = c("Euploid", "Aneuploid"), name = "") +
    ylab("Normalized counts (zeros excluded)") +
    xlab("Embryonic stage") +
    theme(plot.title = element_text(face = "italic")) +
    ggtitle(gene_symbol) 
  
  return(p1)
}

dge_gdf15 <- plot_dge(emtab_sce, results, "GDF15") + scale_y_log10()
dge_zfp42 <- plot_dge(emtab_sce, results, "ZFP42")

dge_aneuploidy_to_supplement <- setorder(dge_aneuploidy, p.value) %>%
  select(c("gene_id", "gene_symbol", "estimate", "std.error", "p.value", "q.value", "AME", "AME_SE", "AME_P")) %>%
  setnames(c("ensembl_id", "symbol", "beta", "beta_se", "beta_p", "beta_q", "ame", "ame_se", "ame_p"))

fwrite(dge_aneuploidy_to_supplement, file = here("results/dge_results.txt"), 
       quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

dge_interaction_to_supplement <- setorder(cell_type_specific_anova, p.value, gene_symbol) %>%
  select(c("gene_id", "gene_symbol", "term", "statistic", "df", "p.value", "q.value")) %>%
  setnames(c("ensembl_id", "symbol", "term", "chi_sq", "df", "anova_p", "anova_q"))

fwrite(dge_interaction_to_supplement, file = here("results/dge_celltype_interaction_results.txt"), 
       quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

###

devtools::session_info()

# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 3.6.1 (2019-07-05)   
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2019-12-30
# 
# ─ Packages ───────────────────────────────────────────────────────────────────
# package              * version   date       lib source
# assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
# backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
# beeswarm               0.2.3     2016-04-25 [1] CRAN (R 3.6.1)
# Biobase              * 2.46.0    2019-10-29 [1] Bioconductor
# BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
# BiocParallel         * 1.20.0    2019-10-30 [1] Bioconductor
# Biostrings             2.54.0    2019-10-29 [1] Bioconductor
# bitops                 1.0-6     2013-08-17 [1] CRAN (R 3.6.1)
# boot                   1.3-22    2019-04-02 [2] CRAN (R 3.6.1)
# broom                * 0.5.2     2019-04-07 [1] CRAN (R 3.6.1)
# callr                  3.3.2     2019-09-22 [1] CRAN (R 3.6.1)
# cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
# cluster                2.1.0     2019-06-19 [2] CRAN (R 3.6.1)
# colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
# cowplot              * 1.0.0     2019-07-11 [1] CRAN (R 3.6.1)
# crayon                 1.3.4     2017-09-16 [2] CRAN (R 3.6.1)
# data.table           * 1.12.8    2019-12-09 [1] CRAN (R 3.6.1)
# DelayedArray         * 0.12.1    2019-12-17 [1] Bioconductor
# desc                   1.2.0     2018-05-01 [1] CRAN (R 3.6.1)
# devtools               2.2.1     2019-09-24 [1] CRAN (R 3.6.1)
# digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)
# dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
# ellipsis               0.3.0     2019-09-20 [1] CRAN (R 3.6.1)
# fansi                  0.4.0     2018-10-05 [2] CRAN (R 3.6.1)
# farver                 2.0.1     2019-11-13 [1] CRAN (R 3.6.1)
# fastmatch              1.1-0     2017-01-28 [1] CRAN (R 3.6.1)
# fgsea                * 1.12.0    2019-10-29 [1] Bioconductor
# fs                     1.3.1     2019-05-06 [1] CRAN (R 3.6.1)
# generics               0.0.2     2018-11-29 [1] CRAN (R 3.6.1)
# GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
# GenomeInfoDbData       1.2.2     2019-12-07 [1] Bioconductor
# GenomicAlignments      1.22.1    2019-11-12 [1] Bioconductor
# GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
# ggbeeswarm           * 0.6.0     2017-08-07 [1] CRAN (R 3.6.1)
# ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
# ggrepel              * 0.8.1     2019-05-07 [1] CRAN (R 3.6.1)
# glue                   1.3.1     2019-03-12 [2] CRAN (R 3.6.1)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 3.6.1)
# gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
# here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
# IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
# labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
# lattice                0.20-38   2018-11-04 [2] CRAN (R 3.6.1)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
# lifecycle              0.1.0     2019-08-01 [2] CRAN (R 3.6.1)
# lme4                 * 1.1-21    2019-03-05 [1] CRAN (R 3.6.1)
# magrittr               1.5       2014-11-22 [2] CRAN (R 3.6.1)
# margins              * 0.3.23    2018-05-22 [1] CRAN (R 3.6.1)
# MASS                   7.3-51.4  2019-03-31 [2] CRAN (R 3.6.1)
# Matrix               * 1.2-17    2019-03-22 [2] CRAN (R 3.6.1)
# MatrixModels           0.4-1     2015-08-22 [1] CRAN (R 3.6.1)
# matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
# memoise                1.1.0     2017-04-21 [2] CRAN (R 3.6.1)
# minqa                  1.2.4     2014-10-09 [1] CRAN (R 3.6.1)
# moments                0.14      2015-01-05 [1] CRAN (R 3.6.1)
# msigdbr              * 7.0.1     2019-09-04 [1] CRAN (R 3.6.1)
# MultiAssayExperiment * 1.12.0    2019-10-29 [1] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
# nlme                   3.1-140   2019-05-12 [2] CRAN (R 3.6.1)
# nloptr                 1.2.1     2018-10-03 [1] CRAN (R 3.6.1)
# pbmcapply            * 1.5.0     2019-07-10 [1] CRAN (R 3.6.1)
# pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
# pkgbuild               1.0.6     2019-10-09 [1] CRAN (R 3.6.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
# pkgload                1.0.2     2018-10-29 [1] CRAN (R 3.6.1)
# plyr                   1.8.5     2019-12-10 [1] CRAN (R 3.6.1)
# prediction             0.3.14    2019-06-17 [1] CRAN (R 3.6.1)
# prettyunits            1.0.2     2015-07-13 [2] CRAN (R 3.6.1)
# processx               3.4.1     2019-07-18 [1] CRAN (R 3.6.1)
# ps                     1.3.0     2018-12-21 [1] CRAN (R 3.6.1)
# purrr                  0.3.3     2019-10-18 [1] CRAN (R 3.6.1)
# quantreg               5.52      2019-11-09 [1] CRAN (R 3.6.1)
# qvalue               * 2.18.0    2019-10-29 [1] Bioconductor
# R6                     2.4.1     2019-11-12 [1] CRAN (R 3.6.1)
# Rcpp                 * 1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
# RCurl                  1.95-4.12 2019-03-04 [1] CRAN (R 3.6.1)
# remotes                2.1.0     2019-06-24 [1] CRAN (R 3.6.1)
# reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.6.1)
# rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
# rprojroot              1.3-2     2018-01-03 [1] CRAN (R 3.6.1)
# Rsamtools              2.2.1     2019-11-06 [1] Bioconductor
# rtracklayer          * 1.46.0    2019-10-29 [1] Bioconductor
# S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
# scales                 1.1.0     2019-11-18 [1] CRAN (R 3.6.1)
# SCnorm               * 1.8.2     2019-11-22 [1] Bioconductor
# sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
# SingleCellExperiment * 1.8.0     2019-10-29 [1] Bioconductor
# SparseM                1.77      2017-04-23 [1] CRAN (R 3.6.1)
# stringi                1.4.3     2019-03-12 [2] CRAN (R 3.6.1)
# stringr                1.4.0     2019-02-10 [2] CRAN (R 3.6.1)
# SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor
# testthat               2.3.1     2019-12-01 [1] CRAN (R 3.6.1)
# tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
# tidyr                  1.0.0     2019-09-11 [1] CRAN (R 3.6.1)
# tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
# usethis                1.5.1     2019-07-04 [1] CRAN (R 3.6.1)
# vctrs                  0.2.1     2019-12-17 [1] CRAN (R 3.6.1)
# vipor                  0.4.5     2017-03-22 [1] CRAN (R 3.6.1)
# viridis              * 0.5.1     2018-03-29 [1] CRAN (R 3.6.1)
# viridisLite          * 0.3.0     2018-02-01 [2] CRAN (R 3.6.1)
# withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
# XML                    3.98-1.20 2019-06-06 [1] CRAN (R 3.6.1)
# XVector                0.26.0    2019-10-29 [1] Bioconductor
# zeallot                0.1.0     2018-01-28 [2] CRAN (R 3.6.1)
# zlibbioc               1.32.0    2019-10-29 [1] Bioconductor
# 
# [1] /home-net/home-4/rmccoy22@jhu.edu/R/x86_64-pc-linux-gnu-library/3.6/gcc/5.5
# [2] /software/apps/R/3.6.1/gcc/5.5.0/lib64/R/library
