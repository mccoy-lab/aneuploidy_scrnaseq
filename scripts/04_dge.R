list_of_packages <- c("broom", "cowplot", "data.table", "dplyr", "fgsea", "GenomicRanges", 
                      "ggbeeswarm", "ggplot2", "ggrepel", "here", "lme4", "margins", 
                      "msigdbr", "MultiAssayExperiment", "pbmcapply", "qvalue", "rtracklayer", 
                      "SCnorm", "SingleCellExperiment", "viridis")

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
                      PrintProgressPlots = TRUE, NCores = 48, useSpikes = FALSE)
  saveRDS(emtab_sce, file = here("ProcessedData/emtab_sce.rds"))
}

# add chromosome location for each gene
# get reference file from conquer
conquer_ref <- readRDS(here("external_metadata/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds"))$gene_granges %>%
  as.data.table()
rowData(emtab_sce)$seqnames <- conquer_ref$seqnames

nb_dge <- function(sce_object, results_data, gene_index, slope = "fixed") {
  
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
  
  if (slope == "fixed") {  
    
    m1 <- glmer.nb(data = cds_gene, 
                   formula = round(counts + 1) ~ (1 | embryo) + (1 | lineage) + num_estage + is_aneuploid, 
                   nAGQ = 0) 
    
    m1 <- suppressWarnings(tidy(m1)) %>%
      mutate(gene_id = gene_id) %>%
      filter(group == "fixed") %>%
      as.data.table()
    
    return(m1)
    
  } else if (slope == "random") {
    
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
                                                  list(value = nb_dge(emtab_sce, results, x, slope = "random"), error_text = error_text)
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

###

devtools::session_info()

