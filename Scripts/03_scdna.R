list_of_packages <- c("cowplot", "data.table", "dplyr", "here", "httr", 
                      "ggdendro", "lme4", "margins", "readxl", "tidyverse")

# Easily install and load packages
install_and_load_packages <- function(pkg){
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
  if(length(new_packages))
    install.packages(new_packages, dependencies = TRUE)
  sapply(list_of_packages, require, character.only = TRUE)
}

install_and_load_packages(list_of_packages)

# zhu et al. 2018
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-017-0007-6/MediaObjects/41588_2017_7_MOESM3_ESM.xlsx"
# to read from URL, uncomment lines below
GET(url, write_disk(tf <- tempfile(fileext = ".xlsx")))

# read and wrangle the data
dt <- read_excel(tf, na = "NA") %>%
  rename_all(list(~ make.names(.))) %>%
  mutate(original_sample_name = Sample_Name) %>%
  filter(Single.Cell.or.Bulk.cell == "Single cell") %>%
  filter(Developmental.Time %in% c("2-cell", "4-cell", "8-cell", "Morula", "Blastocyst")) %>%
  mutate(Dev.Time.Numeric = as.integer(NA)) %>%
  mutate(Dev.Time.Numeric = replace(Dev.Time.Numeric, Developmental.Time == "2-cell", 1)) %>%
  mutate(Dev.Time.Numeric = replace(Dev.Time.Numeric, Developmental.Time == "4-cell", 2)) %>%
  mutate(Dev.Time.Numeric = replace(Dev.Time.Numeric, Developmental.Time == "8-cell", 3)) %>%
  mutate(Dev.Time.Numeric = replace(Dev.Time.Numeric, Developmental.Time == "Morula", 4)) %>%
  mutate(Dev.Time.Numeric = replace(Dev.Time.Numeric, Developmental.Time == "Blastocyst", 5)) %>%
  mutate(is_aneuploid = Chromosome.segregation.error %in% c("mitosis error", "meiosis error")) %>%
  mutate(embryo = as.character(NA)) %>%
  as.data.table()

dt[, Sample_Name := gsub("−", "-", Sample_Name)]
dt[Sample_Name == "Te03-W-s53", Copy.number.of.varied.chromosomes := "1,1,0,1"]

dt[Developmental.Time == "2-cell", embryo := paste0("2C_", map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][3]))]
dt[Developmental.Time == "4-cell", embryo := paste0("4C_", map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][3]))]
dt[Developmental.Time == "4-cell" & grepl("scBS-E-", Sample_Name), embryo := "4C_E"]
dt[Developmental.Time == "8-cell", embryo := paste0("8C_", map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][3]))]
dt[Developmental.Time == "8-cell" & grepl("scBS-Q-", Sample_Name), embryo := "8C_Q"]
dt[Developmental.Time == "Morula", embryo := map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][1])]
dt[Developmental.Time == "Morula" & grepl("scBS-Morula", Sample_Name), embryo := paste0("Mor0", map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][3]))]

dt[Developmental.Time == "Blastocyst" & grepl("scBS-O-*.*-14s", Sample_Name), embryo := "Blast_O_14"]
dt[Developmental.Time == "Blastocyst" & grepl("scBS-O-*.*-5s", Sample_Name), embryo := "Blast_O_5"]
dt[Developmental.Time == "Blastocyst" & grepl("scBS-L-*.*2-s", Sample_Name), embryo := "Blast_L_2"]
dt[Developmental.Time == "Blastocyst" & grepl("scBS-Q-s", Sample_Name), embryo := "Blast_Q"]

dt[Developmental.Time == "Blastocyst" & is.na(embryo), embryo := paste(
  "Blast",
  map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][2]),
  gsub("Icm|Te|Bst", "", map_chr(Sample_Name, function(s) strsplit(s, "-")[[1]][1])),
  sep = "_")]

dt[, lineage := as.character(NA)]
dt[grepl("ICM|Icm", Sample_Name), lineage := "ICM"]
dt[grepl("TE|Te", Sample_Name), lineage := "TE"]

length(unique(dt[Developmental.Time %in% c("Morula", "Blastocyst")]$embryo))
length(unique(dt[Developmental.Time %in% c("Morula", "Blastocyst") & is_aneuploid == TRUE]$embryo))
length(unique(dt[Developmental.Time %in% c("Morula", "Blastocyst") & Chromosome.segregation.error == "mitosis error"]$embryo))
length(unique(dt[Developmental.Time %in% c("Morula", "Blastocyst") & Chromosome.segregation.error == "meiosis error"]$embryo))

### test association with developmental stage and cell type

dt$lineage <- factor(dt$lineage, ordered = FALSE)
dt$lineage <- relevel(dt$lineage, ref = "TE")

m1 <- glmer(data = dt[lineage %in% c("TE", "ICM")], formula = (is_aneuploid == TRUE) ~ (1 | embryo) + lineage, family = binomial)
summary(margins(m1))

### plot heatmaps

dt <- dt[Developmental.Time %in% c("Morula", "Blastocyst")] # restrict to morula and blastocyst stage embryos

aneuploid_index <- which(!is.na(dt$Chromosome.of.CNVs))

get_all_aneuploidies <- function(row_index) {
  aneuploid_embryo <- dt[row_index,]$embryo
  aneuploid_cell <- dt[row_index,]$Sample_Name
  aneuploid_chroms <- str_split(dt[row_index,]$Chromosome.of.CNVs, ",", simplify = TRUE)
  aneuploid_chroms_cn <- as.numeric(str_split(dt[row_index,]$Copy.number.of.varied.chromosomes, ",", simplify = TRUE))
  return(data.table(embryo = aneuploid_embryo, Sample_Name = aneuploid_cell, chr_num = as.integer(gsub("chr", "", aneuploid_chroms)), ploidy = aneuploid_chroms_cn))
}

aneuploid_results <- do.call(rbind, lapply(aneuploid_index, function(x) get_all_aneuploidies(x)))
aneuploid_results <- aneuploid_results[!is.na(chr_num)]

euploid_results <- data.table(Sample_Name = rep(dt$Sample_Name, each = 22), chr_num = 1:22, ploidy = 2)
euploid_results <- merge(dt[, c("embryo", "Sample_Name")], euploid_results, by = "Sample_Name")

results <- rbind(aneuploid_results, euploid_results)[!duplicated(paste(embryo, Sample_Name, chr_num, sep = "_"))] %>%
  setnames(., "Sample_Name", "cell") %>%
  setorder(., embryo, cell, chr_num)

plot_mca_sig <- function(embryo_id, dt) {
  heatmap_data <- pivot_wider(results[(embryo == embryo_id), c("chr_num", "cell", "ploidy")], names_from = chr_num, values_from = ploidy)
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$cell
  heatmap_matrix[heatmap_matrix == 2] <- -5

  distance.row <- dist(heatmap_matrix, method = "euclidean") # same parameters as honeyBADGER
  cluster.row <- hclust(distance.row, method = "ward.D") # same parameters as honeyBADGER

  dendrogram <- ggplot(segment(dendro_data(cluster.row))) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    coord_flip() +
    scale_y_reverse(expand=c(0.2, 0)) +
    theme_dendro() +
    scale_x_reverse() +
    theme(plot.margin = unit(c(0, 0, 0, 10), "pt")) +
    ylab("Cells") +
    NULL

  sample_order <- rev(cluster.row$labels[cluster.row$order])

  dt_to_plot <- results[embryo == embryo_id]
  dt_to_plot$cell <- factor(dt_to_plot$cell, levels = sample_order)
  dt_to_plot$ploidy_factor <- factor(dt_to_plot$ploidy, levels = c(1, 2, 3))

  dt_to_plot$chr_factor <- factor(dt_to_plot$chr_num, levels = 1:22)
  heatmap <- ggplot(data = dt_to_plot, aes(x = chr_factor, y = cell, fill = ploidy_factor)) +
    geom_tile() +
    theme_bw() +
    scale_fill_manual(name = "Copy Number", values = c("blue", "white", "red"), drop = FALSE) +
    xlab("Chromosome") +
    ylab("") +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(5.5, 5.5, 5.5, -5), "pt"),
          panel.grid = element_blank(),
          legend.position = "none")

  plot_grid(dendrogram, heatmap, align = "h", axis = "b", rel_widths = c(0.3, 1), scale = c(1, 1)) +
    draw_label("Cell", x = 0, y = 0.5, vjust = 1.5, angle = 90, size = 11)
}

legend_dt <- data.table(x = 1:2, y = 1:2, ploidy = c("Monosomy", "Trisomy"))
legend_dt$ploidy <- factor(legend_dt$ploidy, levels = c("Monosomy", "Trisomy"))
legend <- get_legend(ggplot(data = legend_dt,
                            aes(x = x, y = y, fill = ploidy)) +
                       geom_tile() +
                       theme_bw() +
                       scale_fill_manual(name = "", values = c("blue", "red"), drop = FALSE))

embryo_heatmaps <- do.call(list,
                           lapply(unique(dt$embryo),
                                  function(x) try(plot_mca_sig(x))))

embryo_heatmaps[[length(embryo_heatmaps) + 1]] <- legend

# plot_grid(plotlist = embryo_heatmaps,
#           labels = unique(dt$embryo),
#           ncol = 3,
#           scale = 0.8)

# zhou et al. 2019; scDNA-seq

url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1500-0/MediaObjects/41586_2019_1500_MOESM3_ESM.zip"
download.file(url, here("41586_2019_1500_MOESM3_ESM.zip"))
unzip(here("41586_2019_1500_MOESM3_ESM.zip"), files = "Supplementary Table 1 Single cells and embryos information.xlsx")

# to read from URL, uncomment lines below
dt_zhou <- read_xlsx(here("Supplementary Table 1 Single cells and embryos information.xlsx"), sheet = 2, skip = 2, na = "") %>%
  as.data.table()
dt_zhou <- dt_zhou[, c(2:4, 11:16)]
dt_zhou <- dt_zhou[!is.na(Sample) & Sample != "Total"] %>%
  setnames(., c("embryo", "volunteer", "sex", 
                "epi_euploid", "pe_euploid", "te_euploid",
                "epi_aneuploid", "pe_aneuploid", "te_aneuploid"))

dt_zhou <- melt(dt_zhou, id.vars = c("embryo", "volunteer", "sex")) %>%
  uncount(value)

dt_zhou[, c("lineage", "ploidy") := tstrsplit(variable, "_", fixed = TRUE)]
dt_zhou[, c("tri", "hv", "day", "ivc", "e") := tstrsplit(embryo, "_", fixed = TRUE)]
dt_zhou[, num_estage := as.numeric(gsub("D", "", day))]
dt_zhou[, is_aneuploid := (ploidy == "aneuploid")]
dt_zhou[, is_trophectoderm :=  (lineage == "te")]

m2 <- glmer(data = dt_zhou, formula = is_aneuploid ~ (1 | embryo) + num_estage + is_trophectoderm, family = "binomial")
summary(margins(m2))

m3 <- glmer(data = dt_zhou, formula = is_aneuploid ~ (1 | embryo) + num_estage * is_trophectoderm, family = "binomial")
summary(m3)

# zhou et al. 2019; scRNA-seq

unzip(here("41586_2019_1500_MOESM3_ESM.zip"), files = "Supplementary Table 2 Sample Information.xlsx")

dt_zhou_rna <- read_xlsx(here("Supplementary Table 2 Sample Information.xlsx"), sheet = 1, na = "") %>%
  as.data.table()

dt_zhou_rna[, is_aneuploid := (CNV == "Abnormal")]
dt_zhou_rna[, num_estage := as.numeric(gsub("D", "", Day))]
dt_zhou_rna[, is_trophectoderm := (Lineage == "TE")]
dt_zhou_rna[, is_te_factor := as.character(NA)]
dt_zhou_rna[is_trophectoderm == TRUE, is_te_factor := "TE"]
dt_zhou_rna[is_trophectoderm == FALSE, is_te_factor := "EPI/PE"]
dt_zhou_rna[, is_aneuploid_numeric := as.numeric(NA)]
dt_zhou_rna[is_aneuploid == TRUE, is_aneuploid_numeric := 1]
dt_zhou_rna[is_aneuploid == FALSE, is_aneuploid_numeric := 0]

m4 <- glmer(data = dt_zhou_rna, formula = is_aneuploid ~ (1 | Ori_Day_Emb) + num_estage + is_trophectoderm, family = "binomial")
summary(margins(m4))

m5 <- glmer(data = dt_zhou_rna, formula = is_aneuploid ~ (1 | Ori_Day_Emb) + num_estage * is_trophectoderm, family = "binomial")
summary(m5)

