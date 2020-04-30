## Contents

`00_get_data.sh`  -  Download Petropoulos et al. (2016) data from conquer, as well as cell type annotations from Stirparo et al. (2018).

`01_data_prep.R`  -  Preprocess expression data and cell type annotations in accordance with `scploid` methods (Griffiths et al., 2017).

`02_detect_aneuploidy.R`  -  Infer aneuploidy in single cells based on signatures of gene expression alteration and allelic imbalance.

`03_scdna.R`  -  Replicate qualitative patterns of mosaic aneuploidy using published calls from Zhu et al. (2018) and Zhou et al. (2019).

`04_dge.R`  -  Perform differential expression analysis to quantify global responses to aneuploidy in single cells of human embryos.

`05_simulation.R`  -  Perform simulations of single-cell gene expression (adapted from Griffiths et al., 2017) and allelic imbalance impacts of aneuploidy to assess performance when combining these signatures.

## References

Griffiths JA, Scialdone A, Marioni JC. 2017. Mosaic autosomal aneuploidies are detectable from single-cell RNAseq data. *BMC Genomics* 18: 904.

Petropoulos S, Edsgärd D, Reinius B, Deng Q, Panula SP, Codeluppi S, Plaza Reyes A, Linnarsson S, Sandberg R, Lanner F. 2016. Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos. *Cell* 165: 1012–1026.

Stirparo GG, Boroviak T, Guo G, Nichols J, Smith A, Bertone P. 2018. Integrated analysis of single-cell embryo data yields a unified transcriptome signature for the human pre-implantation epiblast. *Development* 145: dev158501.

Zhou F, Wang R, Yuan P, Ren Y, Mao Y, Li R, Lian Y, Li J, Wen L, Yan L, et al. 2019. Reconstituting the transcriptome and DNA methylome landscapes of human implantation. *Nature* 572: 660–664.

Zhu P, Guo H, Ren Y, Hou Y, Dong J, Li R, Lian Y, Fan X, Hu B, Gao Y, et al. 2018. Single-cell DNA methylome sequencing of human preimplantation embryos. *Nat Genet* 50: 12–19.
