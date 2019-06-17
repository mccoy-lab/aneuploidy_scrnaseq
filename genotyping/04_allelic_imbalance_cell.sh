#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --job-name=allelic_imbalance_cell
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

ml samtools
ml htslib
ml picard

cd /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous

mkdir ${cell_accession}/tmp

# add read group information
picard AddOrReplaceReadGroups \
  TMP_DIR=${cell_accession}/tmp \
  I=${cell_accession}/${cell_accession}Aligned.sortedByCoord.out.bam \
  O=${cell_accession}/${cell_accession}_reheader.bam \
  RGLB=${cell_accession} \
  RGPL=illumina \
  RGPU=${cell_accession} \
  RGSM=${cell_accession}

# mark duplicates
picard MarkDuplicates \
  TMP_DIR=${cell_accession}/tmp	\
  I=${cell_accession}/${cell_accession}_reheader.bam \
  O=${cell_accession}/${cell_accession}_dedupped.bam \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT \
  M=${cell_accession}/${cell_accession}.metrics

samtools index ${cell_accession}/${cell_accession}_dedupped.bam

# "split'n'trim" the reads and adjust mapping quality
java -Djava.io.tmpdir=/work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${cell_accession}/tmp/ \
  -jar ~/work/progs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T SplitNCigarReads \
  -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -I ${cell_accession}/${cell_accession}_dedupped.bam \
  -o ${cell_accession}/${cell_accession}_split.bam \
  -rf ReassignOneMappingQuality \
  -RMQF 255 \
  -RMQT 60 \
  -U ALLOW_N_CIGAR_READS

# then run the ASE read counter on individual cell BAMs at those heterozygous sites
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php
~/work/progs/gatk-4.0.12.0/gatk \
  ASEReadCounter \
 -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
 -I ${cell_accession}/${cell_accession}_split.bam \
 -V ${embryo}/${embryo}_knownSNP_filtered.vcf.gz \
 -O ${cell_accession}/${cell_accession}.table
