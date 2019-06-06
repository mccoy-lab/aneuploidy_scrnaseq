#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --job-name=allelic_imbalance
#SBATCH --time=3:0:0
#SBATCH --array=1-88%88
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

ml samtools
ml htslib
ml picard

cd /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous

embryo=`sed "${SLURM_ARRAY_TASK_ID}q;d" embryo_list.txt`
accession_list=`cat E-MTAB-3929.sdrf.txt | awk -v var="${embryo}" -F $'\t' '{if ($9 == var) print $31}'`
bam_list=`echo ${accession_list} | tr ' ' '\n' | awk '{print $1"/"$1"Aligned.sortedByCoord.out.bam"}' | tr '\n' ' '`

mkdir ${embryo}

# merge single-cell alignments from the same embryo
samtools merge \
  ${embryo}/${embryo}_merged.bam \
  ${bam_list}

samtools sort \
  -o ${embryo}/${embryo}_sorted.bam \
  -@ 8 \
  ${embryo}/${embryo}_merged.bam

# add read group information
picard AddOrReplaceReadGroups \
  I=${embryo}/${embryo}_sorted.bam \
  O=${embryo}/${embryo}_reheader.bam \
  RGLB=${embryo} \
  RGPL=illumina \
  RGPU=${embryo} \
  RGSM=${embryo}

# mark duplicates
picard MarkDuplicates \
  I=${embryo}/${embryo}_reheader.bam \
  O=${embryo}/${embryo}_dedupped.bam \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT \
  M=${embryo}/${embryo}.metrics

samtools index ${embryo}/${embryo}_merged_reheader.bam

# "split'n'trim" the reads and adjust mapping quality
java -jar ~/work/progs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T SplitNCigarReads \
  -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -I ${embryo}/${embryo}_dedupped.bam \
  -o ${embryo}/${embryo}_split.bam \
  -rf ReassignOneMappingQuality \
  -RMQF 255 \
  -RMQT 60 \
  -U ALLOW_N_CIGAR_READS

~/work/progs/gatk-4.0.12.0/gatk \
  HaplotypeCaller \
  -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -I ${embryo}/${embryo}_split.bam \
  --dont-use-soft-clipped-bases true \
  --standard-min-confidence-threshold-for-calling 20.0 \
  -O ${embryo}/${embryo}.vcf

bgzip ${embryo}/${embryo}.vcf
tabix -p vcf ${embryo}/${embryo}.vcf.gz

~/work/progs/vcfanno_linux64 1KG_AF.toml ${embryo}/${embryo}.vcf.gz > ${embryo}/${embryo}_AF.vcf

bgzip ${embryo}/${embryo}_AF.vcf
tabix -p vcf ${embryo}/${embryo}_AF.vcf.gz

# extract known SNPs from 1KG and limit to heterozygous sites
~/work/progs/gatk-4.0.12.0/gatk \
  SelectVariants \
  -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -V ${embryo}/${embryo}_AF.vcf.gz \
  -O ${embryo}/${embryo}_knownSNP.vcf \
  -L /work-zfs/rmccoy22/resources/reference/gatk/wgs_calling_regions.hg38.interval_list \
  -select-type SNP \
  -select "vc.getGenotype(\""${embryo}"\").isHet()" \
  -select 'vc.hasAttribute("1KG_AF")'

bgzip ${embryo}/${embryo}_knownSNP.vcf
tabix -p vcf ${embryo}/${embryo}_knownSNP.vcf.gz

N=8
(
for cell_accession in ${accession_list}
do
  ((i=i%N)); ((i++==0)) && wait
  (
  # add read group information
  picard AddOrReplaceReadGroups \
    I=${cell_accession}/${cell_accession}Aligned.sortedByCoord.out.bam \
    O=${cell_accession}/${cell_accession}_reheader.bam \
    RGLB=${cell_accession} \
    RGPL=illumina \
    RGPU=${cell_accession} \
    RGSM=${cell_accession}
  
  # mark duplicates
  picard MarkDuplicates \
    I=${cell_accession}/${cell_accession}_reheader.bam \
    O=${cell_accession}/${cell_accession}_dedupped.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=${cell_accession}/${cell_accession}.metrics
  
  samtools index ${cell_accession}/${cell_accession}_dedupped.bam
  
  # "split'n'trim" the reads and adjust mapping quality
  java -jar ~/work/progs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
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
    -V ${embryo}/${embryo}_knownSNP.vcf.gz \
    -O ${cell_accession}/${cell_accession}.table
  ) &
done
)

