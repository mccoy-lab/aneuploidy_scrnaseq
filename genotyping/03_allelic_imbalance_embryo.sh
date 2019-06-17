#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --job-name=allelic_imbalance_embryo
#SBATCH --time=18:0:0
#SBATCH --array=1-88%88
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

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

mkdir ${embryo}/tmp

# add read group information
picard AddOrReplaceReadGroups \
  TMP_DIR=${embryo}/tmp \
  I=${embryo}/${embryo}_sorted.bam \
  O=${embryo}/${embryo}_reheader.bam \
  RGLB=${embryo} \
  RGPL=illumina \
  RGPU=${embryo} \
  RGSM=${embryo}

# mark duplicates
picard MarkDuplicates \
  TMP_DIR=${embryo}/tmp \
  I=${embryo}/${embryo}_reheader.bam \
  O=${embryo}/${embryo}_dedupped.bam \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT \
  M=${embryo}/${embryo}.metrics

samtools index ${embryo}/${embryo}_dedupped.bam

# "split'n'trim" the reads and adjust mapping quality
java -Djava.io.tmpdir=/work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${embryo}/tmp/ \
  -jar ~/work/progs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
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

~/work/progs/gatk-4.0.12.0/gatk \
  VariantFiltration \
  -R /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -V ${embryo}/${embryo}_knownSNP.vcf.gz \
  -O ${embryo}/${embryo}_knownSNP_filtered.vcf \
  -window 35 \
  -cluster 3 \
  --filter-name FS -filter "FS > 30.0" \
  --filter-name QD -filter "QD < 2.0"

bgzip ${embryo}/${embryo}_knownSNP_filtered.vcf
tabix -p vcf ${embryo}/${embryo}_knownSNP_filtered.vcf.gz

for cell_accession in ${accession_list}
do
  sbatch --export=cell_accession=${cell_accession},embryo=${embryo} 04_allelic_imbalance_cell.sh
done
