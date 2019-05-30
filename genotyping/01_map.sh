#!/usr/bin/env bash

# map RNA-seq data in preparation for genotyping

# https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

cd /work-zfs/rmccoy22/resources/reference

/work-zfs/rmccoy22/progs/STAR/bin/Linux_x86_64/STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir /work-zfs/rmccoy22/resources/reference/star_index \
  --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --sjdbGTFfile gencode.v30.annotation.gtf

cd /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous

N=24
(
for cell in ERR*
do
   ((i=i%N)); ((i++==0)) && wait
   /work-zfs/rmccoy22/progs/STAR/bin/Linux_x86_64/STAR \
     --genomeDir /work-zfs/rmccoy22/resources/reference/star_index \
     --readFilesCommand gunzip -c \
     --readFilesIn /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${cell}/${cell}.fastq.gz \
     --outFileNamePrefix /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${cell}/${cell} \
     --outSAMtype BAM SortedByCoordinate & 
done
)
