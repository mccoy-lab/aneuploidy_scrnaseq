#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --job-name=map_scrna
#SBATCH --time=0:30:0
#SBATCH --array=1-1529%100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# map RNA-seq data in preparation for genotyping
# https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

cd /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous

cell=`sed "${SLURM_ARRAY_TASK_ID}q;d" run_list.txt`

/work-zfs/rmccoy22/progs/STAR/bin/Linux_x86_64/STAR \
  --genomeDir /work-zfs/rmccoy22/resources/reference/star_index \
  --readFilesCommand gunzip -c \
  --readFilesIn /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${cell}/${cell}.fastq.gz \
  --outFileNamePrefix /work-zfs/rmccoy22/rmccoy22/mCA/PRJEB11202_Petropolous/${cell}/${cell} \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 8
