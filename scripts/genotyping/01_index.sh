#!/usr/bin/env bash

# generate STAR index

cd /work-zfs/rmccoy22/resources/reference

/work-zfs/rmccoy22/progs/STAR/bin/Linux_x86_64/STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir /work-zfs/rmccoy22/resources/reference/star_index \
  --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --sjdbGTFfile gencode.v30.annotation.gtf
