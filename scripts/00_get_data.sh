#!/bin/bash

mkdir -p ../RawData
if [ ! -f ../RawData/EMTAB3929.rds ]; then
    wget -O ../RawData/EMTAB3929.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/EMTAB3929.rds
fi

mkdir -p ../external_metadata
if [ ! -f ../external_metadata/stirparo2018_tableS4.xlsx ]; then
    wget -O ../external_metadata/stirparo2018_tableS4.xlsx http://www.biologists.com/DEV_Movies/DEV158501/TableS4.xlsx
fi
if [ ! -f ../external_metadata/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds ]; then
    wget -O ../external_metadata/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds
fi

mkdir -p ../ProcessedData
mkdir -p ../results/mosaic_aneuploidy_plots