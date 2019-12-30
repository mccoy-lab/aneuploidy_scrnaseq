#!/usr/bin/env bash

# login to MARCC DTN2 and use aspera to download...
cd /work-zfs/rmccoy22/rmccoy22/mCA/

# Petropolous et al. 2016
enaGroupGet -g read -f fastq -as ~/progs/aspera/aspera_settings.ini PRJEB11202
