#!/usr/bin/env bash

# login to MARCC DTN2 and use aspera to download...
cd /work-zfs/rmccoy22/rmccoy22/mCA/

# Yan et al. 2013
enaGroupGet -g read -f fastq -as ~/progs/aspera/aspera_settings.ini PRJNA153427

# Xue et al. 2013
enaGroupGet -g read -f fastq -as ~/progs/aspera/aspera_settings.ini PRJNA189204

# Blakely et al. 2015 
enaGroupGet -g read -f fastq -as ~/progs/aspera/aspera_settings.ini PRJNA277181

# Petropolous et al. 2016
enaGroupGet -g read -f fastq -as ~/progs/aspera/aspera_settings.ini PRJEB11202
