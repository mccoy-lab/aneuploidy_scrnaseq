# aneuploidy_project


## Author: Margaret R. Starostik


## Data Acquisition

The data set EMTAB3929 (PMID 27062923) was downloaded on Tuesday November 6 13:29 EST 2018 using:

http://imlspenticton.uzh.ch:3838/conquer/

Information downloaded from this website includes the (1) MultiAssay Experiment (EMTAB3929.rds), (2) MultiQC report, (3) scater report, and (4) salmon archive (EMTAB3929_salmo.tar and EMTAB3929 folder).

Supplementary files from Griffiths et al., 2017 were forked on Tuesday, November 6 15:05 EST 2018 from MarioniLab/Aneuploidy2017 on Github: https://github.com/MarioniLab/Aneuploidy2017.

Data used in Griffiths et al., 2017 were obtained using their shell script: sh get_data.sh on Wednesday November 14, 12:33 EST 2018.


## scploid R Package

The package was installed from R using the devtools package and the following commands:

library(devtools)
devtools::install_github("MarioniLab/Aneuploidy2017", subdir = "package")
