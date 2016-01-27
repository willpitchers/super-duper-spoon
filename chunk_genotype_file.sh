#! /bin/bash

## – these '#PBS' lines are for the benefit of the PBS scheduler at our HPC: you'll need to translate these appropriately for your system – ##
#PBS -o /mnt/scratch/pitchers/GWAS/
#PBS -l nodes=1:ppn=1,walltime=01:00:00,mem=16gb
#PBS -M pitchers@msu.edu
## –– ##

nchunks=100           # feel free to change the number of chunks here

# here I set up some variables
genotype=dgrp2.tgeno
nrows=`wc -l dgrp2.tgeno | cut -d ' ' -f 1`
chunklength=`expr $nrows / $nchunks`

# this line splits the dgrp2.tgeno file
split -d -a 4 -l ${chunklength} dgrp2.tgeno DGRP_chunk_

# and this line moves the chunks into their own folder for tidiness
mkdir DGRPchunks && mv DGRP_chunk_* DGRPchunks/
