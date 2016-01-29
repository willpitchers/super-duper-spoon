#! /bin/bash

## – these '#PBS' lines are for the benefit of the PBS scheduler at our HPC: you'll need to translate these appropriately for your system – ##
#PBS -o /mnt/scratch/pitchers/GWAS/
#PBS -l nodes=1:ppn=1,walltime=01:00:00,mem=16gb
#PBS -M pitchers@msu.edu
## –– ##

###### feel free to change the number of chunks here
nchunks=100
###### NB: this number need not divide exactly into the number of variants, `split` function will not throw away data.

# here I set up some variables
genotype=dgrp2.tgeno
nrows=`wc -l dgrp2.tgeno | cut -d ' ' -f 1`
chunklength=`expr $nrows / $nchunks`

# put the header line aside for later
head -1 ${genotype} > headerfile

# this line splits the dgrp2.tgeno file
split -d -a 4 -l ${chunklength} dgrp2.tgeno DGRP_chunk_

# now to fix the headers
for i in seq -w 0001 ${nchunks}
    do echo -e "$(cat headerfile) \n$(cat DGRP_chunk_${i})" > DGRP_chunk_${i}
done

# and this line moves the chunks into their own folder for tidiness
mkdir DGRPchunks && mv DGRP_chunk_* DGRPchunks/
