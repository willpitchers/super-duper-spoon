# Analysis Workflow

## Find a Powerful Computer

  I initially wrote this code while I worked in the [Dworkin Lab](https://www.msu.edu/~idworkin/) at [Michigan State U.](https://www.msu.edu), so I was lucky enough to be able to run my GWAS analyses on the [High Performance Computing Cluster](https://icer.msu.edu/hpcc).

  If you don't have access to a similarly awesome facility, then you might consider running your analyses on [AWS](https://aws.amazon.com) (AWS have been very generous with free compute time for training purposes e.g. for [ANGUS](http://angus.readthedocs.org/en/2014/), so if you're just learning it'd be worth contacting them and asking politely).

  If you have a server of your own then lucky you! I broke the jobs up into chunks that ran in <4hrs to make life easier for the scheduler-robot at the HPCC, and so that jobs could run in parallel. Your mileage will vary depending on the number of cores that you have to run on; you may find that the analysis runs for some days. The good news is that I ran on nodes with 16GB of memory, so if you can wait you could GWAS on a mid-range desktop.

## Get the Right Software

  This is the easy part! As long as you're woking on a linux/unix system (or even OSX) then you already have almost everything you need. You'll want to go install [**R**](https://cran.r-project.org) (I tested this most recently on version 3.2.0) and the packages [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), [car](https://cran.r-project.org/web/packages/car/index.html). If you want to use the plotting script you'll also need [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

## Obtain Data

  1. **Genotype data** – I was unable to find an md5sum for the genotype file that I'm using here. In lieu of getting a hash from the McKay lab, I downloaded the file `dgrp2.tgeno` 3 times to confirm that I got the same hash each time. If *you* download this file then check that it checks out with: `6f4fce6f575dcf392c93ad3d3ff85ccd`

  2. **Phenotype data** – this would be the hard part! Go forth and perform a well-designed experiment.

## Munge Data #######

  I made the decision to fit my phenotype data to the genotype data, since the `dgrp2.tgeno` file is pretty big (~2GB) and therefore slow to work with. The only munging I'll do to the genotype file is to cut it up into chunks to turn one *looooong* job into many jobs of a few hours:
  > `sh chunk_genotype_file.sh`

  Once you've collected your phenotype data I strongly suggest saving it as a flat text file (I'm going to assume `phenotype.csv`) and change the permissions so you can't edit it by mistake. All the necessary manipulation happens within the `UnivarAssocAnalysis.R` or `MultivarAssocAnalysis.R`, but:
  - I've assumed that your DGRP line ID's are formatted in the same way as in the dgrp2.tgeno file, i.e.; "line_21" – "line_913" and called simply "line", so please ensure that this is the case.
  - If your phenotype is univariate and your phenotype file has more than 2 columns – line & phenotype – please edit `UnivarAssocAnalysis.R` at line #29
