
#' ### Association Analysis by Chromosome Chunk ###
#' ### Univariate response phenotype ###


#' In this first code block the script loads the packages required and prints the R version infomation
# ----
R.version
require( WGCNA )
require( dplyr )
require( MASS )


#' This second block reads in the genotype 'chunk' file as specified by the value that replaces "CHUNK_FILE"
#' and sets the filename for the output
# ----
filename <- "DGRP_chunk_0100"#"CHUNK_FILE"
outputfilename <- paste( filename, "_association_output.csv", sep="" )
SNPs <- read.table( filename, header=TRUE, na.strings="-" )


#' Here I'm reading in the phenotype file -- obviously you'll need to edit this part
#' For demonstration I'm just using a set of fake line means generated with th `rnorm` function,
#' but I've tested this script with individual-level datasets too.
# ----
phenotype <- read.csv( "phenotype_file.csv" )
  if ( ncol( phenotype ) >2 ) {
    focal_pheno <- dplyr::select( phenotype, line, take_off_speed )   ##### here "wing_span" is my focal trait: BE SURE TO ENTER YOURS! #####
  }


#' Here I'm extracting just the variant calls from the genotype file and transposing them before
#' combining them with the phenotype data
# ----
snp_start_col <- grep( "line", names( SNPs ) )[1]
SNP_call <- transposeBigData( SNPs[ , snp_start_col:ncol( SNPs ) ] )

  for (i in 1:ncol( SNP_call )) {
      SNP_call[,i] <- factor( SNP_call[,i] )  # this loop ensures that the variant calls are treated as factors for modelling
  }

names( SNP_call ) <- SNPs$id
SNP_call$line <- factor( row.names( SNP_call ))

PandG <- full_join( focal_pheno, SNP_call, by="line" )


#' Defining the Association function
# ----
UV_associations <- function( PandG, MAC=4, Family="gaussian", Covar=FALSE ) {
#   i=1
#   Covar=FALSE

  ## set up objects/variables
  pheno <- as.matrix( PandG[,2] )
  geno <- as.matrix( PandG[, 3:ncol( PandG ) ] )
  n_sites  <- dim( geno )[[2]]

  output_data_frame <- data.frame( matrix( NA, nrow=n_sites, ncol=4 ) )
  names( output_data_frame ) <- c( "estimate", "std_error", "test_stat", "p_value" )

  if( Covar==FALSE ) { my_form <- as.formula( pheno ~  geno[,i] )
                      } else { my_form <- as.formula( paste( "pheno ~  geno[,i] +", Covar ))
                               n_covars <- nrow( attr( terms( my_form ), "factors")) -2
                               covar_out_frame <- data.frame( matrix( NA, nrow=n_sites, ncol= (4*n_covars) )) }

    for ( i in 1:n_sites ) {

    alleles  <- geno[,i]
    minor_allele   <- table( alleles )[[2]]

    if ( minor_allele < MAC ) { output_data_frame[ i, ] <- NA  # NA output if there are too few lines with minor allele
    } else {
          lin_mod <- glm( my_form, family = Family )
          output_data_frame[i ,] <- summary( lin_mod )$coef[ 2, 1:4 ]
        }

    if ( Covar != FALSE & minor_allele > MAC ) {
          covar_out_frame[i ,] <- c( summary( lin_mod )$coef[ 3:( 2+n_covars ), ] )
          covar_names <- sub( "\\$", "_", row.names(summary( lin_mod )$coef)[ 3:( 2+n_covars ) ])
          names( covar_out_frame ) <- paste( rep( covar_names, each=4), "_", c( "estimate", "std_error", "test_stat", "p_value" ), sep="" )
    }
  }
    if ( Covar != FALSE ) { output_data_frame <- cbind( output_data_frame, covar_out_frame ) }
      return( output_data_frame )
}


#' Running the Assocition function
#' you can adjust the minimum acceptable minor allele count (MAC) in the function call,
# ----
output <- UV_associations( PandG )
output <- UV_associations( PandG, MAC=4, Covar="phenotype$take_off_speed" )



# write.csv( output, outputfilename )
