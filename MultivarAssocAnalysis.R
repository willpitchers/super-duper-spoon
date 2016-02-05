
#' ### Association Analysis by Chromosome Chunk ###
#' ### Multivariate response phenotype ###

#' This script execpts to be fed a chunk_ID. It will run a General Linear Model for each variant/SNP,
#' with a single specified response phenotype and an arbitrary number of covariates. The user may also
#' specifiy the desired error distribution; allowing for e.g. response data as counts/frequencies or
#' other data with an non-gaussian error distribution.


#' In this first code block the script loads the packages required and prints the R version infomation
# ---- echo=FALSE
print( paste( "analysis run on", date() ) )
# ----
R.version
require( WGCNA )
require( MASS )
require( dplyr )
require( car )


#' This second block reads in the genotype 'chunk' file as specified by the value that replaces "CHUNK_FILE"
#' and sets the filename for the output
# ----
filename <- "CHUNK_FILE"
outputfilename <- paste( filename, "_association_output.csv", sep="" )
SNPs <- read.table( filename, header=TRUE, na.strings="-" )


#' Here I'm reading in the phenotype file -- obviously you'll need to edit this part.
#' For demonstration I'm just using a set of fake line means generated with th `rnorm` function,
#' but I've tested this script with individual-level datasets too.
# ----
phenotype <- read.csv( "phenotype_file.csv" )


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


#' Joining the Phenotypic and Genotypic data - this makes it easier to feed directly into `glm`
# ----
PandG <- full_join( phenotype, SNP_call, by="line" )


#' Defining the Association function
# ----
MV_associations <- function( PandG, LHS, MAC=4, Family="gaussian", Covar=FALSE ) {

  ## set up objects/variables
  var_start_col <- grep( "[X23LR]{1,2}_\\d+_[A-Z]+", names( PandG ) )[1]
  pheno <- as.matrix( LHS )
  geno <- as.matrix( PandG[, var_start_col:ncol( PandG ) ] )
  n_sites  <- dim( geno )[[2]]

  output_data_frame <- data.frame( matrix( NA, nrow=n_sites, ncol=6 ) )
  names( output_data_frame ) <- c( "wilks", "approx_F", "num_DF", "den_DF", "p_value", "resid" )

  if( Covar==FALSE ) { my_form <- as.formula( pheno ~  geno[,i] )
                      } else { my_form <- as.formula( paste( "pheno ~  geno[,i] +", Covar ))
                               n_covars <- nrow( attr( terms( my_form ), "factors")) -2
                               covar_out_frame <- data.frame( matrix( NA, nrow=n_sites, ncol= (4*n_covars) )) }

    for ( i in 1:n_sites ) {

    alleles  <- geno[,i]
    minor_allele   <- table( alleles )[[2]]

    if ( minor_allele < MAC ) { output_data_frame[ i, ] <- NA  # NA output if there are too few lines with minor allele
    } else {
          lin_mod <- manova( my_form )
          output_data_frame[i ,1:5] <- summary( lin_mod, test="Wilks" )$stats[ 1, 2:6 ]
          output_data_frame[i ,6] <- summary( lin_mod, test="Wilks" )$stats[ 2,1 ]
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
output <- MV_associations( PandG, LHS=select( PandG, wing_span, take_off_speed )  )    # defaults to MAC=4, Family="gaussian"

# here is a more stringent MAC, with a covariate and fitting a different error distribution
# output <- UV_associations( PandG, MAC=7, Covar="phenotype$take_off_speed", Family="quasipoisson" )

write.csv( output, outputfilename )
