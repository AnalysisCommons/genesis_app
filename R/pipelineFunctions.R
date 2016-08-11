library(plyr)
library(gap)



# This function provides a basic check that the phenotype file is a suitable
# for seqMeta.  It checks that the required column names are in the phenotype
# data frame. 
#
checkPhenotype <- function(p, outcome, covariates, id.col=NULL, gender.col=NULL) {
  if (!is.null(id.col)) {
    if (anyDuplicated(p[ , id.col])) {
      stop("Duplicated phenotype ids.")
    }
  }
  
  if(!is.null(gender.col)) {
    if (!(gender.col %in% colnames(p))) {
      msg <- paste(gender.col, "not found in phenotype file.", sep=" ")
      stop(msg)
    }
    gtype <- typeof(p[[gender.col]])
    g <- unique(p[[gender.col]])
    if(gtype == "integer") {
      if (!all(g %in% c(0, 1))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }      
    } else if(gtype == "character") {
      if (!all(g %in% c("F", "M"))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }        
    }     
  } else {
    wmsg <- "No column given to identify Males in phenotype file.  No special 
         handling of the X chromosme."
    warning(wmsg)
  }
  missing.covariates <- !(covariates %in% colnames(p))
  if (any(missing.covariates)) {
    msg <- paste("Covariates:", covariates[missing.covariates], "not found in phenotype file.\n", sep=" ")
    print(colnames(p))
    print(covariates %in% colnames(p))
    print(covariates[covariates %in% colnames(p)])
    stop(msg)
  } 
  return(invisible(NULL)) 
}



# This function reduces a data set to only the variables used in a model
# removing subjects with missing data.  Also, it makes the row names of
# the resulting data fram the subject identifier
#
# JAB addition: subsets to complete cases (i.e. no NAs in outcome or covariates)
#
# p: a data frame containing the variables in the model
#
# formula: a character vector which can be coered to an object of class 
#          "formula" (with as.formula): a symbolic description of the model to
#          be fitted. The details of model specification are given under 
#          'Details' in the "lm" help file.
#
# id: (optional) colunm name identifier of the subjects
#
# gender: (optional) colunm name identifier for the gender classification of
#         the subjects.
#
# returns: data frame with only the columns specified in the formula and with
#          the (optional) row names as the subject identifier.
#

reducePheno <- function(pheno.data, 
                        outcome, 
                        covariates = NULL, 
                        id=NULL, 
                        gender=NULL) {
  checkPhenotype(pheno.data, outcome, covariates, id.col=id, gender.col=gender)   
  
  if (!is.null(id)) {
    rownames(pheno.data) <- pheno.data[ ,id]
  }
  if (!is.null(gender)) {
      gtype <- typeof(pheno.data[[gender]])
      if(gtype == "integer") {
            males <-  pheno.data[[ gender ]] == 0 
            pheno.data[[ gender]] <- "F"
            pheno.data[[gender]][males] <- "M"
      }
  }

  all.terms <- unique(c(outcome, covariates, gender))
  pheno.data <- as.data.frame(pheno.data) 
  pheno <- na.omit(pheno.data[, all.terms, drop=F])
  return(pheno)
}

# Calculate MAF
#
# dose: matrix with dosages rows individuals, columns are variants


Maf <- function(dose){
       aaf <- colMeans(dose,na.rm=T)/2
       return(min(1-aaf, aaf))
}


split.by.comma <- function(cur.string){
    cur.string <- gsub('"', '', cur.string)
    out <- unlist(strsplit(cur.string, ","))
    if (length(out) == 0){
        out = NULL
    }
    return(out)
}
    


