#' mean_records function
#'
#' function to calculate the mean genotypic value at each cycle and stage
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return Real matrix with breeding program cycles in rows and product pipeline stages in columns, each cell being the mean genotypic value for that year and stage
#' @details The records object is a list of lists of populations. This function takes those lists and returns the poplation means in a matrix
#' 
#' @examples
#' recordMeans <- mean_records(records)
#' 
#' @export
mean_records <- function(records){
  return(sapply(records[-1], function(popList) sapply(popList, meanG)))
}
#' framePhenoRec function
#'
#' function to make a data.frame to be used as a source of data to analyze the phenotypic \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return A data.frame of phenotypic records with four columns: 1. The id of individuals; 2. The trial type of the phenotype record; 3. The year the observation was recorded; 4. The phenotypic value
#' @details \code{records} is a list of lists of populations and is primarily useful for maintaining the phenotypic observations across years and stages. For analysis, you need just the phenotypes in a matrix with relevant independent values
#' 
#' @examples
#' phenoDF <- framePhenoRec(records)
#' 
#' @export
framePhenoRec <- function(records){
  allPheno <- data.frame()
  for (trialType in 2:length(records)){
    for (year in 1:length(records[[trialType]])){
      pop <- records[[trialType]][[year]]
      thisPheno <- data.frame(id=pop@id, trialType=names(records)[trialType], year=year, pheno=c(pop@pheno))
      allPheno <- rbind(allPheno, thisPheno)
    }
  }
  return(allPheno)
}

#' iidPhenoEval function
#'
#' function to take a data.frame coming from framePhenoRec and analyze it with individuals as a random effect with an IID covariance matrix
#'
#' @param phenoDF A data.frame of phenotypic observations. See \code{framePhenoRec} for details.
#' @param ppp A list of the product pipeline parameters. See \code{runBreedingScheme} for details. Here we need the named real vector with the error variances
#' @return Named real vector of the BLUPs of all individuals in phenoDF (names are the individual ids), with appropriate weights by error variance of the observation
#' @details Given all the phenotypic records calculate the best prediction of the genotypic value for each individual using all its records
#' 
#' @examples
#' phenoDF <- framePhenoRec(records)
#' iidBLUPs <- iidPhenoEval(phenoDF)
#' 
#' @export
iidPhenoEval <- function(phenoDF, ppp){
  require(lme4)
  # Prepare a vector of error variances since they are heterogeneous
  errVarVec <- numeric(nrow(phenoDF))
  for (n in names(ppp$errVars)){
    errVarVec[grep(n, phenoDF$trialType)] <- ppp$errVars[n]
  }
  fm <- lmer(pheno ~ (1 | id), weights=1/errVarVec, data=phenoDF)
  return(ranef(fm)[[1]])
}
