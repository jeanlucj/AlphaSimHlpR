#' framePhenoRec function
#'
#' function to (do something)
#'
#' @param records [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
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
#' function to (do something)
#'
#' @param phenoDF [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
iidPhenoEval <- function(phenoDF){
  require(lme4)
  # Prepare a vector of error variances since they are heterogeneous
  errVarVec <- numeric(nrow(phenoDF))
  for (n in names(errVars)){
    errVarVec[grep(n, phenoDF$trialType)] <- errVars[n]
  }
  fm <- lmer(pheno ~ (1 | id), weights=1/errVarVec, data=phenoDF)
  return(ranef(fm)[[1]])
}
