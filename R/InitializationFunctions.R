#' initializeSimParam function
#'
#' function to (do something)
#'
#' @param meanDD [value]. Default is meanDD=0.2
#' @param varDD [value]. Default is varDD=0.1
#' @param nSNPperChr [value]. Default is nSNPperChr=100
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
initializeSimParam <- function(meanDD=0.2, varDD=0.1, nSNPperChr=100){
  # Create haplotypes for founder population of 200 outbred individuals
  # Default effective population size for runMacs is 100, so a bit low
  founderPop <- runMacs(nInd=200, nChr=18, segSites=1000, inbred=FALSE)
  
  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderPop)
  # Additive and dominance trait architecture
  SP$addTraitAD(nQtlPerChr=200, mean=0, var=genVarCET, meanDD=meanDD, varDD=varDD, useVarA=FALSE)
  # 100 observed SNPs per chromosome
  SP$addSnpChip(nSNPperChr)
  # Default value for error variance in case none is given
  SP$setVarE(varE=errVarSDN)
  return(list(founderPop, SP))
}
