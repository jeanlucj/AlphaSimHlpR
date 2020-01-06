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
framePhenoRec <- function(records, bsp){
  allPheno <- data.frame()
  for (stage in 2:length(records)){
    errVar <- bsp$errVars[names(records)[stage]]; names(errVar) <- NULL
    for (year in 1:length(records[[stage]])){
      pop <- records[[stage]][[year]]
      thisPheno <- data.frame(id=pop@id, stage=names(records)[stage], year=year, errVar=errVar, pheno=c(pop@pheno))
      allPheno <- rbind(allPheno, thisPheno)
    }
  }
  return(allPheno)
}

#' makeGRM function
#'
#' function to make a genomic relationship matrix to be used to analyze the phenotypic \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param SP The AlphaSimR SimParam object
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across years and stages. For GEBV analysis, you need the GRM of these individuals. \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all individuals that have been phenotyped. The GRM also includes the unphenotyped new F1 individuals in records[[1]]
#' 
#' @examples
#' grm <- makeGRM(records)
#' 
#' @export
makeGRM <- function(records, SP){
  require(sommer)
  allPop <- mergePops(c(records[[1]], records[[2]]))
  allPop <- allPop[!duplicated(allPop@id)]
  grm <- A.mat(pullSnpGeno(allPop, simParam=SP) - 1)
  return(grm)
}

#' iidPhenoEval function
#'
#' function to take a data.frame coming from framePhenoRec and analyze it with individuals as a random effect with an IID covariance matrix
#'
#' @param phenoDF A data.frame of phenotypic observations. See \code{framePhenoRec} for details.
#' @return Named real vector of the BLUPs of all individuals in phenoDF (names are the individual ids), with appropriate weights by error variance of the observation
#' @details Given all the phenotypic records calculate the best prediction of the genotypic value for each individual using all its records
#' 
#' @examples
#' phenoDF <- framePhenoRec(records)
#' iidBLUPs <- iidPhenoEval(phenoDF)
#' 
#' @export
iidPhenoEval <- function(phenoDF){
  require(lme4)
  phenoDF$errVar <- 1/phenoDF$errVar # Make into weights
  fm <- lmer(pheno ~ (1 | id), weights=errVar, data=phenoDF)
  return(as.matrix(ranef(fm)[[1]])[,1]) # Make into matrix to get names
}

#' grmPhenoEval function
#'
#' function to take a data.frame coming from framePhenoRec and GRM and analyze them with individuals as a random effect with a GRM covariance matrix
#'
#' @param phenoDF A data.frame of phenotypic observations. See \code{framePhenoRec} for details
#' @param grm A genomic relationship matrix
#' @return Named real vector of the BLUPs of all individuals in phenoDF (names are the individual ids), with appropriate weights by error variance of the observation
#' @details Given all the phenotypic records calculate the GEBV for each individual using all its records
#' 
#' @examples
#' phenoDF <- framePhenoRec(records)
#' grm <- makeGRM(records)
#' grmBLUPs <- grmPhenoEval(phenoDF, grm)
#' 
#' @export
grmPhenoEval <- function(phenoDF, grm){
  require(sommer)
  levels(phenoDF$id) <- levels(as.factor(rownames(grm))) # Enable prediction
  phenoDF$errVar <- 1/phenoDF$errVar # Make into weights
  fm <- mmer(pheno ~ 1,
             random= ~ vs(id, Gu=grm),
             method="EMMA",
             rcov= ~ units,
             weights=errVar,
             data=phenoDF)
  return(fm$U[[1]][[1]])
}
