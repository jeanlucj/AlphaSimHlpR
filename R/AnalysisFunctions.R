
#' mean_records function
#'
#' function to calculate the mean genotypic value at each cycle and stage
#' WARNING: right now, checks are kept in the population records and they are included in the means that are calculated
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return Real matrix with breeding program cycles in rows and product pipeline stages in columns, each cell being the mean genotypic value for that year and stage
#' @details The records object is a list of lists of populations. This function takes those lists and returns the poplation means in a matrix
#' 
#' @examples
#' recordMeans <- mean_records(records)
#' 
#' @export
mean_records <- function(records){
  chkID <- records$bsp$checks@id
  return(sapply(records$records[1 + 1:records$bsp$nStages], function(popList) sapply(popList, function(popMat) return(mean(popMat$genoVal[!(popMat$id %in% chkID)])))))
}

#' elementWise function
#'
#' given a list of arrays (could be matrices), all of the same dimension, and a function, returns an array of the same dimensions as the objects of the list, for which each element is the application of the function to the vector of cells in a position across all arrays in the list.
#'
#' @param arrayList The list of arrays to which the function will be applied
#' @param fnc The function that can be applied to a vector
#' @return Array of the function values
#' 
#' @examples
#' cellMeans <- elementWise(lapply(replicRecords, mean_records))
#' 
#' @export
elementWise <- function(arrayList, fnc=mean){
  arDim <- dim(arrayList[[1]])
  arMat <- sapply(arrayList, c)
  return(array(apply(arMat, 1, fnc), dim=arDim))
}

#' framePhenoRec function
#'
#' function to make a data.frame to be used as a source of data to analyze the phenotypic \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @return A data.frame of phenotypic records with four columns: 1. The id of individuals; 2. The trial type of the phenotype record; 3. The year the observation was recorded; 4. The phenotypic value
#' @details \code{records} is a list of lists of populations and is primarily useful for maintaining the phenotypic observations across years and stages. For analysis, you need just the phenotypes in a matrix with relevant independent values
#' 
#' @examples
#' phenoDF <- framePhenoRec(records, bsp)
#' 
#' @export
framePhenoRec <- function(records, bsp){
  allPheno <- tibble()
  for (stage in 1 + 1:bsp$nStages){
    for (year in 1:length(records[[stage]])){
      thisPheno <- records[[stage]][[year]] %>% dplyr::mutate(year=year)
      allPheno <- bind_rows(allPheno, thisPheno)
    }
  }
  return(allPheno)
}

#' phenoRecFromPop function
#'
#' function to make a tibble to be added to \code{records}
#'
#' @param pop The population from which to extract phenotypic records. Has to have been phenotyped
#' @param bsp The breeding scheme parameter list
#' @param stage At what stage of the breeding scheme this population was phenotyped. Necessary to determine the error variance and degree of replication
#' @param checks Whether this was a population of experimentals or checks. Necessary to determine the degree of replication
#' @return A tibble of phenotypic records with six columns: 1. The id, 2. The id of the mother, 3. The id of the father, 4. The name of the stage, 5. The actual phenotype, 6. The error variance of the phenotype
#' @details The tibbles coming from this function will be incorporated into \code{records} useful for maintaining the phenotypic observations across years and stages. For analysis, you need these phenotypes coupled to their ids and error variances
#' 
#' @examples
#' phenoDF <- phenoRecFromPop(pop, bsp, stage)
#' 
#' @export
phenoRecFromPop <- function(pop, bsp, stage, checks=FALSE){
  # Entries and replicates have different numbers of stages
  nReps <- if_else(checks, bsp$chkReps[stage], bsp$nReps[stage])
  varE <- bsp$gxyVar + (bsp$gxlVar + bsp$gxyxlVar + bsp$errVars[stage] / nReps) / bsp$nLocs[stage]
  # Set this up for the lmer method distinguishing checks from experimentals
  phenoRec <- tibble(id=pop@id, mother=pop@mother, father=pop@father, stage=bsp$stageNames[stage], isChk=if_else(checks, "check", "exptl"), pheno=pheno(pop), genoVal=gv(pop), errVar=varE)
  return(phenoRec)
}

#' makeGRM function
#'
#' function to make a genomic relationship matrix to be used to analyze the phenotypic \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object. Needed to pull the SNP genotypes
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across years and stages. For GEBV analysis, you need the GRM of these individuals. \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all individuals that have been phenotyped. The GRM also includes the unphenotyped new F1 individuals in records[[1]]
#' 
#' @examples
#' grm <- makeGRM(records, bsp, SP)
#' 
#' @export
makeGRM <- function(records, bsp, SP){
  require(sommer)
  
  allPop <- records$F1
  # Only make GRM of individuals that are specified in the TP
  allID <- NULL
  tpc <- bsp$trainingPopCycles[1]
  if (tpc){
    lastGen <- max(records$F1@fixEff)
    tpc <- min(lastGen, tpc)
    allID <- records$F1@id[records$F1@fixEff %in% lastGen + (1 - tpc):0]
  }
  for (stageNum in 1 + 1:bsp$nStages){
    tpc <- bsp$trainingPopCycles[stageNum]
    if (tpc){
      lastGen <- length(records[[stageNum]])
      tpc <- min(lastGen, tpc)
      for (cyc in lastGen + (1 - tpc):0){
        allID <- c(allID, records[[stageNum]][[cyc]]$id)
      }
    }
  }
  allID <- unique(allID)
  if (!is.null(bsp$checks)) allID <- setdiff(allID, bsp$checks@id)
  allID <- allID[order(as.integer(allID))]
  allPop <- allPop[allID]
  
  if (!is.null(bsp$checks)){
    putInChks <- setdiff(bsp$checks@id, allPop@id)
    if (length(putInChks > 0)) allPop <- c(allPop, bsp$checks[putInChks])
  }
  return(A.mat(pullSnpGeno(allPop, simParam=SP) - 1))
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
#' phenoDF <- framePhenoRec(records, bsp)
#' iidBLUPs <- iidPhenoEval(phenoDF)
#' 
#' @export
iidPhenoEval <- function(phenoDF){
  require(lme4)
  phenoDF$errVar <- 1/phenoDF$errVar # Make into weights
  phenoDF <- phenoDF %>% dplyr::mutate(entryChk=if_else(isChk=="check", id, "-1"))
  fm <- lmer(pheno ~ entryChk + (1|id:isChk), weights=errVar, data=phenoDF)
  blup <- as.matrix(ranef(fm)[[1]])[,1] # Make into matrix to get names
  names(blup) <- (names(blup) %>% strsplit(":", fixed=T) %>% unlist %>%
                  matrix(nrow=2))[1,]
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){
    namesBlup <- names(blup)
    blup <- tapply(phenoDF$pheno, phenoDF$id, mean)
    names(blup) <- namesBlup
  }
  return(blup)
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
#' phenoDF <- framePhenoRec(records, bsp)
#' grm <- makeGRM(records, bsp, SP)
#' grmBLUPs <- grmPhenoEval(phenoDF, grm)
#' 
#' @export
 

grmPhenoEval <- function(phenoDF, grm){
#  if("asreml"%in%installed.packages()) {
    suppressMessages(require(asreml)); suppressMessages(require(Matrix)); suppressMessages(require(synbreed))

    phenoDF <- phenoDF %>% complete(id,nesting(stage, year))
 
    grm <- grm[order(as.numeric(rownames(grm))), order(as.numeric(colnames(grm)))]
    phenoDF <- phenoDF[with(phenoDF,order(as.numeric(id), year)),]
  
    grmPD <- nearPD(grm, keepDiag = TRUE) # Compute the nearest positive definite matrix to an approximate one
    G <- matrix(grmPD[[1]]@x, nrow = grmPD[[1]]@Dim[1])
    G <- G + diag(1e-6, nrow(G)) #
    attr(G, "dimnames") <- grmPD[[1]]@Dimnames
    class(G) <- "relationshipMatrix"
  
    phenoDF$id <- factor(phenoDF$id, levels=rownames(G)) # Enable prediction
    phenoDF$wgt <- 1/phenoDF$errVar # Make into weights
    
    Ginv <- write.relationshipMatrix(G, file = NULL, sorting = "ASReml",
                                     type = c("ginv"), digits = 10) # Invert the G matrix and change to a sparse matrix as required by ASReml package
    names(Ginv) <- c("row", "column", "coefficient")
    attr(Ginv, "rowNames") <- rownames(G)
    attr(Ginv, "colNames") <- colnames(G)
    attr(Ginv, "INVERSE") <- TRUE
    
    fm <- asreml(pheno ~ 1,
                 random = ~ vm(id,Ginv),
                 residual = ~ id(units),
                 weights = wgt,
                 data = phenoDF,
                 workspace = 128e06,
                 na.action = na.method(x = "include",y = "include"))
    
    blup <- summary(fm, coef = T)$coef.random[,"solution"]
    names(blup) <- sapply(strsplit(names(blup), split = "_", fixed = T), function(x) (x[2]))
#    blup <- blup[order(match(names(blup),rownames(grm)))]                  
#} else {
#  require(sommer)
#  print("You decide to use sommer package")
#  
#  phenoDF$id <- factor(phenoDF$id, levels = rownames(grm)) # Enable prediction
#  phenoDF$wgt <- 1 / phenoDF$errVar # Make into weights
#  
#  fm <- mmer(pheno ~ 1,
#             random = ~ vs(id, Gu = grm),
#             method = "EMMA",
#             rcov = ~ units,
#             weights = wgt,
#             data = phenoDF,
#             verbose = F,
#             date.warning = F)
#  blup <- fm$U[[1]][[1]]
#}
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){
    namesBlup <- names(blup)
    blup <- tapply(phenoDF$pheno, phenoDF$id, mean)
    names(blup) <- namesBlup
    }
  return(blup)
}

#' selCritIID function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects are IID
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (not used, here for uniformity)
#' @return An IID BLUP of the trait of the candidates
#' @details Accesses all individuals in \code{records} to pick the highest among candidates. If candidates do not have records, a random sample is returned
#' 
#' @examples
#' allPop <- mergePops(records[[2]])
#' candidates <- allPop@id
#' parents <- allPop[selCritIID(records, candidates, bsp, SP)]
#' 
#' @export
selCritIID <- function(records, candidates, bsp, SP){
  phenoDF <- framePhenoRec(records, bsp)
  # Candidates don't have phenotypes so return random vector
  if (!any(candidates %in% phenoDF$id)){ 
    crit <- runif(length(candidates))
  } else{
    crit <- iidPhenoEval(phenoDF)
    crit <- crit[candidates]
  }
  names(crit) <- candidates
  return(crit)
}

#' selCritGRM function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects covary according to a GRM
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples 
#' candidates <- records[[1]][[1]]@id
#' parents <- records[[1]][[1]][selCritGRM(records, candidates, bsp, SP)]
#' 
#' @export
selCritGRM <- function(records, candidates, bsp, SP){
  grm <- makeGRM(records, bsp, SP)
  if (!any(candidates %in% rownames(grm))){ 
    crit <- runif(length(candidates))
  } else{
    phenoDF <- framePhenoRec(records, bsp)
    # Remove individuals with phenotypes but who do not have geno records
    phenoDF <- phenoDF[phenoDF$id %in% rownames(grm),]
    crit <- grmPhenoEval(phenoDF, grm)
    crit <- crit[candidates]
  }
  names(crit) <- candidates
  return(crit)
}
