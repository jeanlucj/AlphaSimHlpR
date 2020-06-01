#' plotRecords function
#'
#' function to make a plot from the data in \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return Makes a plot of the gain over cycles of selection and returns the dataset used to make that plot
#' @details \code{records} is a list of lists of populations and phenotype matrices useful for maintaining the phenotypic observations across years and stages
#' 
#' @examples
#' exptSummary <- plotRecords(records)
#' 
#' @export
plotRecords <- function(replicRecords){
  # Need stageNames and nStages from bsp
  bsp <- replicRecords[[1]]$bsp
  allStgOut <- NULL
  for (i in 1:length(replicRecords)){
    allStgOut <- allStgOut %>% bind_rows(replicRecords[[i]]$records$stageOutputs %>% dplyr::mutate(repNum=paste0("R", i)))
  }
  
  nYears <- length(unique(allStgOut$year))
  stageNames <- unique(allStgOut$stage)
  nStages <- length(stageNames)
  nS2f <- floor((nStages)/2)
  nS2c <- ceiling((nStages)/2)
  stageOrd <- c(2* nS2f:1, 1:nS2c *2-1)
  so <- paste0("S", stageOrd); names(so) <- stageNames
  allStgOut <- allStgOut %>% dplyr::mutate(so=so[stage])
  
  bp <- brkdn.plot(genValMean ~ so + year, data=allStgOut, lwd=2, cex=1, col=order(so), md="std.error", stagger=1/(nYears+1)/3/nStages, xlab="Year", ylab="Mean Genotypic Value", main="")
  legend("topleft", legend=stageNames, col=1:(nStages), lwd = 2, cex=0.5, horiz = T)
  
  bp <- brkdn.plot(gvOfBestCrit ~ so + year, data=allStgOut, lwd=2, cex=1, col=order(so), md="std.error", stagger=1/(nYears+1)/3/nStages, xlab="Year", ylab="Clone Evaluated Best", main="")
  legend("topleft", legend=stageNames, col=1:(nStages), lwd = 2, cex=0.5, horiz = T)
  
  bp <- brkdn.plot(genValSD ~ so + year, data=allStgOut, lwd=2, cex=1, col=order(so), md="std.error", stagger=1/(nYears+1)/3/nStages, xlab="Year", ylab="Genotypic Value Std. Dev.", main="")
  legend("topright", legend=stageNames, col=1:(nStages), lwd = 2, cex=0.5, horiz = T)
  
  # Gain across versus within cycles
  soLastBest <- allStgOut %>% dplyr::filter(stage==last(bsp$stageNames) & cycle >= 0)
  soF1mean <- allStgOut %>% dplyr::filter(stage=="F1" & cycle <= max(soLastBest$cycle)+1)
  soF1mean <- soF1mean %>% dplyr::mutate(acrossCycGain=lead(genValMean) - genValMean)
  soF1mean <- soF1mean %>% dplyr::filter(cycle <= max(soLastBest$cycle))
  soF1mean <- soF1mean %>% dplyr::mutate(withinCycGain=(soLastBest$gvOfBestCrit - genValMean) / (bsp$nStages + 1))
  
  gains <- soF1mean %>% tidyr::pivot_longer(cols=ends_with("Gain"), names_to="gainType", values_to="gain")
  
  bp <- brkdn.plot(gain ~ gainType + cycle, data=gains, lwd=2, cex=1, col=1:2, md="std.error", xlab="Cycle", ylab="Gain per Year", main="")
  legend("topright", legend=c("Across Cycles", "Within Cycles"), col=1:2, lwd = 2, cex=0.5, horiz = T)
  
  return(dplyr::select(allStgOut, repNum, year, stage, cycle, genValMean, genValSD, gvOfBestCrit, highestGV, nContribToPar))
}
