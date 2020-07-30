# GOING FOR LOW on both criteria

#' idx1 function
#'
#' function to return the number of points lower than the query point
#' in the first column of a matrix
#'
#' @param pt Real query point
#' @param nonDom Real matrix 
#' @return Integer number of values in nonDom <= to pt
#' 
#' @details Utility for Pareto front
idx1 <- function(pt, nonDom){
  sum(nonDom[,1] <= pt[1])
}

#' insertRow function
#'
#' Insert a vector (vec) into a matrix (mat) at row idx
#'
#' @param vec Real row vector
#' @param idx Integer insertion point 
#' @param mat Real matrix 
#' @return Augmented matrix
#' 
#' @details Utility for Pareto front
insertRow <- function(vec, idx, mat){
  top <- NULL; if (idx) top <- mat[1:idx,,drop=F]
  return(rbind(top, vec, mat[-(1:idx),,drop=F]))
}

#' testDom function
#'
#' Test if any of the points in nonDom are dominated by pt
#'
#' @param pt Real double
#' @param nonDom Real two-column matrix 
#' @return Logical vector size number of matrix rows
#' 
#' @details Utility for Pareto front
testDom <- function(pt, nonDom){
  idx <- sum(nonDom[,1] < pt[1])
  c(rep(F, idx), nonDom[-(1:idx), 2] > pt[2])
}

#' checkAddRemove function
#'
#' Test if any of the points in nonDom are dominated by pt
#'
#' @param cand Real double that you want to check if it is dominated or not
#' @param nonDom Real two-column matrix, with vectors in rows that are 
#' non-dominated and sorted from low to high for their first column
#' @return A new nonDom matrix, adding cand if it is non-dominated and
#' removing any rows that are dominated by it
#' 
#' @details Utility for Pareto front
checkAddRemove <- function(cand, nonDom){
  idx <- idx1(cand, nonDom)
  # if idx == 0 cand is lower than all nonDom in the first column
  # so it cannot be dominated
  # else, there are some points lower than it, and if any are also lower
  # in the second coordinate, the point is dominated
  dominated <- idx > 0 & any(nonDom[1:idx,2] < cand[2])
  if (!dominated){
    nonDom <- insertRow(cand, idx, nonDom)
    newDom <- testDom(cand, nonDom)
    if (any(newDom)){
      nonDom <- nonDom[!newDom,]
    }
  }
  return(nonDom)
}

#' returnNonDom function
#'
#' Find the non-dominated points in a matrix
#'
#' @param mat Real two-column matrix
#' @param dir1Low Logical: desired direction in down for column 1
#' @param dir2Low Logical: desired direction in down for column 2
#' @return nonDom Real two-column matrix, with vectors in rows that are 
#' non-dominated and sorted from low to high for their first column
#' 
#' @details Utility for Pareto front
#' @export
returnNonDom <- function(mat, dir1Low=T, dir2Low=T){
  # function tries to go low for both columns so reverse columns if F
  if (!dir1Low) mat[,1] <- -mat[,1]
  if (!dir2Low) mat[,2] <- -mat[,2]
  # find the minima of columns 1 and 2
  min1 <- which(mat[,1] == min(mat[,1]))
  if (length(min1) > 1) min1 <- min1[which.min(mat[min1, 2])]
  nonDom <- mat[min1,,drop=F]
  mat <- mat[-min1,,drop=F]
  min2 <- which(mat[,2] == min(mat[,2]))
  if (length(min2) > 1) min2 <- min2[which.min(mat[min2, 1])]
  nonDom <- rbind(nonDom, mat[min2,,drop=F])
  mat <- mat[-min2,,drop=F]
  for (i in 1:nrow(mat)){
    nonDom <- checkAddRemove(mat[i,,drop=F], nonDom)
  }
  if (!dir1Low) nonDom[,1] <- -nonDom[,1]
  if (!dir2Low) nonDom[,2] <- -nonDom[,2]
  return(nonDom)
}
