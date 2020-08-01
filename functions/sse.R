## ---------------------------
##
## Script name: sse.R
##
## Purpose of script: clustering assessment
##
## Author: Georgia Doing
##
## Date Created: 2020-07-31
##
## 
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes: This script contains function(s) to assess clustering solutions 
## using the sum of squared error calculation (SSE). Specifically the ration of
## within groups SSE to total SSE.
##   
##
## ---------------------------



#' Function description
#'
#' @param groupA A dataframe with samples as columns and gene expression as rows
#' @param groupB A dataframe with samples as columns and gene expression as rows
#' @return A list of numeric values: list('error_ratio','error','a_error','b_error')
#' @examples
#' sse_ratio(data.frame('A' = c(1:10)), data.frame('B'=c(2:11)))
#' @export

sse_ratio <- function(groupA, groupB){
  
  a_mean <- rowMeans(groupA)
  b_mean <- rowMeans(groupB)
  mean <- (a_mean + b_mean) / 2
  
  a_error <- sum(apply(groupA, 2, function(x) (a_mean - x)**2))
  b_error <- sum(apply(groupB, 2, function(x) (b_mean - x)**2))
  
  error <- sum(apply(groupA, 2, function(x) (mean - x)**2)) + 
    sum(apply(groupB, 2, function(x) (mean - x)**2))
  
  error_ratio <- (a_error + b_error) / error
  
  return(list(error_ratio = error_ratio, error = error, a_error = a_error, b_error = b_error))
}
