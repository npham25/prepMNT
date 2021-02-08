#' @title Preprocess MNTs data
#'
#' @description     This package preprocesses the MNTs data by removing MNTs which have few available values and imputing the zero values
#'
#' @param x,pctzero,freqCut,uniqueCut
#'
#' @return
#'
#' @examples prepMNT(data_mnt_1, pctzero=75, freqCut=19, uniqueCut=20)
#'
#' @export prepMNT
inst_pkg <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

inst_pkg(c("dplyr","imputeLCMD","caret"))
require(dplyr)
require(imputeLCMD)
require(caret)
    
prepMNT <- function(x, pctZero=75, freqCut = 19, uniqueCut = 20) {
  res <- caret::nearZeroVar(x, freqCut = freqCut, uniqueCut = uniqueCut,saveMetrics = TRUE, names = FALSE,
                     foreach = FALSE, allowParallel = TRUE)
  #Count number of zeroes and NAs
  res$pctZero <- as.vector(apply(x, 2, function(x)  100*sum(x == 0, na.rm=T)/length(na.omit(x))))
  #Exclude column added
  res$exclude <- res$pctZero == 100 |
    (res$freqRatio > freqCut & res$percentUnique < uniqueCut & res$pctZero > pctZero)

  #Remove lipids in x
  keep <- row.names(res)[res$exclude == FALSE]
  x <- x[, colnames(x) %in% keep]
  #Replace 0 by ""
  x[x == 0] <- ""
  x <- as.matrix(sapply(x, function(x) as.numeric(as.character(x))))
  return (QRILC_wrapper(x))
}

QRILC_wrapper <- function(data, ...) {
  result <- data %>% log %>% impute.QRILC(., ...) %>% extract2(1) %>% exp
  return(result)
}


