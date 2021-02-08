require(imputeLCMD)
require(magrittr)

prepMNT <- function(x, pctZeroCut=75, freqCut = 19, uniqueCut = 20) 
  {
  if (is.null(dim(x))) 
    x <- matrix(x, ncol = 1)
  
  x <- as.matrix(apply(x,2,function(x) as.numeric(as.character(x))))
  #the ratio of the most common value to the second most common value
  freqRatio <- apply(x, 2, function(data) {
    t <- table(data[!is.na(data)])
    if (length(t) <= 1) {
      return(0)
    }
    w <- which.max(t)
    return(max(t, na.rm = TRUE)/max(t[-w], na.rm = TRUE))
  })
  #	the percentage of distinct values out of the number of total samples
  lunique <- apply(x, 2, function(data) length(unique(data[!is.na(data)])))
  percentUnique <- 100 * lunique/apply(x, 2, length)
  #Percentage of zeros
  pctZero <- as.vector(apply(x, 2, function(x)  100*sum(x == 0, na.rm=T)/length(na.omit(x))))
  #Exclude MNTs which have 100% zeros or > 75%
  rmv <- pctZero == 100 | 
    (freqRatio > freqCut & percentUnique < uniqueCut & pctZero > pctZeroCut)
  #Remove lipids in x
  x <- x[, !rmv]
  #Replace 0 by ""
  x[x == 0] <- ""
  x <- as.matrix(apply(x,2,function(x) as.numeric(as.character(x))))
  return (QRILC_wrapper(x))
  }

QRILC_wrapper <- function(data, ...) {
  result <- data %>% log %>% imputeLCMD::impute.QRILC(., ...) %>% magrittr::extract2(1) %>% exp
  return(result)
}

`%>%` <- function(lhs, rhs) {
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  kind <- 1L
  env <- parent.frame()
  lazy <- TRUE
  .External2(magrittr_pipe)
}
