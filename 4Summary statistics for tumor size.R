GM_mean <- function(InertCellNum){
  # GM_mean: calculate the geometric mean from an array of cell number
  GM_x = log(InertCellNum)
  X = mean(GM_x)
  return(exp(X))
}

getGMmeanRaw <- function(DataAll){
  # getGMmeanRaw: calculate the absolute cell number GMmean for each sgID without normalization
  GMmeanAll <- c()
  for (sgID in sgIDAll){
    GMmeanAll <- c(GMmeanAll, GM_mean(DataAll$CellNum[DataAll$sgID == sgID]))
  }
  return(GMmeanAll)
}

getGMmean <- function(DataAll){
  # getGMmean:  Calculate LN mean of each sgID relative to Inert from Data Frame
  # Call: GM_mean
  InertGMmean <- GM_mean(DataAll$CellNum[DataAll$sgID %in% Inert])
  
  # calculate the fraction of sgID changes relative to the total in all four samples
  GMmeanAll <- getGMmeanRaw(DataAll)
  GMmeanAll <- GMmeanAll/InertGMmean
  return(GMmeanAll)
}

LN_mean <- function(InertCellNum){
  # LN_mean: calculate the absolute LNmean from an array of cell number
  # Call: Null
  # MLE of mean of data, presuming a LogNormal Distribution.
  # ### old version of LNmean
  LN_x = log(InertCellNum)
  X = mean(LN_x)
  X2 = var(LN_x)
  return (exp(X + 0.5*X2))
}

getLNmeanRaw <- function(DataAll){
  # getLNmeanRaw: calculate the absolute cell number LNmean for each sgID without normalization
  LNMeanAll <- c()
  for (sgID in sgIDAll){
    LNMeanAll <- c(LNMeanAll, LN_mean(DataAll$CellNum[DataAll$sgID == sgID]))
  }
  return(LNMeanAll)
}

getLNMean <- function(DataAll){
  # getLNMean:  Calculate LN mean of each sgID relative to Inert from Data Frame
  # Call: LN_mean
  InertLNMean <- LN_mean(DataAll$CellNum[DataAll$sgID %in% Inert])
  
  # calculate the fraction of sgID changes relative to the total in all four samples
  LNMeanAll <- getLNmeanRaw(DataAll)
  LNMeanAll <- LNMeanAll/InertLNMean
  return(LNMeanAll)
}