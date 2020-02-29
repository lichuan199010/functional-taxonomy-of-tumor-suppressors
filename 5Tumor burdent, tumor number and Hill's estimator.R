getRelTumorNum <- function(DataAll){
  # calculate the fraction of sgID changes relative to the total
  TumorNum <- c()
  for (sgID in sgIDAll){
    TumorNum <- c(TumorNum, sum(DataAll$sgID == sgID))
  }
  Ind <- which(sgIDAll %in% Inert)
  return(TumorNum/sum(TumorNum[Ind]))
}


getRelTumorBurden <- function(DataAll){
  # Calculate the tumor burden relative to Inert
  # calculate the fraction of sgID changes relative to the total
  CellNum <- c()
  for (sgID in sgIDAll){
    CellNum <- c(CellNum, sum(DataAll$CellNum[DataAll$sgID == sgID]))
  }
  Ind <- which(sgIDAll %in% Inert)
  return(CellNum/sum(CellNum[Ind]))
}

getHill <- function(Array, Perc=0.95, Cutoff = 200){
  # Hill's estimator
  # get the kth percentile
  Array <- Array[Array>=Cutoff]
  Upper = quantile(Array, Perc, na.rm = T)
  Below = max(Array[Array < Upper])
  Above = Array[Array >= Upper]
  HillEst = mean(log(Above) - log(Below)) # the larger this number is, the more extreme the tail is.
  return (c(HillEst, length(Above)))
}
