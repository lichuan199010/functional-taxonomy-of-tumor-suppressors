# convert read number to tumor cell number based on the spike-ins

### Make a summary the data ####
rm(list=ls())
sgIDList <- read.csv(file = "SgIDList.txt", sep = "\t", stringsAsFactors = F)

SampleFinal <- read.csv(file = "SampleInfo_All.txt", sep = "\t", stringsAsFactors = F)

SpiInfo <- read.csv(file = "BarcodeCountSummary.txt", sep = "\t", stringsAsFactors = F)
SpiBCSet <- SpiInfo$Barcode
SpiBCCore <- SpiBCSet[1:3]

FakeBC <- c()
Read2Eql <- c()
ReadNum <- c()
CellNum <- c()
sgIDNum <- c()
sgIDNum200 <- c()
FileID <- c()

for (i in 1:dim(SampleFinal)[1]){
  FileNames <- paste0(SampleFinal$UniqueCode[i], "_BarcodeClean.txt")
  Data <- read.csv(file = FileNames, header = T, sep="\t", stringsAsFactors = F)
  Spi <- Data$Count[Data$sgID == "SpiNew" & Data$BC %in% SpiBCCore]
  
  FakeBC <- c(FakeBC, sum(Data$sgID == "SpiNew" & !Data$BC %in% SpiBCSet))
  Correct <- sum(Spi[1:3])/3/500000
  Read2Eql <- c(Read2Eql, 2/Correct)
  
  MouseInfo <- SampleFinal[i, c(5,13:18)]
  Data$CellNum <- Data$Count/Correct
  
  Data <- cbind(Data, MouseInfo)
  ReadNum <- c(ReadNum, sum(Data$Count))
  CellNum <- c(CellNum, sum(Data$CellNum))
  write.table(Data, paste0(SampleFinal$UniqueCode[i], "_BarcodeFinal.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = F,
              col.names = TRUE)
  sgIDNum <- c(sgIDNum, dim(Data)[1])
  sgIDNum200 <- c(sgIDNum200, sum(Data$CellNum >= 200))
  FileID <- c(FileID, SampleFinal$UniqueCode[i])
}

Output <- data.frame(FileID, FakeBC, Read2Eql, ReadNum, CellNum, BCNum=sgIDNum, BCNum200Cells=sgIDNum200)
Output <- cbind(Output, SampleFinal)
write.table(Output, "SampleInfo.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F,
            col.names = TRUE)