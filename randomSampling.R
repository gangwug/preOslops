###### randomly select 97% (519) of all epidermis samples for maximizing the circadian signal
rm(list=ls())
library(ape)
library(dplyr)
options(stringsAsFactors = FALSE)
## On R version 3.6, need to set below
RNGkind(sample.kind = "Rounding")
set.seed(10)
##the 50000 randomly selection of 519 samples
rnumv <- 50000
samplev <- 519
##the expression profile
dataD <- read.csv("./CYCLOPSv3.0.2.1/example.ToRunCYCLOPS.csv", stringsAsFactors = FALSE)
rownames(dataD) <- dataD$Gene.Symbol
##the benchmark correlation matrix
benchD <- read.csv("./CYCLOPSv3.0.2.1/skinBenchmarkMatrixARNTLorder.csv", stringsAsFactors = FALSE)
colnames(benchD)[1] <- colnames(dataD)[1] <- "Gene.Symbol"
rownames(benchD) <- benchD$Gene.Symbol
##the overlapped genes
bothID <- inner_join(data.frame(id = dataD$Gene.Symbol), data.frame(id = benchD$Gene.Symbol), by = "id")
max_index <- max_zstat <- 0
##the randomly sampling step
for (i in 1:rnumv)  {
  tindex <- sample(2:ncol(dataD), samplev)
  clock_dataD <- dataD[,c(1, tindex)]
  corA <- as.matrix(benchD[bothID$id,bothID$id])
  corB <- cor(t(clock_dataD[bothID$id,-1]), method = "spearman")
  colnames(corB) <- rownames(corB) <- bothID$id
  ##calcualate the similarity between matrix with 'ape' package
  simD <- mantel.test(corA, corB, nperm = 1000)
  if (simD$z.stat > max_zstat)  {
    max_zstat <- simD$z.stat
    max_index <- c(i, 1, tindex)
  }
  if (!(i %% 5000))  {print(i)}
}
max_zstat
max_index
outD <- dataD[, max_index[-1] ]
write.csv(outD, file = "./CYCLOPSv3.0.2.1/example.97Max.ToRunCYCLOPS.csv", row.names = FALSE)
