####use Oscope to do selection on the eigen genes and output the eigen clusters used for CYCLOPS
###select groups of possible periodic eigenenes using Oscope
library(Oscope)
indir <- "./CYCLOPSv3.0.2.1/CYCLOPSout/"
infiles <- grep("EigengeneExp.csv", dir(indir, full.names = TRUE), value = TRUE)
outfiles <- gsub("EigengeneExp.csv", "EigengeneExp_OscopeCluster.csv", infiles)
for (ii in 1:length(infiles))  {
  infile <- infiles[ii]
  outfile <- outfiles[ii]
  eigenD <- read.csv(infile)
  OscopeData <- as.matrix(eigenD[,-1])
  dimnames(OscopeData) <- list("r"=eigenD[,1], "c"=colnames(eigenD)[-1] )
  DataInput <- NormForSine(OscopeData)
  ##apply sine model on the eigengenes
  SineRes <- OscopeSine(DataInput, parallel = TRUE)
  KMRes <- try(OscopeKM(SineRes, quan=0.4, maxK = 10), silent = TRUE)     
  outD <- NULL
  if (class(KMRes) == "list")  {
    ##check the number of each element
    cnum = unlist(lapply(KMRes, length))
    cindex = as.numeric(which(cnum > 1))
    if (length(cindex) < length(cnum))  {
      KMRes = KMRes[cindex]
    }
    ##Flag clusters with small within-cluster sine scores and/or small within-cluster phase shifts
    ToRM <- FlagCluster(SineRes, KMRes, DataInput)
    KMResUse <- KMRes[-ToRM$FlagID]
    if (length(KMResUse))  {
      eigenName <- sapply(KMResUse, function(z) {paste(z, collapse = "|")} )
      eigenIndex <- sapply(KMResUse, function(z) { gz <- gsub("eigen_(\\d+)_\\S+", "\\1", z, perl = TRUE)
      return(paste(sort(gz), collapse = "|"))
      } )
      outD <- data.frame(groupName = names(KMResUse), eigenName = eigenName, eigenIndex = eigenIndex)
    }
  }
  ##select those pair eigen genes with defined cutoff
  ##it may need to adjust this cutoff according to specific data
  cutoff <- as.numeric( quantile(SineRes$DiffMat, probs = seq(0, 1, by=0.1) )[2] )
  similarD <- SineRes$DiffMat
  peigenName <- peigenIndex <- NULL
  eigenID <- rownames(similarD)
  rown <- length(eigenID)
  for (i in 1:(rown-1) )  {
    for (j in (i+1):(rown) )  {
      if (similarD[j,i] <= cutoff)  {
        peigenName <- c(peigenName , paste(eigenID[i], eigenID[j], sep = "|") )
        peigenIndex <- c(peigenIndex, paste(sort(c(i,j)), collapse ="|"))
      }
    }
  }
  ##get the selected pair cluster
  pairD <- NULL
  if (length(peigenName))  {
    pairD <- data.frame(groupName = paste("cluster0", 1:length(peigenName), sep=""),
                        eigenName = peigenName, eigenIndex = peigenIndex )
  }
  ##output the Oscope selected and pair cluster
  outD <- rbind(outD, pairD)
  if (length(outD))  {
    ##get the duplicated rows
    dup <- which(duplicated(outD$eigenIndex) == TRUE)
    if (length(dup))  {
      write.csv(outD[-dup,], file = outfile, row.names = FALSE)
    }  else  {
      write.csv(outD, file = outfile, row.names = FALSE)
    }
  }
}
