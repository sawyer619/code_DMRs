### R-code for R-SVM
### use leave-one-out / Nfold or bootstrape to permute data for external CV
### build SVM model and use mean-balanced weight to sort genes on training set
### and recursive elimination of least important genes
### author: Dr. Xin Lu, Research Scientist
###   Biostatistics Department, Harvard School of Public Health

library(e1071)
library(randomForest)

## read in SVM formated data in filename
## the format following the defination of SVMTorch
## the first line contains 2 integer: nSample  nFeature+1 
## followed by a matrix, each row for one sample, with the last column being +/1 1 for class label
ReadSVMdata <- function(filename)
{
  dd <- read.table( filename, header=F, skip=1)
  x <- as.matrix( dd[, 1:(ncol(dd)-1)] ) 
  y <- factor( dd[, ncol(dd)] )
  
  ret <- list(x=x, y=y)
}

## create a decreasing ladder for recursive feature elimination
CreatLadder <- function( Ntotal, pRatio=0.75, Nmin=5 )
{
  x <- vector()
  x[1] <- Ntotal
  for( i in 1:100 )
  {
    pp <- round(x[i] * pRatio)
    if( pp == x[i] )
    {
      pp <- pp-1
    }          
    
    if( pp >= Nmin )
    {
      x[i+1] <- pp
    } else
    {
      break
    }
  }
  
  x
}

## R-SVM core code
## input:
##    x: row matrix of data
##    y: class label: 1 / -1 for 2 classes
##    CVtype: 
##        integer: N fold CV
##        "LOO":    leave-one-out CV
##        "bootstrape": bootstrape CV
##    CVnum:   number of CVs
##        LOO: defined as sample size
##        Nfold and bootstrape:  user defined, default as sample size
## output: a named list
##    Error: a vector of CV error on each level
##    SelFreq: a matrix for the frequency of each gene being selected in each level
##             with each column corresponds to a level of selection
##             and each row for a gene
##          The top important gene in each level are those high-freqent ones
RSVM <- function(x, y, ladder, CVtype, CVnum=0 )
{
  ## check if y is binary response
  Ytype <- names(table(y))
  if( length(Ytype) != 2) 
  {
    print("ERROR!! RSVM can only deal with 2-class problem")
    return(0)
  }
  
  yy <- vector( length=length(y))
  yy[which(y==Ytype[1])] <- 1
  yy[which(y==Ytype[2])] <- -1        
  y <- yy
  
  x <- as.data.frame(x)  # data.frame
  y <- as.factor(y)  # factor
  
  ## check ladder
  if( min(diff(ladder)) >= 0 )
  {
    print("ERROR!! ladder must be monotonously decreasing")
    return(0);
  }
  
  if( ladder[1] != ncol(x) )
  {
    ladder <- c(ncol(x), ladder)
  }
  
  nSample <- nrow(x)
  nGene   <- ncol(x)
  SampInd <- seq(1, nSample)
  
  if( CVtype == "LOO" )
  {
    CVnum <- nSample
  } else
  {
    if( CVnum == 0 )
    {
      CVnum <- nSample
    }
  }
  
  ## vector for test error and number of tests
  ErrVec <- vector( length=length(ladder))
  names(ErrVec) <- paste("Lev_", ladder, sep="")
  ErrVec.train <- vector( length=length(ladder))
  names(ErrVec.train) <- paste("Lev_", ladder, sep="")
  Error.list <- list(NULL)
  length(Error.list) <- CVnum
  Error_train.list <- list(NULL)
  length(Error_train.list) <- CVnum
  
  nTests <- 0
  nTrains <- 0
  
  SelFreq <- matrix( 0, nrow=nGene, ncol=length(ladder))
  colnames(SelFreq) <- paste("Lev_", ladder, sep="")
  rownames(SelFreq) <- colnames(x)
  
  ## for each CV    
  for( i in 1:CVnum )
  {
    
    #cat("CVnum:", i ,"\n")
    ## split data
    if( CVtype == "LOO" )
    {
      TestInd <- i
      TrainInd <- SampInd[ -TestInd]
    } else
    {
      if( CVtype == "bootstrape" )
      {
        TrainInd <- sample(SampInd, nSample, replace=T )
        TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
      } else
      {
        ## Nfold
        TrainInd <- sample(SampInd, nSample*(CVtype-1)/CVtype )
        TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
      }
    }
    
    nTests <- nTests + length(TestInd)
    nTrains <- nTrains + length(TrainInd)
    
    ## in each level, train a SVM model and record test error
    xTrain <- x[TrainInd, ]
    yTrain <- y[TrainInd]
    
    xTest  <- x[TestInd,]
    yTest  <- y[TestInd]
    
    ## index of the genes used in the 
    SelInd <- seq(1, nGene)
    names(SelInd) <- colnames(x)
    for( gLevel in 1:length(ladder) )
    {
      #cat("gLevel:", gLevel ,"\n")
      #cat("SelInd:", str(SelInd) ,"\n")
      
      ## train SVM model and test error
      rfmodel <- randomForest(xTrain[, SelInd], yTrain, ntree=500,
                                keep.forest=T, importance=TRUE)
      if( CVtype == "LOO" )
      {
        svmpred <- predict(rfmodel, matrix(xTest[SelInd], nrow=1) )
        svmpred.train <- predict(rfmodel, xTrain[, SelInd] )
      } else
      {
        #cat("xTest[, SelInd]:", str(xTest[, SelInd]) ,"\n")
        #cat("rfmodel:",str(rfmodel),"\n")
        svmpred <- predict(rfmodel, xTest[, SelInd] )
        svmpred.train <- predict(rfmodel, xTrain[, SelInd] )
      }
      #cat("svmpred:",str(svmpred),"\n")
      #cat("yTest",str(yTest),"\n")
      ErrVec[gLevel] <- ErrVec[gLevel] + sum(svmpred != yTest )
      ErrVec.train[gLevel] <- ErrVec.train[gLevel] + sum(svmpred.train != yTrain )
      #cat("yTest.error:", mean(svmpred != yTest ) ,"\n")
      #cat("yTrain.error:", mean(svmpred.train != yTrain ) ,"\n")
      
      ## weight vector
      Wrf <- importance(rfmodel)
      #cat("Wrf:", str(Wrf) ,"\n")
      W <- Wrf[,"MeanDecreaseGini"]
      W <- sort(W, decreasing = T)
      
      if( gLevel < length(ladder) )
      {
        #SelInd <- SelInd[which(rkW > (ladder[gLevel] - ladder[gLevel+1]))]
        nameSelInd <- names(W)[1:ladder[gLevel+1]]
        
        ## record the genes selected in this ladder
        SelFreq[nameSelInd, gLevel] <- SelFreq[nameSelInd, gLevel] +W[1:ladder[gLevel+1]]
        SelInd <- SelInd[nameSelInd]
      }
    }
    #cat(ErrVec, "\n")
    #cat("ErrVec/nTests:", ErrVec/nTests,"\n")
    Error.list[[i]] <- ErrVec/nTests
    Error_train.list[[i]] <- ErrVec.train/nTrains
    
  }
  
  Error.test <- ErrVec/nTests
  Error.train <- ErrVec.train/nTrains
  bootError <- Error.test*0.632 + Error.train*0.368  # bootstrap error
  ret <- list(ladder=ladder, Error=Error.test, SelFreq=SelFreq, Error_train=Error.train,
              bootError=bootError, Error.list=Error.list,  Error_train.list=Error_train.list)
  return(ret)
  
}

SummaryRSVM <- function( RSVMres )
{
  ERInd <- min( which(RSVMres$Error == min(RSVMres$Error)) )
  MinLevel <- RSVMres$ladder[ERInd]
  FreqVec <- RSVMres$SelFreq[, ERInd]
  
  SelInd <- which( rank(FreqVec) >= (RSVMres$ladder[1]-MinLevel) )
  
  # bootstrap
  ERInd.boot <- min( which(RSVMres$bootError == min(RSVMres$bootError)) )
  MinLevel.boot <- RSVMres$ladder[ERInd.boot]
  FreqVec.boot <- RSVMres$SelFreq[, ERInd.boot]
  SelInd.boot <- which( rank(FreqVec.boot) >= (RSVMres$ladder[1]-MinLevel.boot) )
  
  #    print("MinCV error of", min(RSVMres$Error), "at", MinLevel, "genes" )
  ret.boot <- list( MinER=min(RSVMres$bootError), MinLevel=MinLevel.boot, SelInd=SelInd.boot)
  ret <- list( MinER=min(RSVMres$Error), MinLevel=MinLevel, SelInd=SelInd,bootstrapInfo=ret.boot,
               RSVMres=RSVMres)
  return(ret)
}


# setwd("D:/ChAMP_data_DMRs/data_KIRC")
# data <- read.csv('TrainafterPreProData_6.csv')
# dataSam <- data[sample(1:nrow(data),250,replace=F),]
# x <- as.matrix(dataSam[,-c(1,2)])
# y <- dataSam[,2]
# y[y==0] <- -1
# y <- factor(y)
# 
# ladder <- CreatLadder( ncol(x), pRatio=0.75, Nmin=5 )
# 
# RSVMres <- RSVM(x, y, ladder, 5, CVnum=0 )
# fSummaryRSVM <- SummaryRSVM( RSVMres )

