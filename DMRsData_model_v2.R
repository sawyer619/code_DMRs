# 程序说明：欠采样建模函数集群
# 时间：2016.5.25

RemoveSmallVar <- function(preData, pct=1.0/3){
  # 函数说明：滤除预处理后数据方差较小的位点
  # 
  # Args：
  # preData：数据
  # pct：滤除数据所占百分比
  #
  # Return:
  # varPreData：滤除方差较小位点后的数据
  
  var.data <- apply(preData,1,var)
  #var.data <- unlist(var.data)
  
  newvar <- var.data[!(names(var.data) %in% c("target"))]
  newvar <- sort(newvar)
  
  # 画图
  d.var <- density(newvar)
  imageName1 <- "trainData_densityPlot.pdf"
  main1 <- paste("Density plot of var of trainData (",length(newvar)," probes)",sep="")
  pdf(imageName1)
  #plot(d.var)
  plot(d.var, main=main1, xlab="var of DMRs site")
  dev.off()
  
  dropSite <-  round(length(newvar) * pct)
  nameNeedSite <- names(newvar[dropSite:length(newvar)])
  nameNeedSite <- c(nameNeedSite, "target")
  varPreData <- preData[nameNeedSite,]
  return(varPreData)
  
}

EvaluateSystem <- function(y_test, svmpred){
  # 函数说明：构建多种评价指标，SP,SE,ACC,MCC
  # 
  # Args：
  # svmpred：预测值
  # y_test：真实值
  #
  # return：
  # 评价指标
  
  if (length(svmpred) != length(y_test)) 
    stop("svmpred and y_test must have the same length")
  
  
  commonSY <- data.frame(y_test, svmpred)
  normExam <- commonSY[which(commonSY[,1]== -1), ]
  if(dim(normExam)[1]==0){
    TN <- 0
    FP <- 0
  }
  else{
    tableNorm <- table(normExam[,2])  # 小到大排列，-1为正常，1为癌症
    nameNorm <- names(tableNorm)
    if(length(nameNorm)>1){
      TN <- tableNorm["-1"]
      FP <- tableNorm["1"]
    }
    if(length(nameNorm) == 1){
      if(nameNorm[1]=="-1"){
        TN <- tableNorm[1]
        FP <- 0
      }
      if(nameNorm[1]=="1"){
        TN <- 0
        FP <- tableNorm[1]
      }
      
    }
  }
  
  
  cancerExam <- commonSY[which(commonSY[,1]==1), ]
  if(dim(cancerExam)[1]==0){
    TP <- 0
    FN <- 0
  }
  else{
    tableCancer <-table(cancerExam[,2])  # 大到小排列，-1为正常，1为癌症
    nameCancer <- names(tableCancer)
    if(length(nameCancer)>1){
      TP <- tableCancer["1"]
      FN <- tableCancer["-1"]
    }
    if(length(nameCancer) == 1){
      if(nameCancer[1]=="1")
      {
        TP <- tableCancer[1]
        FN <- 0
      }
      if(nameCancer[1]=="-1")
      {
        TP <- 0
        FN <- tableCancer[1]
      }
      
    }
  }
  
  SP <- as.numeric(TN/(TN + FP))  # 特异度 正确判断非病人的率
  SE <- as.numeric(TP/(TP + FN))  # 灵敏度 正确判断病人的率
  ACC <- as.numeric((TP + TN)/(TP + FN + TN + FP))  # 正确率
  MCC <- as.numeric((TP*TN - FP*FN)/sqrt((TP + FN)*(TP + FP)*(TN + FP)*(TN + FN)))
  tabley <- table(y_test)
  
  return(list(SP=SP, SE=SE, ACC=ACC, MCC=MCC, taby_test=tabley))
}

ConvertSVMdata <- function(dd ){
  # 函数说明：将数据转换为SVM可识别数据
  # 
  # Args：
  # dd：ConvertPreproToModelData函数转换的数据
  #
  # Returns：
  # list类型，x：特征集，y：label
  
  namesOfdd <- colnames(dd) %in% c('target')
  x <- as.matrix( dd[!namesOfdd] )
  target.dd <- dd$target
  #target.dd[which(target.dd==0)] = -1
  y <- factor( target.dd )
  return(list(x=x, y=y))
}

ModelSVM <- function(dd, CV.int=10, CVnum=20){
  # 函数说明：建立模型训练，寻找最具有诊断效应位点
  # 
  # Args：
  # dd：数据
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  
  x <- dd$x
  y <- dd$y
  if(ncol(x) > 10000){
    ladder1 <- CreatLadder( ncol(x), pRatio=0.5, Nmin=10000 )
    ladder2 <- CreatLadder( 10000, pRatio=0.85, Nmin=5 )
    ladder <- c(ladder1, ladder2)
  }
  
  if(ncol(x) <= 10000 & ncol(x) > 1000){
    ladder1 <- CreatLadder( ncol(x), pRatio=0.7, Nmin=1000 )
    ladder2 <- CreatLadder( 1000, pRatio=0.95, Nmin=5 )
    ladder <- c(ladder1, ladder2)
  }
  
  if(ncol(x) <= 1000){
    ladder <- CreatLadder( ncol(x), pRatio=0.95, Nmin=500 )
  }
  
  
  RSVMres <- RSVM(x, y, ladder, CVtype=CV.int, CVnum=CVnum )   #CVnum次数（参数选择地点）
  fSummaryRSVM <- SummaryRSVM( RSVMres )
  return(fSummaryRSVM)
}



UndersampledTumor <- function(tumorData, normalData, underTimes){
  # 函数说明：对癌症样本欠采样后和正常样本组合成训练样本，剩余为测试样本
  #
  # Args：
  # tumorData：癌症样本
  # normalData：正常样本
  # underTimes：欠采样癌症样本为正常样本的几倍
  #
  # return：
  #
  
  len.norma <- dim(normalData)[1]
  len.tumor <- dim(tumorData)[1]
  
  if(underTimes*len.norma >= len.tumor)
    stop("underTimes is too big, tumor examples aren't too many!")
  
  SampInd <- seq(1, len.tumor)
  TrainInd <- sample(SampInd, floor(underTimes*len.norma), replace=F )
  TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
  
  testData <- tumorData[TestInd,]
  trainData <- rbind(tumorData[TrainInd,], normalData)
  
  return(list(testData=testData, trainData=trainData))
}


OneTimeModelSVM_underTumor <- function(x, preData, underTimes, CV.int, CVnum){
  # 函数说明：一次R_svm训练，寻找差异甲基化位点, 采用欠采样癌症样本方式
  #
  # Args：
  # x：第几次重复
  # preData：训练数据
  # underTimes: 欠采样癌症样本为正常样本的几倍
  #
  # return：
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  time3 <- Sys.time()
  
  cat("第几次寻找差异甲基化位点：", x,"\n")
  
  normalData <- subset(preData, target==-1)
  tumorData <- subset(preData, target==1)
  traintest <- UndersampledTumor(tumorData, normalData, underTimes)  # 欠采样癌症样本
  
  # R_svm建模
  trainData <- ConvertSVMdata(traintest$trainData)
  fSummaryRSVM <- ModelSVM(trainData, CV.int, CVnum)
  SelInd <- fSummaryRSVM$SelInd
  siteName <- colnames(trainData$x)
  siteName <- siteName[SelInd]
  
  # 测试集验证
  testData <- ConvertSVMdata(traintest$testData)
  svmres <- svm(trainData$x[, SelInd], trainData$y, scale=F, type="C-classification",
                kernel="linear" )
  svmpred <- predict(svmres, testData$x[, SelInd] )
  svmpred <- as.numeric(as.vector(svmpred))
  #str(svmpred)
  #message('svmpred:',svmpred)
  y_test <- as.numeric(as.vector(testData$y))
  #str(y_test)
  #message('testData.y:',y_test)
  
  # 评价指标
  evaluateSystem <- EvaluateSystem(y_test, svmpred)
  accuracyRate <- sum(svmpred == y_test )/dim(testData$x)[1]  # 正确率
  
  time4 <- Sys.time()
  cat("一次建模时间：",time4-time3, "\n")
  
  return(list(siteName=siteName, SummaryRSVM=fSummaryRSVM, evaluate=evaluateSystem, 
              accuracyRate= accuracyRate))
  
}

OversampledNorm <- function(tumorData, normalData, overTimes){
  # 函数说明：对正常样本过采样后和癌症样本组合成训练样本
  #
  # Args：
  # tumorData：癌症样本
  # normalData：正常样本
  # overTimes：过采样正常样本为癌症样本的几分之一
  #
  # return：
  #
  
  len.norma <- dim(normalData)[1]
  len.tumor <- dim(tumorData)[1]
  
  SampInd <- seq(1, len.norma)
  nOver <- floor(len.tumor/overTimes) - len.norma
  if(nOver > 0){
    TrainInd <- sample(SampInd, nOver, replace=T )  # 过采样正常样本为癌症样本的一半
    #TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
    
    #testData <- normalData[TestInd,]
    normalData.new <- rbind(normalData, normalData[TrainInd,])
    normalData <- normalData.new
    rm(normalData.new)
    
  }
  trainData <- rbind(tumorData, normalData)
  #return(list(testData=testData, trainData=trainData))
  return(trainData)
}

OneTimeModelSVM_overNorm <- function(x, preData, overTimes, CV.int, CVnum){
  # 函数说明：一次R_svm训练，寻找差异甲基化位点，采用过采样正常样本方式
  #
  # Args：
  # x：第几次重复
  # preData：训练数据
  # overTimes：过采样正常样本为癌症样本的几分之一
  #
  # return：
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  time3 <- Sys.time()
  
  cat("第几次寻找差异甲基化位点：", x,"\n")
  
  normalData <- subset(preData, target==-1)
  tumorData <- subset(preData, target==1)
  trainData <- OversampledNorm(tumorData, normalData, overTimes)  # 欠采样癌症样本
  rm(normalData)
  rm(tumorData)
  cat("过采样完成:", Sys.time()-time3,"\n")
  # R_svm建模
  trainData <- ConvertSVMdata(trainData)
  fSummaryRSVM <- ModelSVM(trainData, CV.int, CVnum)
  SelInd <- fSummaryRSVM$SelInd
  siteName <- colnames(trainData$x)
  siteName <- siteName[SelInd]
  cat("建模完成:", Sys.time()-time3,"\n")
  
  # 测试集验证
  testData <- preData
  rm(preData)
  testData <- ConvertSVMdata(testData)
  svmres <- svm(trainData$x[, SelInd], trainData$y, scale=F, type="C-classification",
                kernel="linear" )
  svmpred <- predict(svmres, testData$x[, SelInd] )
  svmpred <- as.numeric(as.vector(svmpred))
  #str(svmpred)
  #message('svmpred:',svmpred)
  y_test <- as.numeric(as.vector(testData$y))
  #str(y_test)
  #message('testData.y:',y_test)
  cat("测试完成:", Sys.time()-time3,"\n")
  
  # 评价指标
  evaluateSystem <- EvaluateSystem(y_test, svmpred)
  accuracyRate <- sum(svmpred == y_test )/dim(testData$x)[1]  # 正确率
  
  time4 <- Sys.time()
  cat("一次建模时间：",time4-time3, "\n")
  
  return(list(siteName=siteName, SummaryRSVM=fSummaryRSVM, evaluate=evaluateSystem, 
              accuracyRate= accuracyRate))
  
}



VerifyDMRsGene <- function(x, preData, geneDmrs){
  # 函数说明：交叉验证，验证最终所得差异甲基化位点对分类影响
  # 
  # x：L1层10折所得列号
  # preData：训练数据
  # geneDmrs: 差异甲基化位点
  # 
  # return：评价指标
  
  geneDmrs <- c(geneDmrs, "target")
  data <- preData[,geneDmrs]
  data_trainL2 <- data[-x,]
  data_testL2 <- data[x,]
  trainData <- ConvertSVMdata(data_trainL2)
  testData <- ConvertSVMdata(data_testL2)
  
  svmres <- svm(trainData$x, trainData$y, scale=F, type="C-classification",
                kernel="linear" )
  svmpred <- predict(svmres, testData$x )
  svmpred <- as.numeric(as.vector(svmpred))
  #str(svmpred)
  #message('svmpred:',svmpred)
  y_test <- as.numeric(as.vector(testData$y))
  #str(y_test)
  #message('testData.y:',y_test)
  
  # 评价指标
  evaluateSystem <- EvaluateSystem(y_test, svmpred)
  accuracyRate <- sum(svmpred == y_test )/dim(testData$x)[1]  # 正确率
  
  return(list(evaluate=evaluateSystem, acc=accuracyRate))
}


createFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
  # 函数说明：交叉验证
  # 
  # Args：
  # y-label
  # 
  # return：
  # index
  
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}