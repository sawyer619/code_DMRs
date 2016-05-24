#-----------------------------------------------#
# 程序说明：建模识别最具有诊断效应的位点
# 时间：2016.4.21
#-----------------------------------------------#

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
  if(ncol(x) > 5000){
    ladder1 <- CreatLadder( ncol(x), pRatio=0.5, Nmin=5000 )
    ladder2 <- CreatLadder( 5000, pRatio=0.7, Nmin=1000 )
    ladder3 <- CreatLadder( 1000, pRatio=0.85, Nmin=10 )
    ladder <- c(ladder1, ladder2, ladder3)
  }
  
  if(ncol(x) <= 5000 & ncol(x) > 1000){
    ladder1 <- CreatLadder( ncol(x), pRatio=0.7, Nmin=1000 )
    ladder2 <- CreatLadder( 1000, pRatio=0.85, Nmin=10 )
    ladder <- c(ladder1, ladder2)
  }
  
  if(ncol(x) <= 1000){
    ladder <- CreatLadder( ncol(x), pRatio=0.85, Nmin=10 )
  }
  
  
  RSVMres <- RSVM(x, y, ladder, CVtype=CV.int, CVnum=CVnum )   #CVnum次数（参数选择地点）
  fSummaryRSVM <- SummaryRSVM( RSVMres )
  return(fSummaryRSVM)
}


##--确定特征数权重--##

SplitDataL1 <- function(x, data, CV.int=10, CVnum=20){
  # 函数说明：L1层划分后，进一步进行第2次划分，并确定进行几次L2层交叉验证
  # 
  # Args：
  # x：L1层10折所得列号
  # data：训练数据
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  # 
  # Returns：
  # list(accuracyRate=accuracyRate,SelInd=SelInd)
  # accuracyRate：L1层测试集正确率，SelInd：最具诊断效果甲基位点位置
  time1 <- Sys.time()
  
  data_trainL2 <- data[-x,]
  data_testL2 <- data[x,]
  trainData <- ConvertSVMdata(data_trainL2)
  testData <- ConvertSVMdata(data_testL2)
  #message('yTrain:',table(trainData$y))
  #message('yTest:',table(testData$y))
  
  fSummaryRSVM <- ModelSVM(trainData, CV.int, CVnum)
  SelInd <- fSummaryRSVM$SelInd
  
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
  
  
  time2 <- Sys.time()
  cat("交叉验证中一次建模时间：",time2-time1, "\n")
  
  return(list(accuracyRate=accuracyRate, fSummaryRSVM=fSummaryRSVM, evaluate=evaluateSystem))
  
}

OneTimeModelSVM <- function(x, data, CV.int, CVnum){
  # 函数说明：一次R_svm训练，寻找差异甲基化位点
  #
  # Args：
  # x：第几次重复
  # data：训练数据
  #
  # return：
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  time3 <- Sys.time()
  
  cat("第几次寻找差异甲基化位点：", x,"\n")
  trainData <- ConvertSVMdata(data)
  fSummaryRSVM <- ModelSVM(trainData, CV.int, CVnum)
  SelInd <- fSummaryRSVM$SelInd
  siteName <- colnames(trainData$x)
  siteName <- siteName[SelInd]
  
  time4 <- Sys.time()
  cat("一次建模时间：",time4-time3, "\n")
  
  return(list(siteName=siteName, SummaryRSVM=fSummaryRSVM))
  
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
