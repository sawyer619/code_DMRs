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


OneSVMmodel <- function(x, data, CV.int=10, CVnum=20){
  # 函数说明：对一个数据集划分训练集与测试集，训练并预测
  #
  # Args：
  # x：data数据集对应某些列号，x为交叉验证所得
  # data：数据集
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  #
  # Returns：
  # list(accuracyRate=accuracyRate,fSummaryRSVM=fSummaryRSVM)
  # accuracyRate：测试集准确率，fSummaryRSVM：最具诊断效果位点位置
  
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
  #message('svmpred:',svmpred)
  y_test <- as.numeric(as.vector(testData$y))
  #message('testData.y:',y_test)
  accuracyRate <- sum(svmpred == y_test )/dim(testData$x)[1]
  return(list(accuracyRate=accuracyRate, fSummaryRSVM=fSummaryRSVM))
}


OneSVML2 <- function(x, data_train, CV.L2=10, CV.int=10, CVnum=20){
  # 函数说明：nest cross vacation第2层交叉验证
  # 
  # Args：
  # x：第几次重复L2交叉验证
  # data_train：经过L1层划分后的训练集
  # CV.L2；划分L2层数据为多少折
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  # 
  # Returns：
  # 多个这种类型list(accuracyRate=accuracyRate, fSummaryRSVM=fSummaryRSVM)
  
  message('L2第几次重复交叉验证：', x)
  folds_L2 <- createFolds(data_train$target, k=CV.L2) #
  
  start.L2 <- Sys.time()
  cv_results_L2 <- lapply(folds_L2,OneSVMmodel,data_train, CV.int, CVnum)
  end.L2 <- Sys.time()
  cat("交叉验证建模时间：",end.L2 - start.L2)
  
  acRate <- lapply(cv_results_L2,function(x){ return(x$accuracyRate) })
  acRate <- unlist(acRate)
  ERInd <-  which(acRate == max(acRate))
  newERInd <- lapply(ERInd, function(x){return(length(cv_results_L2[[x]]$fSummaryRSVM$SelInd))})
  newERInd <- unlist(newERInd)
  #message('newERInd:',newERInd)
  maxERInd <-  max(which(newERInd == max(newERInd)))
  #message('maxERInd:',maxERInd)
  return(cv_results_L2[[maxERInd]])
}

OneSVML2_V2 <- function(x, data_train, CV.L2=10, CV.int=10, CVnum=20){
  # 函数说明：nest cross vacation第2层交叉验证
  # 
  # Args：
  # x：第几次重复L2交叉验证
  # data_train：经过L1层划分后的训练集
  # CV.L2；划分L2层数据为多少折
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  # 
  # Returns：
  # 多个这种类型list(accuracyRate=accuracyRate, fSummaryRSVM=fSummaryRSVM)
  
  message('L2第几次重复交叉验证：', x)
  folds_L2 <- createFolds(data_train$target, k=CV.L2) #
  test.x <- folds_L2[[1]]
  
  start.L2 <- Sys.time()
  cv_results_L2 <- OneSVMmodel(test.x, data_train, CV.int, CVnum)
  end.L2 <- Sys.time()
  cat("一次建模时间：",end.L2 - start.L2)
  
  return(cv_results_L2)
}
##--确定特征数权重--##

SplitDataL1 <- function(x, preData.new, selWeight.DMRsite=200, times.crossL2=20, CV.L2=10, 
                        CV.int=10, CVnum=20){
  # 函数说明：L1层划分后，进一步进行第2次划分，并确定进行几次L2层交叉验证
  # 
  # Args：
  # x：L1层10折所得列号
  # selWeight.DMRsite：提取前面多少重要位点
  # times.crossL2：L2层进行重复多少次10折交叉验证
  # CV.L2；划分L2层数据为多少折
  # CV.int：选择几折交叉验证，默认为10折
  # CVnum：number of CVs
  # 
  # Returns：
  # list(accuracyRate=accuracyRate,SelInd=SelInd)
  # accuracyRate：L1层测试集正确率，SelInd：最具诊断效果甲基位点位置
  
  
  data_train <- preData.new[-x,]
  data_test <- preData.new[x,]
  timesCross <- 1:times.crossL2  # 进行多少次L2层训练（参数选择地点）
  cvResultsL2_20 <- lapply(timesCross, OneSVML2_V2, data_train, CV.L2, CV.int, CVnum)  # 
  str(cvResultsL2_20)
  siteWeight <- lapply(cvResultsL2_20,function(x){return(x$fSummaryRSVM$SelInd)})
  siteWeight <- unlist(siteWeight)
  tableSiteWeight <- table(siteWeight)
  tableSiteWeight <- sort(tableSiteWeight, decreasing = T)  #权重从大到小排序
  #save(tableSiteWeight,file='tableSiteWeight.Rdata')
  siteEaxmple <- as.numeric(names(tableSiteWeight))
  message('length of siteEaxmple:',length(siteEaxmple))
  if (length(siteEaxmple) < selWeight.DMRsite)
    SelInd <- siteEaxmple
  else
    SelInd <- siteEaxmple[1:selWeight.DMRsite]  #挑选多少个重要位点特征（参数选择地点）
  
  trainData <- ConvertSVMdata(data_train)
  testData <- ConvertSVMdata(data_test)
  svmres <- svm( trainData$x[, SelInd],  trainData$y, scale=F, type="C-classification", kernel="linear" )
  svmpred <- predict(svmres, testData$x[, SelInd] )
  #message('svmpred:',svmpred)
  y_test <- as.numeric(as.vector(testData$y))
  #message('testData.y:',y_test)
  accuracyRate <- sum(svmpred == y_test )/dim(testData$x)[1]
  siteName <- colnames(trainData$x)
  return(list(accuracyRate=accuracyRate,SelInd=SelInd,siteName=siteName[SelInd]))
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
