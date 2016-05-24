
# 整个程序参数
args <- commandArgs(TRUE)  # linux下数据所在文件夹下
fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
# fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  # 数据所在文件夹


start <- Sys.time()


source('R_svm.R')
source('DMRsData_model.R')
#source('DMRsData_model.R', encoding = 'UTF-8')

cat("导入预处理后的数据\n")

# 步骤三：建立模型寻找最具有诊断功能的差异甲基化位点
# csv方式读入数据
# fileTrainData <- paste(fileData,'trainData.csv', sep = "/")
# #fileTrainData <- paste(fileData,'trainData_test.csv', sep = "/")
# preData.t <- read.csv(fileTrainData, stringsAsFactors =F)
# rownames(preData.t) <- preData.t[ ,1]
# trainData <- preData.t[,-1 ]
# rm(preData.t)

# Rdata方式读入数据
fileTrainData <- paste(fileData,'trainData.Rdata', sep = "/")
load(fileTrainData)

preData.dropSmallVar <- RemoveSmallVar(trainData, pct=1.0/3)  # 滤除1/3方差较小位点
rm(trainData)
preData <- t(preData.dropSmallVar)
preData <- as.data.frame(preData)
rm(preData.dropSmallVar)

# 参数设置：

CV.L1 <- 10  # nest cross vacation第1层选几折交叉验证，默认：10
CV.int <- 10  # R_svm交叉验证，默认：10
CVnum <- 20  # R_svm交叉验证次数，默认：20

# 建模寻找差异甲基化位点
cat("交叉验证寻找最佳差异甲基化位点，开始\n")
folds_L1 <- createFolds(preData$target,k=CV.L1)
cv_results_L1 <- lapply(folds_L1,SplitDataL1, preData, CV.int, CVnum)
save(cv_results_L1,file = paste(fileData, 'cv_results_L1.Rdata', sep = "/"))

# 寻找最佳差异甲基化位点特征数
accuracyRate <- lapply(cv_results_L1,function(x){return(x$accuracyRate)})
accuracyRate <- unlist(accuracyRate)
MinLevel <- lapply(cv_results_L1,function(x){return(x$fSummaryRSVM$MinLevel)})
MinLevel <- unlist(MinLevel)
accLevel <- data.frame(accuracyRate,MinLevel)
orderAcc <- accLevel[ order(-accLevel[,1], -accLevel[,2]),]
greatLevel <- orderAcc[1,2]
cat("最佳差异甲基化位点数目：", greatLevel, "\n")

# 对全部数据R_svm建模，寻找特征数
cat("多次重复寻找差异甲基化位点\n")
timesCross <- 1:10
times_results <- lapply(timesCross, OneTimeModelSVM, preData, CV.int, CVnum)
save(times_results,file = paste(fileData, 'times_results.Rdata', sep = "/"))

siteWeight <- lapply(times_results,function(x){return(x$siteName)})
siteWeight <- unlist(siteWeight)
tableSiteWeight <- table(siteWeight)
tableSiteWeight <- sort(tableSiteWeight, decreasing = T)  #权重从大到小排序
#save(tableSiteWeight,file='tableSiteWeight.Rdata')
siteEaxmple <- names(tableSiteWeight)
message('length of siteEaxmple:',length(siteEaxmple))
if(length(siteEaxmple) <= greatLevel)
  DMRsGene_site <- siteEaxmple
if(length(siteEaxmple) > greatLevel)
  DMRsGene_site <- siteEaxmple[1:greatLevel]  #挑选多少个重要位点特征（参数选择地点）
save(DMRsGene_site, file = paste(fileData, 'DMRsGene_site.Rdata', sep = "/"))
write.csv(DMRsGene_site, file = paste(fileData, 'DMRsGene_site.csv', sep = "/"), 
          row.names = F)

end <- Sys.time()
cat("程序运行时间：", end - start, "\n")