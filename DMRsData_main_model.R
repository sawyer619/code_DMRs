
# 整个程序参数
args <- commandArgs(TRUE)  # linux下数据所在文件夹下
fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
# fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  # 数据所在文件夹


start <- Sys.time()

#setwd("F:/ChAMP_data_DMRs/data_KIRC")

#library(caret)

source('R_svm.R')
source('DMRsData_model.R', encoding = 'UTF-8')

# 步骤三：建立模型寻找最具有诊断功能的差异甲基化位点
# csv方式读入数据
# fileTrainData <- paste(fileData,'trainData.csv', sep = "/")
# #fileTrainData <- paste(fileData,'trainData_test.csv', sep = "/")
# preData.t <- read.csv(fileTrainData, stringsAsFactors =F)
# rownames(preData.t) <- preData.t[ ,1]
# TrainData <- preData.t[,-1 ]
# rm(preData.t)

# Rdata方式读入数据
fileTrainData <- paste(fileData,'trainData.Rdata', sep = "/")
load(fileTrainData)

preData.dropSmallVar <- RemoveSmallVar(TrainData, pct=1.0/3)  # 滤除1/3方差较小位点
rm(TrainData)
preData <- t(preData.dropSmallVar)
preData <- as.data.frame(preData)
rm(preData.dropSmallVar)

# 参数设置：

CV.L1 <- 10  # nest cross vacation第1层选几折交叉验证，默认：10
CV.L2 <- 10  # nest cross vacation第2层选几折交叉验证，默认：10
selWeight.DMRsite <- 300  # 选择前多少个最重要甲基化位点位置，默认：200
times.crossL2 <- 10  # nest cross vacation第2层重复多少次，默认：20
CV.int <- 10  # R_svm交叉验证，默认：10
CVnum <- 10  # R_svm交叉验证次数，默认：10

# 建模寻找差异甲基化位点
folds_L1 <- createFolds(preData$target,k=CV.L1)
cv_results_L1 <- lapply(folds_L1,SplitDataL1, preData, selWeight.DMRsite,times.crossL2,
                        CV.L2, CV.int, CVnum)
save(cv_results_L1,file = paste(fileData, 'cv_results_L1.Rdata', sep = "/"))

# 确定最终版差异甲基化位点
siteWeight <- lapply(cv_results_L1,function(x){return(x$siteName)})
siteWeight <- unlist(siteWeight)
tableSiteWeight <- table(siteWeight)
tableSiteWeight <- sort(tableSiteWeight, decreasing = T)  #权重从大到小排序
#save(tableSiteWeight,file='tableSiteWeight.Rdata')
DMRsGene_site <- names(tableSiteWeight[which(tableSiteWeight >= 4)])
save(DMRsGene_site, file = paste(fileData, 'DMRsGene_site.Rdata', sep = "/"))
write.csv(DMRsGene_site, file = paste(fileData, 'DMRsGene_site.csv', sep = "/"), 
          row.names = F)

end <- Sys.time()
end - start