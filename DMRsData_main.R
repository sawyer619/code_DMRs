# 整个程序参数
args <- commandArgs(TRUE)  # linux下数据所在文件夹下
fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
# fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  # 数据所在文件夹

start <- Sys.time()

#setwd("C:/DMRs_data/Liver")

library(ChAMP)
library(stringr)

source('DMRsData_preprocess.R', encoding = 'UTF-8')
source('R_svm.R')
source('DMRsData_model.R', encoding = 'UTF-8')

# 步骤一：从信息文件中获取每个甲基化样本信息，写入KIRC_all_lung_test_set.csv
fileInfo <- FindFileInfo(fileData)
fileSummyLung <- 'summary_lung.csv'  #标准数据形式
filelLungTestSet <- SureEachIdatLabel(fileInfo, fileSummyLung, fileData)


# 步骤二：提取多个批次进行数据预处
#filelLungTestSet <-  'F:/ChAMP_data_DMRs/code_DMRs/data_KIRC/KIRC_all_lung_test_set.csv'
# 得到批次名称
namesOfBatch <- HowManyBatchOfData(filelLungTestSet, dropSmallBatch = 1)
standardBatch <- namesOfBatch[1:2]
#standardBatch <- c(6285617025,6285650047)
# fileA <- "/home/scfan/Data/ForYing/code_liver/afterComBatData"
# existBatch <- FindExistCSV(fileA)
# standardBatch <- c(standardBatch,existBatch)

otherBatch <- setdiff(namesOfBatch, standardBatch)

fileWriteInfo <- paste(fileData, sep = "/",
                       "DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")
dealStandardBatchData <- DealSomeBatchData(NULL, standardBatch, filelLungTestSet,fileWriteInfo, fileData)
save(dealStandardBatchData, file=paste0(fileData,"afterComBatData/dealStandardBatchData.Rdata"))
dealAllBatchData <- lapply(otherBatch, DealSomeBatchData,standardBatch,filelLungTestSet, fileWriteInfo, fileData)
save(dealAllBatchData, file=paste0(fileData,"/afterComBatData/dealAllBatchData.Rdata"))

# 求纠正批次效应后，剩下共同的位点
# load("afterComBatData/dealStandardBatchData.Rdata")
# load("afterComBatData/dealAllBatchData.Rdata")
standardExample <- dealStandardBatchData$col.name  # 标准批次样本名
com.siteOfAllData <- lapply(dealAllBatchData, function(x){return(x$row.name)})
com.site <- Reduce(intersect, com.siteOfAllData) 
# 将标准批次写入csv文件中
#infoStandData <- ComSiteExample_stand(dealStandardBatchData, com.site)
infoStandData <- ComSiteExample_v1(dealStandardBatchData, standardExample, com.site)
save(infoStandData, file=paste0(fileData,"finalDataOfBatch/infoStandData.Rdata"))
infoFinaData <- lapply(dealAllBatchData, ComSiteExample_v1, standardExample, com.site)
save(infoFinaData, file=paste0(fileData,"finalDataOfBatch/infoFinaData.Rdata"))

# 合并数据，行为位点，列为样本
# load("finalDataOfBatch/infoStandData.Rdata")
# load("finalDataOfBatch/infoFinaData.Rdata")

#base <- "finalDataOfBatch"
base <- strsplit(infoStandData[[1]], split="/")[[1]][1]
fileBase <- paste(fileData, base, sep = "/")
# filePreproData <- UnionData(base, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)
# 
# # 将其转换为SVM识别数据，并加上target
# #filePreproData <- "Data_prepro.csv"
# fileTrainData <- ConvertPreproToModelData(filePreproData, fileData)  #转化文件使其模型能够识别

fileTrainData <- UnionAndConvertDataToModel(fileBase, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)


# 步骤三：建立模型寻找最具有诊断功能的差异甲基化位点
#fileTrainData <- paste(fileData,'trainData.csv', sep = "/")
#fileTrainData <- paste(fileData,'trainData_test.csv', sep = "/")
preData.t <- read.csv(fileTrainData, stringsAsFactors =F)
rownames(preData.t) <- preData.t[ ,1]
preData.t <- preData.t[,-1 ]
preData.dropSmallVar <- RemoveSmallVar(preData.t, pct=1.0/3)  # 滤除1/3方差较小位点
rm(preData.t)
preData <- t(preData.dropSmallVar)
preData <- as.data.frame(preData)
rm(preData.dropSmallVar)

# 参数设置：

CV.L1 <- 5  # nest cross vacation第1层选几折交叉验证，默认：10
CV.L2 <- 5  # nest cross vacation第2层选几折交叉验证，默认：10
selWeight.DMRsite <- 200  # 选择前多少个最重要甲基化位点位置，默认：200
times.crossL2 <- 2  # nest cross vacation第2层重复多少次，默认：20
CV.int <- 10  # R_svm交叉验证，默认：10
CVnum <- 2  # R_svm交叉验证次数，默认：20

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

