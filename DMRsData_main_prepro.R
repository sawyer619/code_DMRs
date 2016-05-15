#--------------------------------------#
#程序说明：整个甲基化KIRC处理流程
#时间：2016.5.7
#--------------------------------------#

# 整个程序参数
args <- commandArgs(TRUE)  # linux下数据所在文件夹下
fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
# fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  # 数据所在文件夹

start <- Sys.time()

library(ChAMP)
library(stringr)

#source('DMRsData_preprocess.R', encoding = 'UTF-8')
source('DMRsData_preprocess.R')

# 步骤一：从信息文件中获取每个甲基化样本信息，写入KIRC_all_lung_test_set.csv
fileInfo <- FindFileInfo(fileData)
fileSummyLung <- 'summary_lung.csv'  #标准数据形式
filelLungTestSet <- SureEachIdatLabel(fileInfo, fileSummyLung, fileData)


# 步骤二：提取多个批次进行数据预处
#filelLungTestSet <-  'F:/ChAMP_data_DMRs/code_DMRs/data_KIRC/KIRC_all_lung_test_set.csv'
# 得到批次名称
namesOfBatch <- HowManyBatchOfData(filelLungTestSet, dropSmallBatch = 1)
standardBatch <- namesOfBatch[1:6]
otherBatch <- setdiff(namesOfBatch, standardBatch)

fileWriteInfo <- paste(fileData, sep = "/","DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")
dealStandardBatchData <- DealSomeBatchData(NULL, standardBatch, filelLungTestSet,fileWriteInfo, fileData)
save(dealStandardBatchData, file=paste0(fileData,"/afterComBatData/dealStandardBatchData.Rdata"))
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
save(infoStandData, file=paste0(fileData,"/finalDataOfBatch/infoStandData.Rdata"))
infoFinaData <- lapply(dealAllBatchData, ComSiteExample_v1, standardExample, com.site)
save(infoFinaData, file=paste0(fileData,"/finalDataOfBatch/infoFinaData.Rdata"))

# 合并数据，行为位点，列为样本
# load("finalDataOfBatch/infoStandData.Rdata")
# load("finalDataOfBatch/infoFinaData.Rdata")

#base <- "finalDataOfBatch"
base <- strsplit(infoStandData[[1]], split="/")[[1]][1]
fileBase <- paste(fileData, base, sep = "/")
# filePreproData <- UnionData(fileBase, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)
# 
# # 将其转换为SVM识别数据，并加上target
# #filePreproData <- "Data_prepro.csv"
# fileTrainData <- ConvertPreproToModelData(filePreproData, fileData)  #转化文件使其模型能够识别

fileTrainData <- UnionAndConvertDataToModel(fileBase, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)


end <- Sys.time()
end - start

