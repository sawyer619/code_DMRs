#  程序说明：欠采样癌症样本建立模型，找出差异甲基化位点
#  时间：2016.5.25
#source('DMRsData_main_model_v2.R', encoding = 'UTF-8')

# 整个程序参数
# args <- commandArgs(TRUE)  # linux下数据所在文件夹下
# fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  # 数据所在文件夹


start <- Sys.time()


source('R_svm.R')
# source('DMRsData_model_v2.R')
source('DMRsData_model_v2.R', encoding = 'UTF-8')

cat("导入预处理后的数据\n")

# 步骤三：建立模型寻找最具有诊断功能的差异甲基化位点
# # csv方式读入数据
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


# 参数：
underTimes <- 1
CV.int <- "bootstrape"  # R_svm交叉验证，默认：10  "bootstrape"
CVnum <- 20  # R_svm交叉验证次数，默认：20
timesCross <- 1:12  # 重复建模次数
freqDMRs <- 9  # 选择频次>= freqDMRs 的差异甲基化位点
CV.varify <- 10

# 欠采样癌症数据，R_svm建模，寻找特征数
cat("多次重复寻找差异甲基化位点\n")
times_results_v2 <- lapply(timesCross, OneTimeModelSVM, preData, underTimes, CV.int, CVnum)
save(times_results_v2,file = paste(fileData, 'times_results_v2.Rdata', sep = "/"))

siteWeight <- lapply(times_results_v2,function(x){return(x$siteName)})

cat("交叉验证，验证差异甲基化位点有效性\n")
folds_L1 <- createFolds(preData$target,k=CV.varify)  # 交叉验证得出的差异甲基化位点识别率
# 求timesCross次重复后共同的差异甲基化位点
com_gene <- Reduce(intersect, siteWeight)  # 求交集
if(length(com_gene) == 0){
  cat("欠采样建模，重复多次，交集为空\n")
}else{
  cat("交集，CV\n")
  com_CVverify <- lapply(folds_L1, VerifyDMRsGene, preData, com_gene)
  save(com_CVverify, file = paste(fileData, 'com_CVverify.Rdata', sep = "/"))
}
save(com_gene, file = paste(fileData, 'com_gene.Rdata', sep = "/"))
write.csv(com_gene, file = paste(fileData, 'com_gene.csv', sep = "/"), 
          row.names = F)


# 按频率次数提取前面频率高的位点
siteWeight <- unlist(siteWeight)
tableSiteWeight <- table(siteWeight)
tableSiteWeight <- sort(tableSiteWeight, decreasing = T)  #权重从大到小排序
#save(tableSiteWeight,file='tableSiteWeight.Rdata')
siteEaxmple <- names(tableSiteWeight)
message('length of siteEaxmple:',length(siteEaxmple))
bigFreq <- tableSiteWeight[which(tableSiteWeight >= freqDMRs)]
bigFreq_gene <- names(bigFreq)
if(length(bigFreq_gene) == 0){
  cat("取频次前", freqDMRs, "，差异甲基化位点为空\n")
}else{
  cat("频次前", freqDMRs,"次, CV\n")
  bigFriq_CVverify <- lapply(folds_L1, VerifyDMRsGene, preData, bigFreq_gene)
  save(bigFriq_CVverify, file = paste(fileData, 'bigFriq_CVverify.Rdata', sep = "/"))
}
save(bigFreq_gene, file = paste(fileData, 'bigFreq_gene.Rdata', sep = "/"))
write.csv(bigFreq_gene, file = paste(fileData, 'bigFreq_gene.csv', sep = "/"), 
          row.names = F)

end <- Sys.time()
cat("程序运行时间：", end - start, "\n")

