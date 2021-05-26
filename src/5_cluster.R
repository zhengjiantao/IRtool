library(flexclust)
library(SNFtool)
library(PINSPlus)
library(ConsensusClusterPlus)
library(SC3)
library(SingleCellExperiment)
library(infotheo)
library(SIMLR)
library(Seurat)

NMI_cal <- function(mylist1,mylist2){
    return(2*(entropy(mylist1)+entropy(mylist2)-entropy(cbind(mylist1,mylist2)))
           /(entropy(mylist1)+entropy(mylist2)))
}

ID = "AD_GSE95587"
setwd(paste("/home/fangzy/zhengjt",ID,sep="/"))
datalists = c("IR_gene_read_counts.txt", "IR_transcript_read_counts.txt")
label = read.csv(paste(ID,"design.txt",sep="_"), header=FALSE, sep="\t", row.names=1)
colnames(label) <- c("cell_type")

res_ARI = data.frame()
res_RI = data.frame()
res_NMI = data.frame()

k = 2

for (d in datalists){
    datalist = paste(ID,d,sep="_")
    print(datalist)
    data = read.csv(datalist,header=TRUE,sep="\t",row.names="ID")
    data = log2(as.matrix(data) + 1)
    colnames(data) = paste("sample",seq(1,ncol(data)),sep = "")
    rownames(data) = paste("gene",seq(1,nrow(data)),sep = "")
    set.seed(123)

    # 1.SC3
    print("SC3")
    sce <- SingleCellExperiment(
        assays = list(
            counts = as.matrix(data),
            logcounts = as.matrix(data)
        ), 
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce <- sc3(sce, ks = k, gene_filter = FALSE, biology = FALSE)
    res = as.data.frame(colData(sce)[,1])
    tab = table(res[,1], label[,1])
    ARI = randIndex(tab, correct = TRUE, original = TRUE)
    SC3_ARI = ARI[[1]]
    SC3_RI = ARI[[2]]
    SC3_NMI = NMI_cal(res[,1],label[,1])
    
    # 2.CC
    # print("CC")
    # dc = sweep(data,1, apply(data,1,median))
    # dc = as.matrix(dc)
    # rcc_SC3 = ConsensusClusterPlus(dc, maxK = 10, plot='png', 
    #                                title = "CC")
    # rm(dc)
    # tab <- table(rcc_SC3[k][[1]]$consensusClass, label[,1])
    # ARI = randIndex(tab, correct = TRUE , original = TRUE)
    # CC_ARI = ARI[[1]]
    # CC_RI = ARI[[2]]
    # CC_NMI = NMI_cal(rcc_SC3[k][[1]]$consensusClass,label[,1])
    
    #3. SNF 竖着为样本量
    # data 是行是基因，列是样本
    # print("SNF")
    # ## First, set all the parameters:
    # K = 20;  	# number of neighbors, usually (10~30)
    # alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
    # T = 20; 	# Number of Iterations, usually (10~20)
    # SNF_data = standardNormalization(t(data))
    # Dist = (dist2(as.matrix(SNF_data),as.matrix(SNF_data)))^(1/2)
    # ## next, construct similarity graphs
    # W = affinityMatrix(Dist, K, alpha)
    # SNF_result = spectralClustering(W, k)
    # #ARI
    # tab = table(SNF_result, label[,1])
    # ARI = randIndex(tab,correct = TRUE , original = TRUE)
    # SNF_ARI = ARI[[1]]
    # SNF_RI = ARI[[2]]
    # SNF_NMI = NMI_cal(SNF_result, label[,1])
    # rm(SNF_data)
    
    # # 4. SIMLR
    # print("SIMLR")
    # cluster_SIMLR <- SIMLR(X = data, c = k)
    
    # tab <- table(cluster_SIMLR$y$cluster, label[, 1])
    # ARI <- randIndex(tab, correct = TRUE, original = TRUE)
    # SIMLR_ARI <- ARI[[1]]
    # SIMLR_RI <- ARI[[2]]
    # SIMLR_NMI <- NMI_cal(cluster_SIMLR$y$cluster, label[, 1])

    tmp_ARI = data.frame(dataset=datalist, SC3=SC3_ARI) #, CC = CC_ARI, SIMLR = SIMLR_ARI, SNF = SNF_ARI)
    res_ARI = rbind(tmp_ARI, res_ARI)
    tmp_RI = data.frame(dataset=datalist, SC3=SC3_RI) #, CC = CC_RI, SIMLR = SIMLR_RI, SNF = SNF_RI)
    res_RI = rbind(tmp_RI, res_RI)
    tmp_NMI = data.frame(dataset=datalist, SC3=SC3_NMI) #, CC = CC_NMI, SIMLR = SIMLR_NMI, SNF = SNF_NMI)
    res_NMI = rbind(tmp_NMI, res_NMI)
}
write.csv(res_RI, paste(ID, "res_RI.csv", sep = "_"))
write.csv(res_ARI, paste(ID, "res_ARI.csv", sep = "_"))
write.csv(res_NMI, paste(ID, "res_NMI.csv", sep = "_"))
