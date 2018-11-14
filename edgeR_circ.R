library("readr")
source("./R_function.R")
options(stringsAsFactors=F)

library(Biobase)
########cut off init
p_value <- 0.05 # p value of FDR
lfc <- lfc  #logFC
args <-commandArgs(T)
#args[1] <- "matrix/11_7//filter_input_count.mat"

########load data
countData1 <- as.data.frame(read_delim(args[1],delim = "\t"))
dim(countData1)
rownames(countData1) <- countData1$id
countData=countData1[,-1]

###########pheno data loading
colData <- get_pheno(colnames(countData),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal")
countData <- countData[,as.character(rownames(colData))]
colData$Type <- as.factor(colData$Type)
type_level <- levels(as.factor(colData$Type))
comb <- combn(type_level,2)

#########filter low-abandance circRNA ; the step has been done in node1:filter_circ
# countData <- countData[which(rowSums(countData > 0) >= 2),] #a Count>0 in at least 2 samples
# dim(countData)

################type specific
tumor_circ_tmp <- countData[,grep("T",colnames(countData))]
normal_circ_tmp <- countData[,grep("N",colnames(countData))]
tumor_uniq_circ_mat <- countData[which(rowSums(normal_circ_tmp) == 0 & rowSums(tumor_circ_tmp > 0)>= 2),]
normal_uniq_circ_mat <- countData[which(rowSums(tumor_circ_tmp) == 0 & rowSums(normal_circ_tmp > 0)>= 2),]
tumor_circ <- as.character(rownames(tumor_uniq_circ_mat))
normal_circ <- as.character(rownames(normal_uniq_circ_mat))
print(paste("tumor specific circRNAs:",dim(tumor_uniq_circ_mat)[1]))
print(paste("normal specific circRNAs:",dim(normal_uniq_circ_mat)[1]))
countData$id= as.character(rownames(countData))
shared_circ <- filter(countData, !id%in%c(tumor_circ,normal_circ))
rownames(shared_circ) <- shared_circ$id
shared_circ <- shared_circ[,!colnames(shared_circ)%in%"id"]
countData <- countData[,!colnames(countData)%in%"id"]
Circ_norm_edgeR=cpm(countData) ###norm expres_mat
print(paste("all share circRNAs:",dim(shared_circ)[1]))

##############edgeR
source("./R_function.R")

#sharedCirc_edgeR_tmp <- edgeR_test(expre_mat = shared_circ,group_mat = colData,test_method = "LRT" )
sharedCirc_edgeR <-edgeR_test(expre_mat = shared_circ,group_mat = colData)

DE_sharedCirc_res = subset(sharedCirc_edgeR,FDR <= p_value & abs(logFC) >= lfc)
select <- DE_sharedCirc_res[order(DE_sharedCirc_res$logFC, decreasing = TRUE),]
sharedCirc_DE_edgeR=Circ_norm_edgeR[rownames(select),]
#all DE
DE_list <- c(rownames(DE_sharedCirc_res),tumor_circ,normal_circ)

# ################plot
# pdf(file = paste0("DE/",comb[2,1], "_vs_", comb[1,1], "_DE_10_17.edgeR.pdf"))
# ########Vocano plot
# volcano_plot(sharedCirc_edgeR,c("Normal","Tumor"))
# ########heatmap
# heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR)),title_hp = "Heatmap of Different expressed all DE circRNA")
# heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR)),
#               title_hp = "Heatmap of Different expressed all DE circRNA",cluster_rule = "row")
# heatmap_house(sharedCirc_DE_edgeR,get_pheno(colnames(sharedCirc_DE_edgeR)),title_hp = "Heatmap of Different expressed shared_circRNA")
# ########PCA
# PCA_plot(DE_list,Circ_norm_edgeR,colData)
# 
# dev.off()

###output table
DE_sharedCirc_res$id=rownames(DE_sharedCirc_res)
#write.table(DE_sharedCirc_res, file = "DE/DE_sharedCircRNA_edgeR_10_17.xls",sep = "\t", quote = FALSE, row.names = F)
##all DE
final_DE <- data.frame(id=DE_list,type=c(rep("shared_circRNA",dim(DE_sharedCirc_res)[1]),
                                                 rep("tumor_uniq",length(tumor_circ)),rep("normal_uniq",length(normal_circ))))
tumor_circ_tmp <-  Circ_norm_edgeR[,grep("T",colnames(Circ_norm_edgeR))]
normal_circ_tmp <-  Circ_norm_edgeR[,grep("N",colnames(Circ_norm_edgeR))]
uniqCirc_res_df <- function(uniqcirc,uniqcirc_mat,default_FC=10^6){
  res_df<-data.frame(logFC=rep(default_FC,length(uniqcirc)),
                     logCPM=log2(rowMeans(uniqcirc_mat[uniqcirc,])+1),
                     F=rep("NA",length(uniqcirc)),PValue=rep(0,length(uniqcirc)),
                     FDR=rep(0,length(uniqcirc)),id=uniqcirc)
  return(res_df)
}

final_DE_res <- rbind(DE_sharedCirc_res,uniqCirc_res_df(tumor_circ,tumor_circ_tmp),uniqCirc_res_df(normal_circ,normal_circ_tmp,default_FC = -10^6))
final_DE_res <- merge(final_DE,final_DE_res,by="id")
final_DE_res$logFC <- as.numeric(final_DE_res$logFC)
final_DE_res = mutate(final_DE_res,
                      expres_regu = case_when(
                        logFC >= lfc ~ "up_regu",
                        logFC <= -lfc ~ "down_regu"
                      ))

#write.table(final_DE_res, file = "DE/DE_allcircRNA_filter_edgeR_10_17.xls",sep = "\t", quote = FALSE, row.names = F)

####################DE enrichment analysis
DE_circ_gene = filter(circ_feature,id%in%DE_list)
gene_list=unique(DE_circ_gene$symbol)
length(gene_list)

GoKegg(gene_list,"./")

