
############seperate matrix
paired_sample <- read.table("sample_info/paired.lst",col.names = "paired_sample")
unpaired_sample <- read.table("sample_info/unpaired.lst",col.names = "unpaired_sample")

paired_mat <- circ_input_RPM[,as.character(paired_sample$paired_sample)]
unpaired_mat <- circ_input_RPM[,as.character(unpaired_sample$unpaired_sample)]

##########form a function for wilxon-test for expre_matrix

row_wilcox <- function(group_mat,x,test_mode=""){
  group_mat$Type <- as.factor(group_mat$Type)
  typelevel <- levels(as.factor(group_mat$Type))
  comb_2 <- combn(typelevel,2)
  group1 <- as.character(rownames(subset(group_mat,Type==comb_2[1])))
  group2 <- as.character(rownames(subset(group_mat,Type==comb_2[2])))
  if (test_mode=="paired"){
    res_wix0 <- wilcox.test(x[which(rownames(group_mat)%in%group1)],x[which(rownames(group_mat)%in%group2)],paired = T)
  } else {
    res_wix0 <- wilcox.test(x[which(rownames(group_mat)%in%group1)],x[which(rownames(group_mat)%in%group2)])
  }
  
  res_wix <- c(Pvalue=res_wix0$p.value,statistic=res_wix0$statistic)
  return(res_wix)
}

mat_wilcox <- function(expre_mat,group_mat,test_mode){
  expre_mat=paired_mat
  test_mode="paired"
  keep <- rowSums(expre_mat > 0) >= 2 #a Count>0 in at least 3 samples
  expre_mat <- expre_mat[keep,]
  cat(dim(expre_mat))
  group_mat <- get_pheno(colnames(expre_mat))
  expre_mat <- expre_mat[,as.character(rownames(group_mat))]
  res_wix_lst <- apply(expre_mat,1,function(x) row_wilcox(group_mat,x,test_mode))
  res_wix_lst = as.data.frame(t(res_wix_lst))
  res_wix_lst$FDR = p.adjust(res_wix_lst$Pvalue,method = "BH")
  res_wix_lst$BY = p.adjust(res_wix_lst$Pvalue,method = "bonferroni")
  dim(subset(res_wix_lst,Pvalue <=0.05))
  return(res_wix_lst)
}

################################edgeR with paired samples

##############edgeR
source("/data1/yeying/m6a_circ/script/R_function.R")

paired_edgeR <- edgeR_test(expre_mat = paired_mat,design_mode = "paired",test_method = "LRT")
unpaired_edgeR <- edgeR_test(expre_mat = unpaired_mat,test_method = "LRT" )
unpaired_edgeR_tmp <-edgeR_test(expre_mat = unpaired_mat)
################
DE_unpaired_res = subset(unpaired_edgeR_tmp,FDR <= 0.05 & abs(logFC) >= 0.58)
unpaired_norm_edgeR=cpm(unpaired_mat) ###用标准化后的矩阵聚类
select <- DE_unpaired_res[order(DE_unpaired_res$logFC, decreasing = TRUE),]
unpaired_DE_edgeR=unpaired_norm_edgeR[rownames(select),]

volcano_plot(unpaired_edgeR_tmp,c("Normal","Tumor"))
heatmap_house(unpaired_DE_edgeR,get_pheno(colnames(unpaired_mat)),title_hp = "Heatmap of Different expressed  with unpaired samples")
###############
