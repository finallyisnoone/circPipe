######circRNA and its co-linear transcript
#mRNA

library(dplyr)
library(ggplot2)
library(readr)

###########mean expression of circRNA and co-linear
##########data loading
linear_FPKM0 =  read_delim("/data1/yeying/m6a_circ/analysis/co_linear/gencode.rsem.fpkm.mat",delim = "\t")
colnames(linear_FPKM0)[-1] = as.character(paste0("S_",colnames(linear_FPKM0)[-1]))
linear_FPKM = as.matrix(linear_FPKM0[,as.character(colnames(countData_methStatus))[-c(1:2)]])
rownames(linear_FPKM) = as.character(linear_FPKM0$id)
dim(linear_FPKM)
#circ_RPM_mat0 = circ_input_RPM_m6A2#read_delim("/data1/yeying/m6a_circ/analysis/co_linear/circ_input_correct_RPM_feature_9_14.matrix",delim = "\t")
colnames(circ_input_RPM_0) = gsub("_RPM","",colnames(circ_input_RPM_0))
circ_RPM_mat = as.matrix(circ_input_RPM_feature[,grep("S_",colnames(circ_input_RPM_feature))])
rownames(circ_RPM_mat)=as.character(circ_input_RPM_feature$id)
head(circ_RPM_mat)
circ_RPM_mean = data.frame(gene = circ_input_RPM_feature$ensemble_id,
                           circ_name=circ_input_RPM_feature$id,mean_circ=rowMeans(circ_RPM_mat))
linear_FPKM_mean = data.frame(gene=as.character(linear_FPKM0$id),mean_linear=rowMeans(linear_FPKM))
circ_linear_mean = merge(circ_RPM_mean,linear_FPKM_mean,by="gene",all.x = T)
circ_linear_mean_over1 = filter(circ_linear_mean, mean_linear > 1)

save.image("/data1/yeying/m6a_circ/script/.RData") ###for R studio crashing ......

###############cor  : circ and co-linear with mean expression level
cor.test(circ_linear_mean_over1$mean_circ,circ_linear_mean_over1$mean_linear,method = "kendall")
cormethod="spearman"  #"pearson", "kendall", "spearman"
pcutoff=0.05

for (i in c("pearson", "kendall","spearman")) {
  print(i)
  print(cor.test(circ_linear_mean_over1$mean_circ,circ_linear_mean_over1$mean_linear,method = i))
  #plot 
  CoxExpPoint = Cor_plot(circ_linear_mean_over1, "mean_circ","mean_linear", cormethod= i)+
    ggtitle("cor of expression level of circRNAs and co-linear")
  print(CoxExpPoint)
}

###################correlation by samples
####reform df for cor analysis
linear_FPKM_tmp  <-  linear_FPKM0[,c("id",as.character(sample_name))]
colnames(linear_FPKM_tmp)[-1] <- paste0(colnames(linear_FPKM_tmp)[-1],"_linear")
circ_RPM_mat_tmp <- circ_input_RPM_feature[,c("ensemble_id","id","m6Astatus",as.character(sample_name))]
colnames(circ_RPM_mat_tmp) <- c(colnames(circ_RPM_mat_tmp)[c(1:3)],paste0(colnames(circ_RPM_mat_tmp)[-c(1:3)],"_circ"))
names(linear_FPKM_tmp)[1] <- "ensemble_id"
###
circ_linear_cor.mat <- left_join(circ_RPM_mat_tmp,linear_FPKM_tmp,by="ensemble_id")
dim(circ_linear_cor.mat)
rm(linear_FPKM_tmp,circ_RPM_mat_tmp)

############
circ_linear_cor.df <- cor_for_mat(cor_df.mat = circ_linear_cor.mat,tag1 = "_linear",tag2 = "_circ")
#write_delim(cor_df,path = "matrix/circ_linear_RPM_cor_by_sample.matrix",delim = "\t",na = "NA",col_names = TRUE)
circ_linear_cor.df <- na.omit(circ_linear_cor.df) #remove NA row
summary(circ_linear_cor.df)

#####################plot  
#pdf("plot/11_7/crosstalk_circ_co-linear.pdf")
#####hist of cor estimate
ggplot(circ_linear_cor.df,aes(x=cor_value))+theme_classic()+geom_histogram()+
  ggtitle("hist of correlation estimate of circRNA and co-linear")+   
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black", hjust = 0.25),
        axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))

#####ecdf plot of cor estimate
source("/data1/yeying/m6a_circ/script/R_function.R")
##group by circ methylation status
ks_res <- ks.test(filter(circ_linear_cor.df , m6Astatus == "m6A")$cor_value, 
                  filter(circ_linear_cor.df , m6Astatus == "non-m6A")$cor_value)
ECDF_plot(df = circ_linear_cor.df,value_var = "cor_value" , group_var = "m6Astatus" ,test_result = ks_res ,
          plot_title = "The correlation of co-lnear of m6A and non-m6A circRNA")+
  scale_color_manual(values=c("#DE8E43","#364D54")) ### color config manually

##group by circ and co-linear methylation status
mRNA_m6Astatus_mat=read_delim("/data1/yeying/m6a_circ/analysis/co_linear/gencode.rsem.fpkm_m6Astatus.txt",delim ="\t")
mRNA_m6Astatus <- data.frame(ensemble_id=mRNA_m6Astatus_mat$id,m6Astatus_linear=mRNA_m6Astatus_mat$m6Astatus)
cor_df_linear <- merge(circ_linear_cor.df,mRNA_m6Astatus,by="ensemble_id")
cor_df_linear <- cor_df_linear %>%
    mutate(meth_type = case_when(
            m6Astatus == "m6A" & m6Astatus_linear == "m6A" ~ "co-meth",
            m6Astatus == "non-m6A" & m6Astatus_linear == "m6A" ~ "linear_meth_uniq",
            m6Astatus == "m6A" & m6Astatus_linear == "non-m6A" ~ "circ_meth_uniq",
            m6Astatus == "non-m6A" & m6Astatus_linear == "non-m6A" ~ "co-nonMeth"))

kw_res <- kruskal.test(cor_value ~ as.factor(meth_type) ,data=cor_df_linear)
ECDF_plot(df = cor_df_linear,value_var = "cor_value" , group_var = "meth_type" ,test_result = kw_res ,
          plot_title = "The correlation of co-linear of m6A and non-m6A circRNA")+scale_colour_jama()
  
# dev.off()
# rm(mRNA_m6Astatus_mat)
# save.image("/data1/yeying/m6a_circ/script/.RData") ### R studio crashing ......

###################peak distribution

density_plot <- function(anno_file){
  circ_pos_anno <- read.delim(anno_file,sep = "\t",header  = F)
  peak_freq=c(as.numeric(as.vector(circ_pos_anno[which(circ_pos_anno[,13]=="5UTR"),14])),
              as.numeric(as.vector(circ_pos_anno[which(circ_pos_anno[,13]=="CDS"),14]))+100,
              as.numeric(as.vector(circ_pos_anno[which(circ_pos_anno[,13]=="3UTR"),14]))+200)
  #plot(density(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A peak density",xlab=NA,main=NA,col="darkblue",xaxt="n")
  #plot(table(peak_freq)/length(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A coding peak density",xlab=NA,main=NA,col="darkblue",xaxt="n",type="l") #xingyang
  lines(table(peak_freq)/length(peak_freq),xlim=c(0,300),lwd=3,xlab=NA,main=NA,col="darkblue",xaxt="n",type="l")
  abline(v=100,lty=2,col="black",lwd=1)
  abline(v=200,lty=2,col="black",lwd=1)
  text(x=c(60,160,260),y=-0.0015,pos=2,labels=c("5'UTR","CDS","3'UTR"),xpd=TRUE,font=2)
}

density_plot("/data1/yeying/m6a_circ/analysis/distribution//all_circ_pos.anno.txt")
density_plot("/data1/yeying/m6a_circ/analysis/distribution//all_circ_peak.anno.txt")
density_plot("/data1/yeying/m6a_circ/analysis/distribution//all_linear_bed12.anno.txt")

##############motif



               