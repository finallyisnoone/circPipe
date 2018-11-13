library(survival)
library(dplyr)
library(survminer)
library(ggplot2)
library(readr)
#setwd("/data/xingyang/m6A_zhengjian")
options(stringsAsFactors=F)

#
#DE_DM_circ_mat0 = read_delim("/data1/yeying/m6a_circ/analysis/survival/DM_DE_input_correct_RPM_m6Astatus_9_14.matrix",delim = "\t")
str(DE_DM_circ_mat)

DE_DM_circ_mat0 = filter(circ_input_RPM_0,id%in%DE_DM)
colnames(DE_DM_circ_mat0) = gsub("_RPM","",colnames(DE_DM_circ_mat0))
DE_DM_circ_mat = as.matrix(DE_DM_circ_mat0[,-c(1,2)])
DE_DM_circ_mat=t(DE_DM_circ_mat)
colnames(DE_DM_circ_mat)=as.character(DE_DM_circ_mat0$id)
DE_DM_circ_mat=cbind(t_sample_id=gsub(pattern = "S_",replacement = "",as.character(rownames(DE_DM_circ_mat))),DE_DM_circ_mat)
colnames(DE_DM_circ_mat)
#####load clinic data
clinc_data=read.csv("survival/m6A_clinc_analysis_2.csv")
os_data0=merge(clinc_data,DE_DM_circ_mat,by="t_sample_id")
dim(os_data0)
os_data = os_data0[which(!is.na(os_data0$sample_id_num)),]
os_data$t_sample_id=as.character(os_data$t_sample_id)

str(os_data)
surv_data=Surv(time = os_data$month,event = os_data$status==1)

i=which(colnames(os_data)=="n_sample_id")+1
n=ncol(os_data)
survival_res=data.frame()
 # i=60
 # k=5
while(i<=n){
  if (length(which(os_data[,i] == 0)) > 40) {
    i=i+1
    next
  }
    p.Val_os_peak=data.frame() 
    os_data[,i]=as.numeric(as.character(os_data[,i]))
    temp=os_data[,i]
    for (k in round(length(temp)*0.25):round(length(temp)*0.75)){
      n_low = k
      n_high = length(temp)-k
      temp[order(temp)[(k+1):length(temp)]]=1
      temp[order(temp)[1:k]]=0
      temp=as.factor(temp)
      #analysis
      gene_surv <- surv_data~temp
      gene_survfit <- survfit(gene_surv)
      gene_survdiff <- tryCatch(survdiff(gene_surv), error = function(e) return(NA))
      gene_coxph=coxph(surv_data~temp)
      summary_gene_survfit <- summary(gene_survfit)
      summary_gene_coxph <- summary(gene_coxph)
      #output
      gene=colnames(os_data)[i]
      HR = summary_gene_coxph$coefficients[1,2]
      HR_2 = (gene_survdiff$obs[1]/gene_survdiff$exp[1])/(gene_survdiff$obs[2]/gene_survdiff$exp[2])
      #####p value
      pvalue_chisq = as.numeric(pchisq(gene_survdiff$chisq, length(gene_survdiff$n) - 1, lower.tail = FALSE))
      p_sctest <- as.numeric(summary_gene_coxph$sctest[[3]])       #Score (logrank) test
      p_waldtest <- as.numeric(summary_gene_coxph$waldtest[[3]])   #Wald test
      p_logtest <- as.numeric(summary_gene_coxph$logtest[[3]])
      up95 = summary_gene_coxph$conf.int[1,4]#exp(log(HR) + qnorm(0.975)*sqrt(1/fit1$exp[2]+1/fit1$exp[1]))
      low95 = summary_gene_coxph$conf.int[1,3]#exp(log(HR) - qnorm(0.975)*sqrt(1/fit1$exp[2]+1/fit1$exp[1]))
      zero_num = length(which(os_data[,i] == 0))
      p.Val_os_peak=rbind(p.Val_os_peak,c(gene,i,p_logtest,p_waldtest,p_sctest,pvalue_chisq,HR,HR_2,n_high,n_low,low95,up95,zero_num))
      }
    colnames(p.Val_os_peak)=c( "gene","col_number","p_logtest","p_waldtest","p_sctest","pvalue_chisq", 
                               "HR","HR_2", "n_high", "n_low","low95","up95","zero_num")
    survival_res = rbind(survival_res,as.data.frame(p.Val_os_peak))
    i=i+1
}

#########significance
sig_survival_res=as.data.frame(t(apply(survival_res,1,function(x){p_min <-min(as.numeric(x[3:6]))
return(c(x,p_min=p_min))
})))%>%filter(p_min <= 0.05 )%>%arrange(p_min,zero_num)

unique(sig_survival_res$gene)
names(sig_survival_res)[1] <- "id"
sig_survival_DE_DM_feature = merge(sig_survival_res,circ_feature,by = "id")
#write_delim(sig_survival_DE_DM_feature,path = "survival/sig_survival_DE_DM_feature_11_7.txt",delim = "\t",col_names = T)

#write.table(sig_survival_res,"/data1/yeying/m6a_circ/pancreatic/survival/sig_survival_res.txt",col.names = T,row.names = F,quote = F,sep="\t")

##########plot

plot_survial <- function(i,k,data_df,title_sur){
  #####i is col_num in sur
   # i=35
   # k=14
   # data_df=os_data
  temp=as.numeric(as.character(data_df[,i]))
  n_low = k
  n_high = length(temp)-k
  temp[order(temp)[(k+1):length(temp)]]=1
  temp[order(temp)[1:k]]=0
  sur_data=data.frame(time = data_df$month,event = data_df$status,count=temp)
  p <- ggsurvplot(survfit(surv_data~temp),data = sur_data,title=title_sur,
                  palette = c("#00468B", "#F24747"),legend.labs=c("low","high"),
                 risk.table=T)
  return(p)
}

#######plot all sig circRNA
#pdf(file = "/data1/yeying/m6a_circ/pancreatic/plot/11_7/survival_res_11_7.pdf")

for (c in 1:dim(sig_survival_res)[1]){
  gene_name=paste(sig_survival_DE_DM_feature$symbol[c],
                  sig_survival_DE_DM_feature$id[c],sep = ":")
  title_sur=paste("survival analysis of",gene_name)
  col_num=as.numeric(sig_survival_res$col_number[c])
  cutoff=as.numeric(sig_survival_res$n_low[c])
  cat(gene_name,"P value = ",sig_survival_DE_DM_feature$p_min[c],"\n",
               "HR =",sig_survival_DE_DM_feature$HR[c])
   p <- plot_survial(col_num,cutoff,os_data,title_sur=title_sur)+
  print(p)
  
}

#dev.off()
####################other risk factor
str(os_data)

###########multi vocano plot
#colnames(result)=c("Peak.id","Clinic.factor","P.value","log2FC")
clinc_data_circ <- os_data0[,1:28]
rownames(clinc_data_circ) = paste0("S_",clinc_data_circ$t_sample_id)
tumor_meth.mat = countData[,grep("T",colnames(countData))]
clinc_data_circ = clinc_data_circ[colnames(tumor_meth.mat),]
####PNI perineural invasion

risk_factors = c("gender","Differentiation","smoke","drink","PNI","Angiogenesis",
                 "node","TNM_3","Metastasis","Relapse","status") #"Metastasis.Relapse"

clinc_data_circ$TNM_3[which(is.na(clinc_data_circ$TNM_3))] <- 1
all_risk_factor = data.frame()
attach(clinc_data_circ)
i=1
for ( i in 1:length(risk_factors) ) {
  pheno.df <- data.frame(Type= as.factor(get(risk_factors[i])))
  rownames(pheno.df) = rownames(clinc_data_circ)
  MethCirc_edgeR <-edgeR_test(expre_mat = tumor_meth.mat,group_mat = pheno.df)
  factor_res.df = data.frame(Peak.id = rownames(MethCirc_edgeR),
                             Clinic.factor =rep(risk_factors[i],dim(MethCirc_edgeR)[1]),
                             P.value = MethCirc_edgeR$PValue , log2FC = MethCirc_edgeR$logFC)
  all_risk_factor <- rbind(all_risk_factor,factor_res.df)
} 

detach(clinc_data_circ)

peak2gene= filter(circ_feature, id%in%rownames(tumor_meth.mat))
peak2gene = data.frame(Peak.id = peak2gene$id,gene = peak2gene$symbol)

clinic_heatmap_for_common_hyper=function(peak_associate_with_clinic,top=10,peak2gene){
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$P.value<0.05),]
  peak_associate_with_clinic=merge(peak_associate_with_clinic,peak2gene,by="Peak.id",all.x = T)
  peak_associate_with_clinic=peak_associate_with_clinic[which(!is.na(peak_associate_with_clinic$gene)),]
  peak_associate_with_clinic=as.data.frame(cbind(peak_associate_with_clinic,type="non"))
  peak_associate_with_clinic=peak_associate_with_clinic[order(peak_associate_with_clinic$P.value),]
  i=1
  Clinic.factor=unique(peak_associate_with_clinic$Clinic.factor)
  n=length(Clinic.factor)
  while(i<=n){
    if(length(peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'Peak.id'])<top){
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'type']='sig'
    } else {
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i])[1:top],'type']='sig'
    }
    i=i+1
  }
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$type=="sig"),]
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  peak_associate_with_clinic$Peak.id=paste(peak_associate_with_clinic$Peak.id,":",peak_associate_with_clinic$gene,sep="")
  # peak_associate_with_clinic$Peak.id=ordered(peak_associate_with_clinic$Peak.id,
  #                                            levels=names(table(peak_associate_with_clinic$Peak.id)[
  #                                              order(table(peak_associate_with_clinic$Peak.id),decreasing = T)]))
  peak_associate_with_clinic=data.frame(cbind(peak_associate_with_clinic,Difftype="Up regulation"))
  peak_associate_with_clinic[which(peak_associate_with_clinic$log2FC<0),"Difftype"]="Down regulation"
  #peak_associate_with_clinic$Clinic.factor=ordered(peak_associate_with_clinic$Clinic.factor,
  #                                                 levels=rev(c("age","gender","drink","smoke","TNM_3","xueguan","node","nerve",
  #                                                              "fenhua_class : 1  vs Rest","fenhua_class : 2  vs Rest",
  #                                                              "fenhua_class : 3  vs Rest")))
  write.table(sort(peak_associate_with_clinic$Peak.id),"survival/peak_ass.temp.txt",row.names = F,
              quote = F)
  p=ggplot(peak_associate_with_clinic)+
    geom_point(aes(y=peak_associate_with_clinic$Clinic.factor,
                   x=peak_associate_with_clinic$Peak.id,size=(-log10(peak_associate_with_clinic$P.value)),
                   color=as.factor(peak_associate_with_clinic$Difftype)))+
    scale_color_manual(values=c("#5707d5","#d50709"),name="log2(Flod change)")+
    scale_size_continuous(name="-log10(P value)")+
    xlab("Clinic factor")+ylab("Tumor specific peaks top10 in each clinic factor")+
    theme_classic()+theme(axis.text.x = element_text(angle=90))
  return(p)
}

ggsave(p,filename = "survival/all_risk_factor_analysis.pdf",width = 20,height = 8)
clinic_heatmap_for_common_hyper(peak_associate_with_clinic = all_risk_factor,
                                top = 10,peak2gene = peak2gene) ##### peak2gene : peakID genesymbol



