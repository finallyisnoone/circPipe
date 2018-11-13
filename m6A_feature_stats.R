
####distribution
library(readr)
library(ggplot2)
source("/data1/yeying/m6a_circ/script/R_function.R")
#circ_feature= read_delim("matrix/uniq_circ_m6Astatus_anno.txt","\t")
circ_feature= read_delim("matrix/11_7/all_filter_circ_m6Astatus_anno.txt","\t")
#circ_name = paste(circ_feature$chr,circ_feature$start,circ_feature$end,circ_feature$strand,sep = "_")
circ_feature$position = circ_feature$id

circ_feature_m6A= subset(circ_feature,m6Astatus=="m6A")
circ_feature_nonm6A= subset(circ_feature,m6Astatus=="non-m6A")
########feature: gene type 
plot_type_pie <- function(feature_mat) {
  df=as.data.frame(table(feature_mat$type1))
  names(df)[1] <- "type"
  df$Freq=as.numeric(as.character(df$Freq))
  p=pie_plot(df)
  return(p)
}

#pdf(file = "plot/11_7/genomic distribution of circRNAs.pdf", width=10.31,height=6.23)
m6A_typepie=plot_type_pie(circ_feature_m6A)
m6A_typepie+ggtitle("Genomic distribution of m6A circRNAs")
nonm6A_typepie=plot_type_pie(circ_feature_nonm6A)
nonm6A_typepie+ggtitle("Genomic distribution of nonm6A circRNAs")
all_circ_typepie=plot_type_pie(as.data.frame(circ_feature))
all_circ_typepie+ggtitle("Genomic distribution of all circRNAs")
#dev.off()

tmp <- xtabs(~m6Astatus+type1,circ_feature)
tmp <- chisq.test(tmp)
tmp$p.value
#######exon number
# circ_feature_m6A= m6A_circ
# circ_feature_nonm6A= non_m6A_circ

within7_freq <- function(circ_feature_df){
  circ_feature_stats <- as.data.frame(table(filter(circ_feature_df,circ_type != "ciRNA")$exon_number))
  names(circ_feature_stats)[1] <- "exon_number"
  #circ_feature_stats$exon_number <- as.character(circ_feature_stats$exon_number )
  circ_feature_stats$Freq <- as.numeric(as.character(circ_feature_stats$Freq))
  circ_feature_stats$exon_number <- as.character(circ_feature_stats$exon_number )
  circ_stats_over7 <- c(paste0("over",circ_feature_stats[7,1]),sum(circ_feature_stats$Freq[8:dim(circ_feature_stats)[1]]))
  circ_feature_stats <- rbind(circ_feature_stats[1:7,],circ_stats_over7)
  circ_feature_stats$Freq <- as.numeric(as.character(circ_feature_stats$Freq))
  return(circ_feature_stats)
}

circ_feature_m6A_stats = within7_freq(circ_feature_m6A)
circ_feature_nonm6A_stats = within7_freq(circ_feature_nonm6A)
circ_feature_m6A_stats$type <- rep("m6A",dim(circ_feature_m6A_stats)[1])
circ_feature_nonm6A_stats$type <- rep("nonm6A",dim(circ_feature_nonm6A_stats)[1])
###########
all_circ_exon <- rbind(circ_feature_m6A_stats,circ_feature_nonm6A_stats)
all_circ_exon$Freq <- as.numeric(as.character(all_circ_exon$Freq))
#####percent
all_circ_exon$Percentage <- c(circ_feature_m6A_stats$Freq/sum(circ_feature_m6A_stats$Freq),
                              circ_feature_nonm6A_stats$Freq/sum(circ_feature_nonm6A_stats$Freq))

all_circ_exon2 = filter(all_circ_exon, exon_number != "over7")
#pdf(file = "plot/11_7/Percentage of exons circRNAs spanning.pdf")
ggplot(all_circ_exon2, aes(x = exon_number, y = Freq,fill=type)) + scale_fill_jama()+theme_classic()+
   theme(legend.position=c(0.9,0.85))+scale_y_continuous(expand = c(0,0))+#scale_x_continuous(expand = c(0,0))+
  geom_bar(stat = "identity",position = "dodge")+ggtitle("the number of exons circRNAs spanning ")

ggplot(all_circ_exon2, aes(x = exon_number, y = Percentage,fill=type)) + scale_fill_jama()+theme_classic()+
  theme(legend.position=c(0.9,0.85))+scale_y_continuous(expand = c(0,0))+
  geom_bar(stat = "identity",position = "dodge")+ggtitle("the number of exons circRNAs spanning ")
#dev.off()


#### the length of circRNA

circ_feature$length = circ_feature$end-circ_feature$start
sum(circ_feature$length > 50000)
len_m6Astatus2 = filter(circ_feature,length > 100 & length < 10000 )

summary(len_m6Astatus2$length)
boxplot(len_m6Astatus2$length)
len_df= as.data.frame(len_m6Astatus2$length)
plot(ecdf(len_m6Astatus2$length),do.point=F,verticals=T)

ggplot()+stat_ecdf(data = len_m6Astatus2,aes(x = length), 
                   geom = "step", size = 1) +scale_color_lancet()+theme_bw()+
  ggtitle("the distrubution of circRNA spanning distance")+ylab("Frequency")+xlab("Length")

############


##########peaks distribution


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

peak_freq=c(as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="5UTR"),14])),
            as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="CDS"),14]))+100,
            as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="3UTR"),14]))+200)

#plot(density(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A peak density",xlab=NA,main=NA,col="darkblue",xaxt="n")
plot(table(peak_freq)/length(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A coding peak density",xlab=NA,main=NA,col="darkblue",xaxt="n",type="l") #xingyang
abline(v=100,lty=2,col="black",lwd=1)
abline(v=200,lty=2,col="black",lwd=1)
text(x=c(50,150,250),y=-0.001,pos=2,labels=c("5' UTR","CDS","3' UTR"),xpd=TRUE,font=2)





