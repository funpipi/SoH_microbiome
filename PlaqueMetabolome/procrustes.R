p <- c("vegan","ade4","reshape2","ggplot2","pheatmap","dplyr","randomForest", 
       "cluster", "clusterSim", "plotly", "corrplot", "pROC")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

#-------------------------------
setwd('~/Documents/SoH_project/PlaqueMetabolome')

filename1<-"PlaqueMS275.txt"
prefix1<-"metabolomics"
metadata.filename1<-"PlaqueMS275_Map.txt"
filename2="SoH261.taxa.genus.Abd.xls"
prefix2<-"microbiome-genus"
metadata.filename2<-"SoH261_Map.txt"
outputPath<-"./Procustes_analysis/"
dir.create(outputPath)
f1<-read.table(filename1,header=TRUE,sep="\t",row.names=1, quote=""); f1<-f1[order(rownames(f1)), ]
f2<-read.table(filename2,header=TRUE,sep="\t",row.names=1); f2<-f2[order(rownames(f2)), ]
m1<-read.table(metadata.filename1,header=T,sep="\t",row.names=1); m1<-m1[order(rownames(m1)), ]
m2<-read.table(metadata.filename2,header=T,sep="\t",row.names=1); m2<-m2[order(rownames(m2)), ]
cat("The order of samples in the metabolome data and metadata are equal:", all.equal(rownames(f1), rownames(m1)))
cat("The order of samples in the microbiome data and metadata are equal:", all.equal(rownames(f2), rownames(m2)))

#--------------------------------------------------
# To identified intersected sample IDs
#--------------------------------------------------
Shared_SampleIDs<-intersect(rownames(f1), rownames(f2))
f1_s<-f1[Shared_SampleIDs,]; 
m1_s<-m1[Shared_SampleIDs,]; 
cat("The order of samples in the metabolome data and metadata are equal:", all.equal(rownames(f1_s), rownames(m1_s)))
f2_s<-f2[Shared_SampleIDs,]; 
m2_s<-m2[Shared_SampleIDs,]; 
cat("The order of samples in the microbiome data and metadata are equal:", all.equal(rownames(f2_s), rownames(m2_s)))
library("graphics")
dt<-table(m2[, c("Timepoint", "Group", "Gender", "Smoking")])
mosaicplot(dt, shade = TRUE, las=2, main = "SoH clinical study")
#--------------------------------------------------
# Distance matrix
#--------------------------------------------------
dm_f1<-vegdist(f1_s, method = "bray")
dm_f2<-vegdist(f2_s, method = "bray")
#--------------------------------------------------
#PCA using ade4 package
#--------------------------------------------------
dm_f1.pco <- dudi.pco(dm_f1,scan=FALSE, nf=4)
summary(dm_f1.pco)
dm_f1.Axis1<-dm_f1.pco$li[,1]
dm_f1.Axis2<-dm_f1.pco$li[,2]
#--------------------------------------------------
dm_f2.pco <- dudi.pco(dm_f2,scan=FALSE, nf=4)
dm_f2<-dm_f2.pco$li[,1]
dm_f2<-dm_f2.pco$li[,2]
###################################################
# procuste analysis
###################################################
pro <- procuste(dm_f1.pco$tab, dm_f2.pco$tab, nf = 4)
pro.test<-procuste.randtest(dm_f1.pco$tab, dm_f2.pco$tab, 1000)
pro.test
pdf(paste(outputPath, prefix1,"-",prefix2,".procruste.all.pdf",sep=""),15,15)
plot(pro,clab=0.1,cpoint=2, add.p=TRUE)
dev.off()

X_ax<-data.frame(SampleID=rownames(pro$scorX), m1_s, prefix=rep(prefix1, nrow(m1_s)), ax1=pro$scorX$ax1, ax2=pro$scorX$ax2) 
Y_ax<-data.frame(SampleID=rownames(pro$scorY), m1_s, prefix=rep(prefix2, nrow(m1_s)), ax1=pro$scorY$ax1, ax2=pro$scorY$ax2)
ax<-rbind(X_ax, Y_ax)
#--------------------------------------------------
p<-ggplot(ax, aes(x=ax1, y=ax2, color=Timepoint, shape=prefix)) + 
  geom_point(aes(x=ax1, y=ax2, group=SampleID), alpha=0.4)+
  geom_path(aes(x=ax1, y=ax2, group=SampleID), alpha=0.4)+
  annotate("text", label=paste("Observed_similirity=",formatC(pro.test$obs,digits=4,format="f"),sep=""), x=max(ax$ax1)*0.7, y=max(ax$ax2)*0.95)+
  annotate("text", label=paste("Monte-Carlo test: P=",formatC(pro.test$pvalue,digits=4,format="f"),sep=""), x=max(ax$ax1)*0.7, y=max(ax$ax2)*0.8)+
  xlab("Axis 1")+ylab("Axis 2")+
  theme_bw()
ggsave(filename=paste(outputPath,"/procrustes-Timepoint.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=6, height=5)
ggplotly(p)

#--------------------------------------------------
p<-ggplot(ax, aes(x=ax1, y=ax2, color=Timepoint, shape=prefix)) + 
  geom_point(aes(x=ax1, y=ax2, group=SampleID), alpha=0.4)+
  geom_path(aes(x=ax1, y=ax2, group=SampleID), alpha=0.4)+ 
  xlab("Axis 1")+ylab("Axis 2")+
  facet_wrap(~ Group)+
  theme_bw()
ggsave(filename=paste(outputPath,"/procrustes-Group-Timepoint.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=6, height=5)
ggplotly(p)
#--------------------------------------------------
len_XY<-rowSums((ax[which(ax$prefix=="metabolomics"), c("ax1", "ax2")] - ax[which(ax$prefix=="microbiome-genus"), c("ax1", "ax2")])^2)
len_XY_m<-data.frame(m1_s, len_XY)

#--------------------------------------------------
p<-ggplot(len_XY_m, aes(x=Timepoint, y=len_XY)) + geom_boxplot() + 
  geom_point(aes(x=Timepoint, y=len_XY, group=Host_ID), alpha=0.4)+
  geom_path(aes(x=Timepoint, y=len_XY, group=Host_ID, color=as.factor(Host_ID)), alpha=0.4)+
  ylab("Dissimilarity between microbiome and metabolomics")+
  facet_wrap(~ Group)+
  theme_bw()
ggsave(filename=paste(outputPath,"/procrustes-dissimilarity-Group-Timepoint.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=9, height=5)
p
ggplotly(p)
#--------------------------------------------------

p<-ggplot(len_XY_m, aes(x=len_XY, y=Bleeding, colour=Group)) + geom_point() + 
   facet_wrap(~ Timepoint, nrow=1) + theme_bw() + 
   xlab("Dissimilarity between microbiome and metabolomics")
ggsave(filename=paste(outputPath,"/len_XY-Bleeding.scatterplot.pdf",sep=""), plot=p, limitsize=TRUE, width=6, height=5)
p
cor.test(len_XY_m$len_XY, len_XY_m$Bleeding)
len_XY_m_DAY0DAY28<-subset(len_XY_m, Timepoint=="DAY0" | Timepoint=="DAY28")
with(len_XY_m_DAY0DAY28, wilcox.test(len_XY~Timepoint))


#--------------------------------------------------
#   Correlation analysis of microbiome and metabolome
#--------------------------------------------------
source("util.R")
#--------------------------------------------------
#   Differential abundant Metabolites from Day0 to Day28 -- unpaired Wilcoxon test
#--------------------------------------------------
m1_EG<-m1[which(m1$Timepoint=="DAY0" | m1$Timepoint=="DAY28"), ]
f1_EG<-f1[which(m1$Timepoint=="DAY0" | m1$Timepoint=="DAY28"), ]
f1_EG<-CleanData(f1_EG)
unpaired_test_res<-BetweenGroup.test(f1_EG, factor(m1_EG[, "Timepoint"]), p.adj.method = "fdr")
subset(unpaired_test_res, Wilcoxon.test_fdr <0.05)

#--------------------------------------------------
#   Differential abundant Metabolites from Day0 to Day28 -- paired Wilcoxon test
#--------------------------------------------------
m1_EG_paired<-m1_EG[which(m1_EG$Host_ID!="H8003"),]
f1_EG_paired<-f1_EG[which(m1_EG$Host_ID!="H8003"),]
paired_test_res<-BetweenGroup.test(f1_EG_paired, factor(m1_EG_paired[, "Timepoint"]), p.adj.method = "fdr")
subset(paired_test_res, Wilcoxon.test_fdr <0.05)
#--------------------------------------------------
#  Heatmap of the differential abundant Metabolites from Day0 to Day28 
#--------------------------------------------------
f1_EG_sig<-f1_EG[, rownames(subset(paired_test_res, Wilcoxon.test_fdr <0.05))]
f1_EG_sig<-f1_EG_sig[order(m1_EG$Timepoint), ]
f1_EG_sig_log<-log.mat(f1_EG_sig)
pheatmap(t(f1_EG_sig_log), cellwidth = 2, cellheight = 1, fontsize_row = 1, frontsize_col=2, cluster_rows=T, cluster_cols=F, filename = paste(outputPath,"/EG.metabo.sig.pheatmap_t.pdf",sep=""))

pheatmap(f1_EG_sig_log, cellwidth =3, cellheight = 2, fontsize_row = 2, frontsize_col=0.5, cluster_rows=F, cluster_cols=T, filename = paste(outputPath,"/EG.metabo.sig.pheatmap.pdf",sep=""))

f1_s_sig<-f1_s[, rownames(subset(paired_test_res, Wilcoxon.test_fdr <0.05))]

#--------------------------------------------------
#   Identified microbiome markers
#--------------------------------------------------
MiG_list<-c('Lep_Leptotrichia.','Fus_Fusobacterium.','Por_Porphyromonas.','Por_Tannerella.','Cam_Campylobacter.','Pas_Aggregatibacter.','Lac_Johnsonella.','Pre_Prevotella.','Vei_Selenomonas.','Lac_Dorea.','Pep_Peptoclostridium.','Lac_Catonella.','TM7_TM7.','Gem_Gemella.','Lac_Lachnoanaerobaculum.','Com_Comamonas.','Act_Actinomycetaceae_Group.','Par_Prevotella.','Act_Actinobaculum.','Spi_Treponema.','Com_Hylemonella.','Aer_Abiotrophia.','Por_Paludibacter.','Lac_Lachnospiraceae_Group.','Rum_Ruminococcaceae_Group.','SR1_SR1.','Det_TG5.','Cor_Atopobium.', 'Mic_Rothia.','Str_Streptococcus.','Bur_Lautropia.','Pas_Haemophilus.','Act_Actinomyces.','Str_Streptococcaceae_Group.','Wil_Williamsia.')
MiG_list_s<-intersect(MiG_list, colnames(f2))
f2_s_sig<-f2_s[, MiG_list_s]
#--------------------------------------------------
#   Corrlation between microbiome and metabolomics markers
#--------------------------------------------------

cor_mat<-cor(f1_s_sig, f2_s_sig, method="spearman")
pheatmap(cor_mat, cellwidth =4, cellheight = 2, fontsize_row = 2, frontsize_col=0.5, cluster_rows=T, cluster_cols=T, filename = paste(outputPath,"/metabo-microbio.corr.pheatmap.pdf",sep=""))
pheatmap(t(cor_mat),  cellwidth =1, cellheight = 4, fontsize_row = 2, frontsize_col=0.5, cluster_rows=T, cluster_cols=T, filename = paste(outputPath,"/metabo-microbio.corr.pheatmap_t.pdf",sep=""))


pdf(paste(outputPath,"/",prefix1,"-",prefix2,"_corrplot.pdf",sep=""))
corrplot.mixed(cor_mat)
dev.off()




#--------------------------------------------------
#   host outilers
#--------------------------------------------------
p<-ggplot(m1[order(m1$Timepoint),], aes(x=Timepoint, y=Bleeding)) + 
  geom_boxplot()+
  #geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,aes(fill=Group),alpha=0.4) +
  geom_point(aes(x=Timepoint, y=Bleeding, group=Host_ID), alpha=0.4)+
  geom_path(aes(x=Timepoint, y=Bleeding, group=Host_ID, color=Host_ID), alpha=0.4)+ 
  facet_wrap(~ Group)+
  theme_bw()
ggsave(filename=paste(outputPath,"/Bleeding.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=8, height=5)
ggplotly(p)
#--------------------------------------------------
Host_IDs_3<-subset(m1, Timepoint=="DAY28" & Bleeding<3)$Host_ID
m1_Host_IDs_3<-m1[which(m1$Host_ID %in% Host_IDs_3),]
p<-ggplot(m1_Host_IDs_3[order(m1_Host_IDs_3$Timepoint),], aes(x=Timepoint, y=Bleeding)) + 
  geom_boxplot()+
  #geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,aes(fill=Group),alpha=0.4) +
  geom_point(aes(x=Timepoint, y=Bleeding, group=Host_ID), alpha=0.4)+
  geom_path(aes(x=Timepoint, y=Bleeding, group=Host_ID, color=Host_ID), alpha=0.4)+ 
  facet_wrap(~ Group)+
  theme_bw()
ggsave(filename=paste(outputPath,"/m1_Host_IDs_Bleeding3_atDay28.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=8, height=5)
ggplotly(p)
#--------------------------------------------------
head(m1)
Host_IDs_45<-subset(m1, Timepoint=="DAY28" & Bleeding>3 & Bleeding<5)$Host_ID
m1_Host_IDs_45<-m1[which(m1$Host_ID %in% Host_IDs_45),]
p<-ggplot(m1_Host_IDs_45[order(m1_Host_IDs_45$Timepoint),], aes(x=Timepoint, y=Bleeding)) + 
  geom_boxplot()+
  #geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,aes(fill=Group),alpha=0.4) +
  geom_point(aes(x=Timepoint, y=Bleeding, group=Host_ID), alpha=0.4)+
  geom_path(aes(x=Timepoint, y=Bleeding, group=Host_ID, color=Host_ID), alpha=0.4)+ 
  facet_wrap(~ Group)+
  theme_bw()
ggsave(filename=paste(outputPath,"/m1_Host_IDs_Bleeding45_atDay28.boxplot.pdf",sep=""), plot=p, limitsize=TRUE, width=8, height=5)
ggplotly(p)

#--------------------------------------------------
pdf(paste(outputPath, prefix1,"-",prefix2,".label.procruste.pdf",sep=""),12,12)
par(mar=c(6,6,6,6))
plot(pro$scorX$ax1, pro$scorX$ax2,pch=21, cex=2, col="white", bg="green", xlab="Axis1", ylab="Axis2",cex.axis=2,cex.lab=1,
     xlim=c((min(pro$scorX$ax1)-0.1*(max(pro$scorX$ax1)-min(pro$scorX$ax1))),(max(pro$scorX$ax1)+0.2*(max(pro$scorX$ax1)-min(pro$scorX$ax1)))), 
     ylim=c((min(pro$scorX$ax2)-0.1*(max(pro$scorX$ax2)-min(pro$scorX$ax2))),(max(pro$scorX$ax2)+0.2*(max(pro$scorX$ax2)-min(pro$scorX$ax2)))))
points(pro$scorY$ax1, pro$scorY$ax2,pch=21, cex=2, col="white", bg="red",cex.axis=2,cex.lab=1,
       xlim=c((min(pro$scorY$ax1)-0.1*(max(pro$scorY$ax1)-min(pro$scorY$ax1))),(max(pro$scorY$ax1)+0.2*(max(pro$scorY$ax1)-min(pro$scorY$ax1)))), 
       ylim=c((min(pro$scorY$ax2)-0.1*(max(pro$scorY$ax2)-min(pro$scorY$ax2))),(max(pro$scorY$ax2)+0.2*(max(pro$scorY$ax2)-min(pro$scorY$ax2)))))
segments(x0=pro$scorX$ax1, y0=pro$scorX$ax2,x1=pro$scorY$ax1, y1=pro$scorY$ax2,ljoin=2,lend=2,col="#0000ff22",lwd=3 )
x.label<-0.5*(pro$scorY$ax1+pro$scorX$ax1)
y.label<-0.5*(pro$scorY$ax2+pro$scorX$ax2)
scatterutil.eti.circ(x.label, y.label, rownames(pro$scorX),1)
legend("topleft",c(filename1,filename2), pch=20, col=c("green","red"), pt.cex=2,cex=1)
mtext(paste("Observed_similirity=",formatC(pro.test$obs,digits=4,format="f"),sep=""),outer=FALSE,line=-2,adj=0.9)
mtext(paste("Monte-Carlo test: P_value=",formatC(pro.test$pvalue,digits=4,format="f"),sep=""),outer=FALSE,line=-3,adj=0.9)     
dev.off()


#--------------------------------------------------

