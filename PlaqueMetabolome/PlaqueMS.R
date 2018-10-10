setwd("~/Documents/SoH_project/PlaqueMetabolome")
data<-read.table("PlaqueMS275_comb.txt",header=T,sep="\t")
#'Sample balance check: to find out which subject missing any timepoint
tab<-table(data[, c(2, 4)])
subj_missing_tp<-names(which(apply(tab, 1, function(x) any(x==0))))
#'To keep subjects only when they have balanced sampling at different timepoints

data_blnc<-data[which(!is.element(data$Subject.ID, subj_missing_tp)),]
write.table(data_blnc, "PlaqueMS252_comb.txt", quote=F, sep="\t", row.names=F)
