library(data.table)
library(dplyr)
bim=fread('MR_pruned_merged.bim')

start=0
for(i in unique(bim$V1)){
	bim_a<-bim[which(bim$V1==i),]
	if(dim(bim_a)[1]!=0){
		if(start==0){
			a<-fread(paste0('/media/leelabsg-storage0/jangho/ukb_imp/ukb_whole_saige_step2_result_LOCO_chr',as.character(i),'.txt'))
			a<-a[which(a$CHR %in% bim_a$V4),]
			start=1
		}else{
			add<-fread(paste0('/media/leelabsg-storage0/jangho/ukb_imp/ukb_whole_saige_step2_result_LOCO_chr',as.character(i),'.txt'))
			add<-add[which(add$CHR %in% bim_a$V4),]
			a<-rbind(a,add)
		}
		print(i)
	}
}

a=a[which(a$p.value>=(0.05/length(unique(a$POS)))),]
tab<-fread('table.txt')
final=left_join(a,tab,by=c('V1'='chr','CHR'='pos'))
write.table(final,file='final_markers_list.txt',row.names=F,quote=F)
