#Import libraries
library(MendelianRandomization)
library(nlmr)
library(data.table)
library(dplyr)
library(seqminer)

#Get the pvalue less than 5e-8 and make the list of them 
MR_res_table=fread('../MR_results_table.csv') #A table to write results on
MR_res_table<-as.data.frame(MR_res_table)

for(var in 303:nrow(MR_res_table)){
	if(MR_res_table$code[var]>=30500){
	system(paste0('wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-',as.character(MR_res_table$code[var]),'-both_sexes-irnt.tsv.bgz'),wait=T)
	system(paste0('wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-',as.character(MR_res_table$code[var]),'-both_sexes.tsv.bgz'),wait=T)
	}else{
system(paste0('wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-',as.character(MR_res_table$code[var]),'-both_sexes-irnt.tsv.bgz'),wait=T)
system(paste0('wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-',as.character(MR_res_table$code[var]),'-both_sexes.tsv.bgz'),wait=T)
	}
	chr=list.files()
	chr<-chr[grepl('.bgz',chr)][1]
	system(paste0('mv ',chr,' cont.gz'),wait=T)
	a<-fread('cont.gz')
	a<-as.data.frame(a)
	if('beta_EUR' %in% colnames(a)){
	a<-a[,c('chr','pos','ref','alt','af_EUR','beta_EUR','se_EUR','pval_EUR')] #Estimates from the European ancestry
	a<-a[complete.cases(a),]
	a<-a[which(a$pval_EUR<5e-8),]   #Significant markers with P<5e-8
	a<-a[which(a$chr!='X'),]

	a$chr<-as.numeric(a$chr)
	write.table(a,file='table.txt',row.names=F,quote=F)
	a<-fread('table.txt')

	bim<-fread('/media/leelabsg-storage0/DATA_1/UKBB/cal/ukb_snp_chr1_v2.bim')
	for(i in 2:22){
	bim_add<-fread(paste0('/media/leelabsg-storage0/DATA_1/UKBB/cal/ukb_snp_chr',as.character(i),'_v2.bim'))
	bim<-rbind(bim,bim_add)
	print(i)
	}
	

	fin<-left_join(a,bim,by=c('chr'='V1','pos'='V4'))
	fin2<-fin[complete.cases(fin),]
	write.table(fin2$V2,file='markers_list.txt',row.names=F,quote=F,col.names=F)

	system('./step2.sh',wait=T)  #Make plink files with only the significant markers with regard to the exposure
	Sys.sleep(40)
	system('./step3.sh',wait=T)  #LD pruning
	Sys.sleep(40)
	system('./step4.sh',wait=T)  #LD pruning
	Sys.sleep(40)
	system('./step5.sh',wait=T)  #Extract only the 34,129 individuals with MRI
	Sys.sleep(40)
	system('Rscript step6.R',wait=T) #Merge files from different chromosomes
	system('Rscript step7.R',wait=T) #Remove markers that directly affect the outcome (delta age)
	
	final=fread('final_markers_list.txt')
		
	#Filter out rare variants for stable estimation of causal estimates
	final<-final[which(final$AF_Allele2>=0.01 & final$af_EUR>=0.01),]
	code=as.numeric(MR_res_table$code[var])
	if(dim(final)[1]>=2){
	
	#Linear MR
	obj<-mr_input(bx=final$beta_EUR,bxse=final$se_EUR, by=final$BETA,byse=final$SE)
	all=mr_allmethods(obj)
	print('Linear')
	print(all@Values[,6]<0.05)
	write.table(all@Values, file='Linear_MR_result.txt',row.names=F,quote=F)
	jpeg('mr_plot.jpg',width=1000,height=1000)
	mr_plot(all)
	dev.off()


	#Nonlinear MR
	pl=openPlink('MR_pruned_merged')
	final=final[which(final$beta_EUR>0),]
	dd<-pl[1:length(pl$fam$V1),which(pl$bim$V2 %in% final$POS)]
	dd[dd==-9]<-NA
	for(i in 1:ncol(dd)){
	dd[is.na(dd[,i]),i]<-mean(dd[,i],na.rm=T)}
	allele_score=rowSums(dd)
	allele_score=as.data.frame(allele_score)
	allele_score$IID<-rownames(allele_score)
	allele_score$IID<-as.numeric(allele_score$IID)
	meta<-fread('/media/leelabsg-storage0/jangho/MR/blood_for_MR.txt')
	meta<-as.data.frame(meta)
	meta<-meta[1:nrow(meta),c(1,which(colnames(meta) %in% paste0('f.',as.character(code),'.0.0')),64)]

	nlmr<-left_join(meta,allele_score,by='IID')
	colnames(nlmr)<-c('IID','X','Y','Z')


	try({mod=piecewise_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)

	print('nonlinear')
	print(mod$p_tests<0.05)
	print(mod$p_heterogeneity<0.05)
	
	print('fracpoly')  
	if(min(nlmr$X)>1){
	mod2=fracpoly_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)
	}else{
		add=1.00001-min(nlmr$X)
		nlmr$X<-add+nlmr$X
		mod2=fracpoly_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)
	}
	print(mod2$p_tests<0.05)
	print(mod2$p_heterogeneity<0.05)
	if(any(all@Values[,6]<0.05)|any(mod$p_tests<0.05)|any(mod$p_heterogeneity<0.05)|any(mod2$p_tests<0.05)|any(mod2$p_heterogeneity<0.05)){
		
		MR_res_table$linear[var]<-any(all@Values[,6]<0.05)
		MR_res_table$nonlinear[var]<-any(mod$p_tests<0.05)|any(mod$p_heterogeneity<0.05)|any(mod2$p_tests<0.05)|any(mod2$p_heterogeneity<0.05)
		system(paste0('mkdir res_', as.character(MR_res_table$code[var])),wait=T)
		#weird 
	       	system(paste0('mv Linear_MR_result.txt	./res_',as.character(MR_res_table$code[var]),'/'),wait=T)
		system(paste0('mv mr_plot.jpg  ./res_',as.character(MR_res_table$code[var]),'/'),wait=T)
	}
	})
	}
	system('rm MR*',wait=T)
	system('rm *.gz',wait=T)
	print('done')
	print(var)
	}
}
write.table(MR_res_table,file='../DONE_blood.txt',row.names=F,quote=F)


