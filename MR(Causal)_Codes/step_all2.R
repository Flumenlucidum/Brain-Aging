#Import libraries
library(MendelianRandomization)
library(nlmr)
library(data.table)
library(dplyr)
library(seqminer)
library(RCIT)

#Get the pvalue less than 5e-8 and make the list of them 
MR_res_table=fread('../MR_result_final.csv') #A table to write results on 
MR_res_table<-as.data.frame(MR_res_table)

for(var in 250:nrow(MR_res_table)){
	code=as.numeric(MR_res_table$code[var])

	file=paste0('instruments_',as.character(code),'.txt')
	
	if(!(file %in% list.files())){
	
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
			a<-a[which(a$pval_EUR<5e-8),]
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

			system('./step2.sh',wait=T)  ##Make plink files with only the significant markers with regard to the exposure
			Sys.sleep(40)
			system('./step3.sh',wait=T)  #LD pruning
			Sys.sleep(40)
			system('./step4.sh',wait=T)  #LD pruning
			Sys.sleep(40)
			system('./step5.sh',wait=T)  #Extract only the 34,129 individuals with MRI
			Sys.sleep(40)
			system('Rscript step6.R',wait=T)  #Merge files from different chromosomes
			system('Rscript step7.R',wait=T)  #Remove markers that directly affect the outcome (delta age)
			
			final=fread('final_markers_list.txt')
			final<-final[which(final$AF_Allele2>=0.01 & final$af_EUR>=0.01),]
			library(dplyr)

			#Nonlinear MR
			try({pl=openPlink('MR_pruned_merged')
			final=final[which(final$beta_EUR>0),]
			dd<-pl[1:length(pl$fam$V1),which(pl$bim$V2 %in% final$POS)]
			dd[dd==-9]<-NA
			for(i in 1:ncol(dd)){
			dd[is.na(dd[,i]),i]<-mean(dd[,i],na.rm=T)}
			
			final=final%>%group_by(POS)%>%slice(which.min(pval_EUR))
			final<-as.data.frame(final)
			
			test2<-final[c('POS','beta_EUR')]
			test1<-as.data.frame(colnames(dd))
			colnames(test1)[1]<-'POS'
			test<-left_join(test1,test2,by='POS')
			score<-rep(0,nrow(dd))

			for(i in 1:nrow(dd)){
				score[i]<-dd[i,]%*%test[,2]
			}
			
			allele_score<-cbind(rownames(dd),score)
			colnames(allele_score)[1]<-'IID'
			allele_score<-as.data.frame(allele_score)
			allele_score$IID<-as.numeric(allele_score$IID)
			allele_score$score<-as.numeric(allele_score$score)
			
			meta<-fread('/media/leelabsg-storage0/jangho/MR/blood_for_MR.txt')
			meta<-as.data.frame(meta)
			meta<-meta[1:nrow(meta),c(1,which(colnames(meta) %in% paste0('f.',as.character(code),'.0.0')),64)]

			nlmr<-left_join(meta,allele_score,by='IID')
			colnames(nlmr)<-c('IID','X','Y','Z')
			write.table(nlmr,file=paste0('test_',as.character(code),'exp.txt'),row.names=F,quote=F)
			cor_XZ=cor.test(nlmr$Z,nlmr$X)
			
			print('cor_XZ')
			print(cor_XZ$p.value<0.05)
			print(cor_XZ$estimate>0)

			if(cor_XZ$p.value<0.05 & cor_XZ$estimate>0){
				d1=RCIT(nlmr$Z,nlmr$Y,nlmr$X)  #Repeat RCIT three times
				d2=RCIT(nlmr$Z,nlmr$Y,nlmr$X)
				d3=RCIT(nlmr$Z,nlmr$Y,nlmr$X)

				if(d1$p>=0.05 & d2$p>=0.05 &d3$p>=0.05){
					mod=piecewise_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)

					print('nonlinear')
					MR_res_table[var,52]<-mod$p_tests[1]
					MR_res_table[var,53]<-mod$p_tests[2]
					MR_res_table[var,54]<-mod$p_heterogeneity[1]
					MR_res_table[var,55]<-mod$p_heterogeneity[2]
					print('fracpoly')
					if(min(nlmr$X)>1){
						mod2=fracpoly_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)
					}else{
						add=1.00001-min(nlmr$X)
						nlmr$X<-add+nlmr$X
						mod2=fracpoly_mr(y=nlmr$Y,x=nlmr$X,g=nlmr$Z)
					}

					MR_res_table[var,56]<-mod2$p_tests[1]
					MR_res_table[var,57]<-mod2$p_tests[2]
					MR_res_table[var,58]<-mod2$p_tests[3]
					MR_res_table[var,59]<-mod2$p_tests[4]
					MR_res_table[var,60]<-mod2$p_heterogeneity[1]
					MR_res_table[var,61]<-mod2$p_heterogeneity[2]
				}
			}

			})
		

		system('rm MR*',wait=T)
		system('rm *.gz',wait=T)
		print('done')
		print(var)
		write.csv(MR_res_table,file='MR_result_current_blood_NONLINEAR_FIN.csv',row.names=F)
		}	

	}
}
