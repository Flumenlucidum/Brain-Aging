#Import Libraries
library(SKAT)
library(data.table)

SetID<-fread('./ukb_whole_chr1.SetID',head=F)
gene_list<-unique(SetID$V1)
n=length(gene_list)%/%100
rem=length(gene_list)%%100

File.Bed<-'./ukb_whole_wes_chr1.bed'
File.Bim<-'./ukb_whole_wes_chr1.bim'
File.Fam<-'./ukb_whole_wes_chr1.fam'
File.Cov<-'./ukb_whole_wes.cov'
FAM_Cov<-Read_Plink_FAM_Cov(File.Fam,File.Cov,Is.binary = F)

#Object file for Null model
obj<-SKAT_Null_Model(Phenotype~COVAR1+COVAR2+COVAR3+COVAR4+COVAR5+COVAR6+COVAR7+COVAR8+COVAR9+COVAR10+COVAR11+COVAR12+cv1+cv2+cv3+cv4, data=FAM_Cov, out_type = 'C')

# If the phenotype is binary one, out_type='D'
#When there is no covariate file, use FAM<-Read_Plink_FAM(File.Fam,Is.binary = F) 
#and object obj<-SKAT_Null_Model(y~1, out_type = 'C')

#Load 100 genes at a time to prevent problems with the memory
for(i in 1:n){

	SetID_part<-SetID[which(SetID$V1 %in% gene_list[(100*(i-1)+1):(100*i)]),]
	write.table(SetID_part,file='ukb_whole_chr1_part.SetID',row.names=F,col.names=F,quote=F)
	File.SetID<-'./ukb_whole_chr1_part.SetID'

#Please set file location and name for SSD file and SSD.info file 
File.SSD<-'./Example1.SSD'
File.Info<-'./Example1.SSD.info'

#Generate and open SSD file for analysis
Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.SetID,File.SSD,File.Info )
SSD.INFO<-Open_SSD(File.SSD,File.Info)
print(SSD.INFO$nSample)
print(SSD.INFO$nSets)

#Analysis
if(i==1){
	out<-SKAT_CommonRare.SSD.All(SSD.INFO,obj)
	out<-out$results
}else{
	out_add<-SKAT_CommonRare.SSD.All(SSD.INFO,obj)
	out_add<-out_add$results
	out<-rbind(out,out_add)
}
}

if(rem!=0){

SetID_part<-SetID[which(SetID$V1 %in% gene_list[(n*100+1):length(gene_list)]),]
	write.table(SetID_part,file='ukb_whole_chr1_part.SetID',row.names=F,col.names=F,quote=F)
	File.SetID<-'./ukb_whole_chr1_part.SetID'

#Please set file location and name for SSD file and SSD.info file 
File.SSD<-'./Example1.SSD'
File.Info<-'./Example1.SSD.info'

#Generate and open SSD file for analysis
Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.SetID,File.SSD,File.Info )
SSD.INFO<-Open_SSD(File.SSD,File.Info)
print(SSD.INFO$nSample)
print(SSD.INFO$nSets)

#Analysis

	out_add<-SKAT_CommonRare.SSD.All(SSD.INFO,obj)
	out_add<-out_add$results
	out<-rbind(out,out_add)
}



write.table(out,file='ukb_whole_wes_result_commonrare_chr1.txt',row.names=F,quote=F)
#close SSD file 
Close_SSD()

