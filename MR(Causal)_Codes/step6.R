system('rm *~', wait=T)
a=list.files()
b=a[which(grepl('MR_pruned',a) & grepl('.bed',a))]
b_not=a[which(grepl('MR_not_pruned',a) & grepl('.bed',a))]
c=as.numeric(substr(b,14,15))
print(c)
c_not=as.numeric(substr(b_not,18,19))
print(c_not)
c_not=c_not[which(!c_not %in% c)]


for(i in c_not){
	system(paste0('mv MR_not_pruned_chr',as.character(i),'.bed', ' MR_pruned_chr',as.character(i),'.bed'))
	system(paste0('mv MR_not_pruned_chr',as.character(i),'.bim', ' MR_pruned_chr',as.character(i),'.bim'))
	system(paste0('mv MR_not_pruned_chr',as.character(i),'.fam', ' MR_pruned_chr',as.character(i),'.fam'))
}
a=list.files()
b=a[which(grepl('MR_pruned',a) & grepl('.bed',a))]
c=as.numeric(substr(b,14,15))
min_num=min(c)
c=c[which(c!=min_num)]
c=c[order(c)]
string=''

for(i in c){
string=paste0(string,paste0('MR_pruned_chr',as.character(i),'.bed'),paste0(' MR_pruned_chr',as.character(i),'.bim'), paste0(' MR_pruned_chr',as.character(i),'.fam','\n'))}

write.table(string,file='mergelist.txt',row.names=F,quote=F,col.names=F)
system(paste0('plink --bfile MR_pruned_chr',as.character(min_num),' --merge-list mergelist.txt --make-bed --out MR_pruned_merged'))
