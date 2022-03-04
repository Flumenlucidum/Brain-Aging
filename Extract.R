library(oro.nifti)
library(neurobase)
library(OpenImageR)
library(dplyr)
library(reticulate)
np=import('numpy')
files <- read.table('Aging_WB_nokin.txt')
#z<-c()
#for (i in 12:31){
	#  z[i-11]<-4*i+3
#}
library(unix)
lim <- rlimit_as(1e20)
print(lim)
rlimit_as(cur = lim$max)
z=c(26:153)
k=1
ID<-rep(0,35914)
#35914
final_array<-array(0,c(128,128,128,200))
print(dim(files)[1])


id_list<-as.numeric(files$V1)
for (j in 1:dim(files)[1]){

	
	IID_num<-id_list[j]
	print(paste0(as.character(j/35914),'%'))

        #if(IID_num %in% id_list){
	try(
	    {t1<-readNIfTI(paste0('/media/leelabsg_storage01/jangho/Aging/',IID_num,'_20252_2/T1/T1_brain_to_MNI.nii.gz'))
		t1_sliced<-t1[1:182,1:218,z]
	    rem=k%%200
	    if(rem==0){
		    rem=200
	    }
	
	for (i in 1:128){
		  t1_resized<-t1_sliced[1:182,1:218,i]
	  t1_resized<-resizeImage(t1_resized,128,128,method = 'nearest')
	  
	  final_array[1:128,1:128,i,rem]<-t1_resized
	
	}
	
	ID[k]<-IID_num
       	if(k%%200==0){
		quot=k/200
		np$save(paste0('/media/leelabsg_storage01/jangho/Aging/final_array_128_full_',as.character(quot),'.npy'),r_to_py(final_array))
		final_array<-array(0,c(128,128,128,200))
	}
	k=k+1
	    })
	
	if(j==dim(files)[1]){
		np$save(paste0('/media/leelabsg_storage01/jangho/Aging/final_array_128_full_last.npy'),r_to_py(final_array))
	}

}


ID<-data.frame(ID)
ID$age<-rep(0,35914)
birth<-read.table('Aging_all_birthyear_race.txt',h=T)
for (i in 1:35914){
	try(ID$age[i]<-2014-birth[which(birth$id==as.numeric(ID[i,1])),2])
}

np$save('/media/leelabsg_storage01/jangho/Aging/ID_128_full.npy',r_to_py(ID))

