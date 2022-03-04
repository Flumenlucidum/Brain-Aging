setwd("C:/Users/main/Desktop/FINAL/brain")
#Import libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(qqman)
library(data.table)
library(calibrate)

#CNN prediction results

a=read.csv('./result/FINAL_result.csv')
#Delta age of the individuals without diseases
mean(abs(a$delta_age[which(a$cv1==1 | a$cv2==1 | a$cv3==1 | a$cv4==1)]))
#Delta age of the individuals with diseases
mean(abs(a$delta_age[which(a$cv1==0 & a$cv2==0 & a$cv3==0 & a$cv4==0)]))
colnames(a)[c(2,3,5,10)]=c('True_Age','Predicted_Age', 'Delta','Corrected_Delta')
colnames(a)
#Plots
ggplot(a,aes(x=True_Age,y=Predicted_Age))+geom_point(color='blue',size=1)+geom_abline(slope=1,intercept = 0)+ggtitle('True Age and Predicted Age')+ theme(plot.title = element_text(hjust = 0.5, face= 'bold',size=20),axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))
ggplot(a,aes(x=True_Age,y=Delta))+geom_hline(yintercept = 0)+geom_point(color='blue',size=1)+geom_abline(slope=1,intercept = 0)+ggtitle('True Age and Delta Age')+ theme(plot.title = element_text(hjust = 0.5, face= 'bold',size=20),axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
ggplot(a,aes(x=True_Age,y=Corrected_Delta))+geom_hline(yintercept = 0)+geom_point(color='blue',size=1)+geom_abline(slope=1,intercept = 0)+ggtitle('True Age and Corrected Delta Age')+ theme(plot.title = element_text(hjust = 0.5, face= 'bold',size=20),axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


#SKAT-CommonRare results
a=read.table('./gene_based_result/DELTA_SKAT_CommonRare.txt',head=T)
head(a)
a<-a[c(1:4)]
a<-a[complete.cases(a),]
head(a)
a<-a[which(a$chr==17),]
c<-read.csv('../Delta_SKAT_CommonRare.csv')
for(i in 1:21){
  print(c$SNP[i] %in% a$Gene.refGene)
  print(i)
}
c$SNP[9]
b<-a
head(a)
colnames(a)[1:4]<-c('CHR','BP','SNP','P')
manhattan(a,main='Delta Age with WES validation(SKAT_CommonRare)')
b<-read.table('./gene_based_result/SKATO_validation_result.txt',head=T)

a=read.csv('./gene_based_result/DELTA_SKAT_CommonRare.txt',sep=' ')
b=read.csv('./validation_ow/DELTA_SKAT_CommonRare_validation.txt',sep=' ')
manhattan(b,col= c('#6699CC','black'),chr = 'chr',bp = 'Start', snp = 'Gene.refGene',p = 'P.value',annotatePval = 0.05/dim(b)[1],suggestiveline = -log10(0.05/dim(b)[1]),annotateTop = F,genomewideline = F,main='validation voxel value')
a=a[which(a$P<0.05/(dim(a)[1])),]
b<-b[which(b$Gene.refGene %in% a$Gene.refGene),]

#Again
setwd("C:/Users/main/Desktop/FINAL/brain")
#Delta age
a<-fread('./gene_based_result/DELTA_SKAT_CommonRare.txt')
a<-a[!grepl(';',a$Gene.refGene,fixed = T),] 
#Fornix
a<-fread('./fornix_thalamus/DELTA_SKAT_CommonRare_fornix.txt')
a<-a[!grepl(';',a$Gene.refGene,fixed = T),] 
#Lower part of the thalamus
a<-fread('./fornix_thalamus/DELTA_SKAT_CommonRare_thalamus.txt')
a<-a[!grepl(';',a$Gene.refGene,fixed = T),]

trace('manhattan',edit=T)
manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray10", "gray60"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))),
                   xlab=xlabel, ylab=expression(-log[10](italic(p))))
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=col[1], ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
      icol=icol+1
    }
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
  
  # Highlight snps from a character vector
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col="green3", pch=20, ...)) 
  }
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), 
           textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 1.2), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 1.2, ...)
    }
  }  
  par(xpd = FALSE)
}

manhattan(a,col= c('#afafaf','#007aff'),chr = 'chr',bp = 'Start', snp = 'Gene.refGene',p = 'P.value',annotatePval = 0.05/dim(a)[1],suggestiveline = -log10(0.05/dim(a)[1]),annotateTop = T,genomewideline = F)

manhattan(a,col= c('#afafaf','#007aff'),chr = 'chr',bp = 'Start', snp = 'Gene.refGene',p = 'P.value',annotatePval = 0.05/dim(a)[1],suggestiveline = -log10(0.05/dim(a)[1]),annotateTop = F,genomewideline = F,main='Fornix Volume')

manhattan(a,col= c('#6699CC','black'),chr = 'chr',bp = 'Start', snp = 'Gene.refGene',p = 'P.value',annotatePval = 0.05/dim(a)[1],suggestiveline = -log10(0.05/dim(a)[1]),annotateTop = T,genomewideline = F)

a<-a[order(a$P.value),]
View(a)
a<-a[which(a$P.value<0.05/dim(a)[1]),]
write.table(a,file='./gene_based_result/0227_revising/delta_fornix_significant.txt',row.names=F,quote=F)
write.table(a,file='./gene_based_result/0227_revising/delta_thalamus_significant.txt',row.names=F,quote=F)

#Significance of SNPs in OMA1, NOSTRIN, and CPNE3 genes
setwd("C:/Users/main/Desktop/FINAL/brain/gene_based_result/omanoscpn/")
a=read.table('ukb_whole_wes_result_commonrare_chr1.txt',head=T)
a=read.table('ukb_whole_wes_result_chr8.txt',head=T)
head(a)
a=a[which(a$P.value!=1),]
a$N.Marker.Rare[a$N.Marker.Rare==1]='Rare'
a$N.Marker.Rare[a$N.Marker.Rare==0]='Common'
a<-a[which(a$SetID>150),]

ggplot(a, aes(x=SetID, y=-log10(P.value))) + geom_point(aes(col=ifelse(SetID>150,'OMA1','DAB1;OMA1'),shape=N.Marker.Rare,size=5)) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(color="Rare variant", x="SNPs", y="-log10(p-value)")  +ggtitle('DAB1 and OMA1') 
ggplot(a, aes(x=SetID, y=-log10(P.value))) + geom_point(aes(col='red',shape=N.Marker.Rare,size=5)) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(shape="Rare variant", x="SNPs", y="-log10(p-value)")  +ggtitle('OMA1') 
ggplot(a, aes(x=SetID, y=-log10(P.value))) + geom_point(aes(col='red',shape=N.Marker.Rare,size=5)) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(shape="Rare variant", x="SNPs", y="-log10(p-value)")  +ggtitle('NOSTRIN') 
ggplot(a, aes(x=SetID, y=-log10(P.value))) + geom_point(aes(col='red',shape=N.Marker.Rare,size=5)) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(shape="Rare variant", x="SNPs", y="-log10(p-value)")  +ggtitle('CPNE3') 

plot(a$SetID,-log10(a$P.value),)


#Histogram of R-squared (Checking imputation quality)
setwd("C:/Users/main/Desktop/FINAL/brain/gene_based_result/omanoscpn/")
a=read.csv('correlation_of_othersnps_with_pval.csv')
a<-a[16:nrow(a),]

hist(a$correlation^2,main='Distribution of R-squared of Each SNP',xlab='R-Squared')
a$chr[a$chr==1]='OMA1'
a$chr[a$chr==2]='NOSTRIN'
a$chr[a$chr==8]='CPNE3'
colnames(a)[1]<-'Gene'
a$R2<-a$correlation^2

ggplot(a, aes(x=R2, fill=Gene)) +geom_histogram(binwidth = 0.08,color=1, alpha=0.5, position = 'dodge')+scale_fill_manual(values = c("red","green","blue"))+ labs(x="R-Squared")+geom_text_repel(data = subset(a, P.val.SKAT.<0.05/25105), aes(label = rsID,y=1,size=3,fontface=2),show.legend = F, box.padding = unit(0.45, "lines"),max.overlaps = Inf)
ggplot(a, aes(x=cor2, fill=Gene)) +geom_histogram(binwidth = 0.08,color=1, alpha=0.5, position = 'dodge')+scale_fill_manual(values = c("red","green","blue"))+ labs(x="R-Squared",y='count')+ggtitle('R-Squared of SNPs')+ theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank())+geom_text_repel(data = subset(a, P.val.SKAT. <0.05/25105), aes(label = rsID,y=5), 
box.padding = unit(0.45, "lines"),max.overlaps = Inf)
ggplot(a, aes(x=R2, fill=Gene)) +geom_histogram(binwidth = 0.08,color=1, alpha=0.5, position = 'dodge')+scale_fill_manual(values = c("red","green","blue"))+ labs(x="R-Squared",y='count')+geom_text_repel(data = subset(a, P.val.SKAT. <0.05/25105), aes(label = rsID,y=1,size=3,fontface=2),show.legend = F, box.padding = unit(0.45, "lines"),max.overlaps = Inf)


#Linear regression (Association)
setwd("C:/Users/main/Desktop/FINAL/brain/MR/")
a<-fread('ALL.txt')
colnames(a)[332]
result<-matrix(0,313,3)

for(i in 1:313){
  result[i,1]<-colnames(a)[i+19]
  sub_data<-select(a,c(3:18,(i+19),333))
  model<-lm(corrected_delta~.,data=sub_data)
  b<-summary(model)
  result[i,2]<-b$coefficients[18,1]
  result[i,3]<-b$coefficients[18,4]
  print(i)
}
colnames(a)[81]
result[65:313,1]
result[6,1]


#1~61 Blood
#65~313 Metabolomics
#for blood
blood<-as.data.frame(result[1:61,1:3])
blood$V3<-as.numeric(blood$V3)
blood<-blood[order(blood$V3),]
for(i in 1:nrow(blood)){
  if(blood[i,3]<(0.05*i)/62){
    print(i)
  }
}

name<-read.csv('blood_name.csv')
blood$V1<-as.numeric(substr(blood$V1,3,7))
blood<-left_join(blood,name,by=c('V1'='code'))
write.table(blood[1:30,1:4],file='regression_0128_blood.txt',sep=';',row.names=F,quote=F)

#for metabolomics 
meta<-as.data.frame(result[65:313,1:3])
meta$V3<-as.numeric(meta$V3)
meta<-meta[order(meta$V3),]
for(i in 1:nrow(meta)){
  if(meta[i,3]<(0.05*i)/249){
    print(i)
  }
}

name<-read.csv('meta_name.csv',head=F)
meta$V1<-as.numeric(substr(meta$V1,2,7))
meta<-left_join(meta,name,by='V1')
write.table(meta[1:51,1:4],file='regression_0128_meta.txt',row.names=F,quote=F,sep=';')

#Assigning a group to each biomarker and visualization
colnames(blood)<-c('code','beta','P','name')
colnames(meta)<-c('code','beta','P','name')
RES<-rbind(blood,meta)
View(RES)
type=read.csv('../saige_single/Nightingale_biomarker_groups_v2.csv')
View(type)
fin<-left_join(RES,type,by=c('code'='field_id'))
View(fin)
fin<-fin[order(fin$code),]
View(fin[which(fin$code>=30000),c('code','name','Group')])
fin$Group[which(fin$code ==30000)]='WBC'
fin$Group[which(fin$code >=30010 & fin$code<= 30070)]='RBC'
fin$Group[which(fin$code >=30080 & fin$code<= 30110)]='Platelet'
fin$Group[which(fin$code >=30120 & fin$code<= 30220)]='WBC'
fin$Group[which(fin$code %in% c(30170,30230))]='RBC'
fin$Group[which(fin$code >=30240 & fin$code<= 30300)]='Reticulocyte'
fin$Group[which(fin$code %in% c(30710,30740))]='Inflammatory biomarkers'
fin$Group[which(fin$code %in% c(30760,30780))]='Cholestrol (blood)'
fin$Group[which(is.na(fin$Group) & fin$code>=30000)]='Other blood biomarkers'
fin$Group[which(is.na(fin$Group) & fin$code<30000)]='Relative lipoprotein lipid concentrations'
fin<-fin[order(fin$Group),]

#Put name of the top two phenotypes in the plot
dd=fin %>% group_by(Group) %>% top_n(2,-log10(P))
dd<-as.data.frame(dd)
dd$name

ggplot(fin, aes(x=Group, y=-log10(P))) + geom_point(aes(col=Group,shape=beta>0),size=5) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(color="Category", x="Phenotypes", y="-log10(p-value)")  +geom_text_repel(data=. %>% mutate(label = ifelse(P < 0.05 & name %in% dd$name, as.character(name), "")), aes(label=label), size=3.5, box.padding = unit(0.7, "lines"))+ geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)+ggtitle('Linear Regression Result') 


#MR results 
#Getting the result of Linear MR 

setwd("C:/Users/main/Desktop/FINAL/brain/MR/Revised0215")
a<-read.csv('../../saige_single/current_MR_result.csv')
a<-a[which(a$Simplemedian!=0 &a$SE!=0),]

a$Egger_estimate<-0
a$Egger_se<-0
a$Egger_p<-0
for(i in 1:nrow(a)){
  idx<-1
  a$Egger_estimate[i]<-a[i,(6*idx+22)]
  a$Egger_se[i]<-a[i,(6*idx+23)]
  a$Egger_p[i]<-a[i,(6*idx+24)]
}

a$WM_estimate<-0
a$WM_se<-0
a$WM_p<-0
for(i in 1:nrow(a)){
  a$WM_estimate[i]<-a[i,10]
  a$WM_se[i]<-a[i,11]
  a$WM_p[i]<-a[i,12]
}
a$IVW_estimate<-0
a$IVW_se<-0
a$IVW_p<-0
for(i in 1:nrow(a)){
  idx=1
  a$IVW_estimate[i]<-a[i,(3*idx+13)]
  a$IVW_se[i]<-a[i,(3*idx+14)]
  a$IVW_p[i]<-a[i,(3*idx+15)]
}

Linear_MR<-a[1:nrow(a),c(2,3,62:70)]
write.csv(Linear_MR,file='./LINEAR_MR_RESULT.csv',row.names=F)
a<-read.csv('LINEAR_MR_RESULT.csv')

#Divide by Group
type=read.csv('../../saige_single/Nightingale_biomarker_groups_v2.csv')
View(type)
fin<-left_join(a,type,by=c('name'='title'))
View(fin)
dim(fin[which(fin$code>=30000),c('code','name','Group')])
fin<-fin[order(fin$code),]
View(fin[which(fin$code>=30000),c('code','name','Group')])
fin$Group[which(fin$code ==30000)]='WBC'
fin$Group[which(fin$code >=30010 & fin$code<= 30070)]='RBC'
fin$Group[which(fin$code >=30080 & fin$code<= 30110)]='Platelet'
fin$Group[which(fin$code >=30120 & fin$code<= 30220)]='WBC'
fin$Group[which(fin$code %in% c(30170,30230))]='RBC'
fin$Group[which(fin$code >=30240 & fin$code<= 30300)]='Reticulocyte'
fin$Group[which(fin$code %in% c(30710,30740))]='Inflammatory biomarkers'
fin$Group[which(fin$code %in% c(30760,30780))]='Cholestrol (blood)'
fin$Group[which(is.na(fin$Group) & fin$code>=30000)]='Other blood biomarkers'
fin$Group[which(is.na(fin$Group) & fin$code<30000)]='Relative lipoprotein lipid concentrations'
fin<-fin[order(fin$Group),]

#Put name of the top two phenotypes in the plot
dd=fin %>% group_by(Group) %>% top_n(2,-log10(WM_p))
dd<-as.data.frame(dd)
dd$name

ggplot(fin, aes(x=Group, y=-log10(IVW_p))) + geom_point(aes(col=Group,shape=IVW_estimate>0),size=5) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(color="Category", x="Phenotypes", y="-log10(p-value)")  +geom_text_repel(data=. %>% mutate(label = ifelse(IVW_p < 0.05 & name %in% dd$name, as.character(name), "")), aes(label=label), size=3.5, box.padding = unit(0.7, "lines"))+ geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)+ggtitle('IVW MR result') 
ggplot(fin, aes(x=Group, y=-log10(WM_p))) + geom_point(aes(col=Group,shape=WM_estimate>0),size=5) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(color="Category", x="Phenotypes", y="-log10(p-value)")  +geom_text_repel(data=. %>% mutate(label = ifelse(WM_p < 0.05 & name %in% dd$name, as.character(name), "")), aes(label=label), size=3.5, box.padding = unit(0.7, "lines"))+ geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)+ggtitle('WM MR result') 
ggplot(fin, aes(x=Group, y=-log10(Egger_p))) + geom_point(aes(col=Group,shape=Egger_estimate>0),size=5) + theme_classic() + theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"),axis.text.x = element_blank(), panel.grid.minor=element_line(colour = "grey", linetype="dashed"), axis.ticks=element_blank()) + labs(color="Category", x="Phenotypes", y="-log10(p-value)")  +geom_text_repel(data=. %>% mutate(label = ifelse(Egger_p < 0.05 & name %in% dd$name, as.character(name), "")), aes(label=label), size=3.5, box.padding = unit(0.7, "lines"))+ geom_hline(yintercept=-log10(0.05), color="red", size=1, alpha=0.5)+ggtitle('Egger MR result') 


#Benjamini-Hochberg procedure (Multiple testing correction)
a<-a[which(a$WM_p<0.05 | a$Egger_p<0.05 | a$IVW_p<0.05),]
View(a)
write.csv(a,file='LINEAR_MR_selected_0.05.csv',row.names=F)
a<-read.csv('LINEAR_MR_selected_0.05.csv')
a<-a[order(a$Egger_p),]
a<-a[order(a$WM_p),]
a<-a[order(a$IVW_p),]

a=a[which(a$code>=30000),]
dim(a)
for(i in 1:28){
  if(a$Egger_p[i]<(0.05*i/62)){
    print(a$name[i])
  }
}

a=a[which(a$code<30000),]
dim(a)
for(i in 1:31){
  if(a$Egger_p[i]<(0.05*i/249)){
    print(a$name[i])
  }
}

#WM
a<-a[order(a$WM_p),]
a=a[which(a$code>=30000),]
dim(a)
for(i in 1:28){
  if(a$WM_p[i]<(0.05*i/62)){
    print(a$name[i])
  }
}

a=a[which(a$code<30000),]
dim(a)
for(i in 1:31){
  if(a$WM_p[i]<(0.05*i/249)){
    print(a$name[i])
  }
}

#IVW
a<-a[order(a$IVW_p),]
a=a[which(a$code>=30000),]
dim(a)
for(i in 1:28){
  if(a$IVW_p[i]<(0.05*i/62)){
    print(a$name[i])
  }
}

a=a[which(a$code<30000),]
dim(a)
for(i in 1:31){
  if(a$IVW_p[i]<(0.05*i/249)){
    print(a$name[i])
  }
}

