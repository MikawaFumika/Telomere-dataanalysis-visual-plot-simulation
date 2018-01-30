setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据")

uuid-barcode
manifest= "gdc_manifest.2017-11-20T13_15_37.583830.txt" #Manifest name 
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1='{
  "filters":{
    "op":"and",
    "content":[
      {
        "op":"in",
        "content":{
          "field":"files.file_id",
          "value":['
Part2=']
                }
},
{
  "op":"=",
  "content":{
  "field":"files.data_type",
  "value":"Gene Expression Quantification"
  }
}
]
},
"format":"tsv",
"fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
"size":"'
Sentence= paste(Part1,id,Part2,manifest_length,'"}',collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)


releaseandbind-data
setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/")
bardatabre=read.table("File_metadatabreast.txt",header=T,check.names=FALSE,sep="\t")
bardatakid=read.table("File_metadatakidney.txt",header=T,check.names=FALSE,sep="\t")
bardatalun=read.table("File_metadatalung.txt",header=T,check.names=FALSE,sep="\t")
dataccle=read.table("CCLE_Expression_2012-09-29.res",fill=TRUE,sep="\t",row.names=F)
a=list.files(path="/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/lung")
dir = paste("/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/lung/",a,sep="") 
n = length(dir) 
datalun = read.table(file = dir[1],header=F,row.names=1)

for (i in 2:n){
  newdata = read.table(file = dir[i],header=F,row.names=1)
  datalun = cbind(datalun,newdata)
}
databre=read.table("/Users/yeyusong/Desktop/在投/Telomere/txt数据/Breast.txt",header=F,row.names=1)
datakid=read.table("/Users/yeyusong/Desktop/在投/Telomere/txt数据/Kidney.txt",header=F,row.names=1)
datalun=read.table("/Users/yeyusong/Desktop/在投/Telomere/txt数据/Lung.txt",header=F,row.names=1)


setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据")
write.table(datakid,"Kidney.txt",row.names=T,col.names=F,sep="\t")
write.table(databre,"Breast.txt",row.names=T,col.names=F,sep="\t")
write.table(datalun,"Lung.txt",row.names=T,col.names=F,sep="\t")


gene1="ENSG00000164362.17" #tert
gene2="ENSG00000254093.7"  #pinx1
gene2="ENSG00000147601.12" #terf1
gene2="ENSG00000128513.13" #pot1
gene1="ENSG00000270141.3"  #terc
gene2="ENSG00000132604.9"  #terf2
gene1="ENSG00000076864.18" #rap1gap1
gene2="ENSG00000132359.12" #rapgap2
gene1="ENSG00000127314.16" #rap1b
gene2="ENSG00000092330.14" #TINF2

sumbre=numeric(0)
for(i in 1:1220)
{sumbre[i]=sum(databre[,i])}
sumkid=numeric(0)
for(i in 1:1158)
{sumkid[i]=sum(datakid[,i])}
sumlun=numeric(0)
for(i in 1:1145)
{sumlun[i]=sum(datalun[,i])}

baridbre=character(0)
baridkid=character(0)
baridlun=character(0)
gen1=numeric(0)
gen2=numeric(0)
a=list.files(path="/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/breast")
n=length(a)
for (i in 1:n)
{
  baridbre[i]=as.character(bardatabre[which(bardatabre$file_name==paste0(a[i],".gz")),"cases.0.samples.0.submitter_id"])
}
gen1=t(databre[gene1,])/sumbre*1000000
gen2=t(databre[gene2,])/sumbre*1000000
supx=10
supy=50
par(mfrow=c(1,3))
plot(gen1[which(as.numeric(substr(baridbre,14,15))<10)]+0.001,gen2[which(as.numeric(substr(baridbre,14,15))<10)]+0.001,log='xy',xlim=c(0.001,100),ylim=c(1,100),type='p',main="Breast Cancer",cex=0.2,col="red",xlab="TERT",ylab="TINF2")
axis(side=1,at=c(0.001,0.01,0.1,1,10,100),labels=c("0",expression(10^"-2"),"0.1","1","10","100"),cex.axis=1)
axis(side=2,at=c(0.1,1,10,100),labels=c("1","10","100"),cex.axis=1)
par(new=TRUE)
plot(gen1[which(as.numeric(substr(baridbre,14,15))>9)]+0.001,gen2[which(as.numeric(substr(baridbre,14,15))>9)]+0.001,log='xy',xlim=c(0.001,100),ylim=c(1,100),xaxt='n',yaxt='n',col="blue",type='p',main="Breast Cancer",cex=0.2,xlab="TERT",ylab="TINF2")

par(mfrow=c(1,3))
cancersamples <- rep(1,length(gen1[which(as.numeric(substr(baridbre,14,15))<10)]))
normalsamples <- rep(2,length(gen1[which(as.numeric(substr(baridbre,14,15))>10)]))
plot(cancersamples,gen1[which(as.numeric(substr(baridbre,14,15))<10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="red",xlab="",ylab="FPKM",main="Breast")
par(new=TRUE)
plot(normalsamples,gen1[which(as.numeric(substr(baridbre,14,15))>10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="blue",xlab="",ylab="FPKM",main="Breast")

par(mfrow=c(1,2))
hist(gen2[which(as.numeric(substr(baridbre,14,15))<10)],xlim=c(0,10))
hist(gen2[which(as.numeric(substr(baridbre,14,15))>10)],xlim=c(0,10))

baridbre=character(0)
baridkid=character(0)
baridlun=character(0)
gen1=numeric(0)
gen2=numeric(0)
b=list.files(path="/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/kidney")
n=length(b)
for (i in 1:n)
{
  baridkid[i]=as.character(bardatakid[which(bardatakid$file_name==paste0(b[i],".gz")),"cases.0.samples.0.submitter_id"])
}
gen1=t(datakid[gene1,])/sumkid*1000000
gen2=t(datakid[gene2,])/sumkid*1000000
par(mfrow=c(1,2))
plot(gen1[which(as.numeric(substr(baridkid,14,15))<10)]+0.001,gen2[which(as.numeric(substr(baridkid,14,15))<10)]+0.001,log='xy',type='p',xlim=c(0.002,100),ylim=c(1,100),xaxt='n',yaxt='n',main="Kidney Cancer",cex=0.2,col="red",xlab="TERT",ylab="TERF1")
axis(side=1,at=c(0.002,0.1,10),labels=c("0",expression(10^"-1"),"10"),cex.axis=1)
axis(side=2,at=c(1,10,100),labels=c("0","10",expression(10^"2")),cex.axis=1)
par(new=TRUE)
plot(gen1[which(as.numeric(substr(baridkid,14,15))>9)]+0.001,gen2[which(as.numeric(substr(baridkid,14,15))>9)]+0.001,log='xy',col="blue",type='p',xlim=c(0.002,100),ylim=c(1,100),xaxt='n',yaxt='n',main="Kidney Cancer",cex=0.2,xlab="TERT",ylab="TERF1")
par(new=TRUE)
plot(0.40,8.8,type='p',pch=24,cex=1.5,col="black",xlim=c(0.002,100),ylim=c(1,100),xaxt='n',yaxt='n',xlab="",ylab="",main="",log='xy',bg="red")
par(new=TRUE)
plot(0.015,10.01,type='p',pch=22,cex=1.5,col="black",xlim=c(0.002,100),ylim=c(1,100),xaxt='n',yaxt='n',xlab="",ylab="",main="",log='xy',bg="blue")

cancersamples <- rep(1,length(gen1[which(as.numeric(substr(baridkid,14,15))<10)]))
normalsamples <- rep(2,length(gen1[which(as.numeric(substr(baridkid,14,15))>10)]))
plot(cancersamples,gen1[which(as.numeric(substr(baridkid,14,15))<10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="red",xlab="",ylab="FPKM",main="Kidney")
par(new=TRUE)
plot(normalsamples,gen1[which(as.numeric(substr(baridkid,14,15))>10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="blue",xlab="",ylab="FPKM",main="Kidney")

baridbre=character(0)
baridkid=character(0)
baridlun=character(0)
gen1=numeric(0)
gen2=numeric(0)
c=list.files(path="/Users/yeyusong/Desktop/在投/Telomere/txt数据/expression_data/lung")
n=length(c)
for (i in 1:n)
{
  baridlun[i]=as.character(bardatalun[which(bardatalun$file_name==paste0(c[i],".gz")),"cases.0.samples.0.submitter_id"])
}
gen1=t(datalun[gene1,])/sumlun*1000000
gen2=t(datalun[gene2,])/sumlun*1000000
plot(gen1[which(as.numeric(substr(baridlun,14,15))<10)]+0.001,gen2[which(as.numeric(substr(baridlun,14,15))<10)]+0.001,log='xy',xaxt='n',yaxt='n',xlim=c(0.002,100),ylim=c(10,1000),type='p',main="Lung Cancer",cex=0.2,col="red",xlab="TERT",ylab="TINF2")
axis(side=1,at=c(0.01,1,100),labels=c(expression(10^"-2"),"1",expression(10^"2")),cex.axis=1)
axis(side=2,at=c(10,100,1000),labels=c("10",expression(10^"2"),expression(10^"3")),cex.axis=1)
par(new=TRUE)
plot(gen1[which(as.numeric(substr(baridlun,14,15))>9)]+0.001,gen2[which(as.numeric(substr(baridlun,14,15))>9)]+0.001,log='xy',xlim=c(0.002,100),ylim=c(10,1000),xaxt='n',yaxt='n',col="blue",type='p',main="Lung Cancer",cex=0.2,xlab="TERT",ylab="TINF2")
par(new=TRUE)
plot(0.9699,65.3,type='p',pch=24,cex=1.5,col="black",xlim=c(0.002,100),ylim=c(10,1000),xaxt='n',yaxt='n',xlab="",ylab="",main="",log='xy',bg="red")
par(new=TRUE)
plot(0.03,86.69,type='p',pch=22,cex=1.5,col="black",xlim=c(0.002,100),ylim=c(10,1000),xaxt='n',yaxt='n',xlab="",ylab="",main="",log='xy',bg="blue")


cancersamples <- rep(1,length(gen1[which(as.numeric(substr(baridlun,14,15))<10)]))
normalsamples <- rep(2,length(gen1[which(as.numeric(substr(baridlun,14,15))>10)]))
plot(cancersamples,gen1[which(as.numeric(substr(baridlun,14,15))<10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="red",xlab="",ylab="FPKM",main="Lung")
par(new=TRUE)
plot(normalsamples,gen1[which(as.numeric(substr(baridlun,14,15))>10)],type='p',xaxt='n',xlim=c(0,3),ylim=c(0,30),cex=0.2,col="blue",xlab="",ylab="FPKM",main="Lung")
