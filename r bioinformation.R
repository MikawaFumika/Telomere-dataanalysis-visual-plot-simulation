setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据")
rnaseq1 <- read.table("rnaseq1.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)
rnacli1 <- read.table("rnacli1.txt",fill=TRUE,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)
rnaseq01 <- read.table("rnaseq01.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)

rnaseq2 <- read.table("rnaseq2.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)
rnacli2 <- read.table("rnacli2.txt",fill=TRUE,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)

rnaseq3 <- read.table("rnaseq3.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)
rnacli3 <- read.table("rnacli3.txt",fill=TRUE,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)

rnaseq4 <- read.table("rnaseq4.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)
rnacli4 <- read.table("rnacli4.txt",fill=TRUE,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)

rnaseq5 <- read.table("rnaseq5.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)
rnacli5 <- read.table("rnacli5.txt",fill=TRUE,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)

rnaseq6 <- read.table("rnaseq6.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)

rnaseq7 <- read.table("rnaseq7.txt",fill=TRUE,header=T,check.names=FALSE,stringsAsFactors = FALSE,row.names=1)

r1g1=as.numeric(r1g1)
rig1=as.numeric(r1g1)
(var(rig1))^0.5
rig1=(rig1-mean(rig1))/(12.25)
hist(rig1)
hist(as.numeric(r1g1))

gene1="TERF2IP"

r1g1 <- rnaseq1[gene1,]
r1g1 <- r1g1[-1]

r2g1 <- rnaseq2[gene1,]
r2g1 <- r2g1[-1]

r3g1 <- rnaseq3[gene1,]
r3g1 <- r3g1[-1]

r4g1 <- rnaseq4[gene1,]
r4g1 <- r4g1[-1]

r5g1 <- rnaseq5[gene1,]
r5g1 <- r5g1[-1]


r6g1 <- rnaseq6[gene1,]
r6g1 <- r6g1[-1]

par(mfrow=c(3,2))
hist(as.numeric(r1g1[abs(r1g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")
hist(as.numeric(r2g1[abs(r2g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")
hist(as.numeric(r3g1[abs(r3g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")
hist(as.numeric(r4g1[abs(r4g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")
hist(as.numeric(r5g1[abs(r5g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")
hist(as.numeric(r6g1[abs(r6g1)>1]),xlim=c(-10,10),breaks=20,main="",xlab="Z-score")

par(mfrow=c(3,2))
hist(as.numeric(r1g1),xlim=c(-10,10),breaks=20)
hist(as.numeric(r2g1),xlim=c(-10,10),breaks=20)
hist(as.numeric(r3g1),xlim=c(-10,10),breaks=20)
hist(as.numeric(r4g1),xlim=c(-10,10),breaks=20)
hist(as.numeric(r5g1),xlim=c(-10,10),breaks=20)
hist(as.numeric(r6g1),xlim=c(-10,10),breaks=20)

r1s=as.numeric(length(r1g1[r1g1<(-1)]))/as.numeric((length(r1g1)))
r1l=as.numeric(length(r1g1[r1g1>1]))/as.numeric((length(r1g1)))
r2s=as.numeric(length(r2g1[r2g1<(-1)]))/as.numeric((length(r2g1)))
r2l=as.numeric(length(r2g1[r2g1>1]))/as.numeric((length(r2g1)))
r3s=as.numeric(length(r3g1[r3g1<(-1)]))/as.numeric((length(r3g1)))
r4s=as.numeric(length(r4g1[r4g1<(-1)]))/as.numeric((length(r4g1)))
r4l=as.numeric(length(r4g1[r4g1>1]))/as.numeric((length(r4g1)))
r5s=as.numeric(length(r5g1[r5g1<(-1)]))/as.numeric((length(r5g1)))
r5l=as.numeric(length(r5g1[r5g1>1]))/as.numeric((length(r5g1)))
r6s=as.numeric(length(r6g1[r6g1<(-1)]))/as.numeric((length(r6g1)))
r6l=as.numeric(length(r6g1[r6g1>1]))/as.numeric((length(r6g1)))


gene1="TERT"
gene2="TERF1"
geexp1 <- rnaseq2[gene1,]
geexp1 <- as.numeric(geexp1[-1])
geexp2 <- rnaseq2[gene2,]
geexp2 <- as.numeric(geexp2[-1])
geneexp1=numeric(0)
geneexp2=numeric(0)
for(i in 1:length(geexp1))
    {if(abs(geexp1[i])>1.5 | abs(geexp2[i])>1.5)
    {geneexp1[i]=geexp1[i];
    geneexp2[i]=geexp2[i]}}
plot(geneexp1,geneexp2,type="p",cex=0.5,col="blue",pch=21,bg="blue",xlim=c(-5,5),ylim=c(-5,5))


mea=numeric(0)
for(i in 1:20000)
{
r <- rnaseq01[i,]
r <- r[-1]
r <- as.numeric(sort(r))
l=r[0.95*length(r)]+1
r=r[r<l]
mea[i]=mean(r)
}
hist(mea,xlim=c(0,5000),breaks=5000)

tert1stage2 <- rnaseq1[rnacli1$SAMPLE_ID[rnacli1$TumorStage == 'Stage III']][17815,]
hist(as.numeric(tert1stage2))


tert1stage1 <- numeric(0)
t1s1l=1
tert1stage2 <- numeric(0)
t1s2l=1
tert1stage3 <- numeric(0)
t1s3l=1
tert1stage4 <- numeric(0)
t1s4l=1
tert="PINX1"
lengt=129
x=rnaseq1
y=rnacli1
for(i in 1:lengt)
{if(y[names(x)[i+2],"TumorStage"]=="Stage I"){tert1stage1[t1s1l]=x[tert,i+2];t1s1l=t1s1l+1}
  else if(y[names(x)[i+2],"TumorStage"]=="Stage II"){tert1stage2[t1s2l]=x[tert,i+2];t1s2l=t1s2l+1}
  else if(y[names(x)[i+2],"TumorStage"]=="Stage III"){tert1stage3[t1s3l]=x[tert,i+2];t1s3l=t1s3l+1}
  else if(y[names(x)[i+2],"TumorStage"]=="Stage IV"){tert1stage4[t1s4l]=x[tert,i+2];t1s4l=t1s4l+1}
}
par(mfrow=c(2,2))
hist(tert1stage1,breaks=20)
hist(tert1stage2,breaks=20)
hist(tert1stage3,breaks=20)
hist(tert1stage4,breaks=20)

tert1stage1 <- numeric(0)
t1s1l=1
tert1stage2 <- numeric(0)
t1s2l=1
tert1stage3 <- numeric(0)
t1s3l=1
tert1stage4 <- numeric(0)
t1s4l=1
tert=16883
lengt=178
x=rnaseq2
y=rnacli2
for(i in 1:lengt)
{if(y[chartr("TCGA","LUSC",names(x)[i+2]),"CLIN_T_STAGE"]=="T1"){tert1stage1[t1s1l]=x[tert,i+2];t1s1l=t1s1l+1}
  else if(y[chartr("TCGA","LUSC",names(x)[i+2]),"CLIN_T_STAGE"]=="T2"){tert1stage2[t1s2l]=x[tert,i+2];t1s2l=t1s2l+1}
  else if(y[chartr("TCGA","LUSC",names(x)[i+2]),"CLIN_T_STAGE"]=="T3"){tert1stage3[t1s3l]=x[tert,i+2];t1s3l=t1s3l+1}
  else if(y[chartr("TCGA","LUSC",names(x)[i+2]),"CLIN_T_STAGE"]=="T4"){tert1stage4[t1s4l]=x[tert,i+2];t1s4l=t1s4l+1}
}
par(mfrow=c(2,2))
hist(tert1stage1,breaks=20)
hist(tert1stage2,breaks=20)
hist(tert1stage3,breaks=20)
hist(tert1stage4,breaks=20)
  