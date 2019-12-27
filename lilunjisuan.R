setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据")

guassianinitial <- read.table("guassianinitialdistri.txt",fill=TRUE,header=F)
uniforminitial <- read.table("uniforminitialdistri.txt",fill=TRUE,header=F)

ceshix=c(1,1)
ceshiy=c(1,2)



library(reshape2)
library(ggplot2)



//风度
library(SDMTools)
mystats <- function(x, na.omit=FALSE){
  if (na.omit)
    x <- x[!is.na(x)]
  m <- mean(x)
  n <- length(x)
  s <- sd(x)
  skew <- sum((x-m)^3/s^3)/n
  kurt <- sum((x-m)^4/s^4)/n - 3
  return(c(n=n, mean=m, stdev=s, skew=skew, kurtosis=kurt))
}


par(mfrow=c(2,1),mar=c(4.2,5,1,1))

attach(mtcars)
layout(matrix(c(1,1,1,2,3,4),2,3,byrow=TRUE))
par(mar=c(3.5,5.2,1,1.5))

//telochange
histh=0.00025

telochange20 <- read.table("telochange20.txt",fill=TRUE,header=F)
telomerechange=vector(mode="numeric",length=0)
for(i in 2:8000)
{telomerelengthchange=(t(telochange20[i,]))
telomerechange=c(telomerechange,telomerelengthchange)}
telomerechange=as.numeric(telomerechange)
histtelochange20=hist(telomerechange,breaks=200,xlim=c(-250,250),ylim=c(1,200000))
hist20density=as.numeric(as.character(histtelochange20$density))
hist20density[46]=hist20density[46]/hist20density[1]*histh
hist20density[1]=histh
hist20mids=as.numeric(as.character(histtelochange20$mids))

telochange33 <- read.table("telochange33.txt",fill=TRUE,header=F)
telomerechange=vector(mode="numeric",length=0)
for(i in 2:5700)
{telomerelengthchange=(t(telochange33[i,]))
telomerechange=c(telomerechange,telomerelengthchange)}
telomerechange=as.numeric(telomerechange)
histtelochange33=hist(telomerechange,breaks=200,xlim=c(-250,250),ylim=c(1,200000))
hist33density=as.numeric(as.character(histtelochange33$density))
hist33density[46]=hist33density[46]/hist33density[1]*histh
hist33density[1]=histh
hist33mids=as.numeric(as.character(histtelochange33$mids))

telochange49 <- read.table("telochange49.txt",fill=TRUE,header=F)
telomerechange=vector(mode="numeric",length=0)
for(i in 2:2800)
{telomerelengthchange=(t(telochange49[i,]))
telomerechange=c(telomerechange,telomerelengthchange)}
telomerechange=as.numeric(telomerechange)
histtelochange49=hist(telomerechange,breaks=200,xlim=c(-250,250),ylim=c(1,200000))
hist49density=as.numeric(as.character(histtelochange49$density))
hist49density[46]=hist49density[46]/hist49density[1]*histh
hist49density[1]=histh
hist49mids=as.numeric(as.character(histtelochange49$mids))

telochange55 <- read.table("telochange55.txt",fill=TRUE,header=F)
telomerechange=vector(mode="numeric",length=0)
for(i in 2:800)
{telomerelengthchange=(t(telochange55[i,]))
telomerechange=c(telomerechange,telomerelengthchange)}
telomerechange=as.numeric(telomerechange)
histtelochange55=hist(telomerechange,breaks=200,xlim=c(-250,250),ylim=c(1,200000))
hist55density=as.numeric(as.character(histtelochange55$density))
hist55density[46]=hist55density[46]/hist55density[1]*histh
hist55density[1]=histh
hist55mids=as.numeric(as.character(histtelochange55$mids))

plot(hist55mids,hist55density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='l',lwd=1.5,col="darkgreen",yaxt="n")
par(new=TRUE)
polygon(c(-228,hist55mids[1:190],700),c(0,hist55density[1:190],0), density = NULL, border = NULL, col = "darkgreen")
#plot(hist55mids,hist55density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='h',cex=0.5,col="darkgreen",yaxt="n")
par(new=TRUE)
plot(hist49mids,hist49density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='l',lwd=1.5,col="red",yaxt="n")
par(new=TRUE)
polygon(c(-228,hist49mids[1:190],700),c(0,hist49density[1:190],0), density = NULL, border = NULL, col = "red")
#plot(hist49mids,hist49density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='h',cex=0.5,col="red",yaxt="n")
par(new=TRUE)
plot(hist33mids,hist33density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='l',lwd=1.5,col="blue",yaxt="n")
par(new=TRUE)
polygon(c(-228,hist33mids[1:190],700),c(0,hist33density[1:190],0), density = NULL, border = NULL, col = "blue")
#plot(hist33mids,hist33density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='h',cex=0.5,col="blue",yaxt="n")
par(new=TRUE)
plot(hist20mids,hist20density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='l',lwd=1.5,col="black",yaxt="n")
par(new=TRUE)
polygon(c(-228,hist20mids[1:190],700),c(0,hist20density[1:190],0), density = NULL, border = NULL, col = rgb(220/256,220/256,220/256))
#plot(hist20mids,hist20density,tck=0.03,las=1,xlim=c(-250,700),mgp=c(0,0.2,0),main="",xlab="",ylab="",ylim=c(0,histh),type='h',cex=0.5,cex.lab=1.3,col="black",yaxt="n")
axis(side=2,at=c(0,0.0001,0.0002,0.00025),tck=0.03,mgp=c(0,0.2,0),las=1,labels=c(expression(0),expression(10^"-4"),expression(2%*%10^"-4"),expression(10^"-1")),cex.axis=1.0)
mtext("Proportion",side=2,line=3,cex=1.4)
mtext(expression(Delta),side=1,line=1.8,cex=1.4,at=220)
mtext("TL",side=1,line=1.8,cex=1.4,at=250)
mtext("A",side=2,line=3.5,cex=2,at=(histh*1),las=1)
legend("topright",legend=c("PD55","PD49","PD33","PD20"),col=c("darkgreen","red","blue","black"),bty="n",cex=1.2,lty=1,lwd=1.5)
deltaL<-t(rbind((hist20mids-2.5),hist20density,hist33density,hist49density,hist55density))
write.table(deltaL,"deltaL",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


skewgaussian <- read.table("skewgaussian.txt",fill=TRUE,header=F)
skewuniform <- read.table("skewuniform.txt",fill=TRUE,header=F)
skewnumbers=133
xskew <- seq(1,length(seq(3,skewnumbers,2)),1)
skews=vector(mode="numeric",length(seq(3,skewnumbers,2)))
for(i in seq(3,skewnumbers,2))
{
skew=(t(skewgaussian[i,]))
skew=as.numeric(skew)
skew <- skew[1:32]
telolength=seq(500,16000,500)
z <-vector(mode="numeric",length=0)
for(j in 1:32)
  {xdistrisample=vector(mode="numeric",length=(skew[j]))
  xdistrisample[] <- telolength[j]
  z=c(z,xdistrisample)}
x=mystats(z)
skews[(i-1)/2]=x[4]
}
plot(xskew,skews,type='l',lty=1,ylim=c(0,1.5),lwd=2,tck=0.03,las=1,xlim=c(0,70),mgp=c(0,0.2,0),col="blue",cex=0.3,xlab="",ylab="")
par(new=TRUE)
skewnumbers=127
xskew <- seq(1,length(seq(3,skewnumbers,2)),1)
skews=vector(mode="numeric",length(seq(3,skewnumbers,2)))
for(i in seq(3,skewnumbers,2))
{
  skew=(t(skewuniform[i,]))
  skew=as.numeric(skew)
  skew <- skew[1:32]
  telolength=seq(500,16000,500)
  z <-vector(mode="numeric",length=0)
  for(j in 1:32)
  {xdistrisample=vector(mode="numeric",length=(skew[j]))
  xdistrisample[] <- telolength[j]
  z=c(z,xdistrisample)}
  x=mystats(z)
  skews[(i-1)/2]=x[4]
}
plot(xskew,skews,ylim=c(0,1.5),type='l',lty=2,lwd=3,tck=0.03,las=1,xlim=c(0,70),col="red",mgp=c(0,0.2,0),cex=0.3,xlab="",ylab="")
experimentalskewx=c(10,33,49,55)
experimentalskew=c(0.37,0.64,1.02,1.28)
par(new=TRUE)
plot(experimentalskewx,experimentalskew,tck=0.03,las=1,ylim=c(0,1.5),xlim=c(0,70),mgp=c(0,0.2,0),cex=1.0,pch=24,xlab="",ylab="",cex.lab=1.3)
par(new=TRUE)
plot(experimentalskewx,experimentalskew,tck=0.03,las=1,ylim=c(0,1.5),xlim=c(0,70),mgp=c(0,0.2,0),cex=1.0,type='l',lwd=1.0,xlab="",ylab="",cex.lab=1.3)
mtext("Skewness",side=2,line=1.5,cex=1.4)
mtext("PD",side=1,line=2,cex=1.4)
mtext("B",side=2,line=2.5,cex=2,at=(1.5*1),las=1)
legend("topleft",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",pch=c(1,1,24),lty=c(1,2,1),cex=1.2,pt.cex=c(0.01,0.01,1))


//pinjunchangdu
meantelomerelength=(t(guassianinitial[1*5,]))
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[meantelomerelength>0]
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
pds=seq(1,length(meantelomerelength),1)
plot(pds,meantelomerelength,tck=0.03,las=1,mgp=c(0, 0.2, 0),type='l',lty=1,cex=0.2,col="blue",lwd=2,xlim=c(0,70),ylim=c(0,10000),xlab="",ylab="")
par(new=TRUE)
meantelomerelength=(t(uniforminitial[1*5,]))
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[meantelomerelength>0]
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
pds=seq(1,length(meantelomerelength),1)
plot(pds,meantelomerelength,mgp=c(0, 0.2, 0),tck=0.03,las=1,type='l',lty=2,cex=0.2,col="red",lwd=3,xlim=c(0,70),ylim=c(0,10000),xlab="",ylab="")
x=c(10,33,49,55)
xtelo=seq(250,500*31+250,500)
meanteloexpriment10=sum(xtelo*z1)/sum(z1)
meanteloexpriment33=sum(xtelo*z3)/sum(z3)
meanteloexpriment49=sum(xtelo*z5)/sum(z5)
meanteloexpriment55=sum(xtelo*z7)/sum(z7)
par(new=TRUE)
plot(x,c(meanteloexpriment10,meanteloexpriment33,meanteloexpriment49,meanteloexpriment55),mgp=c(0, 0.2, 0),tck=0.03,las=1,type='l',lwd=1,,cex=1.2,xlim=c(0,70),ylim=c(0,10000),xlab="",ylab="")
par(new=TRUE)
plot(x,c(meanteloexpriment10,meanteloexpriment33,meanteloexpriment49,meanteloexpriment55),type='p',mgp=c(0, 0.2, 0),tck=0.03,las=1,pch=24,cex=1.2,xlim=c(0,70),ylim=c(0,10000),xlab="",ylab="")
mtext("Mean TL",side=2,line=2.5,cex=1.4)
mtext("PD",side=1,line=2,cex=1.4)
mtext("C",side=2,line=3.2,cex=2,at=(10000*1),las=1)
legend("topright",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",pch=c(1,1,24),lty=c(1,2,1),cex=1.2,pt.cex=c(0.01,0.01,1))

//pingjunchangdu he population
meantelomerelength=(t(guassianinitial[1*5,]))
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[meantelomerelength>0]
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
meantelomerelengthy = (t(guassianinitial[1*3,]))
meantelomerelengthy=as.numeric(meantelomerelengthy)
meantelomerelengthy<-meantelomerelengthy[meantelomerelengthy>0]
meantelomerelengthy<-meantelomerelengthy[!is.na(meantelomerelengthy)]
meantelomerelengthy=meantelomerelengthy[1:length(meantelomerelength)]
meantelomerelengthy=log10(meantelomerelengthy)
plot(meantelomerelength,meantelomerelengthy,tck=0.03,las=1,mgp=c(0, 0.2, 0),type='l',lty=1,cex=0.2,col="blue",lwd=2,xlim=c(2000,8000),ylim=c(0,15),xlab="",ylab="",yaxt="n")
par(new=TRUE)
meantelomerelength=(t(uniforminitial[1*5,]))
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[meantelomerelength>0]
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
meantelomerelengthy = (t(uniforminitial[1*3,]))
meantelomerelengthy=as.numeric(meantelomerelengthy)
meantelomerelengthy<-meantelomerelengthy[meantelomerelengthy>0]
meantelomerelengthy<-meantelomerelengthy[!is.na(meantelomerelengthy)]
meantelomerelengthy=meantelomerelengthy[1:length(meantelomerelength)]
meantelomerelengthy=log10(meantelomerelengthy)
plot(meantelomerelength,meantelomerelengthy,mgp=c(0, 0.2, 0),tck=0.03,las=1,type='l',lty=2,cex=0.2,col="red",lwd=3,xlim=c(2000,8000),ylim=c(0,15),xlab="",ylab="",yaxt="n")
axis(side=2,las=1,at=c(2,4,6,8,10,12,14),labels=c(expression(10^"2"),expression(10^"4"),expression(10^"6"),expression(10^"8"),expression(10^"10"),expression(10^"12"),expression(10^"14")),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
mtext("Population",side=2,line=1.6,cex=1.4)
mtext("Mean TL",side=1,line=2,cex=1.4)
mtext("D",side=2,line=2.5,cex=2,at=(15*1),las=1)
legend("topright",legend=c("Gaussian","Uniform"),col=c("blue","red"),bty="n",pch=c(1,1),lty=c(1,2),cex=1.2,pt.cex=c(0.01,0.01))




//potential
library(rgl)
library(gplots)
mkrchange <- read.table("mchangekrchange.txt",fill=TRUE,header=F)
mchange <- seq(0,10,1)
for(i in 1:9)
{mchange=c(mchange,seq(0,10,1))}
rchange <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
krchange <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
for(i in 1:9)
{krchange=c(krchange,(rchange-0.01*i))}
potential=vector(mode="numeric",110)
for(i in 1:110)
{
  celld=(t(mkrchange[i*5-3,]))
  celld=celld[!is.na(celld)]
  potential[i]=length(celld)
}
Mpotential=matrix(nrow=10,ncol=11)
rownames(Mpotential)<-c("0.1","0.09","0.08","0.07","0.06","0.05","0.04","0.03","0.02","0.01")
colnames(Mpotential)<-c("0","1","2","3","4","5","6","7","8","9","10")
for(i in 1:10)
{for(j in 1:11)
{
  celld=(t(mkrchange[(((i-1)*55+j*5)-3),]))
  celld=celld[!is.na(celld)]
  Mpotential[i,j]=length(celld)
}
}
rm(p)
data_mpopulation <- melt(Mpotential)
p <- ggplot(data_mpopulation,aes(x=Var2,y=Var1))
p <- p + geom_tile(aes(fill=value))
p <- p + scale_fill_gradient(low = "white", high = "red")
p <- p + xlab("M") + ylab("kr") + theme_bw() + theme(panel.grid = element_blank()) + theme(legend.key=element_blank())
p + theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p + labs(title="Cell replicative potential (PD)")



par(mar = c(7,18,4,2),oma=c(0.2,0.2,0.2,0.2),mex=0.5)           # 设定边缘
image(x=1:nrow(Mpotential),y=1:ncol(Mpotential),z=Mpotential,axes=FALSE,xlab="",ylab="",col=heat.colors(10),,main="")  
heatmap(Mpotential, Rowv=NA, Colv=NA, col=cm.colors(256), revC=FALSE, scale='column')
heatmap(Mpotential, Rowv=NA, Colv=NA, col=redgreen, revC=FALSE, scale="column", margins=c(5,10))
heatmap.2(Mpotential,col="bluered",Rowv=F,Colv=F,trace="none",cexCol=1,cexRow=1,xlab="M",ylab="kr",cex.lab=4.0,keysize=2.0,key.title="",key.xlab="Replicative Potential",key.ylab="counts",lhei=c(1,3),lwid=c(1,3))



//siwangshuliang
for(i in 1:15)
{meantelomerelength=(t(mchange[(i*5-3),])
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[meantelomerelength>0]
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
pds=seq(1,length(meantelomerelength),1)
plot(pds,meantelomerelength,type='l',col=rgb((i*15)/255,0,0),cex=0.3,xlim=c(0,300),ylim=c(0,10000),xlab="PD",ylab="death cell number")
par(new=TRUE)}

//duanduanlizhanbi
Mshorttelo=matrix(nrow=10,ncol=11)
rownames(Mshorttelo)<-c("0.1","0.09","0.08","0.07","0.06","0.05","0.04","0.03","0.02","0.01")
colnames(Mshorttelo)<-c("0","1","2","3","4","5","6","7","8","9","10")
for(i in 1:10)
{for(j in 1:11)
{
  celld=(t(mkrchange[(((i-1)*55+j*5)-2),]))
  celld=as.numeric(celld[!is.na(celld)])
  Mshorttelo[i,j]=as.numeric(mkrchange[(((i-1)*55+j*5)-1),min(which(celld==max(celld)))])
}
}
rm(p)
data_mpopulation <- melt(Mshorttelo)
p <- ggplot(data_mpopulation,aes(x=Var2,y=Var1))
p <- p + geom_tile(aes(fill=value))
p <- p + scale_fill_gradient(low = "white", high = "red")
p <- p + xlab("M") + ylab("kr") + theme_bw() + theme(panel.grid = element_blank()) + theme(legend.key=element_blank())
p + theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p + labs(title="D Short telomere proportion")

heatmap.2(Mshorttelo,col="bluered",Rowv=F,Colv=F,trace="none",cexCol=1,cexRow=1,xlab="M",ylab="kr",cex.lab=4.0,keysize=2.0,key.title="",key.xlab="Proportion of short telomere",key.ylab="counts",lhei=c(1,3),lwid=c(1,3))

par(mfrow=c(1,2),mar=c(4.2,4.2,1,1))
for(i in 1:6)
{meantelomerelength=(t(mkrchange[(10*i-6),]))
popshorttelo=(t(mkrchange[(10*i-6),]))
meantelomerelength=as.numeric(meantelomerelength)
meantelomerelength <- meantelomerelength[!is.na(meantelomerelength)]
meantelomerelength <- meantelomerelength[1:(which(popshorttelo==max(popshorttelo)))]
pds=seq(1,length(meantelomerelength),1)
plot(pds,meantelomerelength,type='l',mgp=c(0,0.2,0),tck=0.02,las=1,col=rgb((i*40)/255,0,0),cex=0.3,lwd=1.5,xlim=c(0,100),ylim=c(0,0.75),xlab="",ylab="",yaxt="n")
par(new=TRUE)
}
mtext("Short telomere(%)",side=2,line=2.0,cex=1.4)
mtext("PD",side=1,line=1.3,cex=1.4)
axis(side=2,las=1,at=c(0,0.25,0.5,0.75),labels=c("0%","25%","50%","75%"),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
legend("topleft",legend=c("M=10","M=8","M=6","M=4","M=2","M=0"),cex=1.0,x.intersp=0.2,y.intersp=0.9,col=c(rgb(240/255,0,0),rgb(200/255,0,0),rgb(160/255,0,0),rgb(120/255,0,0),rgb(80/255,0,0),rgb(40/255,0,0)),bty="n",lty=1,lwd=1.5)
mtext("A",side=2,line=2.5,cex=2,at=(0.75*1),las=1)

mvar=c(0,1,2,3,4,5,6,7,8,9,10)
mkrshorttelo=(t(Mshorttelo[1,]))
mkrshorttelo=as.numeric(mkrshorttelo)
plot(mvar,mkrshorttelo,type='p',pch=21,col=rgb(0/255,0,0),bg=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
plot(mvar,mkrshorttelo,type='l',pch=21,col=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
mkrshorttelo=(t(Mshorttelo[4,]))
mkrshorttelo=as.numeric(mkrshorttelo)
plot(mvar,mkrshorttelo,type='p',pch=22,col=rgb(80/255,0,0),bg=rgb(80/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
plot(mvar,mkrshorttelo,type='l',pch=22,col=rgb(80/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
mkrshorttelo=(t(Mshorttelo[7,]))
mkrshorttelo=as.numeric(mkrshorttelo)
plot(mvar,mkrshorttelo,type='p',pch=23,col=rgb(160/255,0,0),bg=rgb(160/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
plot(mvar,mkrshorttelo,type='l',pch=23,col=rgb(160/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
mkrshorttelo=(t(Mshorttelo[10,]))
mkrshorttelo=as.numeric(mkrshorttelo)
plot(mvar,mkrshorttelo,type='p',pch=24,col=rgb(240/255,0,0),bg=rgb(240/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
par(new=TRUE)
plot(mvar,mkrshorttelo,type='l',pch=24,col=rgb(240/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(0.15,0.3),xlab="",ylab="",yaxt="n")
legend("topleft",legend=c("kr=0.01","kr=0.04","kr=0.07","kr=0.1"),y.intersp=1,cex=1,bty="n",bg=c(rgb(240/255,0,0),rgb(160/255,0,0),rgb(80/255,0,0),rgb(0/255,0,0)),col=c(rgb(240/255,0,0),rgb(160/255,0,0),rgb(80/255,0,0),rgb(0/255,0,0)),lty=1,lwd=1.5,pch=c(17,18,15,16))
mtext("Accmulative short telomere(%)",side=2,line=2.0,cex=1.1)
mtext("Telomerase concentration",side=1,line=1.3,cex=1.3)
axis(side=2,las=1,at=c(0.15,0.2,0.25,0.3),labels=c("15%","20%","25%","30%"),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
mtext("B",side=2,line=2.5,cex=2,at=(0.3*1),las=1)

//population & potential
mkrchange <- read.table("mchangekrchange.txt",fill=TRUE,header=F)

par(mfrow=c(2,2),mar=c(4.2,4.2,0.5,1.5))
for(i in 1:6)
{population=(t(mkrchange[(i*10-7),]))
population=as.numeric(population)
population <- population[population>0]
population <- population[!is.na(population)]
population=log10(population)
pds=seq(1,length(population),1)
plot(pds,population,type='l',col=rgb((i*50-50)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,100),ylim=c(3,15),xlab="",ylab="",yaxt="n")
par(new=TRUE)}
mtext("Population",side=2,line=2,cex=1.6)
mtext("PD",side=1,line=2,cex=1.6)
mtext("A",side=2,line=2.7,cex=2,at=(15*1),las=1)
axis(side=2,las=1,at=c(3,5,10,15),labels=c(expression(10^"3"),expression(10^"5"),expression(10^"10"),expression(10^"15")),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
legend("topleft",legend=c("M=10","M=8","M=6","M=4","M=2","M=0"),cex=1.0,x.intersp=0.2,y.intersp=0.5,bty="n",col=c(rgb((6*50-50)/255,0,0),rgb((5*50-50)/255,0,0),rgb((4*50-50)/255,0,0),rgb((3*50-50)/255,0,0),rgb((2*50-50)/255,0,0),rgb((1*50-50)/255,0,0)),lty=1,lwd=1.5)


for(i in 1:5)
{population=(t(mkrchange[(i*110-2),]))
population=as.numeric(population)
population <- population[population>0]
population <- population[!is.na(population)]
population=log10(population)
pds=seq(1,length(population),1)
plot(pds,population,type='l',mgp=c(0,0.2,0),tck=0.03,las=1,col=rgb((i*50)/255,0,0),lwd=1.5,xlim=c(0,400),ylim=c(3,30),xlab="",ylab="",yaxt="n")
par(new=TRUE)
}
mtext("Population",side=2,line=2,cex=1.6)
mtext("PD",side=1,line=2,cex=1.6)
mtext("B",side=2,line=2.7,cex=2,at=(30*1),las=1)
axis(side=2,las=1,at=c(3,5,10,15,20,25,30),mgp=c(0,0.2,0),tck=0.03,las=1,labels=c(expression(10^"3"),expression(10^"5"),expression(10^"10"),expression(10^"15"),expression(10^"20"),expression(10^"25"),expression(10^"30")),cex.axis=1.0)
legend(250,31,legend=c("kr=0.01","kr=0.03","kr=0.05","kr=0.07","kr=0.09"),cex=1.0,x.intersp=0.2,y.intersp=0.5,bty="n",col=c(rgb((5*50-50)/255,0,0),rgb((4*50-50)/255,0,0),rgb((3*50-50)/255,0,0),rgb((2*50-50)/255,0,0),rgb((1*50-50)/255,0,0)),lty=1,lwd=1.5)



//potential
potential0.09=vector(mode="numeric",11)
for(i in 1:11)
{potential=(t(mkrchange[(55+i*5-3),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential0.09[i]=length(potential)}
potential0.03=vector(mode="numeric",11)
for(i in 1:11)
{potential=(t(mkrchange[(385+i*5-3),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential0.03[i]=length(potential)}
potential0.01=vector(mode="numeric",11)
for(i in 1:11)
{potential=(t(mkrchange[(495+i*5-3),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential0.01[i]=length(potential)}
pds=seq(0,10,1)
plot(pds,potential0.09,type='p',pch=22,col=rgb(0/255,0,0),bg=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.09,type='l',pch=22,col=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.03,type='p',pch=23,col=rgb(110/255,0,0),bg=rgb(110/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.03,type='l',pch=23,col=rgb(110/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.01,type='p',pch=24,col=rgb(220/255,0,0),bg=rgb(220/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.01,type='l',pch=24,col=rgb(220/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,10),ylim=c(50,400),xlab="",ylab="")
mtext("Replicative potential",side=2,line=2,cex=1.4)
mtext("Telomerase concentration",side=1,line=2,cex=1.3)
mtext("C",side=2,line=2.7,cex=2,at=(400*1),las=1)
legend("topleft",legend=c("kr=0.01","kr=0.03","kr=0.09"),y.intersp=0.5,cex=1.2,bty="n",bg=c(rgb(220/255,0,0),rgb(110/255,0,0),rgb(0/255,0,0)),col=c(rgb(220/255,0,0),rgb(110/255,0,0),rgb(0/255,0,0)),lty=1,lwd=1.5,pch=c(17,18,15))

potential0=vector(mode="numeric",10)
for(i in 1:10)
{potential=(t(mkrchange[(55*i-55+3),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential0[i]=length(potential)}
potential5=vector(mode="numeric",10)
for(i in 1:10)
{potential=(t(mkrchange[(55*i-55+28),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential5[i]=length(potential)}
potential10=vector(mode="numeric",10)
for(i in 1:10)
{potential=(t(mkrchange[(55*i-55+53),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential10[i]=length(potential)}
pds=seq(0.1,0.01,-0.01)
plot(pds,potential0,type='p',pch=22,col=rgb(0/255,0,0),bg=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0,type='l',pch=22,col=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential5,type='p',pch=23,col=rgb(110/255,0,0),bg=rgb(110/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential5,type='l',pch=23,col=rgb(110/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential10,type='p',pch=24,col=rgb(220/255,0,0),bg=rgb(220/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential10,type='l',pch=24,col=rgb(220/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0.01,0.1),ylim=c(50,400),xlab="",ylab="")
mtext("Replicative potential",side=2,line=2,cex=1.4)
mtext("Inhibitor activity",side=1,line=2,cex=1.3)
mtext("D",side=2,line=2.7,cex=2,at=(400*1),las=1)
legend(0.065,410,legend=c("M=10","M=5","M=0"),bty="n",x.intersp=0.2,y.intersp=0.5,cex=1.2,bg=c(rgb(220/255,0,0),rgb(110/255,0,0),rgb(0/255,0,0)),col=c(rgb(220/255,0,0),rgb(110/255,0,0),rgb(0/255,0,0)),lty=1,lwd=1.5,pch=c(17,18,15))
//7x8


Mpopulation=matrix(nrow=10,ncol=11)
rownames(Mpopulation)<-c("0.1","0.09","0.08","0.07","0.06","0.05","0.04","0.03","0.02","0.01")
colnames(Mpopulation)<-c("0","1","2","3","4","5","6","7","8","9","10")
for(i in 1:10)
{for(j in 1:11)
{
  celld=(t(mkrchange[(((i-1)*55+j*5)-2),]))
  celld=as.numeric(celld[!is.na(celld)])
  Mpopulation[i,j]=log10(max(celld))
}
}
rm(p)
data_mpopulation <- melt(Mpopulation)
p <- ggplot(data_mpopulation,aes(x=Var2,y=Var1))
p <- p + geom_tile(aes(fill=value))
p <- p + scale_fill_gradient(low = "white", high = "red")
p <- p + xlab("M") + ylab("kr") + theme_bw() + theme(panel.grid = element_blank()) + theme(legend.key=element_blank())
p + theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p + labs(title="C: Max cell population (log scale)",cex.lab=2.0)

heatmap.2(Mpopulation,col="bluered",Rowv=F,Colv=F,trace="none",cexCol=1,cexRow=1,xlab="M",ylab="kr",cex.lab=8.0,keysize=2.0,key.title="",key.xlab="Population(log10)",key.ylab="counts",lhei=c(1,3),lwid=c(1,3))



#kr=0
attach(mtcars)
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
par(mar=c(3.5,5.2,1,1))

par(mfrow=c(1,2),mar=c(4.2,4.2,1,1))
kr0 <- read.table("kr0mvar.txt",fill=TRUE,header=F)
kr0figc <- read.table("kr0figc.txt",fill=TRUE,header=F)

for(i in 1:4)
{population=(t(kr0[(i*20-7),]))
population=as.numeric(population)
population <- population[population>0]
population <- population[!is.na(population)]
population=log10(population)
pds=seq(1,length(population),1)
par(new=TRUE)
plot(pds,population,type='l',col=rgb((i*50-50)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),ylim=c(3,200),xlab="",ylab="",yaxt="n")}
par(new=TRUE)
population=(t(kr0figc[(1*10-7),]))
population=as.numeric(population)
population <- population[population>0]
population <- population[!is.na(population)]
population=log10(population)
pds=seq(1,length(population),1)
plot(pds,population,type='l',col=rgb((5*50-50)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),ylim=c(3,200),xlab="",ylab="",yaxt="n")

mtext("Population",side=2,line=2,cex=1.6)
mtext("PD",side=1,line=2,cex=1.6)
mtext("A",side=2,line=3,cex=2,at=(200*1),las=1)
axis(side=2,las=1,at=c(3,40,80,120,160,200),labels=c(expression(10^"3"),expression(10^"40"),expression(10^"80"),expression(10^"120"),expression(10^"160"),expression(10^"200")),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
legend("topleft",cex=1.2,legend=c("M=0.9","M=0.7","M=0.5","M=0.3","M=0.1"),x.intersp=0.2,y.intersp=1,bty="n",col=c(rgb((5*50-50)/255,0,0),rgb((4*50-50)/255,0,0),rgb((3*50-50)/255,0,0),rgb((2*50-50)/255,0,0),rgb((1*50-50)/255,0,0)),lty=1,lwd=1.5)

potential0.09=vector(mode="numeric",5)
for(i in 1:5)
{potential=(t(kr0[(i*20-17),]))
potential=as.numeric(potential)
potential <- potential[potential>0]
potential <- potential[!is.na(potential)]
potential0.09[i]=length(potential)}
pds=c(0,0.2,0.4,0.6,0.8)
plot(pds,potential0.09,type='p',pch=22,col=rgb(0/255,0,0),bg=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,1),ylim=c(50,400),xlab="",ylab="")
par(new=TRUE)
plot(pds,potential0.09,type='l',pch=22,col=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,1),ylim=c(50,400),xlab="",ylab="")
mtext("Replicative potential",side=2,line=2,cex=1.4)
mtext("Telomerase concentration",side=1,line=2,cex=1.4)
mtext("B",side=2,line=3,cex=2,at=(400*1),las=1)

population=(t(kr0figc[(1*10-7),]))
population=as.numeric(population)
population <- population[population>0]
population <- population[!is.na(population)]
population=log10(population)
pds=seq(1,length(population),1)
plot(pds,population,type='l',col=rgb((250)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),ylim=c(3,200),xlab="",ylab="",yaxt="n")
mtext("Population",side=2,line=2,cex=1.6)
mtext("PD",side=1,line=1.5,cex=1.6)
mtext("C",side=2,line=3,cex=2,at=(200*1),las=1)
axis(side=2,las=1,at=c(3,40,80,120,160,200),labels=c(expression(10^"3"),expression(10^"40"),expression(10^"80"),expression(10^"120"),expression(10^"160"),expression(10^"200")),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
legend("topleft",cex=1.2,legend=c("M=1"),x.intersp=0.2,bty="n",col=rgb((5*50)/255,0,0),lty=1,lwd=1.5)





guassianinitial <- read.table("guassianinitialdistri.txt",fill=TRUE,header=F)
uniforminitial <- read.table("uniforminitialdistri.txt",fill=TRUE,header=F)
xtelo=seq(0,500*31,500)
xq=seq(250,(500*500+250),500)
guassianinitial10=t(guassianinitial[7,])
guassianinitial10=as.numeric(guassianinitial10,na.rm=TRUE)/sum(as.numeric(guassianinitial10),na.rm=TRUE)
guassianinitial33=t(guassianinitial[9,])
guassianinitial33=as.numeric(guassianinitial33,na.rm=TRUE)/sum(as.numeric(guassianinitial33),na.rm=TRUE)
guassianinitial49=t(guassianinitial[11,])
guassianinitial49=as.numeric(guassianinitial49,na.rm=TRUE)/sum(as.numeric(guassianinitial49),na.rm=TRUE)
guassianinitial55=t(guassianinitial[13,])
guassianinitial55=as.numeric(guassianinitial55,na.rm=TRUE)/sum(as.numeric(guassianinitial55),na.rm=TRUE)
uniforminitial10=t(uniforminitial[7,])
uniforminitial10=as.numeric(uniforminitial10,na.rm=TRUE)/sum(as.numeric(uniforminitial10),na.rm=TRUE)
uniforminitial33=t(uniforminitial[9,])
uniforminitial33=as.numeric(uniforminitial33,na.rm=TRUE)/sum(as.numeric(uniforminitial33),na.rm=TRUE)
uniforminitial49=t(uniforminitial[11,])
uniforminitial49=as.numeric(uniforminitial49,na.rm=TRUE)/sum(as.numeric(uniforminitial49),na.rm=TRUE)
uniforminitial55=t(uniforminitial[13,])
uniforminitial55=as.numeric(uniforminitial55,na.rm=TRUE)/sum(as.numeric(uniforminitial55),na.rm=TRUE)
z1<-c(2,4,10,17,25,54,72,96,105,140,150,162.5,172,170,168,164,130,125,73,76,52,43,25,21,17,5,13,11,3,5,3,0)
z2=z1/sum(z1)
z3 <- c(23,46,78,88,133,135,112.5,124,162.5,115,135,112.5,77,88.5,72,56,42,19,19,19,9,13,4,5,2,3,0,1,0,2,1,0)
z4=z3/sum(z3)
z5=c(60,132,225,200,217,217,175,140,125,74,50,48,23,16,14,15,14,3,2,2,0,0,2,0,0,0,0,0,0,0,0,0)
z6=z5/sum(z5)
z7 <- c(122,305,410,375,310,265,220,180,120,100,55,55,30,20,15,8,4,4,0,0,4,2,2,2,0,0,0,0,0,0,0,0)
z8=z7/sum(z7)
par(mfrow=c(2,2),mar=c(4.2,4.6,1,1.5))
plot(xtelo,z2,type='l',cex=2,xlab="",tck=0.03,ylab="",cex.lab=1.5,las=1,cex.axis=1,ylim=c(0,0.1),xlim=c(0,16000),lwd=2,col=rgb(100/255,100/255,100/255))
par(new=TRUE)
polygon(xtelo,z2, density = NULL, border = NULL, col = rgb(240/255,240/255,240/255))
par(new=TRUE)
plot(xq,guassianinitial10,type='l',cex=1.5,xlab="",tck=0.03,ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.1),lwd=2,col="blue",xlim=c(0,16000))
mtext("Frequency",side=2,line=2.8,cex=1.4)
par(new=TRUE)
plot(xq,uniforminitial10,type='l',cex=1.5,xlab="",tck=0.03,ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.1),lwd=2,col="red",xlim=c(0,16000))
mtext("TL",side=1,line=2.5,cex=1.4)
text(0,0.09,labels="PD10",cex=1.5,pos=4)
mtext("A",side=2,line=3.2,cex=2,at=(0.1*1),las=1)
legend("topright",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",cex=1.0,lty=1,lwd=1.5)

plot(xtelo,z4,type='l',cex=2,xlab="",ylab="",tck=0.03,cex.lab=1.5,las=1,cex.axis=1,ylim=c(0,0.1),xlim=c(0,16000),lwd=2,col=rgb(100/255,100/255,100/255))
par(new=TRUE)
polygon(c(0,xtelo),c(0,z4), density = NULL, border = NULL, col = rgb(240/255,240/255,240/255))
par(new=TRUE)
plot(xq,guassianinitial33,type='l',cex=1.5,xlab="",tck=0.03,ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.1),lwd=2,col="blue",xlim=c(0,16000))
mtext("Frequency",side=2,line=2.8,cex=1.4)
par(new=TRUE)
plot(xq,uniforminitial33,type='l',cex=1.5,xlab="",tck=0.03,ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.1),lwd=2,col="red",xlim=c(0,16000))
mtext("TL",side=1,line=2.5,cex=1.4)
text(0,0.09,labels="PD33",cex=1.5,pos=4)
mtext("B",side=2,line=3.2,cex=2,at=(0.1*1),las=1)
legend("topright",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",cex=1.0,lty=1,lwd=1.5)

plot(xtelo,z6,type='l',cex=2,xlab="",ylab="",tck=0.03,cex.lab=1.5,las=1,cex.axis=1,ylim=c(0,0.15),xlim=c(0,16000),lwd=2,col=rgb(100/255,100/255,100/255))
par(new=TRUE)
polygon(c(0,xtelo),c(0,z6), density = NULL, border = NULL, col = rgb(240/255,240/255,240/255))
par(new=TRUE)
plot(xq,guassianinitial49,type='l',cex=1.5,tck=0.03,xlab="",ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.15),lwd=2,col="blue",xlim=c(0,16000))
mtext("Frequency",side=2,line=2.8,cex=1.4)
par(new=TRUE)
plot(xq,uniforminitial49,type='l',cex=1.5,tck=0.03,xlab="",ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.15),lwd=2,col="red",xlim=c(0,16000))
mtext("TL",side=1,line=2.5,cex=1.4)
text(0,0.09*1.5,labels="PD49",cex=1.5,pos=4)
mtext("C",side=2,line=3.2,cex=2,at=(0.15*1),las=1)
legend("topright",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",cex=1.0,lty=1,lwd=1.5)

plot(xtelo,z8,type='l',cex=2,xlab="",ylab="",tck=0.03,cex.lab=1.5,las=1,cex.axis=1,ylim=c(0,0.2),xlim=c(0,16000),lwd=2,col=rgb(100/255,100/255,100/255))
par(new=TRUE)
polygon(c(0,xtelo),c(0,z8), density = NULL, border = NULL, col = rgb(240/255,240/255,240/255))
par(new=TRUE)
plot(xq,guassianinitial55,type='l',cex=1.5,tck=0.03,xlab="",ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.2),lwd=2,col="blue",xlim=c(0,16000))
mtext("Frequency",side=2,line=2.8,cex=1.4)
par(new=TRUE)
plot(xq,uniforminitial55,type='l',cex=1.5,tck=0.03,xlab="",ylab="",las=1,cex.lab=1.5,cex.axis=1,ylim=c(0,0.2),lwd=2,col="red",xlim=c(0,16000))
mtext("TL",side=1,line=2.5,cex=1.4)
text(0,0.09*2,labels="PD55",cex=1.5,pos=4)
mtext("D",side=2,line=3.2,cex=2,at=(0.2*1),las=1)
legend("topright",legend=c("Gaussian","Uniform","Experiment"),col=c("blue","red","black"),bty="n",cex=1.0,lty=1,lwd=1.5)



//mutation division control
division1 <- read.table("division1.txt",fill=TRUE,header=F)
division0beta1 <- read.table("division0beta1.txt",fill=TRUE,header=F)
division0beta2 <- read.table("division0beta2.txt",fill=TRUE,header=F)
division0beta3 <- read.table("division0beta3.txt",fill=TRUE,header=F)

division1=(t(division1[5,]))
division0beta1=(t(division0beta1[5,]))
division0beta2=(t(division0beta2[5,]))
division0beta3=(t(division0beta3[5,]))
division1=as.numeric(division1)
division1 <- division1[division1>0]
division1 <- division1[!is.na(division1)]
division1x<-seq(1,length(division1),1)
division0beta1=as.numeric(division0beta1)
division0beta1 <- division0beta1[division0beta1>0]
division0beta1 <- division0beta1[!is.na(division0beta1)]
division0beta1x<-seq(1,length(division0beta1),1)
division0beta2=as.numeric(division0beta2)
division0beta2 <- division0beta2[division0beta2>0]
division0beta2 <- division0beta2[!is.na(division0beta2)]
division0beta2x<-seq(1,length(division0beta2),1)
division0beta3=as.numeric(division0beta3)
division0beta3 <- division0beta3[division0beta3>0]
division0beta3 <- division0beta3[!is.na(division0beta3)]
division0beta3x<-seq(1,length(division0beta3),1)
plot(division1x,division1,type='p',pch=4,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,8000),lwd=0.5,col="blue",xlim=c(0,500))
par(new=TRUE)
plot(division0beta1x,division0beta1,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,8000),lwd=0.5,col="red",xlim=c(0,500))
par(new=TRUE)
plot(division0beta2x,division0beta2,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,8000),lwd=0.5,col="green",xlim=c(0,500))
par(new=TRUE)
plot(division0beta3x,division0beta3,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,8000),lwd=0.5,col="darkgrey",xlim=c(0,500))
mtext("Mean TL",side=2,line=2,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)
legend("bottomright",legend=c("Division rate=1","theta origin","theta*0.001","theta*0.000001"),col=c("blue","red","green","darkgrey"),cex=c(1,1),pch=c(4,24,24,24))

division1=(t(division1[3,]))
division0beta1=(t(division0beta1[3,]))
division0beta2=(t(division0beta2[3,]))
division0beta3=(t(division0beta3[3,]))
division1=as.numeric(division1)
division1 <- division1[division1>0]
division1 <- division1[!is.na(division1)]
division1=log10(division1)
division1x<-seq(1,length(division1),1)
division0beta1=as.numeric(division0beta1)
division0beta1 <- division0beta1[division0beta1>0]
division0beta1 <- division0beta1[!is.na(division0beta1)]
division0beta1=log10(division0beta1)
division0beta1x<-seq(1,length(division0beta1),1)
division0beta2=as.numeric(division0beta2)
division0beta2 <- division0beta2[division0beta2>0]
division0beta2 <- division0beta2[!is.na(division0beta2)]
division0beta2=log10(division0beta2)
division0beta2x<-seq(1,length(division0beta2),1)
division0beta3=as.numeric(division0beta3)
division0beta3 <- division0beta3[division0beta3>0]
division0beta3 <- division0beta3[!is.na(division0beta3)]
division0beta3=log10(division0beta3)
division0beta3x<-seq(1,length(division0beta3),1)
plot(division1x,division1,type='p',pch=4,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,15),lwd=0.5,col="blue",xlim=c(0,500))
par(new=TRUE)
plot(division0beta1x,division0beta1,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,15),lwd=0.5,col="red",xlim=c(0,500))
par(new=TRUE)
plot(division0beta2x,division0beta2,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,15),lwd=0.5,col="green",xlim=c(0,500))
par(new=TRUE)
plot(division0beta3x,division0beta3,type='p',pch=24,cex=0.3,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,15),lwd=0.5,col="darkgrey",xlim=c(0,500))
mtext("Population",side=2,line=2,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)
legend("bottomright",legend=c("Division rate=1","theta origin","theta*0.001","theta*0.000001"),col=c("blue","red","green","darkgrey"),cex=c(1,1),pch=c(4,24,24,24))



//mutation sample
mutsample <- read.table("mutationsample.txt",fill=TRUE,header=F)
i=9
mutsamplereproduce=(t(mutsample[2+i*10,]))
mutsamplepop =(t(mutsample[3+i*10,]))
mutsamplemeantl =(t(mutsample[5+i*10,]))
mutsamplemutpropor1 =(t(mutsample[7+i*10,]))
mutsamplemutpropor2 =(t(mutsample[8+i*10,]))
mutsamplemutpropor31 =(t(mutsample[9+i*10,]))
mutsamplemutpropor32 =(t(mutsample[10+i*10,]))
mutsamplenormal =(t(mutsample[6+i*10,]))

par(mfrow=c(5,1),mar=c(4.2,4.2,1,1))
mutsamplereproduce=as.numeric(mutsamplereproduce)
mutsamplereproduce <- mutsamplereproduce[mutsamplereproduce>0]
mutsamplereproduce <- mutsamplereproduce[!is.na(mutsamplereproduce)]
mutsamplereproducex<-seq(1,length(mutsamplereproduce),1)
plot(mutsamplereproducex,mutsamplereproduce,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,2),lwd=0.5,col="blue",xlim=c(0,1000))
mtext("Reproduce rate",side=2,line=2.5,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)

mutsamplepop=as.numeric(mutsamplepop)
mutsamplepop <- mutsamplepop[mutsamplepop>0]
mutsamplepop <- mutsamplepop[!is.na(mutsamplepop)]
mutsamplepopx<-seq(1,length(mutsamplepop),1)
mutsamplepop=log10(mutsamplepop)
plot(mutsamplepopx,mutsamplepop,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,10),lwd=0.5,col="blue",xlim=c(0,1000))
mtext("population(log10)",side=2,line=2.5,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)

mutsamplemeantl=as.numeric(mutsamplemeantl)
mutsamplemeantl <- mutsamplemeantl[mutsamplemeantl>0]
mutsamplemeantl <- mutsamplemeantl[!is.na(mutsamplemeantl)]
mutsamplemeantlx<-seq(1,length(mutsamplemeantl),1)
plot(mutsamplemeantlx,mutsamplemeantl,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,8000),lwd=0.5,col="blue",xlim=c(0,1000))
mtext("Mean TL",side=2,line=2.5,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)

mutsamplemutpropor1=as.numeric(mutsamplemutpropor1)
mutsamplemutpropor1 <- mutsamplemutpropor1[!is.na(mutsamplemutpropor1)]
mutsamplemutpropor1x<-seq(1,length(mutsamplemutpropor1),1)
plot(mutsamplemutpropor1x,mutsamplemutpropor1,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="orange",xlim=c(0,1000))
mtext("Proportion of MC",side=2,line=2.5,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)
par(new=TRUE)
mutsamplemutpropor31=as.numeric(mutsamplemutpropor31)
mutsamplemutpropor31 <- mutsamplemutpropor31[!is.na(mutsamplemutpropor31)]
mutsamplemutpropor31x<-seq(1,length(mutsamplemutpropor31),1)
plot(mutsamplemutpropor31x,mutsamplemutpropor31,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="red",xlim=c(0,1000))
par(new=TRUE)
mutsamplenormal=as.numeric(mutsamplenormal)
mutsamplenormal <- mutsamplenormal[!is.na(mutsamplenormal)]
mutsamplenormalx<-seq(1,length(mutsamplenormal),1)
plot(mutsamplenormalx,mutsamplenormal,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="green",xlim=c(0,1000))



mutsamplemutpropor2=as.numeric(mutsamplemutpropor2)
mutsamplemutpropor2 <- mutsamplemutpropor2[!is.na(mutsamplemutpropor2)]
mutsamplemutpropor2x<-seq(1,length(mutsamplemutpropor2),1)
plot(mutsamplemutpropor2x,mutsamplemutpropor2,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="orange",xlim=c(0,1000))
mtext("Proportion of MC",side=2,line=2.5,cex=1.5)
mtext("PD",side=1,line=2.5,cex=1.5)
par(new=TRUE)
mutsamplemutpropor32=as.numeric(mutsamplemutpropor32)
mutsamplemutpropor32 <- mutsamplemutpropor32[!is.na(mutsamplemutpropor32)]
mutsamplemutpropor32x<-seq(1,length(mutsamplemutpropor32),1)
plot(mutsamplemutpropor32x,mutsamplemutpropor32,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="red",xlim=c(0,1000))
par(new=TRUE)
plot(mutsamplenormalx,mutsamplenormal,type='p',pch=24,cex=0.2,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),lwd=0.5,col="green",xlim=c(0,1000))



//mutation 
mut0004 <- read.table("mutation0004.txt",fill=TRUE,header=F)

mutall <- read.table("mutation200overall.txt",fill=TRUE,header=F)



mutsample=200
livemut=vector(mode="numeric",mutsample)
phenotype=vector(mode="character",mutsample)
finaltype=vector(mode="character",mutsample)
path=vector(mode="character",mutsample)
phenotype=as.numeric(phenotype)
finaltype=as.numeric(finaltype)
path=as.numeric(path)

mutlabel=9*2000


yz=0.8
for(i in 1:mutsample)
{
  if(as.numeric(mutall[i*10-7+mutlabel,1000])>0 & (as.numeric(mutall[i*10-4+mutlabel,1000])+as.numeric(mutall[i*10-3+mutlabel,1000])+as.numeric(mutall[i*10-2+mutlabel,1000]))<(1-yz))
  {livemut[i]=1}
  if(livemut[i]==1)
  {
    nc=as.numeric(t(mutall[i*10-4+mutlabel,]))
    cc1=as.numeric(t(mutall[i*10-1+mutlabel,]))
    cc2=as.numeric(t(mutall[i*10+mutlabel,]))
    pop=as.numeric(t(mutall[i*10-7+mutlabel,]))
    phenom=as.numeric(t(mutall[i*10-3+mutlabel,]))
    phenokr=as.numeric(t(mutall[i*10-2+mutlabel,]))
    cc=cc1+cc2
    phenostart=min(which(cc>(1-yz)))
    phenoend=min(which(cc>yz))
    sumphenom=sum(pop[phenostart:phenoend]*phenom[phenostart:phenoend])
    sumphenokr=sum(pop[phenostart:phenoend]*phenokr[phenostart:phenoend])
    phenotype[i]=sumphenom/(sumphenom+sumphenokr)
    finalm=as.numeric((mutall[i*10-1+mutlabel,1000]))
    finalkr=as.numeric((mutall[i*10+mutlabel,1000]))
    finaltotal=finalm+finalkr
    finalm=finalm/finaltotal
    finalkr=finalkr/finaltotal
    finaltype[i]=finalm
    if(phenotype[i]>0.5 & finaltype[i]>0.5)
    {path[i]=1}
    else if(phenotype[i]>0.5 & finaltype[i]<0.5)
    {path[i]=2}
    else if(phenotype[i]<0.5 & finaltype[i]>0.5)
    {path[i]=3}
    else
    {path[i]=4}
  }
}

ymutcount=hist(as.numeric(phenotype),col=rgb(245/256,245/256,245/256),breaks=20,main="Precancerous Phenotype",xlab="Proportion of M type(phenotype)",xlim=c(0,1),ylim=c(0,22))
ymutcount=ymutcount$counts
xmutcount=seq(0.05,1,0.05)

par(new=TRUE)
plot(xmutcount,ymutcount3,type='l',cex=0.2,xlab="proportion of phenotype m",ylab="Counts",cex.lab=1.5,cex.axis=1.5,lwd=1.5,col="black",xlim=c(0,1),ylim=c(0,22))

hist(as.numeric(finaltype),breaks=20,main="cancer type",xlab="Proportion of M type(finaltype)",xlim=c(0,1))


sum(livemut)/200
(length(which(path==1))+length(which(path==2)))/sum(livemut)
(length(which(path==3))+length(which(path==4)))/sum(livemut)
length(which(path==1))/(length(which(path==1))+length(which(path==2)))
length(which(path==2))/(length(which(path==1))+length(which(path==2)))
length(which(path==3))/(length(which(path==3))+length(which(path==4)))
length(which(path==4))/(length(which(path==3))+length(which(path==4)))

(length(which(path==1))+length(which(path==3)))/sum(livemut)
(length(which(path==2))+length(which(path==4)))/sum(livemut)


mutsample=200
livemut=vector(mode="numeric",mutsample)
phenotype=vector(mode="character",mutsample)
finaltype=vector(mode="character",mutsample)
path=vector(mode="character",mutsample)
phenotype=as.numeric(phenotype)
finaltype=as.numeric(finaltype)
path=as.numeric(path)

mutlabel=0*2000


yz=0.8
yzmut=0.5
for(i in 1:mutsample)
{
  if(as.numeric(mutall[i*10-7+mutlabel,1000])>0 & (as.numeric(mutall[i*10-4+mutlabel,1000])+as.numeric(mutall[i*10-3+mutlabel,1000])+as.numeric(mutall[i*10-2+mutlabel,1000]))<(1-yz))
  {livemut[i]=1}
  if(livemut[i]==1)
  {
    nc=as.numeric(t(mutall[i*10-4+mutlabel,]))
    cc1=as.numeric(t(mutall[i*10-1+mutlabel,]))
    cc2=as.numeric(t(mutall[i*10+mutlabel,]))
    pop=as.numeric(t(mutall[i*10-7+mutlabel,]))
    phenom=as.numeric(t(mutall[i*10-3+mutlabel,]))
    phenokr=as.numeric(t(mutall[i*10-2+mutlabel,]))
    promutm=cc1+cc2+phenom
    promutkr=cc1+cc2+phenokr
    mstart=min(which(promutm>(yzmut)))
    krstart=min(which(promutkr>(yzmut)))
    if(mstart<krstart)
    {path[i]=1}
    else 
    {path[i]=2}
  }
}

sum(livemut)/200
(length(which(path==1)))/sum(livemut)
(length(which(path==2)))/sum(livemut)

#fig332c
mutprobability=c(0,17,37,69,101,114,151,164,178,187,196)
cmutpro=seq(0,0.01,0.001)
plot(cmutpro,mutprobability,type='p',pch=22,col="black",bg=rgb(0/255,0,0),mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,0.01),ylim=c(0,200),xlab="",ylab="",yaxt="n")
par(new=TRUE)
plot(cmutpro,mutprobability,type='l',pch=22,col="black",mgp=c(0,0.2,0),tck=0.02,las=1,cex=1,xlim=c(0,0.01),ylim=c(0,200),xlab="",ylab="",yaxt="n")
mtext("Number of URP",side=2,line=2,cex=1.4)
mtext("Mutation probability",side=1,line=1.5,cex=1.4)
mtext("C",side=2,line=2.5,cex=2,at=(200*1.0),las=1)
axis(side=2,las=1,at=c(0,50,100,150,200),mgp=c(0,0.2,0),tck=0.03,las=1,labels=c("0","50","100","150","200"),cex.axis=1.0)

#fig332d
a=c(0,1)
b=c(0,1)
plot(a,b,type='l',col=rgb(1,1,1),xlim=c(0.5,4.5),ylim=c(0,100),mgp=c(0,0.2,0),tck=0.02,las=1,xlab="",ylab="",yaxt="n",xaxt="n",axes=FALSE)
axis(side=2,las=1,at=c(0,25,50,75,100),mgp=c(0,0.2,-0.5),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
axis(side=1,las=1,at=c(1,2,3,4),mgp=c(0,0.2,3),tck=0.03,las=1,labels=c("0.001","0.004","0.007","0.01"),cex.axis=1.0)
mtext("Mutation probability",side=1,line=1.5,cex=1.4)
mtext("Samples(%)",side=2,line=1.8,cex=1.4)
par(new=TRUE)
polygon(c(0.6,1.4,1.4,0.6),c(0,0,62,62), density = NULL, border = rgb(255/256,106/256,106/256), col = rgb(255/256,106/256,106/256))
par(new=TRUE)
polygon(c(1.6,2.4,2.4,1.6),c(0,0,64,64), density = NULL, border = rgb(255/256,106/256,106/256), col = rgb(255/256,106/256,106/256))
par(new=TRUE)
polygon(c(2.6,3.4,3.4,2.6),c(0,0,68,68), density = NULL, border = rgb(255/256,106/256,106/256), col = rgb(255/256,106/256,106/256))
par(new=TRUE)
polygon(c(3.6,4.4,4.4,3.6),c(0,0,71,71), density = NULL, border = rgb(255/256,106/256,106/256), col = rgb(255/256,106/256,106/256))
par(new=TRUE)
polygon(c(0.6,1.4,1.4,0.6),c(62,62,100,100), density = NULL, border = rgb(0/256,197/256,205/256), col = rgb(0/256,197/256,205/256))
par(new=TRUE)
polygon(c(1.6,2.4,2.4,1.6),c(64,64,100,100), density = NULL, border = rgb(0/256,197/256,205/256), col = rgb(0/256,197/256,205/256))
par(new=TRUE)
polygon(c(2.6,3.4,3.4,2.6),c(68,68,100,100), density = NULL, border = rgb(0/256,197/256,205/256), col = rgb(0/256,197/256,205/256))
par(new=TRUE)
polygon(c(3.6,4.4,4.4,3.6),c(71,71,100,100), density = NULL, border = rgb(0/256,197/256,205/256), col = rgb(0/256,197/256,205/256))




mutall <- read.table("mutation200overall.txt",fill=TRUE,header=F)

specialmutation=vector(mode="numeric")


//judgement of special mutation
j=1
for(z in 1:2000)
{
i=z*10-10
precancercell=as.numeric(t(mutall[7+i,]))+as.numeric(t(mutall[8+i,]))
cancercell=as.numeric(t(mutall[9+i,]))+as.numeric(t(mutall[10+i,]))
if((max(which(precancercell>0.3))-min(which(precancercell>0.3)))>300)
{specialmutation[j]=z
j=j+1}
}

specipopx=seq(0,1000,1)
i=815
specipop=as.numeric(t(mutall[3+i*10-10,]))
plot(specipopx,specipop,type='l',ylim=c(0,3*10^7),col="darkred",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="",yaxt="n")
par(new=TRUE)
#polygon(c(60,360,360,60),c(-4*10^6,-4*10^6,4*10^7,4*10^7), density = NULL, border = rgb(248/256,248/256,248/256), col = rgb(248/256,248/256,248/256))
#par(new=TRUE)
#polygon(c(360,1200,1200,360),c(-4*10^6,-4*10^6,4*10^7,4*10^7), density = NULL, border = rgb(240/256,240/256,240/256), col = rgb(240/256,240/256,240/256))
#par(new=TRUE)
specipop=as.numeric(t(mutall[3+i*10-10,]))
plot(specipopx,specipop,type='l',ylim=c(0,3*10^7),col="darkred",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="",yaxt="n")
t=0
for(q in c(1,2,0,3))
{ 
  par(new=TRUE)
  i=1+q*350
  specipop=as.numeric(t(mutall[3+i*10-10,]))
  specipop=specipop[specipop>0]
  specipopx=seq(1,length(specipop),1)
  plot(specipopx,specipop,type='l',ylim=c(0,3*10^7),col=rgb((0)/256,(250-t*50)/256,(0)/256),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="",yaxt="n")
t=t+1
} 

mtext("Population",side=2,line=2.7,cex=1.6)
mtext("PD",side=1,line=1.7,cex=1.6)
mtext("A",side=2,line=3,cex=2,at=(3*10^7*1),las=1)
axis(side=2,las=1,at=c(0,1*10^7,2*10^7,3*10^7),labels=c("0",expression(10^"7"),expression(2%*%10^"7"),expression(3%*%10^"7")),mgp=c(0,0.2,0),tck=0.03,las=1,cex.axis=1.0)
#4*6
#legend("topleft",legend=c("Tumorigenesis","Normal"),cex=1.0,x.intersp=0.2,bty="n",col=c("darkred","darkgreen"),lty=1,lwd=1.5)

i=815
specipopx=seq(0,1000,1)
precancercell=as.numeric(t(mutall[7+i*10-10,]))+as.numeric(t(mutall[8+i*10-10,]))
plot(specipopx,precancercell,type='l',ylim=c(0,1),col="orange",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="")
#par(new=TRUE)
#polygon(c(60,360,360,60),c(-4,-4,4,4), density = NULL, border = rgb(248/256,248/256,248/256), col = rgb(248/256,248/256,248/256))
#par(new=TRUE)
#polygon(c(360,1200,1200,360),c(-4,-4,4,4), density = NULL, border = rgb(240/256,240/256,240/256), col = rgb(240/256,240/256,240/256))
par(new=TRUE)
specipopx=seq(0,1000,1)
precancercell=as.numeric(t(mutall[7+i*10-10,]))+as.numeric(t(mutall[8+i*10-10,]))
plot(specipopx,precancercell,type='l',ylim=c(0,1),col="orange",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="")
par(new=TRUE)
cancercell=as.numeric(t(mutall[9+i*10-10,]))+as.numeric(t(mutall[10+i*10-10,]))
plot(specipopx,cancercell,type='l',ylim=c(0,1),col="red",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="")
par(new=TRUE)
normalcell=as.numeric(t(mutall[6+i*10-10,]))
plot(specipopx,normalcell,type='l',ylim=c(0,1),col="green",mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1.5,xlim=c(0,1000),xlab="",ylab="")
mtext("Proportion",side=2,line=1.7,cex=1.6)
mtext("PD",side=1,line=1.7,cex=1.6)
mtext("B",side=2,line=2.5,cex=2,at=(1),las=1)
legend(600,0.8,legend=c("Double mutation cell","Single mutation cell","Normal cell"),cex=1.0,x.intersp=0,y.intersp=2,bty="n",col=c("red","orange","green"),lty=1,lwd=1.5)
