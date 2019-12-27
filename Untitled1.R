setwd("/Users/yeyusong/Desktop/在投/Telomere/txt数据")
unifsample <- read.table("uniforminitialsample1.txt",fill=TRUE,header=F)
unifsamplepopulation=t(unifsample[3,])
x=seq(1,101,1)
plot(x,unifsamplepopulation,log='y')
plot(x,unifsamplepopulation,xlim=c(0,100),ylim=c(1,10^10),yaxt="n",xlab="Generation(PD)",ylab="",log="y",col="darkgrey",type='l',lwd=3,cex.axis=1.5,cex.lab=2.0)
axis(side=2,at=c(1,10^3,10^6,10^9),labels=c(expression(10^"0"),expression(10^"3"),expression(10^"6"),expression(10^"9")),cex.axis=1.5)
mtext("Cell population",side=2,cex=2,line=2.5)

unifsampledistri10=t(unifsample[6,])
unifsampledistri10=as.numeric(unifsampledistri10,na.rm=TRUE)/sum(as.numeric(unifsampledistri10),na.rm=TRUE)
unifsampledistri33=t(unifsample[8,])
unifsampledistri33=as.numeric(unifsampledistri33,na.rm=TRUE)/sum(as.numeric(unifsampledistri33),na.rm=TRUE)
unifsampledistri49=t(unifsample[10,])
unifsampledistri49=as.numeric(unifsampledistri49,na.rm=TRUE)/sum(as.numeric(unifsampledistri49),na.rm=TRUE)
unifsampledistri55=t(unifsample[12,])
unifsampledistri55=as.numeric(unifsampledistri55,na.rm=TRUE)/sum(as.numeric(unifsampledistri55),na.rm=TRUE)

xtelo=seq(500,500*32,500)
xq=seq(500,500*101,500)
z1<-c(2,4,10,17,25,54,72,96,105,140,150,162.5,172,170,168,164,130,125,73,76,52,43,25,21,17,5,13,11,3,5,3,0)
z2=z1/sum(z1)
z3 <- c(23,46,78,88,133,135,112.5,124,162.5,115,135,112.5,77,88.5,72,56,42,19,19,19,9,13,4,5,2,3,0,1,0,2,1,0)
z4=z3/sum(z3)
z5=c(60,132,225,200,217,217,175,140,125,74,50,48,23,16,14,15,14,3,2,2,0,0,2,0,0,0,0,0,0,0,0,0)
z6=z5/sum(z5)
z7 <- c(122,305,410,375,310,265,220,180,120,100,55,55,30,20,15,8,4,4,0,0,4,2,2,2,0,0,0,0,0,0,0,0)
z8=z7/sum(z7)
par(mfrow=c(2,2),mar=c(5,5,2,3))
plot(xtelo,z2,type='s',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),xlim=c(0,16000),lwd=2,col="black")
par(new=TRUE)
plot(xtelo,z2,type='h',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),lwd=2,col="black",xlim=c(0,16000))
par(new=TRUE)
plot(xq,unifsampledistri10,type='l',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),lwd=3,col="red",xlim=c(0,16000))
mtext("Frequency",side=2,line=3,cex=1.6)
text(14000,0.08,labels="A",cex=3.0,pos=4)
plot(xtelo,z4,type='s',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),xlim=c(0,16000),lwd=2,col="black")
par(new=TRUE)
plot(xtelo,z4,type='h',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),lwd=2,col="black",xlim=c(0,16000))
par(new=TRUE)
plot(xq,unifsampledistri33,type='l',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.1),lwd=3,col="red",xlim=c(0,16000))
mtext("Frequency",side=2,line=3,cex=1.6)
text(14000,0.08,labels="B",cex=3.0,pos=4)
plot(xtelo,z6,type='s',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.15),xlim=c(0,16000),lwd=2,col="black")
par(new=TRUE)
plot(xtelo,z6,type='h',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.15),lwd=2,col="black",xlim=c(0,16000))
par(new=TRUE)
plot(xq,unifsampledistri49,type='l',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.15),lwd=3,col="red",xlim=c(0,16000))
mtext("Frequency",side=2,line=3,cex=1.6)
text(14000,0.12,labels="C",cex=3.0,pos=4)
plot(xtelo,z8,type='s',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.2),xlim=c(0,16000),lwd=2,col="black")
par(new=TRUE)
plot(xtelo,z8,type='h',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.2),lwd=2,col="black",xlim=c(0,16000))
par(new=TRUE)
plot(xq,unifsampledistri55,type='l',cex=2,xlab="Telomere length",ylab="",cex.lab=2.0,cex.axis=1.5,ylim=c(0,0.2),lwd=3,col="red",xlim=c(0,16000))
mtext("Frequency",side=2,line=3,cex=1.6)
text(14000,0.16,labels="D",cex=3.0,pos=4)

ezmutsample1 <- read.table("ezmutatingbili1.txt",fill=TRUE,header=F)
ezmutsampleplot1=t(ezmutsample1[4,])
ezmutsample2 <- read.table("ezmutatingbili2.txt",fill=TRUE,header=F)
ezmutsampleplot2=t(ezmutsample2[4,])
ezmutsample3 <- read.table("ezmutatingbili3.txt",fill=TRUE,header=F)
ezmutsampleplot3=t(ezmutsample3[4,])
ezmutsample4 <- read.table("ezmutatingbili4.txt",fill=TRUE,header=F)
ezmutsampleplot4=t(ezmutsample4[4,])

pdmut=seq(1,301,1)
plot(pdmut,ezmutsampleplot1,cex=0.3,col="grey",type='p',ylim=c(0,1),ylab="Mutating proportion",xlab="PD")
par(new=TRUE)
plot(pdmut,ezmutsampleplot2,cex=0.3,col="blue",type='p',ylim=c(0,1),ylab="",xlab="PD")
par(new=TRUE)
plot(pdmut,ezmutsampleplot3,cex=0.3,col="red",type='p',ylim=c(0,1),ylab="",xlab="PD")
par(new=TRUE)
plot(pdmut,ezmutsampleplot4,cex=0.3,col="darkgreen",type='p',ylim=c(0,1),ylab="",xlab="PD")

result1 <- read.table("result1.txt",fill=TRUE,header=F)
beta <- t(result1[2,])
popula <- t(result1[3,])
result=matrix(nrow=2,ncol=501)
result[1,] <- beta
result[2,] <- popula
write.table (result, file ="matrix", sep ="\t", row.names =F, col.names =F, quote =F)

lilun <- read.table("result1.txt",fill=TRUE,header=F)
plive <- t(lilun[2,])
plive=as.numeric(plive)
plive <- plive[plive>0]
plive <- plive[!is.na(plive)]
pdeath=vector(mode="numeric",length(plive))
for(i in 1:length(plive))
{pdeath[i]=1}
for(i in 1:length(plive))
{
  for(j in 1:i)
  {
    if(j==i)
    {
      if(plive[j]>1)
      {pli=0}
      else
      {pli=1-plive[j]}
    }
    else
    {
    if(plive[j]>1)
      {pli=1}
    else
    {pli=plive[j]}
    }
    pdeath[i]=pdeath[i]*pli
  }
}
pdeathx<-seq(1,length(pdeath),1)
plot(pdeathx,pdeath,col="blue")
par(new=TRUE)
sum(pdeathx*pdeath)

pmut0 <- read.table("pmut0.txt",fill=TRUE,header=F)
pmut000015 <- read.table("pmut000015.txt",fill=TRUE,header=F)
pmut000018 <- read.table("pmut000018.txt",fill=TRUE,header=F)
pmut00002 <- read.table("pmut00002.txt",fill=TRUE,header=F)

par(mfrow=c(3,1),mar=c(3,4,1,1))
pmut000018pop=(t(pmut000018[3,]))
pmut000018pop=as.numeric(pmut000018pop)
pmut000018pop <- pmut000018pop[pmut000018pop>0]
pmut000018pop <- pmut000018pop[!is.na(pmut000018pop)]
pmut000018pop=log10(pmut000018pop)
pmut000018popx <-seq(1,length(pmut000018pop),1)
plot(pmut000018popx,pmut000018pop,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(2,8))
polygon(c(25,25,123,123),c(1,9,9,1), density = NULL, border = rgb(210/255,255/255,255/255), col = rgb(210/255,255/255,255/255))
polygon(c(123,123,330,330),c(1,9,9,1), density = NULL, border = rgb(255/255,210/255,210/255), col = rgb(255/255,210/255,210/255))
par(new=TRUE)
plot(pmut000018popx,pmut000018pop,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(2,8))
mtext("Cell population",side=2,line=1.5,cex=1.0)
mtext("PD",side=1,line=1.5,cex=1.0)
par(new=TRUE)
pmut000018pop=(t(pmut0[3,]))
pmut000018pop=as.numeric(pmut000018pop)
pmut000018pop <- pmut000018pop[pmut000018pop>0]
pmut000018pop <- pmut000018pop[!is.na(pmut000018pop)]
pmut000018pop=log10(pmut000018pop)
pmut000018popx <-seq(1,length(pmut000018pop),1)
plot(pmut000018popx,pmut000018pop,type='l',lty=2,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="blue",xlim=c(0,300),ylim=c(2,8))
axis(side=2,las=1,at=c(2,4,6,8),labels=c(expression(10^"2"),expression(10^"4"),expression(10^"6"),expression(10^"8")),mgp=c(0,0.2,0),tck=0.02,las=1,cex.axis=1.0)
legend("bottomright",bg=c(rgb(245/256,245/256,245/256)),legend=c("Mutation","Normal"),col=c("red","blue"),cex=c(1,1),lty=c(1,2))


pmut000018tl=(t(pmut000018[5,]))
pmut000018tl=as.numeric(pmut000018tl)
pmut000018tl <- pmut000018tl[pmut000018tl>0]
pmut000018tl <- pmut000018tl[!is.na(pmut000018tl)]
pmut000018tlx <-seq(1,length(pmut000018tl),1)
plot(pmut000018tlx,pmut000018tl,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(0,8000))
polygon(c(25,25,123,123),c(-1000,9000,9000,-1000), density = NULL, border = rgb(210/255,255/255,255/255), col = rgb(210/255,255/255,255/255))
polygon(c(123,123,330,330),c(-1000,9000,9000,-1000), density = NULL, border = rgb(255/255,210/255,210/255), col = rgb(255/255,210/255,210/255))
par(new=TRUE)
plot(pmut000018tlx,pmut000018tl,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(0,8000))
mtext("Mean TL",side=2,line=2.0,cex=1.0)
mtext("PD",side=1,line=1.5,cex=1.0)
par(new=TRUE)
pmut000018tl=(t(pmut0[5,]))
pmut000018tl=as.numeric(pmut000018tl)
pmut000018tl <- pmut000018tl[pmut000018tl>0]
pmut000018tl <- pmut000018tl[!is.na(pmut000018tl)]
pmut000018tlx <-seq(1,length(pmut000018tl),1)
plot(pmut000018tlx,pmut000018tl,type='l',lty=2,cex=0.3,xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="blue",xlim=c(0,300),ylim=c(0,8000))
legend("bottomright",bg=c(rgb(245/256,245/256,245/256)),legend=c("Mutation","Normal"),col=c("red","blue"),cex=c(1,1),lty=c(1,2))


pmut000018bili=(t(pmut000018[6,]))
pmut000018bili=as.numeric(pmut000018bili)
pmut000018bili <- pmut000018bili[pmut000018bili>0]
pmut000018bili <- pmut000018bili[!is.na(pmut000018bili)]
pmut000018bilix <-seq(1,length(pmut000018bili),1)
plot(pmut000018bilix,pmut000018bili,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(0,1))
polygon(c(25,25,123,123),c(-0.1,1.1,1.1,-0.1), density = NULL, border = rgb(210/255,255/255,255/255), col = rgb(210/255,255/255,255/255))
polygon(c(123,123,330,330),c(-0.1,1.1,1.1,-0.1), density = NULL, border = rgb(255/255,210/255,210/255), col = rgb(255/255,210/255,210/255))
par(new=TRUE)
plot(pmut000018bilix,pmut000018bili,type='l',lty=1,cex=0.3,yaxt="n",xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="red",xlim=c(0,300),ylim=c(0,1))
mtext("Mutated cells",side=2,line=2.0,cex=1.0)
mtext("PD",side=1,line=1.5,cex=1.0)
par(new=TRUE)
pmut000018bili=vector(mode="numeric",length=100)
for(i in 1:100)
{
  pmut000018bili[i]<-0
}
pmut000018bilix <-seq(1,100,1)
plot(pmut000018bilix,pmut000018bili,type='l',lty=2,cex=0.3,xlab="",mgp=c(0,0.2,0),tck=0.02,las=1,ylab="",cex.lab=1.5,cex.axis=1.0,lwd=1.5,col="blue",xlim=c(0,300),ylim=c(0,1))
legend("bottomright",bg=c(rgb(245/256,245/256,245/256)),legend=c("Mutation","Normal"),col=c("red","blue"),cex=c(1,1),lty=c(1,2))
