
###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Rlocal <- "D:/srcldata/R" #Laptop
       }else {Rlocal <- "//igsarmewwssrc/SRCdata/R"} #Network
source(paste(Rlocal,"/Scripts/fxn_drop.unused.factors.R",sep=""))

ree2011 <- read.delim(paste(Rlocal,"/MMSD_virus/smplstatus_20101219_20110418.txt",sep=""))
reeAW <- read.delim(paste(Rlocal,"/MMSD_virus/AtomicWeights.txt",sep=""))

# source(paste(Rlocal,"/Scripts/fxn_categorize.R",sep=""))
# 
# ree2011 <- categorize(Groupdat=reeAW,Datafile=ree2011,
#                    varsdata="Parameter",
#                    varsgroups="Element",
#                    catname="ElementAbbrev")

#ree2011$Sample.Date <- as.factor(ree2011$Sample.Date)
Dates <- levels(ree2011$Sample.Date)

reeAW$ElementAbbrev <- factor(reeAW$ElementAbbrev,as.character(reeAW$ElementAbbrev))

ree2011$Element <- factor(ree2011$Element,as.character(reeAW$ElementAbbrev))
ree2011$NASC_WO_Eu <- ifelse(ree2011$Element == "Eu", NA,ree2011$NASC.normalized.concentration)

#assign index numbers starting at Ce = 1 and ending with Dy = 14 for all samples
ree2011$index <- NA
for(i in 1:14) ree2011$index[which (ree2011$Element==REEaw[i,3])] <- i

pdf(file=(paste(Rlocal,"/MMSD_virus/ree2011_base.pdf",sep="")))
for (i in 1:length(Dates)){
  dfsub <- subset(ree2011, Sample.Date==Dates[i])
  dfsub[,"Station.Name"] <- drop.unused.factors(as.data.frame(dfsub[,"Station.Name"]))
  Stations <- levels(dfsub$Station.Name)
  par (mfrow=c(2,3))
#j=1
for (j in 1:length(Stations)){
 dfsub2 <- subset(dfsub,Station.Name==Stations[j])
  dfsub2[which(dfsub2[,"NASC.normalized.concentration"]==0 | 
    dfsub2[,"NASC.normalized.concentration"]<0),"NASC.normalized.concentration"] <- NA
reeplot <- plot(NASC.normalized.concentration~index,data=dfsub2,
                  main=Stations[j], cex.main=0.6,
                  ylim=c(0.00000001,0.0001),
                  log="y",xaxt="n",pch=1)
  axis(side=1,at=1:14,labels=levels(dfsub2$Element),cex=0.5,las=2)
  if (j == 2) mtext(dfsub2$Sample.Date[1],side=3,line=3)
  m1 <- lm(NASC.normalized.concentration~index,data=dfsub3)
  Gdest <- predict(m1,newdata=data.frame(index=7))#dfsub2[which(dfsub2$Element=="Gd"),"index"])
 par(new="TRUE")
 plot(x=predict(m1,newdata=data.frame(index=1:14)),
      ylim=c(0.00000001,0.0001),ylog=T,col="blue",
      xaxt="n",yaxt="n",ylab="",type="l")

print(reeplot)
}
}
dev.off()

dfsub3 <- dfsub2[-which(dfsub2$Element=="Gd"),]
dfsub3 <- dfsub2[-which(dfsub2$Element=="Eu"),]

m1 <- lm(NASC.normalized.concentration~index,data=dfsub3)
Gdest <- predict(m1,newdata=data.frame(index=7))#dfsub2[which(dfsub2$Element=="Gd"),"index"])
Gdanom <- 
  lines(x=predict(m1))
lines(x=predict(m1,newdata=data.frame(index=1:7)),col=1:7)

points(dfsub2$NASC.normalized.concentration~dfsub2$index)
plot(NASC.normalized.concentration~index,col=1:14,log="y",data=dfsub2)
