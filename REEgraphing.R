
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

# pdf(file=(paste(Rlocal,"/MMSD_virus/ree2011_w_Eu.pdf",sep="")))
# 
# for (i in 1:length(Dates)){
# dfsub <- subset(ree2011,Sample.Date==Dates[i])
# reeplot <- xyplot(NASC.normalized.concentration~Element|Station.Name,data=dfsub,
#                   main=Dates[i], 
#                   scales=list(x=list(rot=90),
#                               y=list(log=TRUE),
#                               alternating=1),
#                   strip=strip.custom(par.strip.text=list(cex=0.65)))
# print(reeplot)
# }
# 
# dev.off()

pdf(file=(paste(Rlocal,"/MMSD_virus/ree2011_w_Eu.pdf",sep="")))

for (i in 1:length(Dates)){
dfsub <- subset(ree2011,Sample.Date==Dates[i])
reeplot <- xyplot(NASC.normalized.concentration~Element|Station.Name,data=dfsub,
                  main=Dates[i], 
                  ylim=c(0.00000001,0.0001),#range(ree2011$NASC.normalized.concentration, na.rm=TRUE),
                  scales=list(x=list(rot=90),
                              y=list(log=TRUE),
                              alternating=1),
                  strip=strip.custom(par.strip.text=list(cex=0.65)),
)
print(reeplot)
}

dev.off()

# Assign index numbers starting at Ce = 1 and ending with Dy = 14 for all samples
ree2011$index <- NA
for(i in 1:14) ree2011$index[which (ree2011$Element==REEaw[i,3])] <- i
  
  
ree2011sub <- subset(ree2011,Element!="La")
ree2011sub <- subset(ree2011sub,Element!="Ce")
ree2011sub <- subset(ree2011sub,Element!="Eu")
#ree2011sub <- subset(ree2011sub,Element!="Gd")
ree2011sub <- subset(ree2011sub,Element!="Lu")


i=2
#for (i in 1:length(Dates))
  GAsub <- subset(ree2011sub, Sample.Date==Dates[i])
  GAsub <- drop.unused.factors(GAsub)
  Stations <- levels(GAsub$Station.Name)
j=1
#for (j in 1:length(Stations)){
  GAsub <- subset(GAsub,Station.Name==Stations[j])
  Gdnorm <- GAsub[which(GAsub$Element=="Gd"),"NASC.normalized.concentration"]
  GAsub <- subset(GAsub,Element!="Gd")
  y <- GAsub$NASC.normalized.concentration
  x <- GAsub$index
  x2 <- x^2
  x3 <- x^3
  poly3 <- lm(y~x+x2+x3)
  ly <- log10(y)
  lpoly3 <- lm(ly~x+x2+x3)
summary(poly3)
if(ree2011$Element== Gdnorm <- 
p3coef <- poly3$coefficients
   Gdpred <- (p3coef[1] + p3coef[2]*7 + p3coef[3]*(7^2) + p3coef[4]*(7^3))
   GDanom <- Gdpred/Gdnorm
  predict(poly3,newdata=7)
   plot(poly3)

summary(lpoly3)

   plot(predict(poly3)~x,log="y",ylim=c((10^(-7)),(10^(-4.5))))
   points(y~x,,col="red")
   points(Gdnorm~c(7),col="blue")
