
############################################################################################

##-----PLOT with water levels?-----------------
library(xts)
library(lubridate)
library(RColorBrewer)
#read files
fil<-list.files("C:/Users/corjanWW/Desktop/RWS_langsdam/Waterstanden/from_waterinfo_rws",full.names = TRUE)
wat.lev<-lapply(fil, read.csv,header=TRUE,sep=";")

wat.lev<-lapply(wat.lev, function(x) x[,c("WAARNEMINGDATUM","REFERENTIE","NUMERIEKEWAARDE")])

#function  to xts
toxts<-function(x){
  d<-paste(x[,1], x[,2],sep=" ")
  d<-x[,1]
  date<-as.POSIXct(d, format = "%d-%m-%Y")# %H:%M:%S
  lev.xts<-xts(x[,3],date)
  return(lev.xts)
}
wat.lev.df<-do.call(rbind,wat.lev)

lev.xts<-toxts(wat.lev.df)
lev.xts<-lev.xts[lev.xts<1000]

#daily water levels
end.pt<-endpoints(lev.xts, "days")
lev.day<-period.apply(lev.xts, end.pt, mean)

dts<-seq(from=2009, to=2021,by=1)
st<- paste0(dts,"-07-01")
en<- c(paste0(dts,"-06-30")[-1],"2022-06-30")
summ.yr<-mapply(function(x,y) window(lev.day, start = as.Date(x), end = as.Date(y)), st,en)
summ.yr.vec<-lapply(summ.yr, function(x) as.numeric(x))
lev.yr.df<-do.call(qpcR:::cbind.na,summ.yr.vec)



##DIscharge - Daily
fil<-list.files("C:/Users/corjanWW/Desktop/RWS_langsdam/Waterstanden/discharge_TielWaal",full.names = TRUE)
Q<-lapply(fil, read.csv,header=TRUE,sep=";")
Q<-lapply(Q, function(x) x[,c("Date","Discharge")])

Q.df<-do.call(rbind,Q)
d<-Q.df[,1]
date<-as.POSIXct(d, format = "%d-%m-%Y")# %H:%M:%S
Q.xts<-xts(Q.df[,2],date)
st<- paste0(dts,"-07-01")
en<- c(paste0(dts,"-06-30")[-1],"2022-06-30")
Q.yr<-mapply(function(x,y) window(Q.xts, start = as.Date(x), end = as.Date(y)), st,en)
Q.yr.vec<-lapply(Q.yr, function(x) as.numeric(x))
Q.yr.df<-do.call(qpcR:::cbind.na,Q.yr.vec)



#Boxplot---
#pdf("./boxplot_Q_lmove_25cmpix.pdf",  width = 6, height =5.5) # The height of the plot in inches 
par(oma=c(1,0,1,1),mfrow=c(2,1),mar=c(1,3,1,1),mgp=c(1.65,0.3,0),tck=-0.01)
#water level boxplot
# layout.show(5)
#water level
#boxplot(lev.yr.df,ylab="Daily Water level Waal (cm +NAP)",pch=20,xaxt="n")

#discharge
boxplot(lev.yr.df,ylab="Daily discharge (m3/s)",pch=20,xaxt="n")

#mid diff- all already annual!
boxplot(t(mid.diff),xaxt="n",ylab=("Midslope lateral movement (m)"),pch=20)
axis(side=1,at=1:4,labels=as.character(c(2016,2018:2020)))
#dev.off()



#---as hydrograph---use dischagre
#easier to create an xts for plotting!
dt.seq<-seq(as.Date("2015-07-01"), as.Date("2016-06-30"), "days")
Q.yr.xts<-xts(Q.yr.df,dt.seq)


colnew= brewer.pal(8,'Accent')[c(1,2,3,5,6)]
cols = brewer.pal(8,'Dark2')[c(5,6)]
#colsg =brewer.pal(9,'Greys')[3:8]
cols2<-c(col1,cols)
dt.st<-
  #colnew<-as.character(yarrr::piratepal(palette = "basel",trans =0.1)[c(1,3,4,5,6,8,10)])
  lwd.seq<-c(1,1,4,1,1,4,4)/1.5
lty.seq<-c(2,2,1,2,2,1,1)
#plot
#pdf("./hydrograph_Q_lmove_25cmpix2.pdf",  width = 10, height =4) # The height of the plot in inches 
par(mar=c(4,4,2,2), mgp=c(1.7,0.57,0), tck=-0.02)
plot(as.zoo(Q.yr.xts[,7:13]),screens=1,col=colnew, lwd=lwd.seq,lty=1,xlab="Months", ylab="Average daily discharge\n (m3/s)")
#legend

legend("topright", legend=c(2015:2021), lwd=lwd.seq, seg.len = 2, bty="n", col=colnew, ncol=2, cex=0.8)
#legend("topright", legend=c(2015:2021), lwd=2, seg.len = 2, bty="n",col=colsnew)
#dev.off()



##Water level plots
#easier to create an xts for plotting!
dt.seq<-seq(as.Date("2015-07-01"), as.Date("2016-06-30"), "days")
lev.yr.xts<-xts(lev.yr.df,dt.seq)

cols = brewer.pal(8,'Accent')[c(1,2,3,5,6)]
#plot
#pdf("./hydrograph_watlev_lmove_25cmpix2.pdf",  width = 10, height =4) # The height of the plot in inches 
par(mar=c(4,4,2,2), mgp=c(1.7,0.57,0), tck=-0.02)
plot(as.zoo(lev.yr.xts[,7:13]),screens=1,col=colnew, lwd=lwd.seq,lty=1,xlab="Months", ylab="Waal water level\n (cm +NAP)",bty="n")
box(lwd=0.5)
abline(h=750, lwd=0.5, lty=4,col='grey66')
#legend
#legend("topright", legend=c(2015:2021), lwd=2, seg.len = 2, bty="n", col=colnew, ncol=3, cex=0.8)
#dev.off()



