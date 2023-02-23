library(raster)
library(sp)
library(rgdal)
library(fields)
library(RColorBrewer)
library(rgeos)
library(graphics)
library(imputeTS)

###Read files###
#run DTM_bank_changes2.R first for needed files.

#-------------read .asc langsdam-----------------------------------
ld<-read.table("C:/Users/corjanWW/Desktop/RWS_langsdam/Bathy/160411 Wamel Toplaag asbuilt.asc",sep=";",dec= ",")
colnames(ld)<-c("X","Y","Z")
coordinates(ld)<- ~X+Y
gridded(ld) <- TRUE
as<- as(ld, "SpatialGridDataFrame")

ld<-raster(ld,"Z")
crs<-CRS("+proj=sterea +lat_0=52.1561605555556 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs ")
proj4string(ld)<-crs
ld

#ld_Wamel<-writeRaster(ld,"LD_Wamel.tif",format="GTiff")
##change resol
ld<- disaggregate(ld, fact=4) #to 25cm
ld

#ld[ld<=0]<-NA
plot(ld)


#distance from ld
ld_dist<-distance(ld)
ld_dist2<-round(ld_dist,2)

#----------------------X-section based on transect-----------------------##
#read transects:
tr<-readOGR("./shp/ext_fc","transects_5m_200")
#to list
tr.list<-list()
for(i in 1:nrow(tr)){
  tr.list[[i]]<-tr[i,]
}


t2<-tnew[7:11] #post-LD DEMS
tdis<-lapply(t2, disaggregate, fact=2)# to 25cm
tdis.st<-stack(tdis)



##-------Extract values to transect points - takes a very long time!-----------
tr.ext<-lapply(tr.list, function(x) extract(tdis.st,x))
tr.ext2<-lapply(tr.ext, data.frame)
ldd<-lapply(tr.list, function(x) unlist(extract(ld_dist2,x)))#add cols dist LD





##########no time? --- START HERE! but first load _25cm.RDATA############################
ldd_copy<-ldd


xsec.data<-list()
for(i in seq_along(tr.ext2)){
  if(is.null(ldd[[i]])==TRUE){
    xsec.data[[i]]<-NULL
  }else{
    xsec<-data.frame(qpcR:::cbind.na(tr.ext2[[i]],ldd[[i]]))
    colnames(xsec)[6]<-"ldd"
    xsec.data[[i]]<-xsec[complete.cases(xsec),]
  }  
}



##Interpolate!!!
xsec.data2<-list()
for(i in seq_along(xsec.data)){
  b<-xsec.data[[i]]
  ldd<-xsec.data[[i]]$ldd
  if(length(ldd) == 0){
    mer<-NULL
  }else{
    ldd<-seq(min(ldd,na.rm=TRUE),max(ldd,na.rm=TRUE),0.01)
    df2<-data.frame(ldd)
    mer<-merge(b,df2,by.x="ldd", by.y="ldd",all.y=TRUE)
    mer[2:6]<-apply(mer[2:6],2,na_interpolation)
    #mer<-mer[complete.cases(mer),]
  }
  
  xsec.data2[[i]]<-mer
}


##MID--------------
mid<-list()
for(i in seq_along(xsec.data2)){
  a<-xsec.data2[[i]]
  if(length(a) == 0){
    b<-NULL
  }else{
    ld<-a$ldd
    ld.mid<-apply(a[,2:6],2, function(x) ld[which.min(abs(x-5.5))])#difference of each elev. to 5, get the closest
    bnk<-apply(a[,2:6],2, function(x) x[which.min(abs(x-5.5))])#difference of each elev. to 5, get the closest
    b<-cbind(bnk,ld.mid)
    colnames(b)<-c("bank_mid","dist2LD")
  }
  mid[[i]]<-b
}


mid.dist<-lapply(mid, function(x) x[,2])
mid.dist2<-do.call(qpcR:::cbind.na,mid.dist)
mid.diff<-apply(mid.dist2,2, function(x) diff(x)*-1)
mid.diff2<-mid.diff[2,]/2 # 2018=2016 as annual
#Lod
mid.diff<-mid.diff-LoD
#mid.diff<-ifelse(abs(mid.diff)>6.5,NA,mid.diff)


####---------Main plot bank retreat---------
#pts.select
pts.select<-mid.diff
pts.select[,-c(40,130,220,310,400,490,590)]<-NA #Selected transects


library(scales)

 #par(mfrow=c(3,1),mar=c(3,3,3,3))
# #PLOT ortho2020
# plotRGB(ortho)
# plot(tr,add=TRUE, col=rf,lty=1,lwd=0.7)

#pdf("./Bank_lateral_tr_br_25cmpix_bl2.pdf",  width = 7, height =8) # The height of the plot in inches 
##plot with ortho
cols = brewer.pal(9,'YlGnBu')
rf <- colorRampPalette(cols)(length(tr))   # make colors

par(mfcol=c(4,1),mar=c(0,3,0,1),mgp=c(2.5,0.7,0),oma=c(3.5,2,2,2),
    tck=-0.015,cex.axis=1.5,cex.lab=1.5)
for(i in 1:nrow(mid.diff)){
  x<-mid.diff[i,]
  t<-pts.select
  if(i %in% 2){
    ylim=c(-11,5)
  }else{
    ylim=c(-5,5)
  }
  
  plot(x,type="p",pch=19,cex=0.85,col=alpha("grey38",0.5), ylim=ylim,xlab="",ylab="",xaxt="n",bty="n",yaxt="n")
  if(i %in% 2){
    points(mid.diff2,pch=19,cex=0.85,col=alpha("grey80",0.5), ylim=c(-11,5),xlab="",ylab="",xaxt="n",bty="n",yaxt="n")
  }
  axis(side=2,las=2)
  box(lwd=0.2)
  abline(h=0,lwd=0.3,lty=2,col="grey59")
  points(t,pch=4,cex=1.85,lwd=2,col="red2")
}
axis(side=1,at =seq(0,ncol(mid.diff),by=50))
mtext(side=2, "Bank midslope lateral movement (m)",outer=TRUE)
mtext(side=1, "Line transect",outer=TRUE, line=2.2)

#dev.off()

# boxplot(t(mid.diff))


######################################################################################################

###--------------------Profile transect-midslope change-----------------------------------------
#Selected transects:
tr.select<-tr.list[c(40,130,220,310,400,490,590)]


#extract elevs to trans.selec
postLD.evel<-tnew[7:11]
library(raster)

elev.tr<-list()
for(i in seq_along(tr.select)){
  elev.tr[[i]]<-lapply(postLD.evel, function(y) extract(y,tr.select[[i]]))
}


elev.tr2<-lapply(elev.tr,function(x) lapply(x, function(y) unlist(y)))
elev.tr2<-lapply(elev.tr2, function(x) do.call(cbind,x))

#elev
elev.tr2<-lapply(elev.tr2, function(x) x-LoD)

##Plot empty
xl1<-c(150,200,200,200,200,380,150)
xl2<-c(280,330,330,330,330,510,280)

#pdf("./t7_xsection3.pdf",  width =11, height =3) # The height of the plot in inches 
par(mfrow=c(1,7), mar=c(3,0,1,0.5),oma=c(1,4,1,1),tck=-0.02)
for(i in seq_along(elev.tr)){
  tr<-elev.tr2[[i]]
  plot(1, type="n", xlab="", ylab="", xlim=c(xl1[[i]],xl2[[i]]), ylim=c(4,10),yaxt="n",bty="n",cex.axis=1.25)
  box(lwd=0.2)
  abline(h=5.5,lty=2)
  if(i %in% 1){
    axis(side=2,las=2,cex.axis=1.5)
  }
  for(j in 1:5){
    lines(tr[,j],pch=19,cex=0.5,col=colnew[[j]])
  }
}
mtext(side=1,"Distance from LD", outer=TRUE,line=-0.5)
mtext(side=2, "Elevation (m +NAP)",outer=TRUE,line=2.5)
#dev.off()




                 