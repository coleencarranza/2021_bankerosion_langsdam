##----------------------OVERBANK yearly analysis - extra!-------------------
getObank<-function(ras,mask){
  m<-mask
  o<- mask(m , mask = m==1 , maskvalue =0)#1=overbank
  r<-mask(ras,o)
  return(r)
}

#
pre.obank<-mapply(function(x,y) getObank(x,y), tnew[1:5],bnk.class.filled[1:5])
post.obank<-mapply(function(x,y) getObank(x,y), tnew[7:11],bnk.class.filled[7:11])



##----------TRANSECTS Plot--------------##
library(raster)
library(scales)
#read transects:
tr<-readOGR("./shp/ext_fc","transects_5m_200")
#to list
tr.list<-list()
for(i in 1:nrow(tr)){
  tr.list[[i]]<-tr[i,]
}

tr.select<-tr.list[c(40,130,220,310,400,490,590)]


##--------------Extract values to transect------------##
#function
toLines<-function(line,rasters){
  tr.lin.yr<-lapply(rasters, function(x) terra::extract(x,line))
  tr.lin.yr<-lapply(tr.lin.yr,unlist)
  tr<-do.call(cbind,tr.lin.yr)
  #tr<-tr[complete.cases(tr),]
  return(tr)
}


pre.tr.lin<-lapply(tr.select,function(x) toLines(x,pre.obank))
post.tr.lin<-lapply(tr.select,function(x) toLines(x,post.obank))




par(mfrow=c(2,7),mar=c(2,0,2,0),oma=c(0,3,0,1))
colyr = colorRampPalette(brewer.pal(9,'YlGnBu'))(20)[c(8,9,12,15,16)]
#PRE
for(i in seq_along(pre.tr.lin)){
  t<-pre.tr.lin[[i]]
  nro<-nrow(t)
  plot(1, type="n", xlab="", ylab="", ylim=c(6,10),xlim=c(1,nro),xaxt="n",yaxt="n")
  if(i== 1)
    axis(side=2)
  for(j in c(1,5)){#1:ncol(t)
    lines(t[,j],col=colyr[j],lwd=1.5)
  }
  
}
legend("topright", lwd=1, legend=c(2009,2013),seg.len = 2,col=colyr[c(1,5)],bty="n")

#POST
for(i in seq_along(post.tr.lin)){
  t<-post.tr.lin[[i]]
  nro<-nrow(t)
  plot(1, type="n", xlab="", ylab="", ylim=c(6,10),xlim=c(1,nro),xaxt="n",yaxt="n")
  if(i== 1)
    axis(side=2)
  for(j in c(1,5)){#1:ncol(t)
    lines(t[,j],col=colyr[j],lwd=1.5)
  }
  
}

legend("topright", lwd=1, legend=c(2015,2020),seg.len = 2,col=colyr[c(1,5)],bty="n")




##-----------------------------------
##Total change
par(mfrow=c(2,1),mar=c(2,0,2,0),oma=c(0,3,0,1))
preLD_obank<-pre.obank[[5]]-pre.obank[[1]]
#less than zero set to NA
preLD_obank[preLD_obank < -0.5] <- NA
preLD_obank[preLD_obank > 0.5] <- NA
plot(preLD_obank,box=FALSE,axes=FALSE,main="preLD")



postLD_obank<-post.obank[[5]]-post.obank[[1]]
#less than zero set to NA
postLD_obank[postLD_obank < -0.5] <- NA
postLD_obank[postLD_obank > 0.5] <- NA
plot(postLD_obank,box=FALSE,axes=FALSE,main="postLD")
