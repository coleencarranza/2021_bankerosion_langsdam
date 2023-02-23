library(raster)
library(sp)
library(rgdal)
library(fields)
library(RColorBrewer)
library(rgeos)
library(graphics)
#read 1_Bank_obank_Kmeans_classif.R first!


#-----------Subtract 2 consecutive rasters---------------------
le<-length(tnew)-1

diff.ras<-list()
for(i in 1:le){
  t2<-tnew[[i+1]]
  t1<-tnew[[i]]
  diff<-t2-t1
  diff.ras[[i]]<-overlay(diff, t1, fun = function(x, y) {x[x == 0 & y == 2.8] <- NA; x})
  
}

lapply(diff.ras,plot)




##-------------Uncertainly removal - LoD after diff-----------------------------
error<-read.csv("./error_dtm_rws.csv", header=TRUE)
st.err<-error$St_dev
#t= 2sigma
t = 1.96

#Uc
Uc<-numeric()
for(i in 1:length(st.err)){
  if(i< length(st.err))
    Uc[i]<-t*sqrt(st.err[i]^2+st.err[i+1]^2)
}

Uc
LoD<-round(mean(Uc,na.rm=TRUE),3)
LoD
#0.025

diff.ras<-lapply(diff.ras, function(x) x-LoD)




###-----------------apply classified image to separate diff in bank and overbank---------------
#clean bnk.class.filled for extra pixels because of focal
bnk.class.filled2<-mapply(function(x,y) mask(x,y), bnk.class.filled2,diff.ras)


for(i in seq_along(bnk.class.filled2)){
  bnk.class.filled2[[i]]<-mask(bnk.class.filled2[[i]],diff.ras[[i]])
}

#mask again to separate bank and obank
obank<-lapply(bnk.class.filled2, function(x) mask(x , mask = x==1 , maskvalue =0))
bank<-lapply(bnk.class.filled2, function(x) mask(x , mask = x==2 , maskvalue =0))


diff.bank<-mapply(function(x,y) mask(x,y), diff.ras,bank)
diff.obank<-mapply(function(x,y) mask(x,y), diff.ras,obank)


ras.obank<-mapply(function(x,y) mask(x,y), tnew[-1],obank)
ras.bank<-mapply(function(x,y) mask(x,y), tnew[-1],bank)


# ##---------------PLOTTING-----------------------
cols = brewer.pal(11,'RdYlBu')
rf <- colorRampPalette(cols)   # make colors
rf=  colorRampPalette(c("darkred","red2", "grey99", "deepskyblue2","blue4"))
par(mar=c(2,2.5,2,2.5),mgp=c(2.5,0.5,0.0),oma=c(1,1,1,3),tck=-0.01)

#color scale
brks<-50
cuts=seq(-9,9,length.out=brks) #set breaks

yr<-substr(ts, 17,20)
yr<-as.numeric(yr)

d<-le-1
for(i in 1:length(diff.bank)){
  r <- rf(length(diff.bank[[i]]))
  plot(diff.bank[[i]],breaks=cuts, col =rf(brks),main=paste(yr[[i+1]],"-",yr[[i]]))
  #  image.plot(diff.bank[[i]],col=r, legend.only=T,zlim=zlim)
}



##--------------------cumulative changes in elevation------------------------------
##Cumulative krib years
diff.kr<-stack(diff.bank[1:4])## index of lower year when diff is calculated
diff.kr<-sum(diff.kr,na.rm=TRUE)
diff.kr[diff.kr==0]<-NA
diff.kr<-focal(diff.kr,w=matrix(1,nrow=3, ncol=3), fun=mean, NAonly=TRUE, na.rm=TRUE)


##Cumulative krib-ld years
diff.krld<-stack(diff.bank[5:6])## index of lower year when diff is calculated
diff.krld<-sum(diff.krld,na.rm = TRUE)
diff.krld[diff.krld==0]<-NA


##Cumulative langsdam years
diff.lang<-stack(diff.bank[7:10])
diff.lang<-sum(diff.lang,na.rm = TRUE)
diff.lang[diff.lang==0]<-NA


##Cumulative changes over 10 years
diff.s<-stack(diff.bank)
diff.sum<-sum(diff.s,na.rm=TRUE)
diff.sum[diff.sum==0]<-NA



###------- volume changes from elevation changes -------------
lang.vol<-diff.lang*0.5*0.5
kr.vol<-diff.kr*0.5*0.5

#im mteric tons
lang.vol.er<-sum(lang.vol[lang.vol<0])/1000
lang.vol.dep<-sum(lang.vol[lang.vol>0])/1000

#
kr.vol.er<-sum(kr.vol[kr.vol<0])/1000
kr.vol.dep<-sum(kr.vol[kr.vol>0])/1000





##Erosion per transect-----------------
#read transects:
tr<-readOGR("./shp/ext_fc","transects_5m_200")
#to list
tr.list<-list()
for(i in 1:nrow(tr)){
  tr.list[[i]]<-tr[i,]
}


#volume pre and post LD per transect
tr.vol.lang<-lapply(tr.list, function(x) extract(lang.vol,x))

tr.vol.lang<-lapply(tr.vol.lang,unlist)
tr.vol.lang.mat<-do.call(qpcR:::cbind.na,tr.vol.lang)
tr

##cumulative changes - colmeans
tr.cumm<-apply(tr.vol.lang.mat,2, sum, na.rm=TRUE)
plot(abs(tr.cumm))
max(abs(tr.cumm))

######################################################
tr.er.lang<-ifelse(tr.vol.lang.mat>=0,NA,tr.vol.lang.mat)
tr.er.lang<-apply(tr.er.lang,2,sum, na.rm=TRUE)
plot(tr.er.lang)


tr.dep.lang<-ifelse(tr.vol.lang.mat<=0,NA,tr.vol.lang.mat)
tr.dep.lang<-apply(tr.dep.lang,2,sum, na.rm=TRUE)
plot(tr.dep.lang)

#net volume per transect
tr.vol.lang<-apply(tr.vol.lang,2, function(x) sum(x, na.rm = TRUE))
plot(tr.vol.lang, type="o")



#####for study area plot------------
library(viridis)
#pdf("./Studyarea_part_elevdiffs.pdf",width=7, height=5)
par(oma=c(1.5,1,1,0),mar=c(0,0,0,1))
layout(matrix(1:2,nrow = 2),heights=2)
#DTM 2013
plot(tif.filled2[[5]],col =viridis(brks),xaxt="n",bty="n",xlab="",ylab="",box=FALSE,axes=FALSE)
text(x=161500,y=432700,"DTM 2013",cex=1.15)
#text(x=159300,y=433200,"pre - LD \n(with groynes)", font=2, cex=1.5)

plot(diff.krld,col =viridis(brks),xaxt="n",bty="n",xlab="",ylab="",box=FALSE,axes=FALSE)
text(x=161300,y=432800,"Bank elevation changes\n after LD installation\n in 2015",cex=1.15)

#dev.off()


#-------------------------PLOT CHANGES IN ELEVATION - per period-------------------------------------
#pdf("./Diff.ras_Wam.pdf",width=8.5, height=9)
par(oma=c(3,1,1,2))
layout(matrix(c(1,5,
                2,5,
                3,5,
                4,5),nrow = 4, ncol = 2,byrow=TRUE ),widths = c(1,0.125))
par(mar=c(1,3,1,0))
#layout(matrix(1:5,nrow = 1, ncol = 5),widths = c(1,1,1,1,0.35))
yr<-yrs[-1]

brks<-50
cuts=seq(-9,9,length.out=brks) #set breaks

image(diff.kr,zlim=c(-3,3), col =viridis(brks),xaxt="n",bty="n",xlab="",ylab="")
box(lwd=0.5)
text(x=161500,y=432700,"2009 - 2013",cex=1.5)
text(x=159300,y=433200,"pre - LD \n(with groynes)", font=2, cex=1.5)


image(diff.krld,zlim=c(-5,5), col =viridis(brks), xaxt="n",bty="n",xlab="",ylab="")
box(lwd=0.5)
text(x=161500,y=432700,"2013 - 2015",cex=1.5)
text(x=159300,y=433200,"LD installation", font=2, cex=1.5)

image(diff.lang,zlim=c(-5,5), col=viridis(brks),xaxt="n",bty="n",xlab="",ylab="")
box(lwd=0.5)
text(x=161500,y=432700,"2015 - 2020",cex=1.5)
text(x=159300,y=433200,"post- LD", font=2, cex=1.5)

image(diff.sum,zlim=c(-5,5), col=viridis(brks),xaxt="n",bty="n",xlab="",ylab="")
axis(side=1)
box(lwd=0.5)
text(x=161500,y=432700,"2009 - 2020",cex=1.5)
text(x=159300,y=433200,"Overall change", font=2, cex=1.5)

# Revert to c(1, 1) layout and adjust legend margins

##Second plot
par(mar=c(25,2,20,0.5))
legend_image <- as.raster(rf(brks),0.65)
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.75, y = seq(0.01,0.99,l=7), labels = seq(-9,9,length.out =  7), cex=1.5)
rasterImage(legend_image, 0.0,0.0,0.45,1)

#dev.off()


##-------------plot only pre and post - LD--------------------------
#dev.off()
#pdf("./Pre_post_totalchange2.pdf",width=7, height=9)

par(oma=c(1.5,1,1,0))
layout(matrix(1:2,nrow = 2),heights=2)

yr<-yrs[-1]
brks<-50

cols<- c("royalblue2","lightgoldenrodyellow","darkorange")
rf <- colorRampPalette(rev(cols)) # make colors

cuts=seq(-9,9,length.out=brks) #set breaks
par(mar=c(0,0,0,0))
plot(diff.kr,zlim=c(-3,1), col=rf(brks),xaxt="n",bty="n",xlab="",ylab="",legend=FALSE,box=FALSE,axes=FALSE)
#box(lwd=0.5)
text(x=161500,y=432700,"2009 - 2013",cex=1.15)
text(x=159300,y=433200,"pre - LD \n(with groynes)", font=2, cex=1.15)

plot(diff.lang,zlim=c(-3,1), col=rf(brks),xaxt="n",bty="n",xlab="",ylab="",legend=FALSE,box=FALSE,axes=FALSE)
#box(lwd=0.5)
text(x=161500,y=432800,"2015 - 2020",cex=1.15)
text(x=159300,y=433200,"post- LD", font=2, cex=1.15)


##horizontal legend
r.range <- c(-5,5)
plot(diff.lang,zlim=c(-5,5), legend.only=TRUE, col=rf(brks), horizontal=TRUE,
     legend.width=2, legend.shrink=0.75,
     axis.args=list(at=seq(-5, 5, by=1),labels=seq(-5, 5, by=1),cex.axis=1),
     legend.args=list(text=expression(paste(Delta,"Elevation (m)")), side=1, font=2, line=2.5, cex=1),
     smallplot=c(0.5,0.75, 0.15,0.2))
#dev.off()




###Extract diff lang to transects - takes a long time!
# ##----------plot yearly changes-------------------
# #pdf("./Diff_ras_Wam_yr_cumu_bnk.pdf",width=15, height=12)
# par(oma=c(3,1,0,0))
# layout(matrix(c(1,5,7,
#                 2,6,8,
#                 3,10,9,
#                 4,10,0,
#                 11,12,13),nrow =5, ncol = 3,byrow=TRUE))
# par(mar=c(1,2,2,3),mgp=c(1,0.35,0),tck=-0.02,cex.axis=0.9,bty="n")
# #layout(matrix(1:5,nrow = 1, ncol = 5),widths = c(1,1,1,1,0.35))
# brks<-50
# cuts=seq(-5,5,length.out=brks) #set breaks
# for(i in seq_along(diff.bank)){
#   image(diff.bank[[i]],zlim=c(-5,5), col =rf(brks),bty="n",xlab="",ylab="")
#   box(lwd=0.5)
#   text(x=161500,y=432700,paste0(yrs[[i+1]]," - ",yrs[[i]]),cex=1.5)
#   if(i==1)
#     title("pre - LD (with groynes)", font=2, cex=1.75)
#   if(i==5)
#     title("LD installation", font=2, cex=1.75)
#   if(i==7)
#     title("post - LD", font=2, cex=1.75)
# }
# #LEGEND
# par(mar=c(5,5,5,6))
# legend_image <- as.raster(rf(brks))
# plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
# text(x=0.75, y = seq(0.01,0.99,l=7), labels = round(seq(-5,5,length.out =  7),1), cex=1.5)
# rasterImage(legend_image, 0.0,0.0,0.45,1)
# 
# 
# ##add cumulative changes
# par(mar=c(1,2,2,3),mgp=c(1,0.35,0),tck=-0.02,cex.axis=0.9,bty="n")
# image(diff.kr,zlim=c(-5,5), col =rf(brks),bty="n",xlab="",ylab="")
# text(x=161300,y=432700,paste0("Cumulative change \n",yrs[[1]]," - ",yrs[[5]]),cex=1.5)
# image(diff.krld,zlim=c(-5,5), col =rf(brks),bty="n",xlab="",ylab="")
# text(x=161300,y=432700,paste0("Cumulative change \n",yrs[[5]]," - ",yrs[[7]]),cex=1.5)
# image(diff.lang,zlim=c(-5,5), col =rf(brks),bty="n",xlab="",ylab="")
# text(x=161300,y=432700,paste0("Cumulative change \n",yrs[[7]]," - ",yrs[[10]]),cex=1.5)
# 
# #dev.off()


##---------------------------PLOT POST LD yearly changes---------------------------
#pdf("./PostLD_Changes2.pdf",width=5, height=7)

par(oma=c(1.5,1,1,0),mar=c(0,0,0,0))
# 
# par(oma=c(1.85,1,1,2),tck=-0.42)
layout(matrix(c(1:4),nrow = 4, ncol =1,byrow=TRUE ),heights=c(2,2,2,2.5))
#par(mar=c(0,0,0,0),mfrow=c(4,1))

rf=  colorRampPalette(c("red2","white","limegreen"))
cols<- c("royalblue2","lightgoldenrodyellow","darkorange")
rf <- colorRampPalette(rev(cols)) # make colors

brks<-20
cuts=seq(-4,4,length.out=brks) #set breaks
for(i in 7:10){
  plot(diff.bank[[i]],zlim=c(-4,4), col =rf(brks),bty="n",xlab="",ylab="",legend=FALSE, axes=FALSE, box=FALSE)
  #box(lwd=0.5)
  text(x=161500,y=433000,paste0(yrs[[i+1]]," - ",yrs[[i]]),cex=1.5)
  #title("post - LD", font=2, cex=1.75)
}

##horizontal legend
r.range <- c(-5,5)
plot(diff.bank[[10]],zlim=c(-4,4), legend.only=TRUE, col=rf(brks), horizontal=TRUE,
     legend.width=2, legend.shrink=0.75,
     axis.args=list(at=seq(-4, 4, by=1),labels=seq(-4,4, by=1),cex.axis=1),
     legend.args=list(text=expression(paste(Delta,"Elevation (m)")), side=1, font=2, line=2.5, cex=1),
     smallplot=c(0.3,0.60, 0.15,0.22))

#dev.off()





###---------postLD bnk area classified---------------------------
postLD<-diff.bank[7:10]

###calculate erosion and deposition----------
bank.ero<-lapply(postLD, function(x) {x[x>0]<-NA;x})#set dep to NA
bank.dep<-lapply(postLD, function(x) {x[x<0]<-NA;x})




##--------------plot volume sediments post LD------------------
##plots
lapply(bank.ero,plot)
lapply(bank.dep,plot)


par(mfrow=c(2,1), mar=c(2,1,4,2))
for(i in seq_along(bank.ero)){
  plot(bank.ero[[i]], main=paste( "Bank Erosion",yr[[i+7]]," - ",yr[[i+6]],sep=" "),
       col=topo.colors(50))
  plot(bank.dep[[i]],main="Bank Deposition",
       col=topo.colors(50))
}

bank.ero.vol<-sapply(bank.ero, function(x) cellStats(x, na.rm=TRUE, stat='sum'))
bank.dep.vol<-sapply(bank.dep, function(x) cellStats(x, na.rm=TRUE, stat='sum'))

#to matrix
bank.ero.dep.vols<-cbind(bank.ero.vol,bank.dep.vol)
#to yearly for 2016-2018
bank.ero.dep.vols[2,]<-bank.ero.dep.vols[2,]/2




#plot volumes annual rates of erosion#
#pdf("./PostLD_dep_ero.pdf",width=5, height=5)
par(mar=c(1,2,1,2),mgp=c(2.2,0.4,0),tck=-0.01, oma=c(2,1,2,1))
plot(abs(bank.ero.dep.vols[,1])/1000, type="o" ,pch=21, bg="darkorange", cex=1.85,ylab="",yaxt="n",xaxt="n",bty="n",
     xaxt="n",ylim=c(10,150))
axis(side=2,las=2)
mtext(side=2, "Annual erosion rate (x1000 m3/yr)",outer=TRUE, line=-0.35)
par(new=TRUE)
plot(bank.ero.dep.vols[,2]/1000,type="o",yaxt="n",pch=21, bg="cornflowerblue", cex=1.85,ylab="", ylim=c(0,30),xaxt='n',bty="n")
axis(side=4,las=2)
mtext(side=4, "Annual deposition rate (x1000 m3/yr)",outer=TRUE, line=-0.35)
axis(side=1,at=1:4, labels=c("2015 - 2016 ","2016 - 2018","2018 - 2019","2019 - 2020"))


legend("topright", legend=c("erosion","deposition"),pch=15, col=c("darkorange","cornflowerblue"),bty="n",pt.cex=1.52)
#dev.off()




##---------------------BARCHART!!!-------------------------------
#pdf("./PostLD_dep_ero_bar.pdf",width=6, height=5)
par(mar=c(3,3,1,2),mgp=c(2.2,0.4,0),tck=-0.01, oma=c(1,1,1,1))
barplot(t(abs(bank.ero.dep.vols)/1000),beside=T,col=c("darkorange","cornflowerblue"),ylim=c(0,140),
        yaxt="n",ylab="Annual sediment volume changes (x1000 m3/yr)",
        names=c("2015 - 2016 ","2016 - 2018","2018 - 2019","2019 - 2020"))
axis(side=2, las=2)
legend("topright", legend=c("erosion","deposition"),pch=15, col=c("darkorange","cornflowerblue"),bty="n",pt.cex=1.52)
#dev.off()

# 
# ###-boxplots per year--------
# vec4plot<-function(x){
#   y<-lapply(x, function(z) z[!is.na(values(z))])
#   yy<-do.call(qpcR:::cbind.na,y)
#   y2018<-yy[,7]/3
#   y.else<-yy[,-7]
#   yfin<-qpcR:::cbind.na(y.else[,1:6],y2018,y2018,y2018,y.else[,7:8])
#   colnames(yfin)<-c(2010:2020)
#   return(yfin)
# }
# 
# all.list<-list(bank.dep,bank.ero,obank.dep,obank.ero)
# all.list.vec<-lapply(all.list, vec4plot)

#------------------------------------------------------