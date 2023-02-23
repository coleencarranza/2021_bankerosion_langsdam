##########-------------------Changes in bathymetry post-LD-------------------
library(raster)
library(sp)
library(rgdal)

#Read bathymetry
fils<-list.files(path="/media/coleen/DDrive/1_HWM/RWS_langsdam/DTM_Wam_old/bathy",pattern= "*.tif", full.names = TRUE )
bathy<-lapply(fils, raster)
bathy<-lapply(bathy, function(x) {proj4string(x)<-CRS("+init=epsg:28992");x})

#subset
poly<- readOGR(dsn="/media/coleen/DDrive/1_HWM/RWS_langsdam/DTM_Wam/shp", layer="poly4fillingNApre2")
e<-extent(poly)
bathy_Wamel<-lapply(bathy, crop, y=e)
bathy_Wamel<-lapply(bathy_Wamel, mask,poly)

#Bathy 2016 in cm
bathy_Wamel[[1]]<-bathy_Wamel[[1]]/100



#bathy changes
max_diff<-bathy_Wamel[[4]]-bathy_Wamel[[1]]



#Get elevation per pixel

#to points

pix<-rasterToPoints(max_diff,spatial = TRUE)

#coords_with_elev

pix_coords<-cbind(pix@coords,pix@data)
colnames(pix_coords)<-c("x","y","elev_diff")


#Export thedf

write.csv(pix_coords,file="./bathy_diff_2020_2016_coords.csv",row.names = FALSE)







#plot!
cols<- c("royalblue2","lightgoldenrodyellow","darkorange")
rf <- colorRampPalette(rev(cols)) # make colors
brks<-5000

#pdf("./AHN_plus_bathy_elevdiffs.pdf",width=7, height=5)
par(mar=c(1,1,1,1),mgp=c(1.25,0.15,0),tck=-0.2)

#BATHY
plot(max_diff,zlim=c(-3,2), col=ocean.pal(brks),xaxt="n",bty="n",xlab="",ylab="",legend=FALSE,box=FALSE,axes=FALSE)
#box(lwd=0.5)
#text(x=161500,y=432800,"2015 - 2020",cex=1.15)
text(x=159300,y=433200,"post- LD", font=2, cex=1.15)
#AHN
plot(diff.lang,zlim=c(-4,2), col=land.pal(brks),xaxt="n",bty="n",xlab="",ylab="",legend=FALSE,box=FALSE,axes=FALSE,add=TRUE)




##horizontal legend -AHN
r.range <- c(-5,5)
plot(diff.lang,zlim=c(-5,2), legend.only=TRUE, col=land.pal(brks), horizontal=TRUE,
     legend.width=2, legend.shrink=0.75,
     axis.args=list(at=seq(-5, 2, by=1),labels=seq(-5,2, by=1),cex.axis=1),
     legend.args=list(text=expression(paste(Delta,"Elevation (m) \n (Topography)")), side=1, font=2, line=2.5, cex=1),
     smallplot=c(0.5,0.75, 0.35,0.4))

plot(diff.lang,zlim=c(-3,2), legend.only=TRUE, col=ocean.pal(brks), horizontal=TRUE,
     legend.width=2, legend.shrink=0.75,
     axis.args=list(at=seq(-5, 2, by=1),labels=seq(-5,2, by=1),cex.axis=1),
     legend.args=list(text=expression(paste(Delta,"Elevation(m) \n  (Bathymetry) ")), side=1, font=2, line=2.5, cex=1),
     smallplot=c(0.5,0.75, 0.85,0.9))

#dev.off()




##Elevation change per year
bathy_year<-list()
x<-length(bathy_Wamel)-1
for(i in 1:x){
        bathy_year[[i]]<-bathy_Wamel[[i+1]]-bathy_Wamel[[i]]
}


#pdf("./Bathy_elevdiffs_year.pdf",width=7, height=7)

par(mfrow=c(3,1), mar=c(3,3,1,1))
lapply(bathy_year,plot,zlim=c(-3,2), col=ocean.pal(brks),xaxt="n",bty="n",xlab="",ylab="")



#dev.off()


#BAthy plot only
#pdf("./Bathy_elevs_year.pdf",width=7, height=7)

par(mfrow=c(4,1), mar=c(3,3,1,1))
lapply(bathy_Wamel,plot,col=ocean.pal(brks),xaxt="n",bty="n",xlab="",ylab="")

#dev.off()





##plot width vs elevation change
sn<-read.csv("/media/coleen/DDrive/1_HWM/RWS_langsdam/DTM_Wam/side_channel_sn/bathy_SN.csv", header=TRUE,stringsAsFactors = FALSE)
# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(sn[, 1:2])


r <- raster(e, ncol=1000, nrow=90)
x <- rasterize(sn[, 1:2], r, sn[,3], fun=mean)
plot(x)


#attempt from spdf to sgdf
sgdf <- as(sn, 'SpatialPixelsDataFrame')

# Convert the SpatialPixelsDataFrame (sgdf) to spatialGridDataFrame
sgdf <- as(sgdf, "SpatialGridDataFrame")



width<-read.csv("/media/coleen/DDrive/1_HWM/RWS_langsdam/DTM_Wam/side_channel_sn/width.csv", header=TRUE,stringsAsFactors = FALSE)


#plot
mn<-width$sn[[1]]
mx<-tail(width$sn,1)

pdf("./side_channel_width.pdf", height = 3, width = 10)
par(mar=c(3,3,1,1), mgp=c(1.75,0.5,0), tck=-0.02)
plot(width$sn,width$width,xaxt="n",type="l", xlab=" X - coordinate (m)",ylab="Side channel width (cm)")
axis(side=1, at=seq(mn,mx,length.out=7), labels = round(seq(min(sn$x),max(sn$x),length.out=7),0))

dev.off()


#########----------------
#Color pallette#

#make palette
ocean.pal <- colorRampPalette(
        c("#000000", "#000209", "#000413", "#00061E", "#000728", "#000932", "#002650", 
          "#00426E", "#005E8C", "#007AAA", "#0096C8", "#22A9C2", "#45BCBB", 
          "#67CFB5", "#8AE2AE", "#ACF6A8", "#BCF8B9", "#CBF9CA", "#DBFBDC", 
          "#EBFDED")
)

land.pal <- colorRampPalette(
        c("#996600", "#B27676", "#C2B0B0", "#E5E5E5", 
          "#FFFFFF")
)


##"#9F7B0D","#336600", "#F3CA89", "#D9A627", 
"#A49019", 

