library(terra)
library(scales)
#----------------------X-section based on transect-----------------------##
#read transects:
tr<-readOGR("./shp/ext_fc","transects_5m_200")
#to list
tr.list<-list()
for(i in 1:nrow(tr)){
  tr.list[[i]]<-tr[i,]
}

##-------------calculate cumulative lateral movement--------------
resp<-stack(diff.bank[c(7,8,10)])
resp<-sum(resp,na.rm = TRUE)
resp[resp==0]<-NA
plot(resp)

resp.bank<-resp
# 
# #total bank area
# sum<-sum(stack(bank[7:10]),na.rm=TRUE)
# sum[sum==0]<-NA
# sum[sum>0]<-1
# plot(sum)
# 
# 
# 
# resp.bank<-mask(resp,sum)
# # resp.bank<-resp
# plot(resp.bank)


# # remove edge effects, buffer in
# buf<-rasterToPolygons(sum, na.rm = TRUE,dissolve = TRUE)
# buf_in<-gBuffer(buf,width =-2)
# sum2<-mask(sum, buf_in)
# 
# par(mfrow=c(2,1))
# plot(sum)
# plot(sum2)
# 
# resp.bank<-mask(resp.bank,sum2)
# plot(resp.bank)

##transects
lines.list<-list()
for(i in 1:nrow(tr@data)){
  lines.list[[i]]<-tr[i,]
}

# par(mfrow=c(2,1))
# plot(resp.lin,type="o",pch=19, ylab="Total post LD bank volume changes (m3)")

##------------ read and compile files:---------------------

init_slope<-terrain(resp.bank, opt="slope", unit="degrees", neighbors=8)#2015
vegd<-raster("Wam_veg_fill_dist.tif")
goy<-raster("groynes_dist2.tif")
sill<-raster("LD_openings_dist.tif")
dstream<-raster("Downstream_dist.tif") #

##Clay and sand##
text<-readOGR(dsn = "./shp", layer = "clay_sand_poly")# clay and sand parts
text.ras<-rasterize(text,resp.bank)
text.ras<-mask(text.ras,resp.bank)



#to list
vars.list<-list(init_slope,vegd,goy,sill,dstream,text.ras)
vars.list<-lapply(vars.list, function(x) resample(x,resp.bank))
vars.list<-lapply(vars.list, function(x) mask(x,resp.bank))


#further mask out lake
fin_mask<-readOGR("./shp","poly_no_bay")


vars.list<-lapply(vars.list, function(x) mask(x,fin_mask))
lapply(vars.list,plot)


#Combine distance to downsrteam and to LD openings - normalize first
norm<-function(x){
  mn<- min(values(x),na.rm = TRUE)
  mx<- max(values(x),na.rm = TRUE)
  denom<-mx - mn
  num <- x- mn
  norm = num / denom
  return(norm)
}
  

LD_open_norm<-norm(vars.list[[4]])
LD_dstream_norm<-norm(vars.list[[5]])
#Product two norms
open_dstream_norm<-LD_dstream_norm * LD_open_norm



#final vars.list
vars.list<-c(vars.list[1:3],open_dstream_norm,vars.list[6])
plot(vars.list[[4]])



#for figure
cols = brewer.pal(9,'Greys')
cols<- c("lightgoldenrodyellow","darkorange","darkred")
rf <- colorRampPalette(cols) # make colors
#pdf("./Vars_spatial3.pdf",width=5, height=9)
par(mfrow=c(5,1),oma=c(1.5,1,1,0),mar=c(0,0,0,1))
#SLOPE
plot(vars.list[[1]],col=rf(30),xaxt="n",bty="n",xlab="",ylab="",box=FALSE,axes=FALSE,legend=TRUE,
     horizontal=TRUE,legend.width=2, legend.shrink=0.75,cex.axis=1,smallplot=c(0.7,0.95, 0.45,0.5))
##OTHERS
for(i in 2:4){
  plot(vars.list[[i]],col=rf(30),xaxt="n",bty="n",xlab="",ylab="",box=FALSE,axes=FALSE,legend=TRUE,
        horizontal=TRUE,legend.width=2, legend.shrink=0.75,cex.axis=1,smallplot=c(0.7,0.95, 0.45,0.5))
}
#SUBSURFACE
plot(vars.list[[5]],col=rf(30)[c(20,5)],xaxt="n",bty="n",xlab="",ylab="",box=FALSE,axes=FALSE,legend=TRUE,
     horizontal=TRUE,legend.width=2, legend.shrink=0.75,cex.axis=1,smallplot=c(0.7,0.95, 0.45,0.5))
#dev.off()



#function!
ras.prep<-function(x){
  x.vals<-values(x)
  y.vals<-values(resp.bank)
  all<- cbind(x.vals,y.vals)
  all<-all[complete.cases(all),]
  return(all)
  
}

ras.vars<-lapply(vars.list, ras.prep)
par(mfrow=c(1,5))
lapply(ras.vars,plot)



###################################################
##Methods 2 ----use transect lines - takes a long time!##
resp.lin<-lapply(lines.list, function(x) terra::extract(resp.bank,x))
r<-lapply(resp.lin,unlist)

#function!
tr.prep<-function(ras){
  x_lin<-lapply(lines.list, function(z) terra::extract(ras,z))
  raslins<-lapply(x_lin,unlist)
  #combine two per list
  all<-mapply(function(x,y) cbind(x,y), raslins,r)
  all<-lapply(all, function(x) x[complete.cases(x),])
  #all.bind<-do.call("rbind",all)
  return(all)
}



tr.vars<-lapply(vars.list, tr.prep)
tr.vars2<-tr.vars




##############no time? --Start here! but LOAD RDATA - controls2#################
#plots#####
#1. all points - overlay rasters
par(mfrow=c(1,5), mar=c(2,2,2,2),oma=c(2,2,0,0))
lapply(ras.vars, plot,pch=19, cex=0.2)



#2. transects
library(viridis)
coltr<-rev(viridis(length(lines.list)))
#coltr = colorRampPalette(brewer.pal(9,'YlGnBu'))(length(lines.list))

#png("Controls_bank2016_2018_2020.png",height=10, width=8,units = "in",res = 600)
set.seed(12345)
par(mfrow=c(3,2),mar=c(3,4,1,1), mgp=c(1.7,0.57,0),tck=-0.015)
xlab<-c("Initial Slope (deg)","Vegetation proximity (m)","Groyne remnamt proximity (m)","Longitudinal dam opening proximity (m)")
ylab<-c("Elevation Change (m)","","","")
for(j in seq_along(tr.vars)){
  pl<-tr.vars[[j]]
  s<-sample(1:length(pl),length(pl)/1)
  pl<-pl[s]
  mx<-max(do.call("rbind",pl)[,1])
  mn<-min(do.call("rbind",pl)[,1])
  #empty plot
 if(j %in% 1:4){
   plot(1, type="n", ylab=expression(paste(Delta,"Elevation (m)")), ylim=c(-5,2), xlim=c(mn, mx), xlab=xlab[[j]],cex.lab=1.5,cex.axis=1.2 )
   abline(h=0, lty=2, col="grey59")
   for(i in seq_along(pl)){
     points(pl[[i]],bg=alpha(coltr[s[[i]]],0.15),col=alpha(coltr[s[[i]]],0.15),pch=21,cex=0.88)
   }
 
 }
}

#add line number before combining
bx<-tr.vars[[5]]
bx<-lapply(bx, function(x) unlist(x))
bx<-lapply(bx, function(x) data.frame(x))

for(i in seq_along(bx)){
  if(nrow(bx[[i]])>0){
  bx[[i]]$lin<-i
  }
}
  
bx<-do.call("rbind",bx)




#BOXPLOT!
boxplot(bx$y~bx$x, col="white", ylim=c(-5,2), names=c("Clay","Sand"),xlab="",ylab=expression(paste(Delta,"Elevation (m)")),
        cex.lab=1.5,cex.axis=1.2 )
set.seed(12345)
# Add data points
levProp <- summary(as.factor(bx$x))/nrow(bx)
bx$jitt<- jitter(bx$x, amount=0.27)
#sample points
s<-sample(1:nrow(bx),nrow(bx)/15)
bxs<-bx[s,]
points(bxs$jitt,bxs$y, pch=20, col=alpha(coltr[bxs$lin],0.5),cex=0.86,lwd=0) 


#dev.off()








################################################
##----------------Hexbins#--------------------------------
library(hextri)
library(scales)
data(airquality)
airquality$o3group<-with(airquality, cut(Ozone, c(0,18,60,Inf)))

with(na.omit(airquality), 
     hextri(Solar.R, Temp, class=o3group, colours=c("skyblue","grey60","goldenrod"), style="size",
            xlab="Solar radiation (langley)", ylab="Temperature (F)")
)


#
library(hexbin)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

ylabs<-c("Slope (degrees)","Proximity to vegetation (m)","Proximity to groyne remnants (m)",
         "Normalized proximity to LD openings\n and downstream direction (-)" )
plot.list<-list()
for(i in 1:4){
  test<-tr.vars[[i]]
  tr.test<-do.call(rbind,test)
  tr.test[,1]<-ifelse(tr.test[,1]==0, NA,tr.test[,1])
  #remove 0.03
  tr.test[,2]<-ifelse(tr.test[,2] >= -0.34 & tr.test[,2]<= -0.29, NA,tr.test[,2])
  #complete.case
  tr.test<-tr.test[complete.cases(tr.test),]
  
  tr.test<-tr.test[!duplicated(tr.test[,1]),]
  tr.test<-tr.test[!duplicated(tr.test[,2]),]
  
  
  # Make the plot
  df<-data.frame(tr.test)
  
  plot.list[[i]]<-ggplot(df, aes(x = x, y = y)) +
        geom_hex(color="black", bins=20) +
    scale_fill_viridis_c(alpha = 0.4, option="magma",direction=-1) +
    guides(fill = guide_colourbar(title = "Count"))+
    scale_x_continuous(name=ylabs[[i]]) +
    scale_y_continuous(name="Lateral movement (m)")+
    theme_bw(base_size = 15)
  
}

#BOXPLOT
data <- data.frame(do.call("rbind",tr.vars[[5]]))
data$x<-ifelse(data$x==1, "Clayey","Sandy")

# Plot
  boxp<-ggplot(data, aes(x=x, y=y, fill=x)) +
  geom_boxplot(outlier.color = "grey38",outlier.size = 0.85) +
  scale_fill_viridis(discrete = TRUE, alpha=0.2,name="Substrate",) +
  labs(x="", y = "Lateral movement (m)")+
  theme_bw(base_size = 15)

boxp




figure <- ggarrange( plot.list[[1]],  plot.list[[2]],  plot.list[[3]],  plot.list[[4]], boxp,
                    labels = c("a", "b", "c","d","e"),hjust= 0.15,ncol = 2, nrow = 3)+
  theme(plot.margin = margin(0.1,0.1,0,0.5, "cm")) 

figure

ggsave('ctrl_multi_magma.png', width = 12, height = 13, dpi = 300)






##########################################################################################################
#-------------Extra! --- multivar stats for controls-------------
#Simple linear regression
par(mfrow=c(2,3), mar=c(2,2,2,2))

lm.vars<-list()

for(i in 1:4){
  test<-tr.vars[[i]]
  tr.test<-do.call(rbind,test)
  tr.test[,1]<-ifelse(tr.test[,1]==0, NA,tr.test[,1])
  #remove 0.03
  tr.test[,2]<-ifelse(tr.test[,2] >= -0.34 & tr.test[,2]<= -0.29, NA,tr.test[,2])
  #complete.case
  tr.test<-tr.test[complete.cases(tr.test),]
  
  tr.test<-tr.test[!duplicated(tr.test[,1]),]
  tr.test<-tr.test[!duplicated(tr.test[,2]),]
  
  
  # Make the plot
  x<-tr.test[,1]
  y<-tr.test[,2]
  new <-data.frame(x)
  p<-predict(lm(y ~ x), new, se.fit = TRUE)
  res = resid(p)
  plot(x, res) 
  abline(0, 0)                  # the horizon
  
}



##PCA####
tr.pc<-lapply(tr.vars, function(x) do.call(rbind,x))
tr.pc.x<-lapply(tr.pc, function(x) x[,1])
tr.pc.x<-do.call(qpcR:::cbind.na,tr.pc.x)

tr.pc.x<-tr.pc.x[complete.cases(tr.pc.x),]

pc <- prcomp(tr.pc.x,
             center = TRUE,
             scale. = TRUE)


print(pc)
screeplot(pc, axisLabSize = 18, titleLabSize = 22)
