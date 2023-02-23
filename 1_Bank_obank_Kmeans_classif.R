library(raster)
library(sp)
library(rgdal)
library(fields)
library(RColorBrewer)
library(rgeos)
library(graphics)


##read filled tifs
ts<-list.files(path = "./dtm/dtm_filled",pattern = "*fill.tif",full.names = TRUE)
tif.filled2<-lapply(ts, raster)
lapply(tif.filled2,plot)

yrs<-as.numeric(substr(ts,17,20)) # add 5 for full


##########################################
# replacing NA's by zero in bank only for subtraction

#import help polygons since pre and post years differ in extent!
poly4NApre<-readOGR(dsn = "./shp", layer = "poly4fillingNApre2")
poly4NApost<-readOGR(dsn = "./shp", layer = "poly4fillingNApost2")


tnew<-list()
for(i in seq_along(tif.filled2)){
  t<-tif.filled2[[i]]
  if(i %in% 1:5){
    plot1<-poly4NApre[poly4NApre$id==1,]
    plot2<-poly4NApre[poly4NApre$id==2,] #inner
  }else{
    plot1<-poly4NApost[poly4NApost$id==1,]
    plot2<-poly4NApost[poly4NApost$id==2,] #inner
  }
  t1<-mask(t,plot1)
  t1[is.na(t1)] <- 2.8 #NA to same elevations
  t1<-mask(t1,plot1)
  t2<-mask(t,plot2)
  tnew[[i]]<-merge(t1,t2)
  plot(tnew[[i]])
}


#------------------------K-means Classify bank and overbank----------------------
#Kmeans function!
kmeans.class<-function(x,cl){
  ## returns the values of the raster dataset and write them in a matrix. 
  image<-x
  nr.all<-getValues(image)
  nr<-nr.all[which(!is.na(nr.all)==TRUE)]
  
  # It is important to set the seed generator because `kmeans` initiates the centers in random locations
  # We want to create 10 clusters, allow 500 iterations, start with 5 random sets using "Lloyd" method
  kmncluster <- kmeans(na.omit(nr), centers = cl, iter.max = 500, nstart = 1, algorithm="Lloyd")
  
  # You can also do it like this
  nr.kmeans<-nr.all
  nr.kmeans[which(!is.na(nr.kmeans)==TRUE)]<-kmncluster$cluster
  img<-image
  values(img)<-nr.kmeans
  plot(img)
  return(img)
}

#dont forget  to set seed
set.seed(99)
bnk.class.filled<-lapply(tnew, function(x) kmeans.class(x,cl=2))


#harmonize classes 
rclss<-function(x){
  m<-c(0,1,0.5,1,2,0.25)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(x, rclmat)
  rc<-rc*4
  return(rc)
}    

lapply(bnk.class.filled, plot)
bnk.class.filled[c(3,6,8:11)]<-lapply(bnk.class.filled[c(3,6,8:11)], rclss)# adjust each time set of rasters are changed!
lapply(bnk.class.filled,plot)

# #export to save for paper!
# fnames<-substr(tif.list,7,nchar(tif.list))
# for(i in seq_along(bnk.class.filled)){
#   writeRaster(bnk.class.filled[[i]],paste0("./dtm/Kmeans/K_",fnames[[i]]),overwrite=TRUE)
# }
# 

#accuracy assessment with validation points
library(caret)
dsn<-"./shp/random_pts"
rdm.pts<-list.files(dsn,"*.shp")
rdm.shp<-lapply(rdm.pts, function(x) readOGR(dsn,substr(x,1,nchar(x)-4)))

accu.Kmeans<-function(x,y){
  pred<-extract(y,x)
  ref<-as.numeric(x$ID)
  CM<-confusionMatrix(as.factor(ref), as.factor(pred))
  CM<-list(CM)
  return(CM)
}

accu<-mapply(function(x,y) accu.Kmeans(x,y), rdm.shp,bnk.class.filled)
accu.over<-lapply(accu,function(x) x$overall)
accu.over<-do.call("rbind",accu.over)
#remove first year, to equal length diff
#write.csv(accu.over,"accu_complete_Kmeans.csv",row.names = FALSE)


#remove first to set same length as DoDs
bnk.class.filled2<-bnk.class.filled[-1]


