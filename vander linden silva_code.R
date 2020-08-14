####Vander Linden, Silva: Dispersals as demographic processes####
####oad libraries####
library(rcarbon) # calibration of 14C dates
library(rgdal) # mapping
library(maptools) # mapping
library(tidyverse) # data manipulation
library(changepoint) # changepoint analysis
library(rworldmap)
library(rnaturalearth) # download georeferenced data
library(sp)
library(raster) # raster creation and manipulation
library(ncdf4)
library(gstat) # geostatistics
library(tmap) # plotting maps
library(sf)
library(sp)
library(rgeos)

####calibration, SPDs and corresponding analyses####
#setup
ncores <- 1
nsim<-1000
realstartBCAD <- -10000
realendBCAD <- -3000
bracket <- 1000
workingstartBP <- abs(realstartBCAD-1950)+bracket
workingendBP <- abs(realendBCAD-1950)-bracket
if (workingendBP<0){ workingendBP <- 0 }
workingstartCRA <- uncalibrate(workingstartBP)$ccCRA
workingendCRA <- uncalibrate(workingendBP)$ccCRA
if (workingendCRA<0){ workingendCRA <- 0 }
#load dates
dates <- read.csv("/home/marc/Dropbox/m_pers/crossdem_2/submission/revision/spd.csv", header=TRUE)
dates<-mutate(dates,DateID=row_number()) # add id number to each entry
# calibrate dates
caldates <- calibrate(x=dates$C14AGE, errors=dates$C14STD,
                      calCurves='intcal13', timeRange=c(workingstartBP,workingendBP), normalised=FALSE, ncores=ncores, calMatrix=TRUE)
# sum dates
bins <- binPrep(sites=dates$SITE, ages=dates$C14AGE, h=100)
allspd <- spd(x=caldates, bins=bins, timeRange=c(workingstartBP,workingendBP), datenormalised=FALSE)

# Fit an exponential model and create a predicted SPD
set.seed(123)
expnull <- modelTest(caldates, errors=dates$C14STD, bins=bins, nsim=10,
                     timeRange=c(9000,6500), model="exponential", ncores=ncores, datenormalised=FALSE)

## Fit a logistic model and create a predicted SPD
time <- seq(max(caldates$metadata$CRA), min(caldates$metadata$CRA), -1) # calBP
fit <- drm(y ~ x, data=data.frame(x=time[501:8451],
                                  y=allspd$grid$PrDens[allspd$grid$calBP >= min(time) & allspd$grid$calBP <= max(time)]), fct=L.3())
pred <- as.numeric(predict(fit,data.frame(x=time)))
predgrid <- data.frame(calBP=time, PrDens=pred)

## Logistic model test
set.seed(123)
logisticmod <- modelTest(caldates, errors=dates$C14STD, bins=bins, nsim=1000, runm=50,
                         timeRange=c(workingstartBP,workingendBP), model="custom", predgrid=predgrid, ncores=ncores, datenormalised=FALSE)

# Generate a smoothed SPD
eurofarm.spd.smoothed = spd(caldates,timeRange=c(9000,6500),bins=DK.bins,runm=100)
# Start values should be adjusted depending on the observed SPD
logFit <- nls(PrDens~SSlogis(calBP, Asym, xmid, scale),data=eurofarm.spd.smoothed$grid,
              control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=7500,scale=-100))
# Generate a data frame containing the fitted values
logFitDens=data.frame(calBP=eurofarm.spd.smoothed$grid$calBP,PrDens=SSlogis(input=eurofarm.spd.smoothed$grid$calBP,
                                                  Asym=coefficients(logFit)[1],xmid=coefficients(logFit)[2],scal=coefficients(logFit)[3]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull <- modelTest(caldates, errors=alldates$C14SD, bins=bins,nsim=nsim,
                     timeRange=c(9000,6500), model="custom",predgrid=logFitDens, runm=100, raw=TRUE)


# Regional permutation test
set.seed(123)
permregions <- permTest(x=caldates, bins=bins, marks=dates$Region, timeRange=c(workingstartBP,workingendBP), runm=50, nsim=nsim)

# Load the logistic model test results
load("/home/marc/Dropbox/m_pers/buffalo_april18/eurofarm_logisticmod.rda")

#Plot spd
tiff("/home/marc/Dropbox/m_pers/crossdem_2/spd.tif",width=1200,height=800,compression="lzw")
ymax <- max(allspd$grid$PrDens)
plot(allspd, cal="BP",xlim=c(9000,6500),xaxt='n')
plot(allspd,runm=200,add=TRUE,type="simple",col="red",lwd=2,lty=2,xaxt='n')
text(x=8900, y=0.95*ymax, labels="All Dates (not-normalised)", font=2, cex=0.8, adj=c(0,0.7))
text(x=8900, y=0.9*ymax, labels=paste("n=",nrow(dates),", sites=",length(unique(dates$SITE)),", bins=",length(unique(bins)),sep=""),
     font=1, cex=0.7, adj=c(0,0.7))
legend(x=8900,y=0.85*ymax, legend=c("SPD"), lty=c("solid","solid"), lwd=c(3,0.5,1), col=c("grey90","grey50"), bty="n", cex=0.7)
axis(side=1, at=seq(4500,9000,500), labels=seq(4500,9000,500), las=2, cex.axis=1)
axis(side=1,at=seq(4500,9000,100),tick=TRUE,labels=FALSE)
dev.off()

# Plot exponential
tiff("/home/marc/Dropbox/m_pers/crossdem_2/exp.tif",compression="lzw",width=1200,height=765,units="px",pointsize=24)
ymax_exp<-max(expnull$result$PrDens)*1.1
plot(expnull, col.obs="darkred", lwd.obs=1, ylim=c(0,ymax_exp), calendar="BP", xlim=c(9000,6500))
lines(expnull$fit$calBP,expnull_v2$fit$PrDens, col="black", lty="dashed", lwd=0.5)
legend(x=8900, y=0.90*ymax_exp, legend=c("SPD (dates not normalised)","Exponential Model", "95% MC envelope","positive deviation",
                                 "negative deviation"),col=c(1, "black","lightgrey",rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),
       lty=c(1,2,1,1,1), lwd=c(0.5,0.5,5,5,5), cex=0.6, bty='n', bg="white", title="")
text(x=8900, y=0.95*ymax_exp, labels="Exponential Fit", font=2, cex=0.7, adj=c(0,0.7))
text(x=8900, y=0.6*ymax_exp, cex=0.6, adj=c(0,0.7), labels=substitute(paste(italic(p), "=", x, sep=""), list(x=round(expnull_v2$pval,4))))
dev.off()

# Plot logistic model
tiff("/home/marc/Dropbox/m_pers/crossdem_2/log.tif",compression="lzw",width=1200,height=765,units="px",pointsize=24)
ymax_log <- max(LogNull$result$PrDens)*1.1
plot(LogNull, col.obs="darkred", lwd.obs=1, ylim=c(0,ymax_log), calendar="BP", xlim=c(9000,6500))
lines(LogNull$fit$calBP, LogNull$fit$PrDens, lty="dotted", lwd=0.5)
legend(x=8900, y=0.90*ymax_log, legend=c("Archaeological 14C", "Fitted Logistic Model",
                                 "95% MC envelope","positive deviation","negative deviation"),
       col=c("darkred", "black","lightgrey",rgb(0.7,0,0,0.2), rgb(0,0,0.7,0.2)), lty=c(1,2,1,1,1),
       lwd=c(1,0.5,5,5,5), cex=0.8, bg="white", title="")
text(x=8900, y=0.88*ymax_exp, labels="Logistic Fit", font=2, cex=0.7, adj=c(0,0.7))
text(x=8900, y=0.5*ymax_log, labels=paste("arch.dates=",nrow(dates),", \nsites=",length(unique(dates$SITE)),", bins=",length(unique(bins)),sep=""), font=1, cex=0.8, adj=c(0,0))
text(x=8900, y=0.4*ymax_log, cex=0.8, adj=c(0,0.7),
     labels=substitute(paste(italic(p),"=", x, " (global, 1000 runs)", sep=""),list(x=round(logisticmod$pval,4))))
dev.off()

# Individual Regions
binh <- 200
check <- dates$Region=="Adriatic"
Addates <- dates[check,]
Adbins <- binPrep(sites=Addates$SITE, ages=Addates$C14AGE, h=binh)
Adengspd <- spd(x=caldates[check], bins=Adbins, timeRange=c(workingstartBP,workingendBP), datenormalised=FALSE)
check <- dates$Region=="Danube"
Dadates <- dates[check,]
Dabins <- binPrep(sites=Dadates$SITE, ages=Dadates$C14AGE, h=binh)
Daengspd <- spd(x=caldates[check], bins=Dabins, timeRange=c(workingstartBP,workingendBP), datenormalised=FALSE)
check <- dates$Region=="Greece"
Grdates <- dates[check,]
Grbins <- binPrep(sites=Grdates$SITE, ages=Grdates$C14AGE, h=binh)
Grspd <- spd(x=caldates[check], bins=Grbins, timeRange=c(workingstartBP,workingendBP), datenormalised=FALSE)


# plot regional permutations
## Plot
xlim <- c(-9000,-6500)
xts <- seq(-8000,2000,1000)
xtl <- seq(-8000,2000,1000)
#dev.new(device=pdf, width=4, height=6.5)
tiff("/home/marc/Dropbox/m_pers/crossdem_2/perm.tif",compression="lzw")
layout(matrix(c(1,2,3,4), 4, 1, byrow=TRUE), widths=4, heights=c(1.5,1.5,1.5,1.9))
par(mar=c(0, 2, 0.2, 1)) #c(bottom, left, top, right)
par(yaxs="i")
par(xaxs="i")
plotymax <- max(permregions$envelope[["1"]][,2], permregions$observed[["1"]]$PrDens)*1.1
tmp <- plot(permregions, focalm="1", calendar="BP", col.obs="red", lwd.obs=1, xlim=c(9000,6500), ylim=c(0,plotymax), xaxt='n',yaxt='n')
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=8900, y=plotymax*0.9, labels="A. Adriatic", font=2, cex=0.9, adj=c(0,0.7))
text(x=8900, y=plotymax*0.8, labels=paste("dates=",nrow(Addates),", sites=",length(unique(Addates$SITE)),", bins=",length(unique(Adbins)),sep=""), font=1, cex=0.8, adj=c(0,0.7))
text(x=8900, y=plotymax*0.7, cex=0.8, adj=c(0,0.7), labels=substitute(paste(italic(p), "=", x, " (global, ", nsim, " runs)", sep=""), list(x=round(min(permregions$pValueList["1"],1-permregions$pValueList["1"]),3),nsim=nsim)))
plotymax <- par("usr")[4]
legend(x=8900, y=plotymax*0.6,legend=c("SPD","95% MC envelope","positive deviation","negative deviation"),col=c("red","lightgrey",rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),lty=c(1,1,1,1),lwd=c(1,5,5,5),cex=0.8, bg="white")
axis(side=2, cex.axis=0.8, at=c(0,1))
par(mar=c(0, 2, 0, 1))
plotymax <- max(permregions$envelope[["2"]][,2], permregions$observed[["2"]]$PrDens)*1.1
plot(permregions, focalm="2", calendar="BP", col.obs="orange", lwd.obs=1, xlim=c(9000,6500), ylim=c(0,plotymax*1.1), xaxt='n',yaxt='n')
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=8900, y=plotymax*0.9, labels="B. Danube", font=2, cex=0.9, adj=c(0,0.7))
text(x=8900, y=plotymax*0.8, labels=paste("dates=",nrow(Dadates),", sites=",length(unique(Dadates$SITE)),", bins=",length(unique(Dabins)),sep=""), font=1, cex=0.8, adj=c(0,0.7))
text(x=8900, y=plotymax*0.7, cex=0.8, adj=c(0,0.7), labels=substitute(paste(italic(p), "=", x, " (global, ", nsim, " runs)", sep=""), list(x=round(min(permregions$pValueList["2"],1-permregions$pValueList["2"]),3),nsim=nsim)))
axis(side=2, cex.axis=0.8, at=c(0,1))
par(mar=c(3, 2, 0, 1))
plotymax <- max(permregions$envelope[["3"]][,2], permregions$observed[["3"]]$PrDens)*1.2
plot(permregions, focalm="3", calendar="BP", col.obs="blue", lwd.obs=1, xlim=c(9000,6500), ylim=c(0,plotymax), xaxt='n',yaxt='n')
axis(side=1, at=seq(-10000,2000,1000), labels=seq(-10000,2000,1000), las=2, cex.axis=0.8)
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=8900, y=plotymax*0.9, labels="C. Greece", font=2, cex=0.9, adj=c(0,0.7))
text(x=8900, y=plotymax*0.8, labels=paste("dates=",nrow(Grdates),", sites=",length(unique(Grdates$SITE)),", bins=",length(unique(Grbins)),sep=""), font=1, cex=0.8, adj=c(0,0.7))
text(x=8900, y=plotymax*0.7, cex=0.8, adj=c(0,0.7), labels=substitute(paste(italic(p), "=", x, " (global, ", nsim, " runs)", sep=""), list(x=round(min(permregions$pValueList["3"],1-permregions$pValueList["3"]),3),nsim=nsim)))
par(mar=c(3, 2, 0, 1))
axis(side=2, cex.axis=0.8, at=c(0,1))
par(yaxs="r")
dev.off()

#point-based permutation test
# Create a data.frame of site locations extracting spatial coordinates
sites <- unique(data.frame(id=dates$SiteID,lat=dates$LAT,lon=dates$LON))
rownames(sites) <- sites$id
sites <- sites[,-1]
# Convert to a SpatialPoints class object
sp::coordinates(sites) <- c("lon","lat")
sp::proj4string(sites) <- sp::CRS("+proj=longlat +datum=WGS84")
#Compute distance matrix
d <- sp::spDists(sites,sites,longlat=TRUE)
#Compute spatial weights
w <- spweights(d,h=100)
breaks <- seq(9000,7000,-500) #500 year blocks
timeRange <- c(9000,7000) #set the timerange of analysis in calBP, older date first
farm_spatial <- sptest(calDates=caldates, bins=bins,timeRange=timeRange, locations=sites,
                      permute="locations",nsim=1000, breaks=breaks, spatialweights=w,ncores=3) 
#plot
base <- getMap(resolution="low") #extract basemap
xrange <- bbox(sites)[1,]
yrange <- bbox(sites)[2,]
## Spatial Permutation Test for Transition 1
par(mar=c(1,1,4,1),mfrow=c(1,2))
plot(base,col="antiquewhite3", border="antiquewhite3",xlim=xrange, ylim=yrange,main="9-8.5 to 8.5-8 kBP \n (Test Results)")
plot(eurospatial,index=1, option="test", add=TRUE, legend=TRUE, legSize=0.7, location="topleft")
## Geometric Growth Rate for Transition 1
plot(base,col="antiquewhite3", border="antiquewhite3", xlim=xrange, ylim=yrange, main="9-8.5 to 8.5-8 kBP \n (Growth Rate)")
plot(eurospatial,index=1, option="raw", add=TRUE, breakRange=c(-0.005,0.005), legend=TRUE,legSize=0.7, location="topleft")

## Spatial Permutation Test for Transition 2
par(mar=c(1,1,4,1),mfrow=c(1,2))
plot(base,col="antiquewhite3", border="antiquewhite3",xlim=xrange, ylim=yrange,main="8.5-8 to 8-7.5 kBP \n (Test Results)")
plot(eurospatial,index=2, option="test", add=TRUE, legend=TRUE, legSize=0.7, location="topleft")
## Geometric Growth Rate for Transition 2
plot(base,col="antiquewhite3", border="antiquewhite3", xlim=xrange, ylim=yrange, main="8.5-8 to 8-7.5 kBP \n (Growth Rate)")
plot(eurospatial,index=2, option="raw", add=TRUE, breakRange=c(-0.005,0.005), legend=TRUE,legSize=0.7, location="topleft")

## Spatial Permutation Test for Transition 3
par(mar=c(1,1,4,1),mfrow=c(1,2))
plot(base,col="antiquewhite3", border="antiquewhite3",xlim=xrange, ylim=yrange,main="8-7.5 to 7 kBP \n (Test Results)")
plot(eurospatial,index=3, option="test", add=TRUE, legend=TRUE, legSize=0.7, location="topleft")
## Geometric Growth Rate for Transition 3
plot(base,col="antiquewhite3", border="antiquewhite3", xlim=xrange, ylim=yrange, main="8-7.5 to 7 kBP \n (Growth Rate)")
plot(eurospatial,index=3, option="raw", add=TRUE, breakRange=c(-0.005,0.005), legend=TRUE,legSize=0.7, location="topleft")



####changepoint analysis####
dens<-eurofarm.spd.smoothed$grid[,2] # extract density from a CalDates object
cpt<-cpt.mean(dens,method="BinSeg")
plot(cpt)
#recreate result of analysis as tibble for plotting in ggplot2
# $time will be used as x axis and needs to be calculated using timeRange used to do calibration in rcarbon
cpt.df<-tibble(time=c(9500,(9500-cpt@cpts[1]),
                      (9500-cpt@cpts[1]),(9500-cpt@cpts[2]),
                      (9500-cpt@cpts[2]),6500),
               cpt=c(cpt@param.est$mean[1],cpt@param.est$mean[1],
                     cpt@param.est$mean[2],cpt@param.est$mean[2],
                     cpt@param.est$mean[3], cpt@param.est$mean[3]))

dens<-runMean(allspd$grid[,2],200)
ggspd<-tibble(time=allspd$grid[,1],dens=dens)
ggspd$time<-1950-ggspd$time
ggspd$dens<-as.numeric(ggspd$dens)
ggspd<-ggspd%>%drop_na(dens)
ggspd<-gather(ggspd,'dens',key="type",value="density")
ggspd$dens<-as.factor(ggspd$density)
g2<-ggplot()+
  geom_segment(aes(x=cpt.df$time[1],y=cpt.df$cpt[1],xend=cpt.df$time[2],yend=cpt.df$cpt[2]),colour="red")+
  geom_segment(aes(x=cpt.df$time[3],y=cpt.df$cpt[3],xend=cpt.df$time[4],yend=cpt.df$cpt[4]),colour="red")+
  geom_segment(aes(x=cpt.df$time[5],y=cpt.df$cpt[5],xend=cpt.df$time[6],yend=cpt.df$cpt[6]),colour="red")+
  geom_line(data=ggspd,aes(x=time,y=density,colour=type),show.legend=FALSE)
g2
g2+
  xlim(9000,6500)+
  scale_x_reverse()+
  scale_color_manual(values=c("grey","black"))+
  scale_alpha_manual(values=c(0.5,1))+
  xlab("Time (calBP)")+
  ylab("Posterior Density")+
  theme_bw()

####kriging####
####load data####
dates.kriging<-filter(dates,CULTURE!="Mesolithic"&CULTURE!="Epipalaeolithic"&CULTURE!="Palaeolithic") # Neolithic only dates
dates.kriging<-filter(dates.kriging,LABNR!="OxA-8610"&LABNR!="OxA-8504") # remove outliers
dates.kriging.highqual<-filter(dates.kriging,SAMPLE=="acorn"|SAMPLE=="bone"|SAMPLE=="tooth"|SAMPLE=="antler"|SAMPLE=="antler tool"|
                             SAMPLE=="bone tool"|SAMPLE=="burnt bone"|SAMPLE=="burnt chaff"|SAMPLE=="cereals"|SAMPLE=="charred bones"|
                             SAMPLE=="charred cereals"|SAMPLE=="charred grain"|SAMPLE=="charred grains"|SAMPLE=="charred seed"|
                             SAMPLE=="charred seeds"|SAMPLE=="charred textile remains"|SAMPLE=="grain"|SAMPLE=="gain (whole)"|
                             SAMPLE=="grains"|SAMPLE=="hazelnut shell"|SAMPLE=="horn"|SAMPLE=="pea"|SAMPLE=="seed"|SAMPLE=="seeds") # for kriging on high-quality samples only
kriging.cal<-calibrate(x=dates.kriging$C14AGE,errors=dates.kriging$C14STD,calCurves='intcal13', normalised=FALSE,verbose=FALSE) # calibrate 14C dates
dates.kriging$median<-medCal(kriging.cal) # add median of calDates to dataframe
kriging.highqual.cal<-calibrate(x=dates.kriging.highqual$C14AGE,errors=dates.kriging.highqual$C14STD,
                                calCurves='intcal13', normalised=FALSE,verbose=FALSE) # high quality only
dates.kriging.highqual$median<-medCal(kriging.highqual.cal) # high quality only
####converting to a spatial dataframe####
coordinates(dates.kriging)<-~LON+LAT
proj4string(dates.kriging) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
coordinates(dates.kriging.highqual)<-~LON+LAT # high quality only
proj4string(dates.kriging.highqual) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # high quality only
etrs<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
dates.kriging.v2<-spTransform(dates.kriging,CRS=etrs)
dates.kriging.highqual.v2<-spTransform(dates.kriging.highqual,CRS=etrs) # high quality only
####create raster mask####
land<-ne_download(scale = 50,category = 'physical',type='land')
land<-spTransform(land,CRS="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
land.etrs<-spTransform(land,CRS=etrs)
box <- as(raster::extent(dates.kriging.v2@bbox), "SpatialPolygons")
proj4string(box) <- etrs
box<-spTransform(box,CRS=etrs)
mask<-gIntersection(box,land.etrs,byid=FALSE,drop_lower_td=FALSE)
grid<-makegrid(box,cellsize = 5000)
coordinates(grid)<-~x1+x2
proj4string(grid)<-etrs
grid<-grid[mask, ]
####sample 14C dates####
r<-raster(box,resolution=50000)
rp<-as(r,'SpatialPolygons')
proj4string(rp)<-etrs
rp<-spTransform(rp,CRS=etrs)
sample<-over(dates.kriging.v2,rp)
dates.kriging.v2$idpol<-sample
dates.kriging.v3<-as.data.frame(dates.kriging.v2)
max.kriging<-group_by(dates.kriging.v3,idpol)%>%
  slice(which.max(median))
max.map<-tibble(lon=max.kriging$LON,lat=max.kriging$LAT)
coordinates(max.kriging)<-~LON+LAT
proj4string(max.kriging)<-etrs
sample.highqual<-over(dates.kriging.highqual.v2,rp)#high quality only
dates.kriging.highqual.v2$idpol<-sample.highqual#high quality only
dates.kriging.highqual.v3<-as.data.frame(dates.kriging.highqual.v2) #high quality only
max.kriging.highqual<-group_by(dates.kriging.highqual.v3,idpol)%>%
  slice(which.max(median)) #high quality only
max.map.highqual<-tibble(lon=max.kriging.highqual$LON,lat=max.kriging.highqual$LAT)#high quality only
coordinates(max.kriging.highqual)<-~LON+LAT#high quality only
proj4string(max.kriging.highqual)<-etrs#high quality only
####kriging####
v.kriging<-variogram(median~1,max.kriging) # define variogram
v.fit<-fit.variogram(v.kriging,vgm(c("Exp","Mat","Ste","Sph")),fit.kappa=TRUE) # fit variogram to data
krg<-krige(median~1,max.kriging,grid,v.fit) # kriging
base.krg<-grid<-raster(box,resolution = 5000,
                       crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
r.krg<-rasterize(krg,base.krg,field=krg$var1.pred)

####kriging high quality####
v.kriging.highqual<-variogram(median~1,max.kriging.highqual) # define variogram
v.fit.highqual<-fit.variogram(v.kriging.highqual,vgm(c("Exp","Mat","Ste","Sph")),fit.kappa=TRUE) # fit variogram to data
krg.highqual<-krige(median~1,max.kriging.highqual,grid,v.fit.highqual) # kriging
r.krg.highqual<-rasterize(krg.highqual,base.krg,field=krg.highqual$var1.pred)

####plot results####
st.mask<-st_as_sf(mask)
points.kriging<-st_as_sf(max.map,coords=c("lon","lat"),crs=etrs)
krig.map<-tm_shape(r.krg)+
  tm_raster(style="equal",n=10,title="Predicted values \n (calBP)",palette="YlOrRd")+
  tm_legend(outside=TRUE)+
  tm_shape(st.mask)+
  tm_borders(col="black")+
  tm_shape(points.kriging)+
  tm_dots(shape=3,size=0.1)+
  tm_scale_bar(breaks=c(0,100,200,300),position=c(0.02,0.05))
krig.map
tmap_save(krig.map,"/home/marc/Dropbox/m_pers/crossdem_2/krig.tiff",width=1200,height=800)

####plot results highqhual####
st.mask<-st_as_sf(mask)
points.kriging.highqual<-st_as_sf(max.map.highqual,coords=c("lon","lat"),crs=etrs)
krig.highqual.map<-tm_shape(r.krg.highqual)+
  tm_raster(style="equal",n=10,title="Predicted values \n (calBP)",palette="YlOrRd")+
  tm_legend(outside=TRUE)+
  tm_shape(st.mask)+
  tm_borders(col="black")+
  tm_shape(points.kriging.highqual)+
  tm_dots(shape=3,size=0.1)+
  tm_scale_bar(breaks=c(0,100,200,300),position=c(0.02,0.05))
krig.highqual.map
tmap_save(krig.highqual.map,"/home/marc/Dropbox/m_pers/crossdem_2/krig_highqual.tiff",width=1200,height=800)
    