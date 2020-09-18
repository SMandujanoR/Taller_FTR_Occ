# housekeeping


n <- 60   # number of sites
T <- 5    # number of primary periods
J <- 15    # number of secondary periods

site <- read.csv("covariates.csv",h=T)# read site covs
site.num<-(site[,2:6])
site.mean<-colMeans(site.num) # Mean of site
site.sd<-sapply(site[,2:6],sd) # Sd of site
#site.z<- (site.num-site.mean)/site.sd # why diferent ?
site.z<-scale(site[,2:6])# scale just numerical values


colnames(site.z)<-c("aspect.s","elevation.s", "slope.s", "canopy.s", "edge.s")
site<-cbind(site, site.z) # paste scaled values to site matrix

years <- data.frame(matrix(rep(2007:2012, each=n), n, T))
years <- data.frame(lapply(years, as.factor))

splist<-read.csv(file="splist.csv")
### select one species number from the list put in a and do not loop
#a<- 4 # change species number acording to the list splist.csv


 #for (a in 1:13){  ## loop to do all species, skip to do one by one species
 #####################  

  
  spname<-as.character(splist[a,]) # extract from table splist species(a)
  

  filesp<-paste(splist[a,1], ".csv", sep = "") # make file name for species(a)
  filespcov<-paste(splist[a,1], "_cov.csv", sep = "") # make file name for covariates of species(a) 
  obs_occasions<-read.csv(filesp,h=T) # read observations for species(a)
  cov_occations<-read.csv(filespcov,h=T) # read covariates of observations for species(a)
  
  occasions<- (cov_occations[,-1])
  y <- obs_occasions[,-1]
  
  umf <- unmarkedMultFrame(y=y,
                           siteCovs = data.frame(site=site),
                           obsCovs=list(occasion=occasions),
                           yearlySiteCovs=list(year=years),
                           numPrimary=T)
  
  #umf                                    # look at data
  #summary(umf)                           # summarize
  nobserv <- sum(rowSums(y, na.rm=T)) # get the number of observations
#   graphtitle <- paste(spname, "n=", nobserv, sep = " ")# make graph title
#   print(plot(umf, main=graphtitle, las=1, cex=1.5))  # plot umf
  
  #obsCovs(umf, matrices = TRUE)       # ver numero de dias as matrix
  
  m0 <- colext(~1, ~1, ~1, ~1, umf)  # fit a null
  m1 <- colext(~1, ~1, ~1, ~occasion, umf) # deteccion depende de numero de dias
  m2 <- colext(~1, ~year-1, ~1, ~1, umf) # year-dependent colonization
  m3 <- colext(~1, ~1, ~year-1, ~1, umf) # year-dependent extinction
  m4 <- colext(~1, ~year-1, ~year-1, ~1, umf) # year-dependent colonization and extinction
  m5 <- colext(~site.elevation.s,~1,~1,~1,umf) # ocupancy dependent of elevation   
  m6 <- colext(~1, ~year-1+I(year)^2,~1,~1,umf) # colonization dependent of time 2  
  m7 <- colext(~1,~1, ~year-1+I(year)^2,~1,umf) # extinction dependent of time 2 
  m8 <- colext(~site.canopy.s,~1,~1,~1,umf) # ocupancy dependent of canopy height  
  m9 <- colext(~1,~site.canopy.s,~1,~1,umf) # colonization dependent of canopy height  
  m10 <-colext(~1,~1,~site.canopy.s,~1,umf) # extinction dependent of canopy height  
  m11 <-colext(~site.elevation.s+site.canopy.s,~1,~1,~1,umf) # ocupancy dependent of elevation and canopy height  
  m12 <-colext(~site.aspect.s,~1,~1,~1,umf) # ocupancy dependent of Aspect
  m13 <-colext(~1,~year-1+I(year)^2,~year-1+I(year)^2,~1,umf) # colonization and extincion dependent of Yr^2
  m14 <-colext(~1,~year-1,~year-1,~occasion,umf) # extinction dependent dependent of Yr and detection of No. days
  m15 <-colext(~site.slope.s,~1,~1,~1,umf) # ocupancy dependent of Slope
  m16 <-colext(~1,~site.ForestType,~1,~1,umf) # colonization dependent of forest type
  m17 <-colext(~1,~1,~site.ForestType,~1,umf) # extinction dependent of forest type
  m18 <-colext(~site.slope.s+site.aspect.s,~1,~1,~1,umf) # ocupancy dependent of Slope
  m19 <-colext(~1,~1,~1,~site.canopy.s,umf) # detecion dependent of Canopy Height
  m20 <-colext(~site.elevation.s+I(site.elevation.s^2),~1,~1,~1,umf) # ocupancy dependent of elevation+Elev2  
  m21 <-colext(~site.edge.s,~1,~1,~1,umf) # ocupancy dependent of Edge
  m22 <-colext(~1,~site.edge.s,~1,~1,umf) # colonization dependent of Edge
  m23 <-colext(~1,~1,~site.edge.s,~1,umf) # extinction dependent of Edge
  m24 <-colext(~1,~1,~1,~site.edge.s,umf) # detection dependent of Edge
  m25 <-colext(~site.edge.s+I(site.edge.s)^2,~1,~1,~1,umf) # ocupancy dependent of Edge2
  m26 <- colext(~site.canopy.s,~year,~1,~1,umf) # ocupancy dependent of canopy height colonization Yr 
  m27 <- colext(~site.canopy.s,~1,~year,~1,umf) # colonization dependent of canopy height extinction Yr  
  m28 <- colext(~site.canopy.s+site.edge.s,~1,~1,~1,umf) # ocupancy dependent of canopy height colonization Yr
  m29 <- colext(~1,~year-1,~year-1,~site.canopy.s,umf) # detection dependent of canopy height colonization extinction Yr   
  m30 <- colext(~site.canopy.s,~1,~1,~site.edge.s,umf) # ocupancy dependent of canopy height detection dependent of Edge 
  m31 <- colext(~site.elevation.s+site.canopy.s+site.ForestType,~1,~1,~1,umf) # ocupancy dependent of elevation, canopy and forest type
  m32 <- colext(~site.slope.s+site.canopy.s+site.ForestType,~1,~1,~1,umf) # ocupancy dependent of Slope, canopy and forest type
  m33 <- colext(~site.aspect.s+site.canopy.s+site.ForestType,~1,~1,~1,umf) # ocupancy dependent of Aspect, canopy and forest type
  m34 <- colext(~1,~site.edge.s,~site.edge.s,~site.canopy.s,umf) #
  m35 <- colext(~site.elevation.s,~site.edge.s,~site.edge.s,~site.canopy.s,umf) #
  m36 <- colext(~1,~site.edge.s,~site.edge.s,~1,umf) # colonization extinction dependent of Edge
  m37 <- colext(~site.slope.s+site.canopy.s,~1,~1,~1,umf) # ocupancy dependent of Slope an canopy 
  m38 <- colext(~1,~year-1,~year-1,~year-1,umf) # colonization extinction and detection dependent of yr 
    
  
  models <- fitList(
    'psi(.)gam(.)eps(.)p(.)' = m0,
    'psi(.)gam(.)eps(.)p(ObsCovs)' = m1,
    'psi(.)gam(yr)eps(.)p(.)' = m2,
    'psi(.)gam(.)eps(yr)p(.)' = m3,
    'psi(.)gam(yr)eps(yr)p(.)' = m4,
    'psi(Elev)gam(.)eps(.)p(.)' = m5,
    'psi(.)gam(yr^2)eps(.)p(.)' = m6,
    'psi(.)gam(.)eps(yr^2)p(.)' = m7,
    'psi(Canopy)gam(.)eps(.)p(.)' = m8,
    'psi(.)gam(Canopy)eps(.)p(.)' = m9,
    'psi(.)gam(.)eps(Canopy)p(.)' = m10,
    'psi(Elev+Canopy)gam(.)eps(.)p(.)' = m11,
    'psi(Aspect)gam(.)eps(.)p(.)' = m12,
    'psi(.)gam(yr^2)eps(yr^2)p(.)' = m13,
    'psi(.)gam(Yr)eps(Yr)p(ObsCovs)' = m14,
    'psi(Slope)gam(.)eps(.)p(.)' = m15,
    'psi(.)gam(FType)eps(.)p(.)' = m16,
    'psi(.)gam(.)eps(FType)p(.)' = m17,
    'psi(Slope+Aspect)gam(.)eps(.)p(.)' = m18,
    'psi(.)gam(.)eps(.)p(Canopy)' = m19,
    'psi(Elevation^2)gam(.)eps(.)p(.)' = m20,
    'psi(Edge)gam(.)eps(.)p(.)' = m21,
    'psi(.)gam(Edge)eps(.)p(.)' = m22,
    'psi(.)gam(.)eps(Edge)p(.)' = m23,
    'psi(.)gam(.)eps(.)p(Edge)' = m24,
    'psi(Edge^2)gam(.)eps(.)p(.)' = m25,
    'psi(Canopy)gam(Yr)eps(.)p(.)' = m26,
    'psi(Canopy)gam(.)eps(Yr)p(.)' = m27,
    'psi(Canopy+Edge)gam(.)eps(.)p(.)' = m28,
    'psi(.)gam(Yr)eps(Yr)p(Canopy)' = m29,
    'psi(Canopy)gam(.)eps(.)p(Edge)' = m30,
    'psi(Elev+Canopy+FType)gam(.)eps(.)p(.)' = m31,
    'psi(Slope+Canopy+FType)gam(.)eps(.)p(.)' = m32,
    'psi(Aspect+Canopy+FType)gam(.)eps(.)p(.)' = m33,
    'psi(.)gam(Edge)eps(Edge)p(Canopy)' = m34,
    'psi(Elev)gam(Edge)eps(Edge)p(Canopy)' = m35,
    'psi(.)gam(Edge)eps(Edge)p(.)' = m36,
    'psi(Slope+Canopy)gam(.)eps(.)p(.)' = m37,
    'psi(.)gam(yr)eps(yr)p(yr)' = m38)

  ms <- modSel(models)
  ms
#######}# End Function
  
  
## Do stuff for model coeficients
#   coef(ms)
#   toExport <- as(ms, "data.frame")
  
  #This part store some models coeficients in a table (mat_models) to compare on screen
  ms_AIC_models<-as.data.frame(ms@ Full[1], row.names = NULL) #store model name
  ma_nPars<-as.data.frame(ms@Full$nPars) #store parameter number
  ms_AIC_values<- as.data.frame(ms@Full$AIC) #store AIC values
  ms_AIC_delta<- as.data.frame(ms@Full$delta) #store AIC delta values
  ms_AIC_AICwt<- as.data.frame(ms@Full$AICwt) #store AIC wt values
  ms_AIC_cumultw<-as.data.frame(ms@Full$cumltvWt) #store model name
  ms_m<-as.data.frame(row.names(ms_AIC_models)) #store m number
  ms_formula<- as.data.frame(ms@Full$formula) #store model formula
  mat_models <- cbind(ms_AIC_models, ma_nPars, ms_AIC_values, ms_AIC_delta, ms_AIC_AICwt, ms_AIC_cumultw) #paste in matrix
  colnames(mat_models)<-c("model", "nPars",'AIC', "delta", "AICwt", "cumltvWt") # change row names

  ##Print los 7 primeros modelos
  print(a)
  print(spname)
  print (mat_models[c(1:7),])

 #} # End Loop all species
  




  ###
  ### PART 3: Do some analysis of the results
  ###      ---  Bootstrap ---
  ###
  
  ## Does this model fit worth a darn?
  # Function returning three fit-statistics.
  # NOTE: na.rm=TRUE !!!!!!
  fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2,na.rm=TRUE)
    chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
  }
  
  ### Take some time
  model_par_boot <- parboot(m36, fitstats, nsim=100, report=1)# put the model number "best"
  model_par_boot # loot at Pr (probabilidad)
  plot(model_par_boot)
  
##
## Model does not fit!
## Three choices:
## 1. We can expand this model or tinker with components of it (NB, ZIP)
##   [takes a long time to run these models]
## 2. We can seek out an implausible data-generating model that fits.
##   [might satisfy referee to have good p-value]
## 3. We can proceed. 
##
##

  
###################################################
  ### code chunk number 12: unmarked.Rnw:207-220
  ### Only Chi square  ---  Bootstrap ---
  ###################################################
  chisq <- function(fm) {
    umf <- getData(fm)
    y <- getY(umf)
    y[y>1] <- 1
    sr <- fm@sitesRemoved
    if(length(sr)>0)
      y <- y[-sr,,drop=FALSE]
    fv <- fitted(fm, na.rm=TRUE)
    y[is.na(fv)] <- NA
    sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
  }
  
  ### Take some time
  model_par_boot_chi <- parboot(m8, statistic=chisq, nsim=100, report=1)
  model_par_boot_chi
  
  
  
#   #returns the estimate of PAO (proportion of sites occupied) at 95% confidence interval.
#   ### better Sugested aproach  WinBugs 
#   re <- ranef(m12)
#   EBUP <- bup(re, stat="mode")
#   CI <- confint(re, level=0.95)
#   rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / 130)
#   
  
  
  
  ###### for eac species only!!!! not for loop.
  #plot the lower AIC model in this case m19 for Pecary
 
  
  plot(1:5,m19@projected.mean[2,],ylim=c(0,1),t='b')
  plot(1:5,m36@projected.mean[1,],ylim=c(0,1),t='b')
  #lines(1:5,m19@smoothed.mean[2,],lty=2)
  #lines(1:5,m19@smoothed.mean[2,],lty=2)
  
  # we predict over range of values on that standarized scale
  # We must fix "canopy, aspect, slope, edge" at some arbitrary value (let's use the mean)
  newCanopy <- data.frame(site.canopy.s=seq(-2.4698, 1.7088, by=0.07)) #add site.aspect.s=0 if the model has aspect
  newAspect <- data.frame(site.aspect.s=seq(-1.3369, 1.7538, by=0.05), site.edge.s=0) #add site.aspect.s=0 if the model has aspect
  newSlope <- data.frame(site.slope.s=seq(-1.5129, 3.3447, by=0.05))
  newEdge <- data.frame(site.edge.s=seq(-1.7333, 2.5911, by=0.05))
  #colnames(newCovs)<- c("site.canopy.s", "date", "time")
  E.p <- predict(m19, type="det", newdata=newCanopy, appendData=TRUE)# type puede ser:det psi, col, ext
  head(E.p)
  
  # Plot it
  plot(Predicted ~ site.canopy.s, E.p, type="l", ylim=c(0,1),
       xlab="Canopy (standardized)",
       ylab="Expected Occupancy")
  lines(lower ~ site.canopy.s, E.p, type="l", col=gray(0.5))
  lines(upper ~ site.canopy.s, E.p, type="l", col=gray(0.5))
  
  # Plot it again, but this time convert the x-axis back to original scale
  plot(Predicted ~ site.canopy.s, E.p, type="l", ylim=c(0,1),
       xlab="Canopy (m)",
       ylab="Expected Occupancy",
       xaxt="n")
  xticks <- -1:2
  xlabs <- xticks*site.sd[4] + site.mean[4] #Use the mean and sd of original value to change label name
  axis(1, at=xticks, labels=round(xlabs, 1))
  lines(lower ~ site.canopy.s, E.p, type="l", col=gray(0.5))
  lines(upper ~ site.canopy.s, E.p, type="l", col=gray(0.5))
  
# what is primo elevation for the willow tit?
# put on calculus hat.....
#   quadratic response: logistic response?
#      y = a + b*x + c*x2
#   differentiate and set to 0:
#      dy/dx = b + 2*c*x = 0
#   solve
#      xopt = -b/(2*c)
#
  
  # Expected occupancy over range of "Aspect" Model No. 12 for Tapir
  newEdge<- data.frame(site.edge.s=seq(-2.8484, 263.6055, by=4.5))
  newEdgeCanopy<-cbind(newEdge,newCanopy)
  #E.p <-   predict(m8, type="det", newdata=newCovs, appendData=TRUE)# puede ser:det psi, col, ext
  E.psi <- predict(m12, type="psi", newdata=newAspect, appendData=TRUE) 
  head(E.psi)
  
  # Plot predictions with 95% CI
  plot(Predicted ~ site.aspect.s, E.psi, type="l", ylim=c(0,1),
       xlab="Aspect (standardized)",
       ylab="Expected occupancy")
  lines(lower ~ site.aspect.s, E.psi, type="l", col=gray(0.5))
  lines(upper ~ site.aspect.s, E.psi, type="l", col=gray(0.5))
  
  # Plot it again, but this time convert the x-axis back to original scale
  plot(Predicted ~ site.aspect.s, E.psi, type="l", ylim=c(0,1),
       xlab="Aspect",
       ylab="Expected Ocupancy",
       xaxt="n")
  xticks <- -1:2
  xlabs <- xticks*site.sd[1] + site.mean[1] #Use the mean and sd of original value to change label name
  axis(1, at=xticks, labels=round(xlabs, 1))
  lines(lower ~ site.aspect.s, E.psi, type="l", col=gray(0.5))
  lines(upper ~ site.aspect.s, E.psi, type="l", col=gray(0.5))
  

##########################
# Require Pkg(raster)
##########################

aspect_map <- raster("D:\\TEAM\\mapcovs\\aspect.txt")
slope_map <- raster("D:\\TEAM\\mapcovs\\slope.txt")
plot(aspect_map)
plot(slope_map)

# To evaluate the predictions under the model we 
# first we have to get the landscape covariates on the same scale as
# the model covariates. i.e., standardized.
#
aspect_map.sd<- cellStats(aspect_map, 'sd')
aspect_map.mean<- cellStats(aspect_map, 'mean')
slope_map.sd<- cellStats(slope_map, 'sd')
slope_map.mean<- cellStats(slope_map, 'mean')
g.aspect_map<-(aspect_map - aspect_map.mean)/aspect_map.sd #map Scaled by sd
g.slope_map<-(slope_map - slope_map.mean)/slope_map.sd #map Scaled by sd

st <- stack(g.slope_map, g.aspect_map)
sum(st)

#plot(g.slope_map)# el mismo mapa con unidades scaled
#### calculator del modelo
# s <- calc(g.aspect_map, fun=function(x){ x[x < 4] <- NA; return(x)} ) )
# as.matrix(s)
# 


### 
### PART 3c:
### SPATIAL ANALYSIS/PREDICTION
###
### Now lets make some really cool spatial predictions......
###
###
###

landscape <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/Swiss_landscape.csv?attredirects=0&d=1")

##VVB data
landscape1<-read.csv("D:\\TEAM\\mapcovs\\slope_asvector_xy.csv")
landscape2<-read.csv("D:\\TEAM\\mapcovs\\aspect_asvector_xy.csv")
landscape2<-cbind(landscape2[,],landscape1[,3])
colnames(landscape2)<-(c("FID","PID","aspect","X","Y","slope"))#invierte X y Y
head(landscape2)
#head(landscape)
## Note integer coordinates - row/column ids



###VBdata
gslope<- landscape2[,"slope"]   # median elevation of quadrat
gaspect<-landscape2[,"aspect"]
ggrid<-landscape2[,c("X","Y")]



# lets plot these variables to see how they look
#
# two options: (1) use my simpleton spatial.plot function
#              (2) stuff the data into a matrix and use image()
#
# grab utility functions including spatial.plot



###VB data
par(mar=c(3,3,3,5),mfrow=c(2,1))
spatial.plot(ggrid,gslope)
spatial.plot(ggrid,gaspect)
# this is cool

# this is even cooler (could also use "raster" package but thats too easy):



###VB data
xmin <- min(ggrid[, 1])
xmax <- max(ggrid[, 1])
ymin <- min(ggrid[, 2])
ymax <- max(ggrid[, 2])
z <- matrix(NA, nrow = length(xmin:xmax), ncol = length(ymin:ymax))


##
## stuff attribute into matrix where matrix element [i,j] = coordinate [x,y] of
## the attribute to be plotted
## In this case we show an example of making an image plot of the ELEVATION attribute:
##



####
####
####  END OF BASIC IMAGE MAKING STUFF
####  Lets repeat that process but using the model predictions of E[N]
####
#



# To evaluate the predictions under the model we 
# first we have to get the landscape covariates on the same scale as
# the model covariates. i.e., standardized.

###VB data
slope.mean<- mean(gslope)
slope.sd <- sd(gslope)
gslope<- (gelev - slope.mean)/slope.sd
aspect.mean<- mean(gaspect)
aspect.sd <- sd(gaspect)
gaspect<- (gaspect - aspect.mean)/aspect.sd


#raster
aspect_map.sd<- cellStats(aspect_map, 'sd')
aspect_map.mean<- cellStats(aspect_map, 'mean')
slope_map.sd<- cellStats(slope_map, 'sd')
slope_map.mean<- cellStats(slope_map, 'mean')
g.aspect_map<-(aspect_map - aspect_map.mean)/aspect_map.sd #map Scaled by sd
g.slope_map<-(slope_map - slope_map.mean)/slope_map.sd #map Scaled by sd

# # remember length = 0 is saturation sampling because length = 1/L
# #
# new<- data.frame(elev=gelev,forest=gforest,length=0)
# pred<-predict(fm10,type="state",newdata=new,appendData=TRUE)
###
### this destroys your computer -> bug in unmarked
###
### instead we have to do this the old-fashioned way:
### look at col names to figure order of parameters
###
par(mfrow=c(1,1))
betavec<-coef(m12)[1:5]
Xg<-cbind(rep(1,length(gaspect)),gaspect,1,1,1)
pred<-exp(Xg%*%(betavec)) #valor de predicho de Psi 

#png("predN.png",width=6,height=4, units="in", res=400)
##
## stuff predictions into our raster
##


par(mar=c(3,3,3,5),mfrow=c(2,1))
spatial.plot(ggrid,gaspect) #Drape on ggrid
spatial.plot(ggrid,pred)


