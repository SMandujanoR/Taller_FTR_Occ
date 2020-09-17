# Several functions require the code developed by Jorge Ahumada which is appended below. 
# Daniel Spitale, Nov 2014.
# update Jan 2015: changed "Sampling.Period" with "Sampling.Event" in the variable names
# update Feb 2015: new functions and tidy
# update Dec 2015: PNAB

# the function joins in a single command three functions which adjust the data
fix.dta <- function(dtaframe) {
    datafix <- f.fix.data(dtaframe)
    datacorrect <- f.correct.DF(datafix)
    dataord <- f.order.data(datacorrect)
}


# the function returns the matrix camera by species by events (independent events set with f.separate.events function)
event.sp <- function(dtaframe, year, thresh) {
    # thresh= minutes
    require(reshape)
    hr <- f.separate.events(dtaframe, thresh)
    sel <- subset(hr, select = c(Sampling.Event, Sampling.Unit.Name, Photo.Date, bin, grp))
    del <- unique(sel)
    dta <- rename(del, c(bin = "value"))
    yrsel <- dta[dta$Sampling.Event == year, ]
    events <- cast(yrsel, Sampling.Unit.Name ~ value, length)
}


# the function returns the matrix sampling effort by sampling unit (camera-days) per year
cam.days <- function(dtaframe, year) {
    yr <- dtaframe[dtaframe$Sampling.Event == year, ]
    yr$ndays <- as.numeric(difftime(yr$End.Date, yr$Start.Date))
    selvar <- subset(yr, select = c(Sampling.Unit.Name, Start.Date, End.Date, ndays))
    cam.days <- unique(selvar)
}


# the function returns the matrix species by hours to be used for activity profile (it uses all the available years)
events.hours <- function(dtaframe) {
    require(reshape)
    require(chron)
    dtaframe$ore <- trunc(dtaframe$Photo.Time, "hours")
    selvar <- subset(dtaframe, select = c(Sampling.Unit.Name, Photo.Date, bin, ore))
    ev.hr <- unique(selvar)
    ev.hr <- rename(ev.hr, c(ore = "value"))
    ev.sp <- cast(ev.hr, value ~ bin, length)
    colnames(ev.sp) <- sub(" ", ".", colnames(ev.sp))
    return(ev.sp)
}

# the function returns the observed species accumulation curve for a given year
acc.curve <- function(dtaframe, year) {
    require(reshape)
    require(vegan)
    yr <- dtaframe[dtaframe$Sampling.Event == year, ]
    mat <- f.matrix.creator(yr)
    pr <- melt(mat)
    colnames(pr) <- c("Sampling.Unit.Name", "Date", "value", "species")
    ev.sp <- cast(na.omit(pr), Sampling.Unit.Name + Date ~ species, sum)
    ac <- specaccum(ev.sp[, -c(1:2)], method = "random", permutations = 100)
    mt <- data.frame(ac$sites, ac$richness, ac$sd)
    colnames(mt) <- c("Camera.trap.days", "species", "sd")
    return(mt)
}

# the function prepares data for multiseason analyses
# adjMult works on list AFTER the function shrink

adjMult<-function(matr) {

    dlist <- llply(matr, .fun = dim)
    rw <- sapply(dlist, "[[", 2)
    mrw <- max(rw)

    adc <- function(yr) {
        madj <- cbind(yr, matrix(NA, nrow = nrow(yr), ncol = mrw - ncol(yr)))
    }

    lad <- llply(matr, .fun = adc)

    for (i in 1:length(lad)) {
        colnames(lad[[i]]) <- paste(names(lad)[i], rep(1:ncol(lad[[i]])), sep = ".")
    }

    addcam <- function(lad) {
        dlad <- as.data.frame(lad)
        dlad$cam <- rownames(dlad)
        dlad
    }

    spcra <- llply(lad, .fun = addcam)

    joined <- join_all(spcra, by = "cam", type = "full", match = "all")
    rownames(joined) <- joined$cam
    joined$cam <- NULL

    joined
}

# adjMult2 works on list BEFORE the function shrink
adjMult2 <- function(matr) {

    dlist <- llply(matr, .fun = dim)
    rw <- sapply(dlist, "[[", 2)
    mrw <- max(rw)

    adc <- function(yr) {
        madj <- cbind(yr, matrix(NA, nrow = nrow(yr), ncol = mrw - ncol(yr)))
        rown <- rownames(yr)
        nd <- as.Date(colnames(yr))
        coln <- seq.Date(nd[1], by = "day", length.out = ncol(madj))
        dimnames(madj) <- list(rown, coln)
        madj
    }

    lad <- llply(matr, .fun = adc)

    addcam <- function(lad) {
        dlad <- as.data.frame(lad)
        dlad$cam <- rownames(dlad)
        dlad
    }

    spcra <- llply(lad, .fun = addcam)

    joined <- spcra[[1]]
    for (i in 2:length(spcra)) {
        joined <- merge(joined, spcra[[i]], by.x = "cam", by.y = "cam", all = T)
    }

    rownames(joined) <- joined$cam
    joined$cam <- NULL

    joined
}

# the function return a shrinked matrix collapsed using nday; if necessary, an X number of columns filled with NA are added
# to adjust the size of the shrinked matrix; be careful that nday makes sense for the size of the matrix, so that
# not many columns of NA are added 
shrink <- function(matrice, nday) {
    dy <- nday
    while (dy < ncol(matrice)) {
        dy <- dy + nday
    }
    addcol <- dy - ncol(matrice)
    if (addcol != 0) {
        matNA <- matrix(NA, nrow = nrow(matrice), ncol = addcol)
        matrice <- data.frame(matrice, matNA)
    }

    period <- ncol(matrice)/nday
    newday <- rep(1:period, each = nday)

    shr <- function(vec) {
        nav <- is.na(vec)
        dom <- all(nav == T)
        if (dom == T) {
            y <- NA
        } else {
            abb <- sum(vec, na.rm = T)
            y <- ifelse(abb == 0, 0, 1)
        }
        return(y)
    }

    matday <- data.frame(newday, t(matrice))
    shrmat <- t(aggregate(matday[, -1], list(matday$newday), shr))

    return(shrmat[-1, ])
}

# the function returns occupancy and detection probability in the 0-1 interval; to do that it sets covs =0; 
# in practice, the function uses only the intercept of the model
coef01 <- function(modello) {
    otb <- coef(modello, type = "state")
    psi.int <- otb[1]
    SE.psi.int <- SE(modello[1])
    SE.psi.int[1]

    dtb <- coef(modello, type = "det")
    p.int <- dtb[1]
    SE.p.int <- SE(modello[2])
    SE.p.int[1]

    occ01 <- exp(psi.int)/(1 + exp(psi.int))
    occSE <- sqrt(((exp(psi.int)/((1 + exp(psi.int))^2))^2) * SE.psi.int[1]^2)

    detp01 <- exp(p.int)/(1 + exp(p.int))
    detpSE <- sqrt(((exp(p.int)/((1 + exp(p.int))^2))^2) * SE.p.int[1]^2)
    ris <- c(occ01, occSE, detp01, detpSE)
    names(ris) <- c("psi", "SEpsi", "p", "SEp")
    return(ris)
}


# "covs.list" returns a list of covariables ready to be used in the function "unmarkedFrameOccu" (in the argument "obsCovs");
# it deletes the camera not matching with the cameras in the species matrix (check name consistence between them!)
# it standardizes the numerical covariables leaving unaffected the nominal (if any)
# the function should be used in analyzing species richness in the occupancy framework
covs.list <- function(matcov, matsp) {
    require(vegan)
    delcam <- which(matcov$Sampling.Unit.Name %in% colnames(matsp))  # matcov must have the variable 'Sampling.Unit.Name'
    covnum <- matcov[delcam, sapply(matcov, is.numeric)]
    covsy <- decostand(covnum, method = "standardize")

    if (ncol(covnum) != ncol(matcov)) {
        nom.var <- matcov[, sapply(matcov, is.factor)]
        covsy <- data.frame(covsy, nom.var[delcam, ])
    }

    creamat <- function(vet) {
        vc <- rep(vet, nrow(matsp))
        mat <- matrix(vc, nrow = nrow(matsp), byrow = T)
        colnames(mat) <- colnames(matsp)
        return(mat)
    }

    lista <- lapply(covsy, creamat)
    return(lista)
}


# the function returns the naive occupancy for the list of species (matrix.list as returned by f.matrix.creator)
naive <- function(matrix.list) {

    naive.f <- function(matrice) {
        somma <- rowSums(matrice, na.rm = T)
        res <- sum(somma > 0)/length(somma)
        return(res)
    }

    naive.list <- lapply(matrix.list, naive.f)
    n.table <- data.frame(species = names(naive.list), naive = unlist(naive.list))
    o.table <- n.table[order(n.table$species), ]
    rownames(o.table) <- NULL
    return(o.table)
}



#############################################################################
# the following code is by Jorge Ahumada
#############################################################################

#script to process raw TEAM files and get them ready for analysis
f.readin.fix.data<-function(){
	require(chron)
	data<-read.csv(file.choose(),h=T,skip=62)
	data<-f.fix.data(data)
	#make sure date info makes sense
	data
	}


#function to create binary matrices for all species at a site and sampling period. Matrix has a 1 if the species was seen in a day a 0 if not seen and NA if not sampled
#The function requires data from one sampling event and will return a list composed of 0,1 matrices, one matrix for each species.

#THIS FUNCTION WORKS WITH NEW TEAM DATA ONLY - do not use with legacy TEAM data
# this works one year at a time. Separate data in different years first
f.matrix.creator<-function(data){
	#results object
	res<-list()
	
	#get the dimensions of the matrix
	
	#list if sanpling units
	cams<-unique(data$Sampling.Unit.Name)
	cams<-sort(cams)
	rows<-length(cams)
	#start and end dates of sampling periods
	min<-min(data$Start.Date)
	max<-max(data$End.Date)
	cols<-max-min+1
	
	#sampling period
	date.header<-seq(from=min,to=max, by=1)
	mat<-matrix(NA,rows,cols,dimnames=list(cams,as.character(date.header)))
	
	#for all cameras, determine the open and close date and mark in the matrix
	start.dates<-tapply(as.character(data$Start.Date),data$Sampling.Unit.Name,unique)
	end.dates<-tapply(as.character(data$End.Date),data$Sampling.Unit.Name,unique)
	
	#outline the sampling periods for each camera j
	for(j in 1:length(start.dates)){
		#for each camera beginning and end of sampling
		low<-which(date.header==start.dates[j])
		hi<-which(date.header==end.dates[j])
		indx<-seq(from=low,to=hi)
		mat[j,indx]<-0
		}
		mat.template<-mat
				#get the species
		species<-unique(data$bin)
		#construct the matrix for each species i
		for(i in 1:length(species)){
			indx<-which(data$bin==species[i])
			#dates and cameras when/where the species was photographed
			dates<-data$Photo.Date[indx]
			cameras<-data$Sampling.Unit.Name[indx]
			dates.cameras<-data.frame(dates,cameras)
			#unique combination of dates and cameras 
			dates.cameras<-unique(dates.cameras)
			#fill in the matrix
			for(j in 1:length(dates.cameras[,1])){
				col<-which(date.header==dates.cameras[j,1])
				row<-which(cams==dates.cameras[j,2])
				mat[row,col]<-1
				}
			mat.nas<-is.na(mat)
			sum.nas<-apply(mat.nas,2,sum)
			indx.nas<-which(sum.nas==rows)
			if(length(indx.nas)>0){
			mat<-mat[,-indx.nas]
			}
	
			res<-c(res,list(mat))
			#return the matrix to its original form
			mat<-mat.template
			}
			
		names(res)<-species
		#res<-lapply(res,f.dum)
		res
	
	}
	

f.check.NA.breaks<-function(vector){
	notna<-which(!is.na(vector))
	if(min(notna)+length(notna)-1==max(notna)) print("ok")
	else print("aggh")
	}
	
f.start.minus.end<-function(data){
	data$End.Date-data$Start.Date
	}
f.start<-function(data){
	data$Start.Date
	}
f.end<-function(data){
	data$End.Date
	}

f.picture.dates<-function(data){
	data$Photo.Date
	}
f.picture.span <-function(data){
	max(data)-min(data)
	}
f.picture.min<-function(data){
	min(data)
	}
f.picture.max<-function(data){
	max(data)
	}


#function to fix the start/stop time of a camera if it is incorrectly entered	
f.start.stop.date.fixer<-function(data){
	
	cam.start.date<-by(data,data$Sampling.Unit.Name,f.start)
	cam.start.date<-lapply(cam.start.date,unique)
	cam.end.date<-by(data,data$Sampling.Unit.Name,f.end)
	cam.end.date<-lapply(cam.end.date,unique)
	
	#cam.span<-(by(data,data$Sampling.Unit.Name,f.start.minus.end))
	#cam.span<-lapply(cam.span,unique)
	
	pic.span<-by(data,data$Sampling.Unit.Name,f.picture.dates)
	min.pic<-lapply(pic.span,f.picture.min)
	max.pic<-lapply(pic.span,f.picture.max)
	#pic.span<-lapply(pic.span,f.picture.span)
	
	indx<-which(as.numeric(cam.start.date)-as.numeric(min.pic)>0 |as.numeric(cam.end.date)-as.numeric(max.pic)<0)
	#figure out which camera has the problem
	#indx<-which(as.numeric(cam.span)-as.numeric(pic.span)<=0)
	if(length(indx)){
	{cam.id<-names(pic.span)[indx]
	print("There are problems with the following cameras:")
	print(cam.id)}
	#for(i in 1:length(indx)){
	#	index<-which(data$Sampling.Unit.Name==cam.id[i])
	#	data$Start.Date[index]<-min.pic[[indx[i]]]
	#	data$End.Date[index]<-max.pic[[indx[i]]]
	#	
	#	}
		}
		else
			print("No problems detected..")
	#data
	}	


#function to convert a list of sampling matrices generated by f.matrix.creator into a data frame that can be used by the unmarked package	
f.convert.to.unmarked<-function(list){
require(unmarked)
nspecies<-length(list)
nrows<-dim(list[[1]])[1]
ncols<-dim(list[[1]])[2]
oldmat<-list()

for(i in 1:nspecies){
	mat<-rbind(oldmat,list[[i]])
	oldmat<-mat
	}	
y<-as.matrix(mat[,-ncols])
rownames(y)<-NULL
colnames(y)<-NULL
species<-gl(n=nspecies,k=nrows,labels=names(list))
siteCovs<-as.data.frame(species)
unmarkedFrameOccu(y=y,siteCovs=siteCovs)	
	}
	
	f.correct.DF<-function(DF){
ind <- sapply(DF, is.factor)
DF[ind] <- lapply(DF[ind], "[", drop=TRUE)
DF
	}

f.fix.data <- function(data){
#This function converts the dates and times in the data into date and time
  #objects that can be manipulated later
  #remove the Sampling date records
  indx<-which(as.character(data$Photo.Date)=="")
  #sdrecs<-data[indx,]
  #data<-data[-indx,]
	data$Photo.Date[indx]<-NA
  data$Photo.Time[indx]<-NA
  data$Photo.Date<-as.Date(data$Photo.Date)
	data$Photo.Time<-times(data$Photo.Time)
	#this line stores the date and time info for each photo into a single object
  time.date<-chron(dates=as.character(data$Photo.Date),times=data$Photo.Time,
                   format=c(dates="y-m-d",times="h:m:s"))
	#split the date from the time for the Camera start date and time
  qwe<-strsplit(as.character(data$Camera.Start.Date.and.Time)," ",fixed=T)
	#extract just the date
  qwe2<-lapply(qwe,function(x){x[1]})
  qwe2<-factor(as.character(qwe2))
	#Put again as a date object in a new column(below)
  qwe2<-as.Date(qwe2)
  #extract times
  qwe3<-lapply(qwe,function(x){x[2]})
	#convert them to a character
  qwe3<-as.character(qwe3)
  #store in a chron object
	qwe3<-chron(times=qwe3,format=c(times="h:m:s"))
  #combine Start date and time into a chron object and replace original
  data$Camera.Start.Date.and.Time<-chron(dates=as.character(qwe2),times=qwe3,
                                         format=c(dates="y-m-d",times="h:m:s"))
  #attach the start date of each camera trap to the data frame
  data<-data.frame(data,Start.Date=qwe2)
	
	#Now do the same but for the End date and time of each camera trap
  qwe<-strsplit(as.character(data$Camera.End.Date.and.Time)," ",fixed=T)
	qwe2<-lapply(qwe,function(x){x[1]})
	qwe2<-factor(as.character(qwe2))
	qwe2<-as.Date(qwe2)
  qwe3<-lapply(qwe,function(x){x[2]})
	qwe3<-as.character(qwe3)
	qwe3<-chron(times=qwe3,format=c(times="h:m:s"))
	data$Camera.End.Date.and.Time<-chron(dates=as.character(qwe2),times=qwe3,
	                                     format=c(dates="y-m-d",times="h:m:s"))
  data<-data.frame(data,End.Date=qwe2)
	#create new variable with binomial - genus species
  bin<-paste(data$Genus,data$Species)
	data<-data.frame(data,bin=bin,td.photo=time.date)
  #sdrecs<-data.frame(sdrecs,End.Date=NA,bin=NA,td.photo=NA)
	#data<-rbind(data,sdrecs)
  data
	}
f.dum<-function(data){
	dum<-apply(data,1,sum,na.rm=T)
	dum<-ifelse(dum>0,0,1)
	data<-data.frame(data,dum=dum)
	data
	}
f.extract.rare.sp<-function(raredata,alldata){
	spnum<-length(raredata[,1])
	oldindx<-numeric()
	for(i in 1:spnum){
		
		indx<-c(oldindx,which(alldata@siteCovs$species==as.character(raredata[i,1])))
		oldindx<-indx
		}
	sp<-siteCovs(alldata)[indx,]
	sp<-factor(sp)
	newufframe<-unmarkedFrameOccu(y=getY(alldata)[indx,],siteCovs=data.frame(species=sp))
	newufframe
	}

#Separate independent photographic events for a species in a given camera trap and date. thresh gives the threshold for considering events separate
#thresh is in minutes
f.separate<-function(data,thresh){
	
	#diff(data$td.photo)
	l<-length(data)
	interval<-diff(data)#abs(c(data[2:l],NA)-data)
	interval<-interval*1440 #convert to minutes
	interval<-as.numeric(interval)
  ev<-1;res<-numeric()
	cond<-interval>thresh #about 5 minutes in between events
	for(i in 1:(l-1)){
		if(!cond[i]) ev<-ev
		else ev<-ev+1
		res<-c(res,ev)
		
		}
	c(1,res)
	}
#test function; not usually used	
f.test.sep<-function(cond){
	l<-length(cond)
	#interval<-c(data$Photo.Time[2:l],NA)-data$Photo.Time
	ev<-1;res<-numeric()
	#cond<-interval>5
	for(i in 1:(l-1)){
		if(!cond[i]) ev<-ev
		else ev<-ev+1
		res<-c(res,ev)
		
		}
	c(1,res)

	
	}
	
#Order the data by Sampling unit name and photo raw name. This will order images chronologically
f.order.data<-function(data){
	indx<-order(data$Sampling.Event,data$Sampling.Unit.Name,data$td.photo)
	data<-data[indx,]
	data
	}
#function to separate independent events, extract from the list and paste together with the data set.
#This function removes records that are NOT images.. e.g. Sampling Date records
f.separate.events<-function(data,thresh){
	
	#e.data<-by(data$td.photo,data$Sampling.Unit.Name,f.separate,thresh)
  indx<-which(is.na(data$td.photo))
  if(length(indx)>0)
    data<-data[-indx,]
  e.data<-f.separate(data$td.photo,thresh)
#e.data<-data.frame(grp=unlist(e.data))
data.frame(data,grp=paste(data$Sampling.Event,".",data$Sampling.Unit.Name,".",e.data,sep=""))

	}
#Simulation to explore the effect of changing the threshold on the number
# of independent events. Thresh range is given as a sequence in mins
f.sim.thres<-function(data,threshRange){
  qwe<-data[,-42]
  res<-data.frame(thresh=threshRange,n.events=NA)
  for(i in threshRange){
  qwe<-f.separate.events(qwe,threshRange[i])
  res[i,2]<-length(unique(qwe[,42]))
  qwe<-qwe[,-42]
  qwe<-f.correct.DF(qwe)
  }
  plot(res[,1],res[,2],xlab="Threshold (min)",ylab="Number of events",type='b')
  res
  }
  
#convert Farenheit to Celsius

f.FtoC<-function(temp) {
  round(5/9*(temp-32))
}

#extract temperatures for a given species and graph
f.extract.temp<-function(data,species) {
  qwe<-data[data$bin==species,]
  qwe<-f.correct.DF(qwe)
  res<-as.numeric(by(qwe$Temp,qwe$grp,mean))
  res
}

f.graph.temp<-function(species,spname,nbins) {
  par(lwd=2)
  truehist(species,xlim=c(15,35),ylim=c(0,0.25),xlab=expression(paste("Temperature (",degree,"C)",sep="")),col="blue",nbins=nbins)
  abline(v=28.6,lty=2,lwd=3)
  #abline(v=26.7,lwd=2)
  title(spname)
}

#create code to separate variables for a single event
f.events<-function(data){
  temp<-mean(data$Temperature)
  site<-unique(as.character(data$Site.Name))
  date<-min(as.character(data$Photo.Date))
  time<-min(as.character(data$Photo.Time))
  sp<-unique(as.character(data$bin))
  sun<-unique(as.character(data$Sampling.Unit.Name))
  lat<-unique(data$Latitude)
  lon<-unique(data$Longitude)
  sap<-unique(data$Sampling.Event)
  mp<-unique(data$Moon.Phase)
c(site,date,time,sap,temp,sp,sun,lat,lon,mp)
  }

#code to put together dataframe with each redcord being an event
f.events.dataframe<-function(data){
  require(chron)
  qwe<-by(data,data$grp,f.events)
  qwe<-as.data.frame(do.call("rbind",qwe))
  names(qwe)<-c("Site.Name","Date","Time","Sampling.Event","Temperature","bin","Sampling.Unit.Name","Latitude","Longitude","Moon.Phase")
  qwe$Date<-as.Date(chron(dates=as.character(qwe$Date),format=c(dates="y-m-d")))
  qwe$Time<-as.POSIXct(as.character(qwe$Time),format="%H:%M:%S")
  qwe$Temperature<-as.numeric(as.character(qwe$Temperature))
  qwe$Latitude<-as.numeric(as.character(qwe$Latitude))
  qwe$Longitude<-as.numeric(as.character(qwe$Longitude))
  qwe$Moon.Phase<-as.numeric(as.character(qwe$Moon.Phase))
  qwe
}  
  
#Code to create temperature event dataframes for a list of species.
#puts them all in a list

f.create.events.splist<-function(splist,fulldata){
  results<-list()
  for(i in 1:length(splist)){
    #Extract the data
    sp<-fulldata[fulldata$bin==splist[i],]
    sp<-f.correct.DF(sp)
    if(dim(sp)[1]<2){
      print(paste("Species ",splist[i]," has insufficient data..",sep=""))
      results<-c(results,list(NULL))
      next
    }
    #Order the data in chronological order
    sp<-f.order.data(sp)
    #Create independent observation events list with a threshold of 5 min
    sp<-f.separate.events(sp,5) 
   #Create a simplified data frame with just the events
    sp<-f.events.dataframe(sp)
    results<-c(results,list(sp))
    print(paste("Species ",splist[i]," processed..",sep=""))
  }
  names(results)<-splist
  results
  
}
f.print.graphs<-function(data){
  path="/Users/jorge/Analyses/TempTV/graphs2/"
for(i in 1:length(data)) {
  newp<-paste(path,names(data)[i],".pdf",sep="")
  pdf(newp)
  qplot(Temperature,data=data[[i]],geom="histogram",binwidth=1,main=names(data)[i])
  ggsave(newp)
}  
  
}
# funcion para asignar camaras faltantes que no tomaron fotos de animales.
# SuName,startDate y endDate deben estar entre comillas. Los demas argumentos no.
# startDate y endDate estan en format yyyy-mm-dd
f.assign.missing<-function(SuName,SuPeriod,startDate,endDate,data){
  rows<-dim(data)[1]
  
  #agregar Sampling Unit Name
  data[rows+1,3]<-SuName
  #agregar StartDate
  data[rows+1,38]<-as.Date(startDate)
  #agregar End Date
  data[rows+1,39]<-as.Date(endDate)
  #agregar sampling unit period
  data[rows+1,6]<-SuPeriod
  data
  
}

f.minusBirds<-function(data){
  indx<-which(data$Class=="AVES")
  data<-data[-indx,]
  data<-f.correct.DF(data)
}

#code to shrink the matrix by half
f.shrink.matrix.half<-function(matrix){
  #if number of columns in the matrix is even
  if(!ncol(matrix)%%2){
    #figure out how many columns
    nc<-ncol(matrix)/2  
    #disagregate into individual matrices
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=nc)
    old.cols<-seq(1,ncol(matrix),2)
    for(i in 1:nc){
      #sum the rows for the column sections
      sum.rows<-apply(matrix[,old.cols[i]:(old.cols[i]+1)],1,sum,na.rm=T)
      #convert to 0s and 1s
      new.matrix[,i]<-ifelse(sum.rows>=1,1,0)
    }  
    new.matrix
  }
  #if the number of columns is not even
  else{
    #store the first column in col1  
    col1<-matrix[,1]
    #convert the matrix to an even matrix
    matrix<-matrix[,-1]
    nc<-ncol(matrix)/2  
    #disagregate into individual matrices
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=nc)
    old.cols<-seq(1,ncol(matrix),2)
    for(i in 1:nc){
      sum.rows<-apply(matrix[,old.cols[i]:(old.cols[i]+1)],1,sum,na.rm=T)
      new.matrix[,i]<-ifelse(sum.rows>=1,1,0)
    }
    cbind(col1,new.matrix)  
  }
}

#collapseData<-f.shrink.matrix(testData)
#umf2<-unmarkedFrameOccu(y=collapsedData)
#umf2
#fmcoll<-occu(~1 ~1,umf2)
#summary(fmcoll)
#plogis(coef(mod,"det"))
#plogis(coef(mod,"state"))
#plogis(coef(mod2,"det"))
#plogis(coef(mod2,"state"))

#code to shrink the matrix to about 15 columns
f.shrink.matrix<-function(matrix){
  
  #disagregate into individual matrices
  nc<-length(seq(1,ncol(matrix),9))  
  new.matrix<-matrix(NA,nr=nrow(matrix),nc=nc)
  rownames(new.matrix)<-rownames(matrix)  
  old.cols<-seq(1,ncol(matrix),9)
  for(i in 1:nc){
    #sum the rows for the column sections
    sum.rows<-apply(matrix[,old.cols[i]:(old.cols[i]+1)],1,sum)
    #convert to 0s and 1s
    new.matrix[,i]<-ifelse(sum.rows>=1,1,0)
  }  
  new.matrix
}

f.plot.jag.res<-function(jags,species.name,model.name){
  #naive occupancy
  naive.occ<-apply(tmp,2,sum,na.rm=T)/apply(tmp,2,function(x){length(which(!is.na(x)))})
  #plot results
  mat<-jags$BUGSoutput$summary
  #mean occupancy from model
  m.occ<-jags$BUGSoutput$mean$psi
  #median occupancy from model
  med.occ<-jags$BUGSoutput$median$psi
  #extract 95% confidence limits for occupancy
  psiCol<-which(rownames(mat)=="psi[1]")
  lo95ci<-mat[psiCol:(psiCol+3),3]
  hi95ci<-mat[psiCol:(psiCol+3),7]
  par(mfrow=c(2,1))
  plot(2009:2012,m.occ,ylim=c(0,1),type='b',xlab="year",ylab="occupancy")
  lines(2009:2012,as.numeric(lo95ci),lty=2)
  lines(2009:2012,as.numeric(hi95ci),lty=2)
  lines(2009:2012,naive.occ,col='red')
  #lines(2007:2011,med.occ,lwd=3)
  title(paste(species.name,"\n",model.name))
  
  #graph lambda
  lambdaCol<-which(rownames(mat)=="lambda[1]")
  m.lambda<-mat[lambdaCol:(lambdaCol+2),1]
  
  lo95ci<-mat[lambdaCol:(lambdaCol+2),3]
  hi95ci<-mat[lambdaCol:(lambdaCol+2),7]
  plot(1:3,m.lambda,ylim=range(lo95ci,hi95ci),type='b',xlab="year interval",ylab="lambda")
  lines(1:3,as.numeric(lo95ci),lty=2)
  lines(1:3,as.numeric(hi95ci),lty=2)
  abline(h=1,lty=3)
  title(paste(species.name,"\n",model.name))
}


