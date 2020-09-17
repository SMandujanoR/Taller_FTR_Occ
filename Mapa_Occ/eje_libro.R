########Modelos de ocupaci?n para una especie una temporada con enfoque de M?xima verosimilitud#####
##Modificado script de Gabriel Andrade
#### Paqueter?a requerida ####
library(vegan)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(ggplot2)
library(gridExtra)
library(boot)
library(corrplot)
library(usdm)
library(ggcorrplot)
library(raster)

setwd("C:/Users/HP Notebook/Desktop/e_libro")

mim<-read.csv("historial.csv",header = T)

y <- as.matrix(mim)

### covariables de sitio
site_cov <- read.csv("cov.csv", header = T)
str(site_cov)
cova<-site_cov[,4:7]

#covariables de observaci?n
obscov<-read.csv("dia.csv",header = T)
date.orig<-as.matrix(obscov)


#An?lisis de correlaci?n
CorelacionSite<-round(cor(cova[,c(1,2,3,4)]),1)
CorelacionSite
ggcorrplot(CorelacionSite, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)

#An?lisis de colinearidad
cov.num2 <-cova[,c(1,2,3,4,5)]
cormat <- cor(cov.num2)
(corrplot(cormat, outline = T, tl.col = "black", mar = c(2,0,1,1.5), title = "Colinearity matrix"))

v<-vif(cov.num2)
v

#estandarizar covariables
msavi.orig<-cova$msavi
aspect.orig<-cova$orientacio
#elevacion.orig<-cova$elevacion
pend.orig<-cova$pendiente

msavi<-decostand(msavi.orig,method="standardize")
orientacion<-decostand(aspect.orig,method="standardize")
#elevacion<-decostand(elevacion.orig, method = "standardize")
pendiente<-decostand(pend.orig, method = "standardize")


#unmarked format
umf<-unmarkedFrameOccu(y=y,siteCovs = data.frame(msavi, orientacion,pendiente),obsCovs = list(date=date.orig))

plot(umf)
summary(umf)

#2.1.1 Modelando la detecci?n de Mimus graysoni
#dejamos la ocupaci?n con el modelo global

summary(fm1 <- occu(~date~msavi+pendiente+orientacion, data=umf))
summary(fm2 <- occu(~1~msavi+pendiente+orientacion, data=umf,control=list(maxit=1000, trace=TRUE, REPORT=1)))

#,control=list(maxit=1000, trace=TRUE, REPORT=1)

fmsoc <- fitList("p(date)psi(msavi+pendiente+orientacion)" =fm1,
                 "p(.)psi(msavi+pendiente+orientacion)"=fm2)

(msoc <- modSel(fmsoc))

#AICc para muestras peque?as
aictab(list(fm1, fm2), modnames = c("p(date)psi(NDVI+ori+pend)", "p(.)psi(NDVI+ori+pen)"), second.ord = TRUE, nobs = 54, sort = TRUE, c.hat = 1)

#ocupaci?n
summary(fm3<-occu(~1~msavi+pendiente+orientacion, data= umf))
summary(fm4<-occu(~1~msavi+orientacion, data= umf))
summary(fm5<-occu(~1~msavi+pendiente, data= umf))
summary(fm6<-occu(~1~orientacion+pendiente, data= umf))
summary(fm7<-occu(~1~msavi, data= umf))
summary(fm8<-occu(~1~pendiente, data= umf))
summary(fm9<-occu(~1~orientacion, data= umf))
summary(fm10<-occu(~1~1, data= umf))


fmsoc1 <- fitList("p(.)psi(msavi+pendiente+orientacion)" =fm3,
                  "p(.)psi(msavi+orientacion)" =fm4,
                  "p(.)psi(msavi+pendiente)" =fm5,
                  "p(.)psi(orientacion+pendiente)" =fm6,
                  "p(.)psi(msavi)"=fm7,
                  "p(.)psi(pendiente)" =fm8,
                  "p(.)psi(orientacion)" =fm9,
                  "p(.)psi(.)" =fm10)


(msoc1 <- modSel(fmsoc1))

#AICc 
aictab(list(fm3,fm4, fm5,fm6, fm7,fm8, fm9, fm10), modnames = c("p(.)psi(msavi+orientacion+pendiente)","p(.)psi(msavi+ori)", "p(.)psi(msavi+pendiente)","p(.)psi(orientacion+pendiente)", "p(.)psi(msavi)", "p(.)psi(pend)", "p(.)psi(ori)","p(.)psi(.)"),  second.ord = TRUE, nobs = 54, sort = TRUE, c.hat = 1)

# Guardar resultados

#save(fm7, file = "modelo_fm7.RData")
#save(msavi, file = "msavi_sd.RData")


#prueba de bondad de ajuste modelo

#p(.)psi(elevacion)
#windows()
#gof.boot<-mb.gof.test(fm7,nsim = 1000) 
#c-hat= 0.4

###intervalo de confianza variables

#confint(fm7, type="det")
#confint(fm7, type="state")

#####################################################
#mejor modelo fm7 p(.)psi(msavi)

###obtener la detecci?n
#backTransform(fm7, type="det")

#Newdat2 <- data.frame(msavi=seq(min(msavi),max(msavi)))

#pre.msavi.fm7 <-predict(fm7, newdata=Newdat2,type="state", appendData=TRUE)

#Graficar predicciones ocupaci?n

#gr?fico ocu msavi fm7

#(msaviplotocu<- ggplot(pre.msavi.fm7,aes((x=msavi),y=Predicted))+ylim(0,1)+ 
#    labs(x="Elevacion",y="Probabilidad de ocupaci?n")+
#    geom_line(data=pre.msavi.fm7)+
#    geom_ribbon(data = pre.msavi.fm7,aes(ymin=lower,ymax=upper),alpha=0.3)+
#    theme_classic())


################MAPAS########################################
#M?todo uno

#msavi

#scaled center:0.4591
#scaled scale: 0.20.03

#(betas<-coef(fm7, type="state"))

#psi(Int) psi(msavi) 
#-2.351697   5.503494 

#logit.psi <- -2.351697 + 5.503494*MSAVI_e
#psi <- exp(logit.psi) / (1 + exp(logit.psi))
#print(spplot(psi, col.regions=terrain.colors(50)))

#####################################################
####Mapa predicci?n de ocupaci?n, SE e IC95%
#para que tarde menos divid? el area en 8 pedazos y estandariz? el raster

###PEDAZO 1

#setwd("C:/Users/HP Notebook/Desktop/e_libro/1")
#MSAVI = raster("./msavi.tif")
#MSAVI_e<- (MSAVI-0.45)/.20

##stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/1msavi", x,".tif"),overwrite=TRUE)})

###PEDAZO 2

#setwd("C:/Users/HP Notebook/Desktop/e_libro/2")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/2msavi", x,".tif"),overwrite=TRUE)})

###PEDAZO 3

#setwd("C:/Users/HP Notebook/Desktop/e_libro/3")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/3msavi", x,".tif"),overwrite=TRUE)})

####PEDAZO 4

#setwd("C:/Users/HP Notebook/Desktop/e_libro/4")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/4msavi", x,".tif"),overwrite=TRUE)})

###PEDAZO 5 Tehuac?n

#setwd("C:/Users/HP Notebook/Desktop/e_libro/5")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/5msavi", x,".tif"),overwrite=TRUE)})

#pedazo 6 Tehuac?n

#setwd("C:/Users/HP Notebook/Desktop/e_libro/6")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/6msavi", x,".tif"),overwrite=TRUE)})

###PEDAZO 7 Tehuac?n

#setwd("C:/Users/HP Notebook/Desktop/e_libro/7")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/7msavi", x,".tif"),overwrite=TRUE)})

###PEDAZO 8 Tehuac?n

#setwd("C:/Users/HP Notebook/Desktop/e_libro/8")

#MSAVI2 = raster("./msavi.tif")
#MSAVI_e<- (MSAVI2-0.45)/.20

#stackCapas = stack(MSAVI_e)

#E.psi <- predict(fm7, type="state", newdata=stackCapas)

#plot(E.psi, axes=FALSE, col=terrain.colors(100))

#lapply(names(E.psi), function(x){writeRaster(E.psi[[x]], paste0("C:/Users/HP Notebook/Desktop/8msavi", x,".tif"),overwrite=TRUE)})


###MAPA FINAL modelo cap?tulo: aqu? ya un? los 8 pedazos en un solo raster para la predicci?n, uno para SE, y otros paro los intervalos de confianza

prediccion = raster("./mod_predic.tif")
SE = raster("./mod_SE.tif")
IC_superior = raster("./mod_upper.tif")
IC_inferior = raster("./mod_lower.tif")

windows()
par(mfrow=c(2,2))
plot(prediccion, axes=F, col=terrain.colors(100), main="Predicci?n")
plot(SE, axes=F, col=terrain.colors(100), main="EE")
plot(IC_inferior, axes=F, col=terrain.colors(100), main="Inferior")
plot(IC_superior, axes=F, col=terrain.colors(100), main="Superior")

