#####################################################
# 
#   Modelos de ocupación una especie/una temporada
#
#####################################################

## Librerias

library(vegan)
library(corrplot)
library(usdm)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(ggplot2)


## Covariables

covar <- read.csv("CovariablesOcc.csv", sep=",")
cov.site <- covar[,-1]

# Extracción covariables númericas/factor

cov.num<-cov.site[,sapply(cov.site,is.numeric)]
cov.std<-decostand(cov.num,method="standardize")
CAM<-cov.site[,sapply(cov.site,is.character)]

# Correlación

cormat <- cor(cov.std)
corrplot(cormat, method = c("circle"), type = "upper", outline = T, tl.col = "black", tl.cex = 0.8, tl.srt = 45,  mar = c(2,0,1,1.5), title = "Matriz de correlación")

# Factor inflación varianza
 
usdm::vif(cov.std)
no_corr <- vifstep(cov.std, th=4) 
no_corr

# Unión covariables

covs<-data.frame(CAM, cov.std)

# Generación marco de datos 

Urci <- read.csv("DetHist_UrciClase.csv", sep=",")
Urcihist <- Urci[,-1]
UmUr <- unmarkedFrameOccu(y=Urcihist, 
                          siteCovs=covs, 
                          obsCovs= NULL)
summary(UmUr)
str(UmUr)

## Proceso observacional

summary(fmUr1 <- occu(~1~1, start=c(1,1), data=UmUr ))

summary(fmUr2 <- occu(~CAM~1, data=UmUr))
summary(fmUr3 <- occu(~CVERT~1, data=UmUr))
summary(fmUr4 <- occu(~DenVeg~1, data=UmUr))

# Lista modelos

fms <- fitList ("p(.)psi(.)"            = fmUr1,
                "p(CAM)psi(.)"          = fmUr2,
                "p(CVERT)psi(.)"        = fmUr3,
                "p(DenVeg)psi(.)"        = fmUr4)
(ms <- modSel(fms))

## Proceso ecológico

summary(fmUr5  <- occu(~DenVeg~ DBEB, data=UmUr))
summary(fmUr6  <- occu(~DenVeg~ DCUL, data=UmUr))
summary(fmUr7 <- occu(~DenVeg~ DPOB, data=UmUr))
summary(fmUr8 <- occu(~DenVeg~ CVERT, data=UmUr))
summary(fmUr9  <- occu(~DenVeg~ Slope, data=UmUr)) 
summary(fmUr10  <- occu(~DenVeg~ MSAVI, data=UmUr))

# Lista modelos

fms <- fitList ("p(DenVeg)psi(DBEB)"                = fmUr5,
                "p(DenVeg)psi(DCUL)"                = fmUr6,
                "p(DenVeg)psi(DPOB)"                = fmUr7,
                "p(DenVeg)psi(CVERT)"               = fmUr8,
                "p(DenVeg)psi(Slope)"               = fmUr9,
                "p(DenVeg)psi(MSAVI)"               = fmUr10
)

(ms <- modSel(fms))

# Obtención valor probabilidad de detección/ocupación

backTransform(linearComb(fmUr6, coefficients=c(1,0), type="state"))


backTransform(linearComb(fmUr6, coefficients=c(1,0), type="det"))

## Importancia de las covariables

confint(fmUr6, type="state")
confint(fmUr6, type="det")

## Prueba de bondad de ajuste

# Cargar archivo Rdata

load("mb_UrciClaseOccu.RData")

#save(mb_Urci6clase, file="mb_UrciClaseOccu.Rdata")

#mb_Urci6clase <- mb.gof.test(fmUr6, nsim = 500, plot.hist = TRUE) 

# Resultado

mb_Urci6clase

# Gráfico GOF-test

par(mar = c(6.50,
            5.00,
            6.50,
            3.00))
hist(mb_Urci6clase$t.star, xlab=expression(paste("Xi"^"2")), ylab="Frecuencia", col="lightgrey",
     font.lab=2, cex.lab=0.9, main="Prueba bondad de ajuste modelo Ocupación") 
abline(v=mb_Urci6clase$chi.square, lty=2, lwd=3, col="red")


## Gráficos de prediccción de las covariables

# Detección

newData <- data.frame(DenVeg=seq(min(covs$DenVeg),max(covs$DenVeg), length=100))
statepredict <-predict(fmUr6, type="det", newdata=newData, appendData=TRUE)
DenVegplot <- ggplot(statepredict,aes(x=DenVeg,y=Predicted))+ylim(0,1)+
  labs(x="Densidad de Vegetación",y="p")+
  geom_ribbon(data = statepredict,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#cfe0e3")+
  geom_line(data = statepredict, colour="#2aa3bb", size=1.2)+
  theme_classic()
DenVegplot

# Ocupación

DCUL<-data.frame(DCUL=seq(min(covs$DCUL),max(covs$DCUL),length=100)) 
PredDCUL <-predict(fmUr6,type="state",newdata=DCUL,appendData=TRUE) 
DCULplot <- ggplot(PredDCUL,aes(x=DCUL,y=Predicted))+ylim(0,1)+
  labs(x="Distancia a cultivos",y=expression(psi))+geom_ribbon(data = PredDCUL,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#cfe0e3")+
  geom_line(data = PredDCUL, colour="#2aa3bb", size=1.2)+
  theme_classic()
DCULplot



