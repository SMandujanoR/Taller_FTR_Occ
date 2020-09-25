#####################################################
# 
#   Modelos Royle-Nichols 
#   Karen Lorena Velásquez C
#   Instituto de Ecología A. C.
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

summary(fmUrRN1 <- occuRN(~1~1, start=c(1,1), data=UmUr ))

summary(fmUrRN2 <- occuRN(~CAM~1, data=UmUr))
summary(fmUrRN3 <- occuRN(~CVERT~1, data=UmUr))
summary(fmUrRN4 <- occuRN(~DenVeg~1, data=UmUr))

# Lista modelos

fms <- fitList ("p(.)psi(.)"            = fmUrRN1,
                "p(CAM)psi(.)"          = fmUrRN2,
                "p(CVERT)psi(.)"        = fmUrRN3,
                "p(DenVeg)psi(.)"       = fmUrRN4)
(ms <- modSel(fms))

## Proceso ecológico

summary(fmUrRN5  <- occuRN(~DenVeg~ DBEB, data=UmUr))
summary(fmUrRN6  <- occuRN(~DenVeg~ DCUL, data=UmUr))
summary(fmUrRN7 <- occuRN(~DenVeg~ DPOB, data=UmUr))
summary(fmUrRN8 <- occuRN(~DenVeg~ CVERT, data=UmUr))
summary(fmUrRN9  <- occuRN(~DenVeg~ Slope, data=UmUr)) 
summary(fmUrRN10  <- occuRN(~DenVeg~ MSAVI, data=UmUr))

# Lista modelos

fms <- fitList ("p(DenVeg)psi(DBEB)"                = fmUrRN5,
                "p(DenVeg)psi(DCUL)"                = fmUrRN6,
                "p(DenVeg)psi(DPOB)"                = fmUrRN7,
                "p(DenVeg)psi(CVERT)"               = fmUrRN8,
                "p(DenVeg)psi(Slope)"               = fmUrRN9,
                "p(DenVeg)psi(MSAVI)"               = fmUrRN10
)

(ms <- modSel(fms))

# Obtención valor probabilidad de detección/ocupación

backTransform(linearComb(fmUrRN6, coefficients=c(1,0), type="state"))


backTransform(linearComb(fmUrRN6, coefficients=c(1,0), type="det"))

## Importancia de las covariables

confint(fmUrRN6, type="state")
confint(fmUrRN6, type="det")

## Prueba de bondad de ajuste

# Cargar archivo Rdata

#load("mb_UrciRN6Clase.RData")

#save(mb_UrciRN6clase, file="mb_UrciRNClase.Rdata")

mb_UrciRN6clase <- mb.gof.test(fmUrRN6, nsim = 500, plot.hist = TRUE) 

# Resultado

mb_UrciRN6clase

# Gráfico GOF-test

par(mar = c(6.50,
            5.00,
            6.50,
            3.00))
hist(mb_UrciRN6clase$t.star, xlab=expression(paste("Xi"^"2")), ylab="Frecuencia", col="lightgrey",
     font.lab=2, cex.lab=0.9, main="Prueba bondad de ajuste modelo Ocupación") 
abline(v=mb_UrciRN6clase$chi.square, lty=2, lwd=3, col="red")


## Gráficos de prediccción de las covariables

# Detección

newData <- data.frame(DenVeg=seq(min(covs$DenVeg),max(covs$DenVeg), length=100))
statepredict <-predict(fmUrRN6, type="det", newdata=newData, appendData=TRUE)
DenVegplot <- ggplot(statepredict,aes(x=DenVeg,y=Predicted))+ylim(0,1)+
  labs(x="Densidad de Vegetación",y="p")+
  geom_ribbon(data = statepredict,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#cfe0e3")+
  geom_line(data = statepredict, colour="#2aa3bb", size=1.2)+
  theme_classic()
DenVegplot

# Ocupación

DCUL<-data.frame(DCUL=seq(min(covs$DCUL),max(covs$DCUL),length=100)) 
PredDCUL <-predict(fmUrRN6,type="state",newdata=DCUL,appendData=TRUE) 
DCULplot <- ggplot(PredDCUL,aes(x=DCUL,y=Predicted))+ylim(0,15)+
  labs(x="Distancia a cultivos",y=expression(lambda))+geom_ribbon(data = PredDCUL,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#cfe0e3")+
  geom_line(data = PredDCUL, colour="#2aa3bb", size=1.2)+
  theme_classic()
DCULplot



