#####################################
# ANÁLISIS OCUPACIÓN A NIVEL DE ESPECIE
# basado en el Capítulo 6 de Rovero y Spitale (2016)
#####################################

# PASO 1: cargar paquetería
library(chron) 
library(reshape)
library(ggplot2)
library(vegan)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(plyr)
library(R2jags)

# PASO 2: cargar recursos adicionales y datos generados en Wild.ID 

source("TEAM library 1.7.R")

team_data <- read.csv(file = "teamexample.csv", sep = ",", h = T, stringsAsFactors = F)
View(team_data)

# PASO 3: adicionar Clase, Orden, Familia empleando la base IUCN

iucn.full <- read.csv("IUCN.csv", sep = ",", h = T)

iucn<-iucn.full[,c("Class","Order","Family","Genus","Species")]

team <- merge(iucn, team_data, all.y = T)
View(team)

# PASO 4: se emplea función para darle el formato adecuado

fd <- fix.dta(team)
View(fd)

yr2009 <- fd[fd$Sampling.Event =="2009.01" & fd$Class=="MAMMALIA",]
View(yr2009)

# PASO 5: cargar datos de covariables

cov <- read.table("covariates.txt", header = TRUE)
View(cov)

workingcam <- which(cov$Sampling.Unit.Name %in% unique(yr2009$Sampling.Unit.Name)) # se eliminan cámaras que no funcionaron

cov.or <- cov[workingcam, ] # solo cámaras 2009

cov.num <- cov.or[,sapply(cov.or, is.numeric)]

cov.std <- decostand(cov.num, method = "standardize")

cov.fac <- cov.or[,sapply(cov.or, is.factor)]  # extrea factores

covs <- data.frame(cov.fac, cov.std)
View(covs)

# PASO 6: crear matrices para cada especie

mat.udz.09 <- f.matrix.creator(yr2009)
names(mat.udz.09) 

naivetable <- naive(mat.udz.09) 
View(naivetable)

# -----------------------------------
# Análisis con el paquete "unmarked"
# ejemplo de análisis con una especie: Sanje mangabey (Cercocebus sanjei) 

Cs <- shrink(mat.udz.09[["Cercocebus sanjei"]],5)
View(Cs)

(umCs <- unmarkedFrameOccu(y = Cs, siteCovs = covs))

# modelos analizados
m0 <- occu(~1 ~1, umCs)
d1 <- occu(~edge ~1, umCs)
d2 <- occu(~border ~1, umCs)
d3 <- occu(~edge + border ~1, umCs)
o1 <- occu(~1 ~border, umCs)
o2 <- occu(~1 ~habitat, umCs)
o3 <- occu(~1 ~habitat + border, umCs)
m1 <- occu(~edge ~border, umCs)         
m2 <- occu(~border ~border, umCs)
m3 <- occu(~edge + border ~border, umCs)  
m4 <- occu(~edge ~habitat, umCs)
m5 <- occu(~border ~habitat, umCs)
m6 <- occu(~edge + border ~habitat, umCs)
m7 <- occu(~edge + border ~habitat + border, umCs)

# a manera de ejemplo se examina el modelo "m1"

m1

backTransform(linearComb(m1, coefficients = c(1,0), type = "det")) 

backTransform(linearComb(m1, coefficients = c(1,0), type = "state"))

# análisis de selección de modelo

dlist <- fitList(Nullo = m0, d1 = d1, d2 = d2, d3 = d3, o1 = o1, o2 = o2, o3 = o3, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5, m6 = m6, m7 = m7)

(selmod <- modSel(dlist, nullmod = "Nullo"))

# predecimos la ocupación de hábitat del modelo "m6"

newhab <- data.frame(habitat = c("Deciduous", "Montane"))

(pred <- predict(m6, type = "state", newdata = newhab, appendData = T))

ggplot(pred,aes(x = habitat, y = Predicted)) +
  geom_point(size = 4) +
  ylab("Predicted Psi Cercocebus sanjei") +
  theme_bw() +
  geom_errorbar(aes(ymin = Predicted-SE, ymax = Predicted + SE), width = 0.2)

# se estandarizan las covariables para crear un raster empleando la media y sd de las covariables medidas en cada cámara

map <- read.table("covs100x100.txt", h = T)  
View(map)

mapst <- data.frame(x = map$x, y = map$y, habitat = map$habitat, edg = (map$edge - mean(cov.or$edge)) / sd(cov.or$edge), border = (map$border - mean(cov.or$border)) / sd(cov.or$border), river = (map$river - mean(cov.or$river)) / sd(cov.or$river))
View(mapst)

# se genera el mapa de probabilidad de ocupación

predmap <- predict(o2, type = "state", newdata = mapst, appendData = T) 

levelplot(predmap$Predicted ~ x + y, map, aspect = "iso", xlab = "Easting (m)", ylab = "Northing (m)", main = "Cercocebus sanjei", col.regions = topo.colors(90))

# ----------------------------------
# FIN SCRIPT

rm(list = ls())
dev.off()

