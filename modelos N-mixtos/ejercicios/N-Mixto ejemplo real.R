###################################################
# CAPITULO 6: MODELOS N-MIXTOS
# Ejemplo 6.9: Estimación de la abundancia. Datos reales de una especie de ave: "Great tit". Pp.254-278.
###################################################

# Leer datos:
tits <- read.table("SwissTits_mhb_2004_2013.csv", header = T, sep = ";")
#show(tits)
head(tits)
#str(tits)
names(tits)

# Para ver las especies
table(tits$latname) # nombre latino
table(tits$name) # nombre en inglés

# -------------------------------------------------
# Para este ejemplo se trabajará con el "Great tit" y las covariables del año 2013. Además, se eliminarán los sitios no muestreados durante ese año.

tits <- tits[tits$latname == "PARMAJ",] # datos del great tit
show(tits)
NA.sites <- which(apply(is.na(tits[,41:43]), 1, sum) == 3) # sitios no muestreados
(tits <- tits[-NA.sites,]) # eliminación de esos sitios

y <- cbind(tits$count131,tits$count132, tits$count133) # conteos
date <- cbind(tits$date131,tits$date132, tits$date133) # día del muestreo
dur <- cbind(tits$dur131,tits$dur132, tits$dur133) # duración del muestreo

# gráficamente
matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey date (1 = April 1)", ylab = "Count of Great Tits")

# ---------------------------------
# Ajuste de modelos
library(unmarked)

# Primer modelo:
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y, siteCovs=data.frame(elev = scale(tits[,"elev"]), forest = scale(tits[,"forest"]), iLength = 1/tits[,"rlength"]), obsCovs = list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)  # Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) # Summarize max counts

###############################
# Se generan modelos

# Modelo 1
###fm1 <- pcount(~ (elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength, umf, control=list(trace=TRUE, REPORT=5))
###summary(fm1)
###save(data, file = "fm1.RData")
load("fm1.RData")

# Modelo 2
###fm2 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength, starts = coef(fm1)[1:31], umf, control=list(trace=TRUE, REPORT=5))
###summary(fm2)
###save(data, file = "fm2.RData")
load("fm2.RData")

# Modelo 3
###fm3 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2) ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength, starts = coef(fm2)[-c(23, 27, 31)], umf, control=list(trace=TRUE, REPORT=5))
###summary(fm3)
###save(data, file = "fm3.RData")
load("fm3.RData")

# Modelo 4
fm4 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2) - elev:I(date^2) - I(date^2):dur ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength, starts = coef(fm3)[-c(21, 28)], umf, control=list(trace=TRUE, REPORT=5))
summary(fm4)
###save(fm4, file = "fm4.RData")

# Modelo 5 Binomial-Poisson
fm5 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2) - elev:I(date^2) - I(date^2):dur ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength - I(elev^2):forest - I(elev^2):I(forest^2), starts = coef(fm4)[-c(9:10)], umf, control=list(trace=TRUE, REPORT=5))
summary(fm5)
###save(fm5, file = "fm5.RData")

# Negative binomial (NB) mixture
fm5NB <- pcount(fm5@formula, starts = c(coef(fm5),0), umf, control=list(trace=TRUE, REPORT=5), mixture = "NB")
summary(fm5NB)
###save(fm5NB, file = "fm5NB.RData")

# Zero-inflated Poisson (ZIP) mixture
fm5ZIP <- pcount(fm5@formula, starts = c(coef(fm5),0), umf, control=list(trace=TRUE, REPORT=5), mixture = "ZIP")
summary(fm5ZIP)
###save(fm5ZIP, file = "fm5ZIP.RData")

# -----
cbind(fm5@AIC, fm5NB@AIC, fm5ZIP@AIC)

round(cbind(rbind("Poisson" = exp(coef(fm5)[1]), "NegBin" = exp(coef(fm5NB)[1]), "ZIP" = exp(coef(fm5ZIP)[1])), rbind(plogis(coef(fm5)[15:17]), plogis(coef(fm5NB)[15:17]), plogis(coef(fm5ZIP)[15:17]))),2)

# ----
# Análisis GoF
library(AICcmodavg)

# ojo para fines de rapidez, solo se similan 10 iteraciones
system.time(gof.P <- Nmix.gof.test(fm5, nsim=10)) 
gof.P
###save(gof.P, file = "gof.P.RData")

system.time(gof.NB <- Nmix.gof.test(fm5NB, nsim=10)) 
gof.NB
###save(gof.NB, file = "gof.NB.RData")

system.time(gof.ZIP <- Nmix.gof.test(fm5ZIP, nsim=10))
gof.ZIP
###save(gof.ZIP, file = "gof.ZIP.RData")

# gráficamente
print(cbind(y, fitted(fm5), residuals(fm5)), 2) 
print(cbind(y, fitted(fm5NB), residuals(fm5NB)), 2) 
print(cbind(y, fitted(fm5ZIP), residuals(fm5ZIP)), 2) 

library(AHMbook)
plot_Nmix_resi(fm5, fm5NB, fm5ZIP)

# Se calcula el RMSE para los 3 modelos
(RMSEP <- sqrt(mean((y - fitted(fm5))^2, na.rm = TRUE))) 
(RMSENB <- sqrt(mean((y - fitted(fm5NB))^2, na.rm = TRUE))) 
(RMSEZIP <- sqrt(mean((y - fitted(fm5ZIP))^2, na.rm = TRUE)))

# Para analizar la estructura espacial
map.Nmix.resi(fm5)
map.Nmix.resi(fm5NB)
map.Nmix.resi(fm5ZIP)

# --------------------------
# Análisis de resultados

lamNewData <- data.frame(elev = (seq(200, 2250,,200) - mean(tits[,"elev"]))/ sd(tits[,"elev"]), forest = 0, iLength = 1/5.1)

pNewData <- data.frame(elev = 0, time = factor("2", levels = c("1", "2", "3")), dur = 0, date = (seq(1,90,,200) - mean(date, na.rm = T))/sd(date, na.rm = T))

predict(fm5, type="state", newdata=lamNewData) 
predict(fm5, type="det", newdata=pNewData) 
predict(fm5NB, type="state", newdata=lamNewData) 
predict(fm5NB, type="det", newdata=pNewData) 
predict(fm5ZIP, type="state", newdata=lamNewData) 
predict(fm5ZIP, type="det", newdata=pNewData)

predictSE(fm5ZIP, newdata=lamNewData, print.matrix = TRUE, type="response", parm.type = "lambda", c.hat = 2.47)

# Predicciones de lambda y p
modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 3.82)

modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 3.82) 

modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 3.82)

modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 3.82) 

# modelo NegBin model
modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 1.79)

modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 1.79) 

modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 1.79)

modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 1.79) 

# modelo ZIP 
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 2.47) 
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 2.47)

# 
rlength <- seq(1, 30, 0.01) 
newData <- data.frame(elev=0, forest=0, iLength=1/rlength)
pred <- predictSE(fm5ZIP, parm.type="lambda", newdata=newData, c.hat = 2.47)

par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.3)
plot(rlength, pred[[1]], type = "l", lwd = 3, col = "blue", frame = F, xlab = "Transect length (km)", ylab = "Exposed population (lambda)", ylim = c(0, 16), axes = F) 
axis(1, at = seq(2,30,2)) ; axis(2)
abline(v = c(1.2, 5.135, 9.4), lwd = 1, lty =2)
matlines(rlength, cbind(pred[[1]]-pred[[2]], pred[[1]]+pred[[2]]), type = "l", lty = 1, lwd = 2, col = "gray")
sat.pred <- predictSE(fm5ZIP, parm.type="lambda", newdata= data.frame(elev=0, forest=0, iLength=0), c.hat = 2.47)
abline(h = sat.pred$fit, lwd = 2, lty = 2, col = "red")

# 
print(cbind("Route length" = rlength, "Exp. pop" = pred[[1]], "Rel. exp. pop" = pred[[1]] / sat.pred$fit), 3)

# Análisis de predicción y estandarización de covariables
ep.orig <- seq(200, 2250, length.out=100) 
(elev.mean <- mean(tits[,"elev"])) 
(elev.sd <- sd(tits[,"elev"]))
ep <- (ep.orig - elev.mean) / elev.sd 
fp.orig <- seq(0, 100, length.out=100) 
(forest.mean <- mean(tits[,"forest"])) 
(forest.sd <- sd(tits[,"forest"]))
fp <- (fp.orig - forest.mean) / forest.sd 
date.orig <- seq(1, 90, length.out=100) 
(date.mean <- mean(date, na.rm = TRUE)) 
(date.sd <- sd(date, na.rm = TRUE))
datep <- (date.orig - date.mean) / date.sd 
dur.orig <- seq(90, 420, length.out=100) 
(dur.mean <- mean(dur, na.rm = TRUE)) 
(dur.sd <- sd(dur, na.rm = TRUE))
durp <- (dur.orig - dur.mean) / dur.sd

# Predicciones en el gradiente de una covariable
newData1 <- data.frame(elev=ep, forest=0, iLength=0, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))

pred1 <- predictSE(fm5ZIP, newdata=newData1, c.hat = 2.47)
pred1[[1]]

pred2 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData1, parm.type = "detect", type = "response", c.hat = 2.47)
pred2[[1]]

newData3 <- data.frame(elev=0, forest=fp, iLength=0, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))
pred3 <- predictSE(fm5ZIP, newdata=newData3, c.hat = 2.47)
pred3[[1]]

newData4 <- data.frame(elev=0, forest=0, iLength=0, date=datep, dur=0, time = factor("2", levels = c("1", "2", "3")))
pred4 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData4, parm.type = "detect", type = "response", c.hat = 2.47)
pred4[[1]]

newData5 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=durp, time = factor("2", levels = c("1", "2", "3")))
pred5 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData5, parm.type = "detect", type = "response", c.hat = 2.47)

newData6 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=0, time = c("1", "2", "3"))
pred6 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData6, parm.type = "detect", type = "response", c.hat = 2.47)

# Gráfico des la predicciones
par(mfrow = c(3,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)

plot(ep.orig, pred1[[1]], type = "l", lwd = 2, col = "blue", xlab = "Elevation (m)", ylab = "Expected abundance", las = 1, ylim = c(0,50), frame = F)
matlines(ep.orig, cbind(pred1[[1]]-pred1[[2]], pred1[[1]]+pred1[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")

plot(fp.orig, pred3[[1]], type = "l", lwd = 2, col = "blue", xlab = "Forest cover (%)", ylab = "Expected abundance", las = 1, ylim = c(0, 18), frame = F)
matlines(fp.orig, cbind(pred3[[1]]-pred3[[2]], pred3[[1]]+pred3[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")

plot(ep.orig, pred2[[1]], type = "l", lwd = 2, col = "blue", xlab = "Elevation (m)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F)
matlines(ep.orig, cbind(pred2[[1]]-pred2[[2]], pred2[[1]]+pred2[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")

plot(date.orig, pred4[[1]], type = "l", lwd = 2, col = "blue", xlab = "Survey date (1 = April 1)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F) 
matlines(date.orig, cbind(pred4[[1]]-pred4[[2]], pred4[[1]]+pred4[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")

plot(dur.orig, pred5[[1]], type = "l", lwd = 2, col = "blue", xlab = "Survey duration (min)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F) 
matlines(dur.orig, cbind(pred5[[1]]-pred5[[2]], pred5[[1]]+pred5[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")

barplot(pred6[[1]], names.arg = c("1", "2", "3"), ylim = c(0,1), ylab = "Expected detection", xlab = "Survey Number")
segments(c(0.7,1.9, 3.1), pred6[[1]]-pred6[[2]], c(0.7,1.9, 3.1), pred6[[1]]+ pred6[[2]], lwd = 2)

# Predicciones en un gradiente de dos covariables 

for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(x=0, y=0, elev=ep[i], forest=fp[j], iLength=0) pred.matrix1[i,j] <- predict(fm5ZIP,type="state", newdata=newData)[1,1]
  } }

pred.matrix2 <- array(NA, dim = c(100, 100)) for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], date=datep[j], dur=0, time = factor("2", levels = c("1", "2", "3")))
    pred.matrix2[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1] }
}

pred.matrix3 <- array(NA, dim = c(100, 100)) for(i in 1:100){
  for(j in 1:100){ 
    newData <- data.frame(elev=ep[i], date=0, dur=durp[j], time = factor("2", levels = c("1", "2", "3")))
  pred.matrix3[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1]
  } 
}

# Gráficos
par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5) mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

image(x=ep.orig, y=fp.orig, z=pred.matrix1, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Forest cover (%)")
contour(x=ep.orig, y=fp.orig, z=pred.matrix1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250)) axis(2, at = seq(0, 100, by = 10))
box()
points(tits$elev, tits$forest, pch="+", cex=1.5)
image(x=ep.orig, y=date.orig, z=pred.matrix2, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Date (1 = April 1)")

contour(x=ep.orig, y=date.orig, z=pred.matrix2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(10, 120, by = 10))
box()
matpoints(tits$elev, date, pch="+", cex=1.5)
image(x=ep.orig, y=dur.orig, z=pred.matrix3, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Duration (min)")
contour(x=ep.orig, y=dur.orig, z=pred.matrix3, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250)) axis(2, at = seq(90, 420, by = 20))
box()
matpoints(tits$elev, dur, pch="+", cex=1.5)

# --------------------------
# FIN SCRIPT
rm(list=ls())
dev.off()
