##########################################
# Exploración y análisis de datos empleando modelos lineales: R, nlme, lme4 y quantreg.
# Autores: Salvador Mandujano y Odalis Morteo-Montiel
##########################################

# Paquetes
library(nlme) # para analizar LM mixtos
library(lme4) # para GLM mixtos
library(mgcv) # para modelos GAMs
library(quantreg) # para regresiones cuantiles
library(psych) # para opciones gráficas
library(akima) # para gráfico de interpolación

# ------------------------------------
# Datos
datos <- read.csv("datos.csv")
attach(datos)

# ------------------------------------
# Exploración y visualización de datos

### Colinearidad
ejem <- data.frame(Conteo, Densidad, HSI, Lluvia)
pairs(ejem)
cor.plot(cor(ejem))

### Prueba de normalidad de los datos
Modelo_LM <- lm(Densidad ~ HSI)
par(mfrow = c(2,2))
plot(Modelo_LM)

### Análisis de datos extremos "outliers"
modelo.inf <- lm(Densidad ~ HSI)
influence.measures(modelo.inf)$is.inf

summary.aov(lm(Densidad[-7] ~ HSI[-7]))

### Exploración visual de diferentes modelos

par(mfrow = c(2, 3), mar = c(3,4,3,3))

# modelo lineal:
plot(Conteo, Densidad, frame.plot = F, las=1, xlab = "Número de fotos", ylab = "Densidad (ind/km2)", main = "Modelo lineal simple")
abline(lm(Densidad ~ Conteo), col = "red") 

# modelo lowess:
plot(Conteo, Densidad, frame.plot = F, las=1, xlab = "Número de fotos", ylab = "Densidad (ind/km2)", main = "Modelo lowess")
lines(lowess(Conteo, Densidad), col = "red") 

# modelo loess:
plot(Conteo, Densidad, frame.plot = F, las=1, xlab = "Número de fotos", ylab = "Densidad (ind/km2)", main = "Modelo loess")
mod1 <- loess(Densidad ~ Conteo) 
x <- 0:120
y <- predict(mod1, data.frame(Conteo = x))
lines(x, y, col = "red")

# modelo polinomial:
plot(Conteo, Densidad, frame.plot = F, las=1, xlab = "Número de fotos", ylab = "Densidad (ind/km2)", main = "Modelo polinomial")
mod3 <- lm(Densidad ~ Conteo + I(Conteo^2))
x <- 0:120
y <- predict(mod3, list(Conteo = x))
lines(x, y, col = "red") 

# modelo GAM:
plot(Conteo, Densidad, frame.plot = F, las=1, xlab = "Número de fotos", ylab = "Densidad (ind/km2)", main = "Modelo GAM")
mod2 <- gam(Densidad ~ s(Conteo))
x <- 0:120
y <- predict(mod2, list(Conteo = x))
lines(x, y, col = "red")

# ------------------------------------
# Modelos lineales (LM) y mixtos (LMM)

### Regresión lineal simple
Modelo_LM <- lm(Densidad ~ HSI)

summary(Modelo_LM)

par(mfrow = c(1,1))
plot(HSI, Densidad, frame.plot = F, las=1, pch = 16, col = "skyblue", cex =2, ylim = c(0,8), xlab = "Calidad del hábitat (HSI)", ylab = "Densidad (venados/km2)")
abline(lm(Densidad ~ HSI), col = "red", lwd = 3)
text(HSI, Densidad, labels = Sitio, pos = 3, cex = 0.6)

# ------------------------------------
### LM transformación logNorm
Modelo_logN <- lm(log(Conteo) ~ Lluvia)
summary(Modelo_logN)
plot(Lluvia, log(Conteo), bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Precipitación (mm)", ylab = "Log(Número de venados contados)", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7, xlim = c(300,1400), ylim = c(2,5))
abline(lm(log(Conteo) ~ Lluvia), col = "red", lwd =3)
text(Lluvia, log(Conteo), labels = Sitio, pos = 3, cex = 0.6)

# ------------------------------------
### LM ANCOVA
Modelo_Anco <- lm(Densidad ~ HSI * TVeg)
summary(Modelo_Anco)
plot(HSI, Densidad, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Calidad de hábitat (HSI)", ylab = "Densidad (venados/km2)", pch = 16, col = "skyblue", cex = 2, main = "", cex.main = 0.7, ylim = c(0,8), xlim = c(0,0.9))
text(HSI, Densidad, labels = TVeg, pos = 3, cex = 0.6)
# el siguiente bucle es para representar las líneas:
Punto <- as.factor(Punto)
for (i in 1:18) {
  J <- Punto == i  
  x1 <- HSI[J]
  y1 <- Modelo_Anco$fitted[J]
  Ord <- order(x1)
  lines(x1[Ord], y1[Ord], col = "red", lwd = 2.5) 
}

# ------------------------------------
### LM regresión múltiple
Modelo_Mul <- lm(Densidad ~ HSI + Lluvia)
summary(Modelo_Mul)
graf_inter <- interp(Lluvia, HSI, Densidad, duplicate = T)
filled.contour(graf_inter, col = topo.colors(15), cex.lab = 1.0, key.title = title(main = "Densidad"), ylab = "Calidad de hábitat (HSI)", xlab = "Precipitación (mm)")
text(Lluvia, HSI, labels = Sitio, cex = 0.9, pos = 1, col = "black")

# ------------------------------------
### LM regresión polinomial
Modelo_Polin <- lm(Densidad ~ Lluvia + I(Lluvia^2) + I(Lluvia^3))
summary(Modelo_Polin)
# predicción de intervalos y ajuste polinomial:
prediccion <- data.frame(Lluvia = seq(300, 1400, 50))
predict(Modelo_Polin, interval = "pred", newdata = prediccion)
pc <- predict(Modelo_Polin, newdata = prediccion, interval = "conf")
# grafico:
plot(Lluvia, Densidad, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Precipitación (mm)", ylab = "Densidad de venados (ins/km2)", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7, xlim = c(300,1400), ylim = c(0,10))
matlines(prediccion$Lluvia, pc, lty= c(1,3,3), col = "red", lwd = c(3,1,1))

# ------------------------------------
### LMM con intercepta aleatoria
Modelo_lmm <- lme(Densidad ~ HSI, random = ~1 | Punto)
summary(Modelo_lmm)
random.effects(Modelo_lmm)
# predicciones componente fijo:
F0a <- fitted(Modelo_lmm, level = 0) 
# predicciones componente aleatorio:
F1a <- fitted(Modelo_lmm, level = 1) 
# gráfico:
plot(HSI, F0a, lwd = 3, type= "l", ylim = c(0,8), las = 1, col = "red", frame.plot = F, cex.main = 0.7, xlab = "Calidad hábitat HSI)", ylab= "Venados (ind/km2)") # línea gruesa: efecto fijo
points(HSI, Densidad, pch = 16, col = "skyblue", cex = 2)
for (i in 1:18) {
  x1a <- HSI[Punto == i]
  y1a <- F1a[Punto == i]
  Ka <- order(x1a)
  lines(sort(x1a), y1a[Ka], col = "red", lty = 2, lwd = 1.5) } # efectos de cada localidad
text(HSI, Densidad, TVeg, cex = 0.7)

# ------------------------------------
### LMM con intercepta y pendiente aleatoria
Modelo_3b <- lme(Densidad ~ HSI, random= ~HSI | Punto, method= "REML")
summary(Modelo_3b)
# predicciones componente fijo:
F0b <- fitted(Modelo_3b, level= 0) 
# predicciones componente aleatorio: 
F1b <- fitted(Modelo_3b, level= 1) 

plot(HSI, F0b, lwd = 3, col = "red", type = "l", las = 1, frame.plot = F, ylim = c(0, 8), cex.main = 0.7, xlab = "Calidad habitat (HSI)", ylab = "Venados (ind/km2)") # línea gruesa: efecto fijo
points(HSI, Densidad, pch = 16, col = "skyblue", cex = 2)

for (i in 1:18) {
  x1b <- HSI[Punto == i]
  y1b <- F1b[Punto == i]
  Kb <- order(x1b)
  lines(sort(x1b), y1b[Kb], lty = 2, col = "red", lwd = 1.5) } # efectos de cada localidad
text(HSI, Densidad, TVeg, cex = 0.7)

# -------------------------------------
# Modelos lineales generalizados (GLM) y mixtos (GLMM)

### GLM regresión Poisson
Modelo_Pois <- glm(Conteo ~ Lluvia, family = "poisson")
summary(Modelo_Pois)
par(mfrow = c(1,1))

plot(Lluvia, Conteo, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Precipitación (mm)", ylab = "Número de venados contados", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7, xlim = c(300,1400), ylim = c(0,120))
par(new = TRUE)
curve(exp(Modelo_Pois$coefficients[1] + (Modelo_Pois$coefficients[2] * x)), 300, 1400, col = "red", ylim = c(0,120), xlab = "", ylab = "", lwd = 3, axes = F)
text(Lluvia, Conteo, labels = Sitio, pos = 3, cex = 0.6)

# -------------------------------------
### GLM prueba-t tipo Poisson
Vacas <- as.factor(Vacas)
Modelo_tSt <- glm(Conteo ~ Vacas-1, family = "poisson")
summary(Modelo_tSt)
boxplot(Conteo ~ Vacas, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Presencia de ganado", ylab = "Número de venados contados", col = "skyblue", cex =2)

# -------------------------------------
### GLM ANOVA de una vía tipo Poisson
Modelo_aovPois <- glm(Conteo ~ TVeg - 1, family = "poisson")
summary(Modelo_aovPois)
plot(Conteo ~ TVeg, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Tipos de vegetación", ylab = "Número de venados contados", col = "skyblue", cex =2)

# -------------------------------------
### GLMM tipo Poisson
Modelo_glmm <- glmer(Conteo ~ HSI + Lluvia + (1 | Sitio), family = "poisson")
summary(Modelo_glmm)

# rescalar variables
(HSI.scal <-  scale(HSI, center= TRUE, scale=TRUE))
(Lluvia.scal <-  scale(Lluvia, center= TRUE, scale=TRUE))

Modelo_glmm2 <- glmer(Conteo ~ HSI.scal + Lluvia.scal + (1 | Sitio), family = "poisson")
summary(Modelo_glmm2)

Modelo_glmm3 <- glmer(Conteo ~ HSI.scal + (1 | Sitio), family = "poisson")
summary(Modelo_glmm3)

Modelo_glmm4 <- glmer(Conteo ~ Lluvia.scal + (1 | Sitio), family = "poisson")
summary(Modelo_glmm4)

AIC(Modelo_glmm, Modelo_glmm2, Modelo_glmm3, Modelo_glmm4)

# -------------------------------------
### GLM con distribución Gamma
par(mfrow = c(2,1))

### Simulación
# datos originales
hist(Densidad, xlim = c(0,10), xlab = "Densidad (venados/km2)", main = "datos originales (n = 18)", nclas = 10, col = "lightblue", border = "white", ylab = "Frecuencia")
abline(v = median(Densidad), lwd = 3, col = "red") 

# simulación de distribución gamma:
rate <- mean(Densidad) / var(Densidad)
shape <- rate * mean(Densidad)
sim_gamma <- rgamma(n = 10000, shape = shape, rate = rate)
hist (sim_gamma, main = "datos simulados (n = 10,000)", col = "lightblue", border = "white", nclass = 80, xlim = c(0, 10), xlab = "Densidad (venados/km2)", ylab = "Frecuencia")
abline(v = median(sim_gamma), lwd = 3, col = "red") 

### Modelo GAMMA
Modelo_Gam <- glm(Densidad + 0.001 ~ HSI, family = Gamma(link = "inverse"))
summary(Modelo_Gam)
par(mfrow = c(1,1))
plot(Densidad ~ HSI, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Calidad habitat (HSI)", ylab = "Venados (ind/km2)", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7, ylim = c(0,8))
y2 <- predict(Modelo_Gam, data.frame(HSI = seq(min(HSI), max(HSI), length = 20)), se = TRUE)
lines(seq(min(HSI), max(HSI), length = 20), 1/(y2$fit), lwd = 3, col = "red")
text(HSI, Densidad, labels = Sitio, pos = 3, cex = 0.6)

# -------------------------------------
### GLM regresión logística o binomial
Modelo_Bin <- glm(datos$Vacas ~ Densidad, family = binomial)
summary(Modelo_Bin)
# creamos una secuencia de datos en "x"
x1 <- seq(0, 20, 0.1) 
# predecimos valores de "y" a partir del modelo
y1 <- predict(Modelo_Bin, list(Densidad = x1), type = "response")

plot(Densidad, datos$Vacas, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Densidad (venados/km2)", ylab = "Presencia de vacas", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7)
lines(x1, y1, col = "red", lwd = 3)

# -------------------------------------
### Regresión cuantil
Modelo_Rquan <- rq(Densidad ~ Lluvia, tau = c(0.05,0.50,0.95))
summary(Modelo_Rquan)
plot(Lluvia, Densidad, bty = "l", frame.plot = F, las = 1, cex.lab = 0.8, xlab = "Precipitación (mm)", ylab = "Densidad de venados (ins/km2)", pch = 16, col = "skyblue", cex =2, main = "", cex.main = 0.7, xlim = c(300,1400), ylim = c(0,10))
abline(rq(Densidad ~ Lluvia, tau = 0.50), col = "red", lwd = 1, lty = 2) 
# modelo cuartiles
taus <- c(0.95)
for( i in 1:length(taus)){
  abline(rq(Densidad ~ Lluvia, tau = taus[i]), lwd = 3, col = "red")
}

# -------------------------------------
## Selección de modelos: criterio Akaike
# poner aquí los nombres de los modelos que se quieren evaluar

modelos <- list(Modelo_x, Modelo_y, Modelo_z)

aic <- unlist(lapply(modelos, AIC))

modelos_n <- c("Modelo_x", "Modelo_y", "Modelo_z")

resultados <- data.frame(modelos_n, aic)

resultados[order(aic),]

# -----------------------------
# FIN SCRIPT

rm(list = ls()) # elimina datos
dev.off() # elimina gráficos

