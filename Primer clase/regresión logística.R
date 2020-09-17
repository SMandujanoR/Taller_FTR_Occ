# Regresión logística o binomial
# Ejemplo de incidiencia (0 y 1) de una especie de mono en fragmentos de selva.

# ---------------------
# leer datos archivos
datos <- read.table("incidencia.txt", header = T)
View(datos)
names(datos)

# ----------------------------------
# modelos para cada variable por separado
modelo_area <- glm(incidencia ~ area, data = datos, family = binomial)
summary(modelo_area)

# gráfico
par(mfrow = c(1,1))
plot(datos$area, datos$incidencia, pch = 21, bg = "skyblue", col = "black", cex = 1.5, frame.plot = F, main ="Modelo área", xlab =  "Tamaño del fragmento", ylab = "Probabilidad de incidencia")

# creamos una secuencia de datos en "x"
x1 <- seq(0, 9, 0.01) # 900 datos

# predecimos valores de "y" a partir del modelo
y1 <- predict(modelo_area, list(area = x1), type = "response")
lines(x1, y1, col = "red", lwd = 4)

# modelo aislamiento
modelo_aislamiento <- glm(incidencia ~ aislamiento, data = datos, family =  binomial)
summary(modelo_aislamiento)

plot(datos$aislamiento, datos$incidencia, pch = 21, bg = "skyblue", col = "black", cex = 1.5, frame.plot = F, main ="Modelo aislamiento", xlab =  "Tamaño del fragmento", ylab = "Probabilidad de incidencia")

x2 <- seq(0,10,0.01)
y2 <- predict(modelo_aislamiento, list(aislamiento = x2), type = "response")
lines(x2, y2, col = "red", lwd = 4)

# -----------------------------
# modelos con ambas variables
modelo_2 <- glm(incidencia ~ area + aislamiento, data = datos, family = "binomial")
summary(modelo_2)

# ----------------------
# gráfico 2D
library(akima)

par(mfrow = c(1, 1))

occ <- interp(datos$area, datos$aislamiento, datos$incidencia, duplicate = T)

filled.contour(occ, col = topo.colors(20), cex.lab = 1, cex.main = 0.9, key.title = title(main= "Inciden"), xlab = "Tamaño del fragmento", ylab = "Aislamiento del fragmento", main = "Interpolación en el \n espacio ecológico")

# ---------------------------------
# FIN SCRIPT

rm(list = ls()) # elimina datos
dev.off() # elimina gráficos
