#####################################
# CAPITULO 4: SIMULACION DE DATOS
# seccion 4.2 Generación/simulación de datos de conteo # Kéry & Royle (2016)
#####################################
# PASO 1: tamaño de muestra y covariables
#######################################
  M <- 267     # número de sitios
  J <- 3       # número de réplicas temporales
  
  set.seed(24) # crea numeros aleatorios replicables
  altitud <- runif(n = M, -1, 1) # altitud estandarizada de cada sitio
  head(altitud, 10) # muestra solo los primeros 10 datos
  bosque <- runif(n = M, -1, 1) # cobertura estandarizada de cada sitio
  head(bosque, 10)
  viento <- array(runif(n = M * J, -1, 1), dim = c(M,J)) # velocidad viento
  head(viento, 10)
  
#####################################
  # PASO 2: simulación del proceso ecológico
#######################################
  mean.lambda <- 2 # abundancia media esperada de la especie
  beta0 <- log(mean.lambda) # lo mismo en escala Log
  beta1 <- -2 # efecto pendiente de "altitud"
  beta2 <- 2 # efecto pendiente de "bosque"
  beta3 <- 1 # efecto de interacción de "altitud" y "bosque"
  
  log.lambda <- beta0 + beta1 * altitud + beta2 * bosque + beta3 * altitud * bosque
  lambda <- exp(log.lambda) # Link: transformación inversa
  head(lambda,10)
  
  par(mfrow = c(1,2), cex.main = 1.5)
  curve(exp(beta0 + beta1 * x), -1, 1, col = "red", frame.plot = F, ylim = c(0,40), xlab = "Altitud", ylab = "Lambda", lwd = 3)
  par(new = TRUE)
  plot(altitud, lambda, ylab = "", xlab = "", ylim = c(0,40), axes = FALSE) 
  text(0.9,36,"a)", cex = 1)
  
  curve(exp(beta0 + beta2 * x), -1, 1,col = "red",frame.plot = F,ylim = c(0,40), xlab = "Cobertura arborea", ylab = "Lambda",lwd = 3)
  par(new = TRUE)
  plot(bosque, lambda, ylab = "", xlab = "", ylim = c(0,40), axes = FALSE) 
  text(0.9,36,"b)", cex = 1)

# -----------------------------
# interacción "altitud" x "Cobertura arbórea"

cov1 <- seq(-1,1,,100) # para "altitud"
cov2 <- seq(-1,1,,100) # para "Cobertura arbórea"
lambda.matrix <- array(NA, dim = c(100,100)) # matriz para predicciones
head(lambda.matrix, 10)

for (i in 1:100) {
  for (j in 1:100) {
    lambda.matrix[i,j] <- exp(beta0 + beta1 * cov1[i] + beta2 * cov2[j] + beta3 * cov1[i] * cov2[j])
  }
}
head(lambda.matrix,10)

par(mfrow = c(1,1), cex.main = 1.5)
mapPalette<- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = cov1, y = cov2, z = lambda.matrix, col = mapPalette(100), xlab = "Altitud", ylab = "Cobertura arborea", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = lambda.matrix, add = T, lwd = 2.0, labcex = 1.5)
matpoints(altitud, bosque, pch = 1, cex = 0.7)

# incorporando la estocasticidad
N <- rpois(n = M, lambda = lambda) # abundancia realizada
sum(N) # tamaño poblacinal total en los M sitios
table(N) # distribución frecuencia de la abundancia

#######################################
# PASO 3: Simulación del proceso observacional
#######################################

mean.detection <- 0.3 # probabilidad detección esperada media
alpha0 <- qlogis(mean.detection) # lo mismo en escala Logit
alpha1 <- 1 # efecto de la "altitud"
alpha2 <- -3 # efecto de la "viento"
alpha3 <- 0 # efectos interacción "altitud" x "viento"

logit.p <- alpha0 + alpha1 * altitud + alpha2 * viento + alpha3 * altitud * viento
p <- plogis(logit.p) # transformación inversa
head(p,10)
mean(p) # probabilidad media por "individuo"

par(mfrow = c(1,2), cex.main = 1)

matplot(altitud, p, pch = 1, ylim = c(0,1.1), axes = T, frame.plot = F, xlab = "Altitud", ylab = "Prob.deteccion (p)")
par(new = TRUE)
curve(plogis(alpha0 + alpha1 * x),-1,1, col= "red", lwd= 3, frame = F, xlab = "", ylab = "", ylim = c(0,1.1))
text(0.9,1,"a)", cex = 1)

curve(plogis(alpha0 + alpha2 * x),-1,1, col = "red", frame.plot = F, ylim = c(0,1.1), xlab = "Viento", ylab = "Prob.deteccion (p)", lwd = 3)
par(new=TRUE)
matplot(viento, p, pch = 1, ylim = c(0,1.1), xlab = "", ylab = "",axes = FALSE) 
text(0.9,1,"b)", cex = 1)

# ---------------------------------
# interacción "altitud" x "viento"

cov1 <- seq(-1,1,,100) # para "altitud"
cov2 <- seq(-1,1,,100) # para "viento"
p.matrix <- array(NA, dim = c(100,100)) # para predicciones de p
head(p.matrix,10)

for (i in 1:100) {
  for (j in 1:100) {
    p.matrix[i,j] <- plogis(alpha0 + alpha1 * cov1[i] + alpha2 * cov2[j] + alpha3 * cov1[i] * cov2[j])
  }
}
head(p.matrix,10)

par(mfrow = c(1,1), mar = c(5,4,2,2), cex.main = 1.5)
mapPalette<- colorRampPalette(c("grey","yellow","orange","red"))
image(x = cov1, y = cov2, z = p.matrix, col = mapPalette(100), xlab = "Altitud", ylab = "Viento", cex.lab = 1)
contour(x = cov1, y = cov2, z = p.matrix, add = T, lwd = 2.0, labcex = 1.5)
matpoints(altitud, viento, pch = 1, cex = 0.7, col = "black")

# aplicación del proceso observacional para producir conteos replicados para cada sitio
C <- matrix(NA, nrow = M, ncol = J)
for (i in 1:J) {
  C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}
head(C,10)

#######################################
# PASO 4: Resultados de la simulación
######################################
head(cbind("N_real"= N, "1er_conteo"= C[,1], "2er_conteo"= C[,2], "3er_conteo"= C[,3]),7)
table(C)

par(mfrow = c(2,2), cex.main = 1.5)

matplot(altitud, C, pch = 16, frame.plot = F, ylim = c(0,38), xlab = "Altitud", ylab = "Conteos (C)")
text(0.9,36,"a)", cex = 1)

matplot(bosque, C, pch = 16, frame.plot = F, ylim = c(0,38), xlab = "Cobertura arborea", ylab = "Conteos (C)")
text(0.9,36,"b)", cex = 1)

matplot(viento, C, pch = 16, frame.plot = F, ylim = c(0,38), xlab = "Viento", ylab = "Conteos (C)")
text(0.9,36,"c)", cex = 1)

hist(C, breaks = 40, col="black", border = "white", ylim = c(0,460), main = "", xlab = "Conteos (C)", ylab = "Frecuencia")
text(28,450,"d)", cex = 1)

sum(N) # abundancia total real en todos los sitios
sum(apply(C,1,max)) # abundancia total observada en todos los sitios
sum(N > 0) # numero real de sitios ocupados
sum(apply(C,1,max) > 0) # numero observado de sitios ocupados

# --------------------------
# FIN SCRIPT
rm(list=ls())

