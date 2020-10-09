############################### 
# Kéry y Royle (2016)
# Cap.6. Modelos N-mixtos
# Sec. 6.6: Simulación del efecto del número de sitios y ocasiones de muestreo, sobre la estimación de Lambda dependiendo de la probabilidad de detección.

# IMPORTANTE: Esta simulación puede servir de base para el diseño del muestreo.
# Pero el proceso puede tardar mucho dependiendo del número de simulaciones y el procesador de las computadoras.
###############################

# Cargar paquetes:
library(unmarked)
library(AHMbook)

# ---------------------------------
# Parámetros de la simulación:

# Número de similaciones: 
simreps <- 500 

# número de sitios:
nsites <- c(10, 20, 40)

# número de ocasiones:
nreps <- c(7, 15, 30)

# construcción de matriz:
estimates <- array(NA, dim = c(2, simreps, 3, 3))

# probabilidad de detección:
p <- array(runif(n = simreps*3*3, 0.01, 0.99), dim = c(simreps, 3, 3))

# ---------------------------------
# Función para la simulación:
for(s in 1:3){ 
  for(r in 1:3){
    for(i in 1:simreps){
      cat("*** Simrep number", i, "***\n")
      data <- simNmix(nsite = nsites[s], nvisit = nreps[r], mean.lam = 5, mean.p = p[i,s,r], show.plot = F)
      umf <- unmarkedFramePCount(y = data$C) 
      fm <- pcount(~1 ~1, umf)
      estimates[,i,s,r] <- coef(fm)
    }
  }
}

# ---------------------------------
# para guardar y luego leer los resultados:
#save(data, file = "sim.RData")
#load("sim.RData")

# gráficos:

#jpeg(filename = "sim2.jpg", width = 8000, height = 7000, units = "px", res = 1200)

par(mfrow = c(3,3), mar = c(4.5, 4.5, 2, 2), cex.lab = 1.2, cex.axis = 1) 

for(s in 1:3){ 
  for(r in 1:3){
    plot(p[,s,r], exp(estimates[1,,s,r]), xlab = "Prob. detección", ylab = expression(lambda), main = "", ylim = c(0, 90), xlim = c(0, 0.5), frame = F, pch = 1, cex = 0.7, col ="skyblue", las = 1) 
    text(0.35, 50, paste("M = ", nsites[s], ", J = ", nreps[r], sep = ""), cex = 1)
    abline(h = 5, col = "red", lwd = 2) 
    lines(smooth.spline(exp(estimates[1,,s,r])~p[,s,r]), col = "blue", lwd = 2)
  } 
}

#dev.off()

# ----------------
# FIN SCRIPT

rm(list=ls())
dev.off()
