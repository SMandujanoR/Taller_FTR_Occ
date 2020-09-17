###################################################
#
# Modelo de ocupación una especie una temporada
#
####################################################

# Cargar datos 
Urci <- read.csv("DetHist_UrciClase.csv", sep=",")
Urcihist <- Urci[,-1]
y <- Urcihist

# Generación datos en formato para el modelo

str( win.data2 <- list(y = y, M = nrow(y), J = ncol(y)) )

# Modelo bayesiano en JAGS

cat(file = "model1.txt",
    "
model {
# Información a priori
   psi ~ dunif(0, 1)
   p ~ dunif(0, 1)
# Probabilidad
   for (i in 1:M) {              # Bucle sobre los sitios
      z[i] ~ dbern(psi)          # Modelo de estado (proceso ecológico)
      for (j in 1:J) {           # Bucle sobre los muestreos replicados
         y[i,j] ~ dbern(z[i]*p)  # Proceso de observación
      }
   }
}

")

# Valores de inicio

zst <- apply(y, 1, max, na.rm=T) # Observed occurrence as starting values for z
inits <- function() {list(z = zst)}

# Parametros que deseo evaluar

params <- c("psi", "p")

# Definición de Cadenas de Marcov Monte Carlo

ni <- 5000; nt <- 1; nb <- 1000; nc <- 3

# Ejecución modelo

library(jagsUI)

fm2 <- jags(win.data2, inits, params, "model1.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb)

print(fm2, dig = 3)


