#####################################
# MODELO N-MIXTO: EL MÁS SIMPLE DE TODOS
#
# Ni ~ Poisson(lambda) # proceso ecológico
#
# cij|Ni ~ Binomial(Ni,p) # proceso observacional
# 
###################################

# SIMULACIÓN ESPECIE ABUNDANCIA INTERMEDIA

par(mfrow = c(1,1))
set.seed(24)

# datos de entrada 
M <- 8 # número de sitios
J <- 60 # número de visitas
c <- matrix(NA, nrow = M, ncol = J)
c
lambda <- 2.5 # abundancia esperada
p <- 0.4 # probabilidad de detección por individuo

# se genera los datos de abundancia real en los sitios
(N <- rpois(n = M, lambda = lambda))

# se genera los datos de conteos en los sitios en J=2 ocaciones
for(j in 1:J){
  c[,j] <- rbinom(n = M, size = N, prob = p)
}
c

# datos conocidos de la población (valores paramétricos):
table(N)  # distribución de la abundancia real
sum(N)    # abundancia real en los sitios
sum(N > 0)  # número real de sitios ocupados
mean(N)   # abundancia promedio real (lambda)

# datos de las observaciones durante los muestreos:
table(apply(c,1,max)) # distribución de la abundancia observada
sum(apply(c,1,max)) # abundancia total observada en los sitios
sum(apply(c,1,max) > 0) # número de sitios observados ocupados
mean(apply(c,1,max)) # abundancia relativa promedio

head(cbind(N = N, count1 = c[,1], count2 = c[,2], count3 = c[,3], count4 = c[,4], count5 = c[,5], count6 = c[,6], count7 = c[,7], count8 = c[,8]))


###########################
# ANALISIS EN UNMARKED
###########################

library(unmarked)

umf <- unmarkedFramePCount(y = c)
summary(umf)
table(c) # OJO: esta tabla solo es para comparar que sale lo mismo
(fm1 <- pcount(~1 ~1, data = umf))

# se transforman los datos a otra escala
backTransform(fm1, "state")
backTransform(fm1, "det")

######################
# MODELO EN JAGS
######################

library(jagsUI) 

# Bundle data 
(win.data <- list(c = c, M = nrow(c), J = ncol(c)))

# Modelo BUGS
sink("model_jags.txt")
cat("
    model{
    
    # Priors
    lambda ~ dgamma(0.001, 0.001) # OJO: observar la distribución empleada
    p ~ dunif(0, 1)
    
    # Likelihood
    for (i in 1:M) {   
      N[i] ~ dpois(lambda) # modelo del proceso ecologico         
    for (j in 1:J) { 
      c[i,j] ~ dbin(p,N[i]) # modelo del proceso observacional  
        } # j
      } # i
    } # modelo
    ",fill = TRUE)
sink()

# Valores iniciales
Nst <- apply(c,1,max)       
inits <- function(){list(N = Nst)}

# Parámetros monitoreados
params <- c("lambda","p")

# MCMC 
ni <- 25000   
nt <- 20  
nb <- 5000   
nc <- 3

# JAGS y posteriors
fm2 <- jags(win.data, inits, params, "model_jags.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(fm2, dig = 3)
plot(fm2)

# --------------------------
# FIN SCRIPT



