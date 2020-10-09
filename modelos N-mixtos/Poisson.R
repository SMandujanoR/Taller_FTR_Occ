# Ejemplo distribución Poisson

#jpeg(filename= "Poisson.jpg", width= 8000, height= 7000, units= "px", res=1200)

par(mfrow= c(2,2))

# Lambda = 1
plot(0:15, dpois(0:15, 1), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados", 
     main= "Distribución Poisson \n Lambda = 1", 
     cex.main= 0.9, cex.lab= 0.9)

# Lambda = 3
plot(0:15, dpois(0:15, 3), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados", 
     main= "Distribución Poisson \n Lambda = 3", 
     cex.main= 0.9, cex.lab= 0.9)

# Lambda = 5
plot(0:15, dpois(0:15, 6), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados", 
     main= "Distribución Poisson \n Lambda = 5", 
     cex.main= 0.9, cex.lab= 0.9)

# Lambda = 10
plot(0:15, dpois(0:15, 10), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados", 
     main= "Distribución Poisson \n Lambda = 10", 
     cex.main= 0.9, cex.lab= 0.9)

#dev.off()

