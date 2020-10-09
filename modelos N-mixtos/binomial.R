# Distribución Binomial

#jpeg(filename= "binomial.jpg", width= 8000, height= 7000, units= "px", res=1200)

par(mfrow= c(2,2), mar = c(5,5,5,5))

# si p= 0.1
plot(0:15, dbinom(0:15, 20, 0.1), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados (n)", 
     main= "Distribución Binomial \n Prob. detección = 0.1", 
     cex.main= 0.8, cex.lab= 0.8, cex.axis = 0.7)

# si p= 0.3
plot(0:15, dbinom(0:15, 20, 0.3), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados (n)", 
     main= "Distribución Binomial \n Prob. detección = 0.3", 
     cex.main= 0.8, cex.lab= 0.8, cex.axis = 0.7)

# si p= 0.6
plot(0:15, dbinom(0:15, 20, 0.6), type= "h", bty= "l", 
     las=1, col= "lightblue", lwd= 6,
     ylab= "Probabilidad", 
     xlab= "No. observados (n)", 
     main= "Distribución Binomial \n Prob. detección = 0.6", 
     cex.main= 0.8, cex.lab= 0.8, cex.axis = 0.7)

# si p=0.9
plot(0:15, dbinom(0:15, 20, 0.9), type= "h", bty= "l",  
     las=1, col= "lightblue", lwd= 6, 
     ylab= "Probabilidad", 
     xlab= "No. observados (n)", 
     main= "Distribución Binomial \n Prob. detección = 0.9", 
     cex.main= 0.8, cex.lab= 0.8, cex.axis = 0.7)

#dev.off()

# -----------------
# FIN SCRIPT

rm(list = ls())
dev.off()