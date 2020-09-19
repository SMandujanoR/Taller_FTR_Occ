#############################
# Paquete: colapsar
# Creado por: SMandujanoR
# Última modificación: Septiembre 18, 2020

# Función creada para simular gráficamente el efecto de la variación en el número de sitios (M), ocasiones (J), Lambda (lambda), probabilidad de detección (p) y ocupación (psi), colpadando los datos en 1,2,3,4,5,6,10,12,15,20 y 30 días (ocasiones).

# La función simultáneamente genera los gráficos para los modelos N-mixtos, Royle-Nichols y Ocupación simple, empleando el paquete unmarked. No emplea covariables.

# La línea roja indica el valor base conocido del parámetro. La diferencia entre la línea roja y los valores simulados en color azul, indica el sesgo.

# Si se quieren guardar los gráficos en formato jpg, se deb activar las funciones:

#jpeg(filename= "nombre.jpg", width= 7000, height= 6000, units= "px", res=1200)
#dev.off()

####################################################

fn.colapso <- function(M, J, lambda, p, psi) {
  require(unmarked)
 
  # simulación de datos
  c1d <- matrix(NA, nrow= M, ncol= J)
  N <- rpois(n= M, lambda= lambda)
  for(j in 1:J){
    c1d[,j] <- rbinom(n= M, size= N, prob= p)
  }
  
  # MODELO N-MIXTO
  # preparación de datos 
  data_1d <- unmarkedFramePCount(y= c1d)
  mod_1d <- pcount(~1 ~1,data= data_1d) 
  N_1d <- backTransform(mod_1d,"state")
  p_1d <- backTransform(mod_1d,"det")
  
  # rutina para colapsar datos
  dias <- 2
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c2d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_2d <- unmarkedFramePCount(y= c2d)
  mod_2d <- pcount(~1 ~1,data= data_2d) 
  N_2d <- backTransform(mod_2d,"state")
  p_2d <- backTransform(mod_2d,"det")
  # ------------------------------
  dias <- 3
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c3d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_3d <- unmarkedFramePCount(y= c3d)
  mod_3d <- pcount(~1 ~1,data= data_3d)
  N_3d <- backTransform(mod_3d,"state")
  p_3d <- backTransform(mod_3d,"det")
  # ---------------------
  dias <- 4
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c4d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_4d <- unmarkedFramePCount(y= c4d)
  mod_4d <- pcount(~1 ~1,data= data_4d) 
  N_4d <- backTransform(mod_4d,"state")
  p_4d <- backTransform(mod_4d,"det")
  # ---------------------
  dias <- 5
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c5d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_5d <- unmarkedFramePCount(y= c5d)
  mod_5d <- pcount(~1 ~1,data= data_5d) 
  N_5d <- backTransform(mod_5d,"state")
  p_5d <- backTransform(mod_5d,"det")
  # ---------------------
  dias <- 6
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c6d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_6d <- unmarkedFramePCount(y= c6d)
  mod_6d <- pcount(~1 ~1,data= data_6d)
  N_6d <- backTransform(mod_6d,"state")
  p_6d <- backTransform(mod_6d,"det")
  # ---------------------
  dias <- 10
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c10d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_10d <- unmarkedFramePCount(y= c10d)
  mod_10d <- pcount(~1 ~1,data= data_10d)
  N_10d <- backTransform(mod_10d,"state")
  p_10d <- backTransform(mod_10d,"det")
  # ---------------------
  dias <- 12
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c12d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_12d <- unmarkedFramePCount(y= c12d)
  mod_12d <- pcount(~1 ~1,data= data_12d)
  N_12d <- backTransform(mod_12d,"state")
  p_12d <- backTransform(mod_12d,"det")
  # ---------------------
  dias <- 15
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c15d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_15d <- unmarkedFramePCount(y= c15d)
  mod_15d <- pcount(~1 ~1,data= data_15d) 
  N_15d <- backTransform(mod_15d,"state")
  p_15d <- backTransform(mod_15d,"det")
  # ---------------------
  dias <- 20
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c20d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_20d <- unmarkedFramePCount(y= c20d)
  mod_20d <- pcount(~1 ~1,data= data_20d)
  N_20d <- backTransform(mod_20d,"state")
  p_20d <- backTransform(mod_20d,"det")
  # ---------------------
  dias <- 30
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c30d <- tapply(c1d, list(row(c1d), group[col(c1d)]), sum)
  data_30d <- unmarkedFramePCount(y= c30d)
  mod_30d <- pcount(~1 ~1,data= data_30d) 
  N_30d <- backTransform(mod_30d,"state")
  p_30d <- backTransform(mod_30d,"det")
  
  # Para generar gráficos
  
  #jpeg(filename= "Nmix_SpCom.jpg", width= 7000, height= 6000, units= "px", res=1200)
  
  par(mfrow = c(3,2), mar = c(5,5,5,5))
  x <- c(1,2,3,4,5,6,10,12,15,20,30)
  
  lambda_mean <- round(cbind(N_1d@estimate, N_2d@estimate, N_3d@estimate, N_4d@estimate, N_5d@estimate, N_6d@estimate, N_10d@estimate, N_12d@estimate, N_15d@estimate, N_20d@estimate, N_30d@estimate),1)
  plot(x, lambda_mean, type = "b", frame = F, las = 1, xlab = "Colapso días", ylab = expression(lambda), pch = 16, main = "N-mixto", cex.axis = 0.75, col = "skyblue")
  abline(h = lambda, col = "red")
  
  pd <- round(cbind(p_1d@estimate, p_2d@estimate, p_3d@estimate, p_4d@estimate, p_5d@estimate, p_6d@estimate, p_10d@estimate, p_12d@estimate, p_15d@estimate, p_20d@estimate, p_30d@estimate), 3)
  plot(x, pd, type = "b", frame = F, las = 1, ylim = c(0,1), xlab = "Colpaso días", ylab = "P.detección", main = "", cex.axis = 0.75, pch = 16, col = "skyblue")
  abline(h = p, col = "red")
  
  #dev.off()

# MODELO ROYLE-NICHOLS
  Occ1d <- c1d
  Occ1d[Occ1d >1] <- 1
  # ----------------------------
  data_1d <- unmarkedFrameOccu(y= Occ1d)
  mod_1d <- occuRN(~1 ~1,data= data_1d) 
  OccRN_1d <- backTransform(mod_1d,"state")
  OccRN_p_1d <- backTransform(mod_1d,"det")
  # --------------------------------
  dias <- 2
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c2d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ2d <- c2d 
  Occ2d[Occ2d >1] <- 1 
  data_2d <- unmarkedFrameOccu(y= Occ2d)
  mod_2d <- occuRN(~1 ~1,data= data_2d) 
  OccRN_2d <- backTransform(mod_2d,"state")
  OccRN_p_2d <- backTransform(mod_2d,"det")
  # ------------------------------
  dias <- 3
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c3d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ3d <- c3d 
  Occ3d[Occ3d >1] <- 1 
  data_3d <- unmarkedFrameOccu(y= Occ3d)
  mod_3d <- occuRN(~1 ~1,data= data_3d) 
  OccRN_3d <- backTransform(mod_3d,"state")
  OccRN_p_3d <- backTransform(mod_3d,"det")
  # ---------------------
  dias <- 4
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c4d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ4d <- c4d 
  Occ4d[Occ4d >1] <- 1 
  data_4d <- unmarkedFrameOccu(y= Occ4d)
  mod_4d <- occuRN(~1 ~1,data= data_4d) 
  OccRN_4d <- backTransform(mod_4d,"state")
  OccRN_p_4d <- backTransform(mod_4d,"det")
  # ---------------------
  dias <- 5
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c5d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ5d <- c5d 
  Occ5d[Occ5d >1] <- 1 
  data_5d <- unmarkedFrameOccu(y= Occ5d)
  mod_5d <- occuRN(~1 ~1,data= data_5d) 
  OccRN_5d <- backTransform(mod_5d,"state")
  OccRN_p_5d <- backTransform(mod_5d,"det")
  # ---------------------
  dias <- 6
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c6d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ6d <- c6d 
  Occ6d[Occ6d >1] <- 1 
  data_6d <- unmarkedFrameOccu(y= Occ6d)
  mod_6d <- occuRN(~1 ~1,data= data_6d) 
  OccRN_6d <- backTransform(mod_6d,"state")
  OccRN_p_6d <- backTransform(mod_6d,"det")
  # ---------------------
  dias <- 10
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c10d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ10d <- c10d 
  Occ10d[Occ10d >1] <- 1 
  data_10d <- unmarkedFrameOccu(y= Occ10d)
  mod_10d <- occuRN(~1 ~1,data= data_10d) 
  OccRN_10d <- backTransform(mod_10d,"state")
  OccRN_p_10d <- backTransform(mod_10d,"det")
  # ---------------------
  dias <- 12
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c12d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ12d <- c12d 
  Occ12d[Occ12d >1] <- 1 
  data_12d <- unmarkedFrameOccu(y= Occ12d)
  mod_12d <- occuRN(~1 ~1,data= data_12d) 
  OccRN_12d <- backTransform(mod_12d,"state")
  OccRN_p_12d <- backTransform(mod_12d,"det")
  # ---------------------
  dias <- 15
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c15d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ15d <- c15d 
  Occ15d[Occ15d >1] <- 1 
  data_15d <- unmarkedFrameOccu(y= Occ15d)
  mod_15d <- occuRN(~1 ~1,data= data_15d) 
  OccRN_15d <- backTransform(mod_15d,"state")
  OccRN_p_15d <- backTransform(mod_15d,"det")
  # ---------------------
  dias <- 20
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c20d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ20d <- c20d 
  Occ20d[Occ20d >1] <- 1 
  data_20d <- unmarkedFrameOccu(y= Occ20d)
  mod_20d <- occuRN(~1 ~1,data= data_20d) 
  OccRN_20d <- backTransform(mod_20d,"state")
  OccRN_p_20d <- backTransform(mod_20d,"det")
  # ---------------------
  dias <- 30
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c30d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ30d <- c30d 
  Occ30d[Occ30d >1] <- 1 
  data_30d <- unmarkedFrameOccu(y= Occ30d)
  mod_30d <- occuRN(~1 ~1,data= data_30d) 
  OccRN_30d <- backTransform(mod_30d,"state")
  OccRN_p_30d <- backTransform(mod_30d,"det")
  # -----------------------------
  # RESULTADOS FINALES
  #jpeg(filename= "Figura 1 Royle-Nichols.jpg", width= 7000, height= 6000, units= "px", res=1200)
  
  x <- c(1,2,3,4,5,6,10,12,15,20,30)
  
  lambda_mean <- round(cbind(OccRN_1d@estimate, OccRN_2d@estimate, OccRN_3d@estimate, OccRN_4d@estimate, OccRN_5d@estimate, OccRN_6d@estimate, OccRN_10d@estimate, OccRN_12d@estimate, OccRN_15d@estimate, OccRN_20d@estimate, OccRN_30d@estimate),1)
  
  plot(x,lambda_mean, type= "b", frame= F, las= 1, xlab= "Colpaso días", ylab= expression(lambda), pch= 16, main= "Royle-Nichols", cex.axis= 0.75, col = "skyblue")
  abline(h= lambda, col= "red")
  
  pd <- round(cbind(OccRN_p_1d@estimate, OccRN_p_2d@estimate, OccRN_p_3d@estimate, OccRN_p_4d@estimate, OccRN_p_5d@estimate, OccRN_p_6d@estimate, OccRN_p_10d@estimate, OccRN_p_12d@estimate, OccRN_p_15d@estimate, OccRN_p_20d@estimate, OccRN_p_30d@estimate),3)
  plot(x,pd, type= "b", frame= F, las= 1, ylim= c(0,1), xlab= "Colapso días", ylab= "P.detección", main= "", cex.axis= 0.75, pch= 16, col = "skyblue")
  abline(h= p, col= "red")
  
  #dev.off()
  
# MODELO OCUPACIÓN SIMPLE
  c1d <- matrix(NA, nrow= M, ncol= J)
  z <- rbinom(n= M, size= 1, prob= psi) 
  for (j in 1:J){
    c1d[,j] <- rbinom(n= M, size= 1, prob= z*p) 
  }
  Occ1d <- c1d 
  Occ1d[Occ1d >1] <- 1 
  
  # ----------------------------
  data_1d <- unmarkedFrameOccu(y= Occ1d)
  mod_1d <- occu(~1 ~1,data= data_1d) 
  Occ_1d <- backTransform(mod_1d,"state")
  Occ_p_1d <- backTransform(mod_1d,"det")
  # --------------------------------
  dias <- 2
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c2d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ2d <- c2d 
  Occ2d[Occ2d >1] <- 1 
  data_2d <- unmarkedFrameOccu(y= Occ2d)
  mod_2d <- occu(~1 ~1,data= data_2d) 
  Occ_2d <- backTransform(mod_2d,"state")
  Occ_p_2d <- backTransform(mod_2d,"det")
  # ------------------------------
  dias <- 3
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c3d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ3d <- c3d 
  Occ3d[Occ3d >1] <- 1 
  data_3d <- unmarkedFrameOccu(y= Occ3d)
  mod_3d <- occu(~1 ~1,data= data_3d) 
  Occ_3d <- backTransform(mod_3d,"state")
  Occ_p_3d <- backTransform(mod_3d,"det")
  # ---------------------
  dias <- 4
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c4d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ4d <- c4d 
  Occ4d[Occ4d >1] <- 1 
  data_4d <- unmarkedFrameOccu(y= Occ4d)
  mod_4d <- occu(~1 ~1,data= data_4d) 
  Occ_4d <- backTransform(mod_4d,"state")
  Occ_p_4d <- backTransform(mod_4d,"det")
  # ---------------------
  dias <- 5
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c5d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ5d <- c5d 
  Occ5d[Occ5d >1] <- 1 
  data_5d <- unmarkedFrameOccu(y= Occ5d)
  mod_5d <- occu(~1 ~1,data= data_5d) 
  Occ_5d <- backTransform(mod_5d,"state")
  Occ_p_5d <- backTransform(mod_5d,"det")
  # ---------------------
  dias <- 6
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c6d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ6d <- c6d 
  Occ6d[Occ6d >1] <- 1 
  data_6d <- unmarkedFrameOccu(y= Occ6d)
  mod_6d <- occu(~1 ~1,data= data_6d) 
  Occ_6d <- backTransform(mod_6d,"state")
  Occ_p_6d <- backTransform(mod_6d,"det")
  # ---------------------
  dias <- 10
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c10d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ10d <- c10d 
  Occ10d[Occ10d >1] <- 1 
  data_10d <- unmarkedFrameOccu(y= Occ10d)
  mod_10d <- occu(~1 ~1,data= data_10d) 
  Occ_10d <- backTransform(mod_10d,"state")
  Occ_p_10d <- backTransform(mod_10d,"det")
  # ---------------------
  dias <- 12
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c12d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ12d <- c12d 
  Occ12d[Occ12d >1] <- 1 
  data_12d <- unmarkedFrameOccu(y= Occ12d)
  mod_12d <- occu(~1 ~1,data= data_12d) 
  Occ_12d <- backTransform(mod_12d,"state")
  Occ_p_12d <- backTransform(mod_12d,"det")
  # ---------------------
  dias <- 15
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c15d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ15d <- c15d 
  Occ15d[Occ15d >1] <- 1 
  data_15d <- unmarkedFrameOccu(y= Occ15d)
  mod_15d <- occu(~1 ~1,data= data_15d) 
  Occ_15d <- backTransform(mod_15d,"state")
  Occ_p_15d <- backTransform(mod_15d,"det")
  # ---------------------
  dias <- 20
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c20d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ20d <- c20d 
  Occ20d[Occ20d >1] <- 1 
  data_20d <- unmarkedFrameOccu(y= Occ20d)
  mod_20d <- occu(~1 ~1,data= data_20d) 
  Occ_20d <- backTransform(mod_20d,"state")
  Occ_p_20d <- backTransform(mod_20d,"det")
  # ---------------------
  dias <- 30
  clases <- J/dias
  group <- sort(rep(1:clases,dias))
  c30d <- tapply(Occ1d, list(row(Occ1d), group[col(Occ1d)]), sum)
  Occ30d <- c30d 
  Occ30d[Occ30d >1] <- 1 
  data_30d <- unmarkedFrameOccu(y= Occ30d)
  mod_30d <- occu(~1 ~1,data= data_30d) 
  Occ_30d <- backTransform(mod_30d,"state")
  Occ_p_30d <- backTransform(mod_30d,"det")
  # -----------------------------
  # RESULTADOS FINALES
  #jpeg(filename= "Figura 1 Single Occupancy.jpg", width= 7000, height= 6000, units= "px", res=1200)
  x <- c(1,2,3,4,5,6,10,12,15,20,30)
  
  occupancy <- round(cbind(Occ_1d@estimate, Occ_2d@estimate, Occ_3d@estimate, Occ_4d@estimate, Occ_5d@estimate, Occ_6d@estimate, Occ_10d@estimate, Occ_12d@estimate, Occ_15d@estimate, Occ_20d@estimate, Occ_30d@estimate),1)
  plot(x,occupancy, type= "b", frame= F, las= 1, ylim= c(0,1), xlab= "Colpaso días", ylab= expression(psi), pch= 16, main= "Ocupación", cex.axis= 0.75, col = "skyblue")
  abline(h= psi, col= "red")
  
  pd <- round(cbind(Occ_p_1d@estimate, Occ_p_2d@estimate, Occ_p_3d@estimate, Occ_p_4d@estimate, Occ_p_5d@estimate, Occ_p_6d@estimate, Occ_p_10d@estimate, Occ_p_12d@estimate, Occ_p_15d@estimate, Occ_p_20d@estimate, Occ_p_30d@estimate),3)
  plot(x,pd, type= "b", frame= F, las= 1, ylim= c(0,1), xlab= "Colapso días", ylab= "P. detección", main= "", cex.axis= 0.75, pch= 16, col = "skyblue")
  abline(h= p, col= "red")
  
  #dev.off()
}
# --------------------------
# FIN SCRIPT

