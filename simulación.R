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

# cargar el paquete como:
source("colapsar.R")

# definir los escenarios a simular:

# especie común:
fn.colapso(M = 50, J = 90, lambda = 5.0, p = 0.3, psi = 0.3)

# especie rara y poco muestreada:
fn.colapso(10, 30, 0.1, 0.1, 0.1)

fn.colapso(50, 90, 0.1, 0.1, 0.1)

##################################################
# FIN SCRIPT

rm(list=ls())
dev.off()
