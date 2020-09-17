
# Modelo condicional RW para la co-ocurrencia de dos especies----------

# Cargar las librerias

require(wiqid) # Por simplicidad en este caso vamos a trabajar con wiqid, pero también es posible trabajar con RPresence

# Cargar las bases de datos----------------
# En este caso vamos a cargar historias de detección de dos especies de mamíferos el lince y la zorra gris. Recordemos que este modelo asume una especie dominante y una subordinada. En este caso el lince será la especie dominante.

# En los modelos de co-ocurrencia las historias de detección de todas las especies deben tener las mismas dimensiones, es decir el mismo número de sitios y ocasiones de muestreo.

# El funcionamiento de `wiqid` es muy sencillo y solo necesitamos cargar nuestras historias de detección por separado creando un objeto para cada archivo `.csv`.


# Datos de lince:
Lyru <- read.csv("data/Lynx.hist5.csv") # sp dominante

#Datos de zorra gris:
Urci <- read.csv("data/Fox.hist5.csv") # sp subordinada

### Ajustar los modelos-----------------

# Lo primero que tenemos que considerar es que vamos a modelar teniendo en cuenta posibles hipótesis de relación espacial entre las especies seleccionadas, en ese sentido el modelo RW puede contestar las siguientes preguntas:
   #1) ¿La ocupación de una especie depende de la otra?
   #2) ¿La detección de ambas especies difiere si la otra especie esta presente?
   #3) ¿La probabilidad de detección de una especie depende de la detección de la otra especie, cuando ambas especies estuvieron presentes?
  
# Vamos a utilizar la función `occ2sps()` cuya sintaxis es muy simple; en este caso solo necesitamos poner la historia de captura de la especie dominante (Lince) y la especie subordinada (zorra gris) respetando ese orden.

m1 <- occ2sps(Lyru, Urci)
m1

# El objeto `m1` es una lista compleja de la cual nos interesa el objeto *real* que corresponde a los valores de estimados reales (ya transformados) de cada uno de los parámetros estimados. Precisamente estos valores son los que visualizamos en el código y en este caso se asume que psiBA=psiBa es decir que la probabilidad de ocupación de la especie B es igual esté o no presente la especie A, lo cual simplemente interpretamos como PsiBa = PsiB

m1$real # que es el objeto real

# Ahora empecemos a modelar "interacciones"


#psiA, psiBA, psiBa, pA=rA, pB=rBA=rBa
m2 <- occ2sps(Lyru, Urci, model = list(psiBA ~ 1))
m2
# Bajo una hipótesis de segregación ¿que esperamos de los parámetros de psi?

# Continuemos modelando, pero ahora vamos a investigar la relación en la detección condicionada a la presencia de ambas especies. En este caso le indicamos al modelo que *pA* y *pB* son diferentes de *rA* y *rBA* respectivamente. De nuevo bajo una hipótesis de segregación esperamos que p_{A},p_{B}>r_{A},r_{BA}

#psiA, psiBA=psiBa, pA, pB, rA, rBA=rBa
m3 <- occ2sps(Lyru, Urci, model = list(rA~ 1, rBa~ 1))
m3

#psiA, psiBA, psiBa, pA, pB, rA, rBA=rBa
m4<- occ2sps(Lyru, Urci, model = list(psiBA ~ 1, rA~ 1, rBa~ 1))
m4

# ¿ En que se diferencian m3 y m4 ?

# Por último vamos a analizar como la detección o no detección del Lince afecta la detección de la zorra gris . En este caso le indicamos al modelo que *rBA* y *rBa* son diferentes. En ese sentido bajo la hipótesis de segregación esperaríamos que $r_{Ba}>r_{BA}$.

#psiA, psiBA=psiBa, pA, pB, rA, rBA, rBa
m5<- occ2sps(Lyru, Urci, model = list(rA~ 1, rBA~ 1, rBa~ 1))
m5

#psiA, psiBA, psiBa, pA, pB, rA, rBA, rBa
m6<- occ2sps(Lyru, Urci, model = list(psiBA~1, rA~ 1, rBA~ 1, rBa~ 1))
m6

# Por el momento nos quedamos con los modelos construidos y vamos a proceder a escoger el mejor o los mejores modelos mediante el criterio de información de Akaike corregido para muestras pequeñas (AICc).

lista <- AICc(m1,m2,m3,m4,m5,m6)
row.names(lista)<- c("PsiA,PsiBA=PsiBa,pA=rA,pB=rB", "PsiA,PsiBA,PsiBa,pA=rA,pB=rB", "PsiA,PsiBA=PsiBa,pA,rA,pB,rB", "PsiA,PsiBA,PsiBa,pA,rA,pB,rB", "PsiA,PsiBA=PsiBa,pA,rA,pB,rBa,rBA", "PsiA,PsiBA,PsiBa,pA,rA,pB,rBa,rBA")

AICtable(lista, digits=3)

# Según la la regla de `Delta` <2 tanto `m3` como `m5` son igualmente plausibles. ¿Que tienen en común estos dos modelos?. Ambos modelos asumen que p es diferente de r  y según los valores estimados, la probabilidad de detectar una zorra gris cuando no es detectado el Lince (r_{Ba}=0.28) es mayor que cuando el lince es detectado en la misma ocasión de muestreo (r_{BA}=0.18).

# Veamooslo de manera gráfica -------------
best_mod1 <- data.frame(m3$real) # aquí solamente llamo los valores
best_mod2 <- data.frame(m5$real)

library(ggplot2)
library(gridExtra)

m3plot <- ggplot(best_mod1, aes(x = rownames(best_mod1), y = est)) +
  ylim(0,1) + labs( x= "Parámetros", y = "Probabilidad") +
  geom_point(data = best_mod1) + 
  geom_errorbar(data = best_mod1, aes(ymin = lowCI, ymax = uppCI), alpha = 0.8, width = 0.25, colour = "black") +
  theme_classic()

m5plot <- ggplot(best_mod2, aes(x = rownames(best_mod2), y = est))+ 
  ylim(0,1) + labs( x = "Parámetros", y = NULL) +
  geom_point(data = best_mod2) + 
  geom_errorbar(data = best_mod2, aes(ymin = lowCI, ymax = uppCI), alpha = 0.8, width = 0.25, colour = "black") + 
  theme_classic()

grid.arrange(m3plot,m5plot, ncol=2, top= "Valores estimados para los parámetros de los modelos m3 y m5")


# Calcular el indice o factor de interacción----------------
# El  FIE es la relación de cuanto las especies tienden a co-ocurrir en una unidad comparada con la hipótesis nula de que ambas ocurren de manera independiente. Valores de FIE<1 indican que las especies co-ocurren menos de lo esperado (exclusión o evasión), valores de FIE>1 sugieren que ambas especies tienden a co-ocurrir más de lo esperado, y si FIE=1 quiere decir que la ocupación de ambas especies es independiente o aleatoria.

#Parámetros reales
psiA <- m3$real[1,1]
psiBa <- m3$real[2,1]
psiBA <- m3$real[3,1]

#Factor de interacción de especies
FIE <- (psiA*psiBA)/(psiA*(psiA*psiBA+(1-psiA)*psiBa))
FIE

# También podemos calcular rho, el cual sigue el mismo principio que FIE pero para analizar las interacciones en la detección de las especies

rBa <- m5$real[7,1]
rBA <- m5$real[8,1]

rho <- (rBA/(1-rBA))/(rBa/(1-rBa))
rho

# ¿ Que nos indica el valor de Rho?
