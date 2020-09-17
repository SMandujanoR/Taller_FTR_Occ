# Modelo de ocupación para la co-ocurrencia de dos o más especies----------
## Gabriel Andrade-Ponce

# Cargar las librerias
require(ggplot2) # Para las gráficas
require(unmarked) # Para correr los modelos

# Cargar las historias de detección
# Es importante que todas tengan el mismo número de ocasiones y sitios de muestreo
#Datos de zorra gris:
Urci <- read.csv("data/Fox.hist5.csv")
# Datos de lince:
Lyru <- read.csv("data/Lynx.hist5.csv")
# Datos de conejos:
Syfl <- read.csv("data/Rabbit.hist5.csv")

# Ahora creamos una lista que contenga todas las historias

y <- list(as.matrix(Urci),
          as.matrix(Lyru),
          as.matrix(Syfl))
# Indica los nombres de las especies:
names(y) <- c("zorra", "lince", "conejo")

# Crear el objeto de trabajo en unmarked
int <- unmarkedFrameOccuMulti(y = y, # lista de historias de detección
      siteCovs = NULL, # covariables de sitio
      obsCovs = NULL) # covariables de observación
summary(int)
plot(int)

# Un objeto importante es el fDesing, que es la matriz de los parámetros y nos va a servir para generar las formulas del modelo
int@fDesign

# Es necesario crear las formulas para cada proceso, en el proceso ecológico cada formula "~" corresponderá a las columnas del fDesing
formupsi <- c("~1","~1","~1","~0","~0","~0","~0")

# en el proceso observacional cada formula corresponde a cada una de las especies
formup <- c("~1", "~1", "~1")
# En este caso vamos a correr el modelo más sencillo (sin interacciones)
summary(m1 <- occuMulti(formup, formupsi, data = int))

#  Ahora ajustemos modelos con "interacciones"

formupsi2 <- c("~1","~1","~1","~1","~0","~0","~0")
formup2 <- c("~1","~1","~1")
summary(m2 <- occuMulti(formup2, formupsi2, data = int)) 

formupsi3 <- c("~1","~1","~1","~0","~1","~0","~0")
formup3 <- c("~1","~1","~1")
summary(m3 <- occuMulti(formup3, formupsi3, data = int)) 

formupsi4 <- c("~1","~1","~1","~0","~0","~1","~0")
formup4 <- c("~1","~1","~1")
summary(m4 <- occuMulti(formup4, formupsi4, data = int))

formupsi5 <- c("~1","~1","~1","~1","~1","~1","~0")
formup5 <- c("~1","~1","~1")
summary(m5 <- occuMulti(formup5, formupsi5, data = int))

# Escoger un mejor modelo mediante AIC--------------------

# Construimos una lista con todos los modelos candidatos
fl <- fitList("independencia" = m1,
              "Fox+Lynx" = m2,
              "Fox+Rabbit" = m3,
              "Lynx+Rabbit" = m4,
              "interacción pareada" = m5)

modSel(fl) # Función de selección

# Explorar el mejor modelo-------
# Vamos a construir objetos con los valores de probabilidad de ocupación condicional de la zorra gris, para ello vamos a utilizar la función predict de unmarked

# Zorra | conejo presente:
foxplusrabbit <- predict(m5, 'state', species = 'zorra', cond = 'conejo')[1,] 

# Zorra | conejo ausente
foxminusrabbit <- predict(m5, 'state', species = 'zorra', cond = '-conejo')[1,] 

# Zorra | lince presente
foxpluslynx <- predict(m5, 'state', species = 'zorra', cond = 'lince')[1,] 

# Zorra | lince ausente
foxminuslynx <- predict(m5, 'state', species = 'zorra', cond = '-lince')[1,] 

# Ahora unimos todas los objetos en una sola base de datos
predall2 <- rbind(foxplusrabbit, foxminusrabbit, foxpluslynx, foxminuslynx)
interactionfox <- c("conejo", "-conejo", "lince", "-lince")
predall2 <- data.frame(interactionfox, predall2)

# Finalmente podemos crear un gráfico con los valores de predicción

ocufoxinteractplot <- ggplot(predall2, aes(x = predall2$interactionfox, y = Predicted)) + ylim(0,1) + labs(main = "Probabilidad de ocupación condicional de zorra gris", x = "Presencia y ausencia de especies", y = paste(c("zorra gris", expression(psi)))) + geom_point(data = predall2) + geom_errorbar(data = predall2, aes(ymin = lower, ymax = upper), alpha = 0.8, width = 0.25, colour = "black") + theme_classic()
ocufoxinteractplot

# Que tan significativa es la estimación de los parámetros
confint(m5, type = "state")
