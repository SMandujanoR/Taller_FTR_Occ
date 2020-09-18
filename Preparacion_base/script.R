#Script para crear la base de 0 y 1 que se utilizará en los modelos de ocupación.

#Primero debemos habilitar el paquete camtrapR con la siguiente función

library(camtrapR)

#Después debemos leer la base de datos con la actividad de nuestras cámaras y la de registros de especies. 

actividad_cam <- read.csv("actividad de estaciones.csv", header = T)
dataframe_registros <- read.csv("registros_especies.csv", header = T)

#El siguiente paso es obtener una base de datos con los días que opero cada cámara de la siguiente forma 

actividad_cam_2 <- cameraOperation(CTtable = actividad_cam,
                   stationCol = "Station",
                   setupCol = "Fecha_colocacion",
                   retrievalCol = "Fecha_retiro",
                   hasProblems = TRUE,
                   dateFormat = "%d/%m/%Y")

View(actividad_cam_2)

#Por último utilizaremos el objeto actividad_cam_2 y dataframe_registros para obtener la base de datos de 0 y 1 con la siguiente función:

str(dataframe_registros)

hist_captura_conejo <- detectionHistory(recordTable = dataframe_registros, 
camOp = actividad_cam_2,
stationCol = "Station", 
speciesCol = "Species", 
recordDateTimeCol = "DateTimeOriginal", 
recordDateTimeFormat = "%Y-%m-%d", 
species = "Sylvilagus floridanus", 
occasionLength = 1, 
day1 = "station",
includeEffort = FALSE,
timeZone = "America/Mexico_City")

View(hist_captura_conejo)

#Ahora podemos guardar nuestra base en un archivo csv

write.csv(hist_captura_conejo, "historia captura conejo.csv")

#En occasionLength se debe especificar el número de días de cada ocasión de detección, es decir "1" dará la detección del venado en cada día de actividad de la cámara y "10" agrupará las detecciones de venado cada diez días

#Hagamos el ejemplo

hist_captura_conejo10 <- detectionHistory(recordTable = dataframe_registros, 
camOp = actividad_cam_2,
#output = "binary",
stationCol = "Station", 
speciesCol = "Species", 
recordDateTimeCol = "DateTimeOriginal", 
recordDateTimeFormat = "%Y-%m-%d", 
species = "Sylvilagus floridanus", 
occasionLength = 10, 
day1 = "station",
includeEffort = TRUE,
scaleEffort = TRUE, 
timeZone = "America/Mexico_City",
#unmarkedMultFrameInput = TRUE
)

#Vemos que contiene nuestro objeto
str (hist_captura_conejo10)

#Ahora extraemos cada uno de los dataframe

detecciones <- hist_captura_conejo10 [[1]]

duración <- hist_captura_conejo10 [[2]]

#Por último los guardamos en archivo csv

write.csv(detecciones, "detecciones_conejo.csv")
write.csv(duracion, "duracion_estandarizada.csv")

#----------------------------------------------------

#Fin del script 




