# Modelo de ocupación condicional RW para interacción de dos especies: `wiqid` y `RPresence`

library(wiqid)
library(RPresence)

### *Gabriel Andrade-Ponce* {-}

### Cargar las historias de detección

#Datos para RPresence
#Datos de las dos especies
twosphist <- (read.csv("data/lynx_foxRW.csv"))[,-1]
#Colapso
# R presence requiere un objeto donde  por cada sitio y muestreo: 0 = a ninguna especie detectada, 1= especie A detectada, 2=especie B detectada y 3= ambas especies detectadas.
nsites <- nrow(twosphist)/2
nsrvys <- ncol(twosphist)
cov1 = cov2 = NULL
datos = twosphist[1:nsites,] + 2*twosphist[nsites+1:nsites,]

## Construcción de los modelos

### Generar el dataframe de RPresence
pao <- createPao(datos, unitcov = NULL, survcov = NULL, title = "Lynx-fox")
summary(pao)

### Ajustar los modelos

#### Modelo RW en `RPresence`

#Ajustar primer modelo, identico a m1 de wiqid
m1P <- occMod(model = list(psi~SP, p~SP), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA=PsiBa, pA=rA, pB=rB")

#Matriz de diseño de m1P
m1P$dmat$psi

#valores beta de m1P
m1P$beta$psi

# Valores reales de m1P para PsiA
head(m1P$real$psiA)

#Valor de parámetro derivado nu
head(m1P$derived$nu)

#Ajustar el resto de modelos identicos a los creados en wiqid
m2P <- occMod(model = list(psi~SP + INT, p~SP), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA, PsiBa, pA=rA, pB=rB")

m3P <- occMod(model = list(psi~SP, p~SP + INT_o + SP:INT_o), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA=PsiBa, pA, rA, pB, rB")

m4P <- occMod(model = list(psi~SP + INT, p~SP+ INT_o + SP:INT_o), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA, PsiBa, pA, rA, pB, rB")

m5P <- occMod(model = list(psi~SP, p~SP +  INT_d + INT_o+SP:INT_o), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA=PsiBa, pA, rA, pB, rBa, rBA")

m6P <- occMod(model = list(psi ~ SP + INT, p ~ SP + INT_d + INT_o + SP:INT_o), data = pao, type = "so.2sp.1", param = "PsiBA", modname = "PsiA, PsiBA, PsiBa, pA, rA, pB, rBa, rBA")

#Matriz de diseño de m2P para los valores de ocupación
m2P$dmat$psi


#Tabla de AICc en RPresence
models <- list(m1P,m2P,m3P,m4P,m5P,m6P)
results <- createAicTable(models,use.aicc = TRUE )
summary(results)

#Valor derivado nu en RPresence
m5P$derived$nu[1,1]

#Valor derivado rho en Rpresence
m5P$derived$rho[1,1]

#Valores beta de la ocupación en el modelo m3P
m3P$beta$psi

#Valores beta de la detección en el modelo m3P
m3P$beta$p

