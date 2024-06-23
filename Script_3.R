# Carga de Datos

rm(list=ls()) # Limpiar el entorno actual de trabajo 

library(readr) # Cargar la librer√≠a para leer ficheros CSV y TXT
library(dplyr) # Cargar la librer√≠a para manipulaci√≥n de datos

# Leer los nombres de columnas
nombres_columnas <- read_lines("column_names.txt")

View(nombres_columnas)

# Leer los datos del fichero gene_expression.csv
datos_gen <- read_csv2("gene_expression.csv", col_names = FALSE)

View(datos_gen)

# Leer los datos del fichero con el nombre de las muestras y el tipo de muestras 
datos_muestras <- read_csv2("classes.csv", col_names = FALSE)

View(datos_muestras)


# Asignar los nombres de columnas
colnames(datos_gen) <- nombres_columnas

View(datos_gen)

# Asignar nombres a columna X1 X2
colnames(datos_muestras) <- c("muestra", "clase")

View(datos_muestras)

# Combinar los dataframes datos_gen y datos_muestras
datos_gen$muestra <- datos_muestras$muestra
datos_gen$clase <- datos_muestras$clase

# Comprobar datos
head(datos_gen)

dim(datos_gen)

str(datos_gen)

View(datos_gen)

# Preparar Datos para Limpieza 

# Crear un dataframe para limpieza excluyendo las √∫ltimas dos columnas
datos_para_limpiar = select(datos_gen, -muestra, -clase)

View(datos_para_limpiar)

# Identificar qu√© columnas no son num√©ricas
no_numerico <- sapply(datos_para_limpiar, function(x) !is.numeric(x))

View(no_numerico)

# Imprimir los nombres de las columnas no num√©ricas
no_numerico <- names(datos_para_limpiar)[no_numerico]
print(no_numerico)

str(no_numerico)

# Verificar los datos almacenados como texto
# datos_para_limpiar$CFB

# Visualizar el tipo de cada una de las columnas no num√©ricas
tipos_no_numericos <- sapply(datos_para_limpiar[no_numerico], class)
print(tipos_no_numericos)

# Transformar n√∫meros almacenados como texto a datos num√©ricos
# asumiendo que se necesita revisar todas las columnas y convertir donde sea aplicable
datos_para_limpiar <- datos_para_limpiar %>%
  mutate(across(where(is.character), ~parse_number(.))) # Convierte texto a num√©ricos

#---> datos_para_limpiar  esta ya con todas columnas numericas

View(datos_para_limpiar)

str(datos_para_limpiar)

# Identificar valores nulos en el conjunto de datos
# Calcular el n√∫mero de valores NA por columna

#-----------------------#

any(is.na(datos_para_limpiar)) #--> con esto al dar false no es necesario el resto

na_por_columna <- sapply(datos_para_limpiar, function(x) sum(is.na(x)))

# Para obtener una visi√≥n m√°s general, podemos calcular el total de NAs en todo el conjunto de datos
total_nas <- sum(na_por_columna)
print(paste("Total de valores NA en el conjunto de datos:", total_nas))

# Para ver las filas que contienen al menos un NA, puedes usar:
filas_con_na <- which(rowSums(is.na(datos_para_limpiar)) > 0)
print(paste("Hay", length(filas_con_na), "filas con al menos un valor NA."))

#--------------------------#

# Identificar valores infinitos en el conjunto de datos
# Calcular el n√∫mero de valores infinitos por columna
infinitos_por_columna <- sapply(datos_para_limpiar, function(x) sum(is.infinite(x)))

# Para obtener una visi√≥n m√°s general, podemos calcular el total de infinitos en todo el conjunto de datos
total_infinitos <- sum(infinitos_por_columna)
print(paste("Total de valores infinitos en el conjunto de datos:", total_infinitos))

# Normalizar los datos (Datos limpios normalizados)
# Escalar cada columna para que tenga media ~0 y desviaci√≥n est√°ndar ~1
datos_normalizados <- as.data.frame(lapply(datos_para_limpiar, scale))

View(datos_normalizados)

anyNA(datos_normalizados)

#--------------------------------#
  ##################################################
 # MDS (No Supervisado - ReducciÛn de dimensiones)#
##################################################

library(stats) # Necesario mara usar MDS (funcion cmdscale)
library(ggplot2) # Gr·ficos
library(plotly) # Gr·ficos 3D

View(datos_para_limpiar) # Datos limpios 

dim(datos_para_limpiar)

View(datos_normalizados)

#datos_limpios <- select(datos_para_limpiar, -MIER3, -ZCCHC12, -RPL22L1)


# Calculo de matriz de distancias euclideas
#matriz_distancias <- dist(datos_para_limpiar, method = "euclidean")

matriz_distancias <- dist(datos_normalizados, method = "euclidean")

#datos_normalizados

# AplicaciÛn MDS de 2 componentes
mds_2D <- cmdscale(matriz_distancias, k = 2, eig = TRUE)


# Preparar dataframe para su visualizaciÛn 2D


datos_2D <- as.data.frame(mds_2D$points)

datos_2D$clase <- datos_muestras$clase

# Gr·fica 2D

ggplot(datos_2D, aes(x = V1, y = V2, color = clase)) +
  geom_point() +
  ggtitle("MDS de 2 Componentes") +
  xlab("Componenete 1") +
  ylab("Componenete 2")


# AplicaciÛn MDS de 3 componentes

mds_3D <- cmdscale(matriz_distancias, k = 3, eig = TRUE)



# Preparar dataframe para su visualizaciÛn 3D


datos_3D <- as.data.frame(mds_3D$points)

colnames(datos_3D) <- c("Componente1", "Componente2", "Componente3")

#colnames(datos_3D) <- c("Componente1", "Componente2", "Componente3")

datos_3D$clase <- datos_muestras$clase

# Gr·fica 3D

grafica3D <- plot_ly(datos_3D, 
                     x = ~Componente1, 
                     y = ~Componente2, 
                     z = ~Componente3, 
                     color = ~clase, 
                     colors = c("#8EF52D", "#F8AA18", "#F82218", "#00B4AB", "#581845"), 
                     type = 'scatter3d', 
                     mode = 'markers') %>%
  layout(scene = list(xaxis = list(title = 'Componente 1'),
                      yaxis = list(title = 'Componente 2'),
                      zaxis = list(title = 'Componente 3')),
         title = 'MDS - 3 Componentes')

# Mostrar el gr·fico
grafica3D
  
#---------------------------------#
  ########################################
 # K-Means (No Supervisado - Clustering)#
########################################

library(factoextra) # Permite la visualizaciÛn de los grupos creados

#datos_normalizados$clase <- datos_muestras$clase

View(datos_normalizados)

datos_normalizados_limpios <- select(datos_normalizados, -MIER3, -ZCCHC12, -RPL22L1)

# PreparaciÛn de datos para K-means

# Convertir clase tumores a factorial

datos_muestras$clase <- as.factor(datos_muestras$clase)

anyNA(datos_normalizados_limpios)

# Ver cuantos niveles hay

levels(datos_muestras$clase)


# Convertir Clase tumores a numÈrica
clase_numerica <- as.numeric(datos_muestras$clase)

View(clase_numerica)
str(clase_numerica)


# Unir columnas de Clase de tumores numÈricos a los datos normalizados
datos_kmeans <- cbind(datos_normalizados_limpios, clase_numerica)

View(datos_kmeans)

anyNA(datos_kmeans)


# C·lculo valor Ûptimo de grupos K (con curva del codo)

set.seed(123) # Semilla aleatorizaciÛn

k_optimo <- function(k) {
  kmeans(datos_kmeans, k, nstart = 15)$tot.withinss
}


# Suma valor k

k_suma <- 1:15

valor_k <- sapply(k_suma, k_optimo)

# Grafica curva del codo

elbow_plot <- qplot(k_suma, valor_k, geom = "line") +
  ggtitle("Curva del codo para eleccion K") +
  xlab("N˙mero de gupos K") +
  ylab("Total Suma Cuadrados")

print(elbow_plot)

# ComprobaciÛn k Ûptimo centrado en k=5 k=6 k=7 junto a la gr·fica

duda_k3 <- 3

duda_k4 <- 4

duda_k5 <- 5

duda_k6 <- 6

duda_k7 <- 7

resultado_kmeans <- kmeans(datos_kmeans, centers = duda_k5, nstart = 25)

# AÒadir el resultado de la agrupacion a los datos iniciales

datos_kmeans$agrupaciones <- as.factor(resultado_kmeans$cluster)

View(datos_kmeans)

#------------------#

# Para hacer las observaciones lo suyo es hacerlo con la reducciÛn de dimensiones previas de la PCA y 
# agrupar asi todos las columnas de tumores y visualizarlas con las agrupaciones--> datos_normalizados_limpios (luego aÒadirle clase_numerica para hacer cluster_plot) 

resultado_pca <- prcomp(datos_normalizados_limpios, center = TRUE)

print(summary(resultado_pca)) # se cogen 2 dimensiones PC1 y PC2

pca_datos <- data.frame(resultado_pca$x[, 1:2], agrupaciones = as.factor(resultado_kmeans$cluster))

#---------------------#


# Gr·fica ggplot

cluster_plot <- ggplot(pca_datos, aes(x = PC1, y = PC2, color = agrupaciones)) + 
  geom_point() + 
  ggtitle(paste("Grupos Identificados con K-Means (K =", duda_k5, ")")) + 
  xlab("Componente Principal 1 (PC1)") + 
  ylab("Componente Principal 2 (PC2)")

print(cluster_plot)

# Se comprueba con duda_k3, duda_k4, duda_k5, duda_k6, duda_k7 tambiÈn para ver cual se ve mejor --> K=5 3 agrupaciones claras 2 con ruido

  ######################################
 # DBSCAN (No supervisado- Clustering)#
######################################

library(dbscan) # Necesario para aplicar DBSCAN

# PreparaciÛn de datos

datos_kmeans2 <- cbind(datos_normalizados_limpios, clase_numerica)

# ConfiguraciÛn DBSCAN

eps_value <- 0.5 # Distancia m·xima de 2 puntos para considerarse vecinos

minPts_value <- 5 # N˙mero mÌnimo de puntos para considerarse unn grupo

# EjecuciÛn DBSCAN

dbscan_result <- dbscan(datos_kmeans2[, -ncol(datos_kmeans2)], eps = eps_value, minPts = minPts_value)

# AÒadir agrupaciones a datos iniciales

datos_kmeans$agrupaciones <- as.factor(dbscan_result$cluster)

# Gr·fica

cluster_plot2 <- ggplot(pca_datos, aes(x = PC1, y = PC2, color = agrupaciones)) + 
  geom_point() + 
  ggtitle(paste("Grupos Identificados con DBSCAN (eps =", eps_value, ", minpts =", minPts_value, ")")) + 
  xlab("Componente Principal 1 (PC1)") + 
  ylab("Componente Principal 2 (PC2)")


print(cluster_plot2)















