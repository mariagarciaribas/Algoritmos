
# Actividad Grupal - Algor?tmos e Inteligencia Artificial

# -- ?scar N??ez L?pez

# -- Mar?a Garc?a Ribas

###################################################
# Procesamiento de datos #
###################################################

rm(list=ls()) # Limpiar el entorno actual de trabajo 

# Cargar las librerías necesarias

library(readr) # Cargar la librería para leer ficheros CSV y TXT
library(dplyr) # Cargar la librería para manipulación de datos
library(e1071)  # Cargar la librería para SVM, SVM Gaussiano y Naive Bayes
library(caret)  # Cargar la librería para la división de los datos y evaluación del modelo
library(randomForest)  # Cargar la librería para Random Forest
library(dplyr)  # Cargar la librería para la manipulación de datos

# Leer los nombres de columnas
nombres_columnas <- read_lines("column_names.txt")

# Leer los datos del fichero gene_expression.csv
datos_gen <- read_csv2("gene_expression.csv", col_names = FALSE)

# Asignar los nombres de columnas
colnames(datos_gen) <- nombres_columnas

# Verificar datos
head(datos_gen)

# Leer los datos del fichero con el nombre de las muestras y el tipo de muestras 
datos_muestras <- read_csv2("classes.csv", col_names = FALSE)

# Asignar nombres de columna
colnames(datos_muestras) <- c("muestra", "clase")

# Combinar los dataframes datos_gen y datos_muestras
datos_gen$muestra <- datos_muestras$muestra
datos_gen$clase <- datos_muestras$clase

# Comprobar datos
head(datos_gen)

# Preparar Datos para Limpieza 

# Crear un dataframe para limpieza excluyendo las últimas dos columnas
datos_para_limpiar = select(datos_gen, -muestra, -clase)

# Transformar los valores almacenados como texto a datos numéricos
datos_para_limpiar <- datos_para_limpiar %>%
  mutate(across(where(is.character), ~parse_number(.))) 

# Contar valores NA en cada columna
colSums(is.na(datos_para_limpiar))

# Eliminar filas con valores NA
datos_limpios <- datos_para_limpiar[complete.cases(datos_para_limpiar),]

# Verificar nuevamente si hay valores NA
any(is.na(datos_limpios))

# Identificación de Datos Infinitos
infinitos_por_columna <- sapply(datos_para_limpiar, function(x) sum(is.infinite(x)))

print(infinitos_por_columna)

# Normalización de Datos
datos_normalizados <- as.data.frame(lapply(datos_para_limpiar, scale))

# Añadir de nuevo las columnas "muestra" y "clase"
datos_normalizados$muestra <- datos_gen$muestra
datos_normalizados$clase <- datos_gen$clase

# Eliminar las columnas con errores 
datos_para_pca <- select(datos_normalizados, -MIER3, -ZCCHC12, -RPL22L1)

# Eliminar la columna "muestra"
datos_pca <- datos_para_pca[ , !(names(datos_para_pca) %in% c("muestra"))]

datos_pca$clase <- as.factor(ifelse(datos_pca$clase == "CGC", 1, 
                                    ifelse(datos_pca$clase == "CHC", 2, 
                                           ifelse(datos_pca$clase == "CFB", 3, 
                                                  ifelse(datos_pca$clase == "AGH", 4, 5)))))

# Separar las características y la variable objetivo
caracteristicas <- datos_pca[ , !(names(datos_pca) %in% c("clase"))]
respuesta <- datos_pca$clase

  ###################################################
 # PCA (No Supervisado - Reducci?n de dimensiones) #
###################################################


# Realizar PCA y reducir a 3 dimensiones

pca <- prcomp(caracteristicas, center = TRUE, scale. = TRUE)
datos_pca_reducidos <- as.data.frame(pca$x[, 1:3])
datos_pca_reducidos$clase <- respuesta

# Mostrar la varianza explicada por los componentes principales

varianza_explicada <- summary(pca)$importance[2, 1:3]
print(paste("Varianza explicada por los 3 primeros componentes:", round(varianza_explicada * 100, 2), "%"))


  ##################################################
 # MDS (No Supervisado - Reducci?n de dimensiones)#
##################################################


View(datos_para_limpiar) # Datos limpios 

dim(datos_para_limpiar)

View(datos_normalizados)

# Calculo de matriz de distancias euclideas

matriz_distancias <- dist(datos_normalizados, method = "euclidean")

# Aplicaci?n MDS de 2 componentes

mds_2D <- cmdscale(matriz_distancias, k = 2, eig = TRUE)

# Preparar dataframe para su visualizaci?n 2D

datos_2D <- as.data.frame(mds_2D$points)

datos_2D$clase <- datos_muestras$clase

# Gr?fica 2D

ggplot(datos_2D, aes(x = V1, y = V2, color = clase)) +
  geom_point() +
  ggtitle("MDS de 2 Componentes") +
  xlab("Componenete 1") +
  ylab("Componenete 2")


# Aplicaci?n MDS de 3 componentes

mds_3D <- cmdscale(matriz_distancias, k = 3, eig = TRUE)

# Preparar dataframe para su visualizaci?n 3D

datos_3D <- as.data.frame(mds_3D$points)

colnames(datos_3D) <- c("Componente1", "Componente2", "Componente3")

datos_3D$clase <- datos_muestras$clase

# Gr?fica 3D

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

# Mostrar el gr?fico

grafica3D


  ########################################
 # K-Means (No Supervisado - Clustering)#
########################################


View(datos_normalizados)

datos_normalizados_limpios <- select(datos_normalizados, -MIER3, -ZCCHC12, -RPL22L1)

# Preparaci?n de datos para K-means

# Convertir clase tumores a factorial

datos_muestras$clase <- as.factor(datos_muestras$clase)

anyNA(datos_normalizados_limpios)

# Ver cuantos niveles hay

levels(datos_muestras$clase)

# Convertir Clase tumores a num?rica

clase_numerica <- as.numeric(datos_muestras$clase)

View(clase_numerica)
str(clase_numerica)

# Unir columnas de Clase de tumores num?ricos a los datos normalizados

datos_kmeans <- cbind(datos_normalizados_limpios, clase_numerica)

View(datos_kmeans)

anyNA(datos_kmeans)

# Eliminar columnas no num?ricas

datos_kmeans <- select(datos_kmeans, -clase, -muestra)

is.numeric(datos_kmeans)

View(datos_kmeans)

# C?lculo valor ?ptimo de grupos K (con curva del codo)

set.seed(123) # Semilla aleatorizaci?n

k_optimo <- function(k) {
  kmeans(datos_kmeans, k, nstart = 15)$tot.withinss
}


# Suma valor k

k_suma <- 1:15

valor_k <- sapply(k_suma, k_optimo)

# Gr?fica curva del codo

elbow_plot <- qplot(k_suma, valor_k, geom = "line") +
  ggtitle("Curva del codo para eleccion K") +
  xlab("N?mero de gupos K") +
  ylab("Total Suma Cuadrados")

print(elbow_plot)

# Comprobaci?n k ?ptimo centrado en k=5 k=6 k=7 junto a la gr?fica (tambi?n se extiende comprobaci?n)

duda_k3 <- 3

duda_k4 <- 4

duda_k5 <- 5

duda_k6 <- 6

duda_k7 <- 7

resultado_kmeans <- kmeans(datos_kmeans, centers = duda_k5, nstart = 25)

# A?adir el resultado de la agrupacion a los datos iniciales

datos_kmeans$agrupaciones <- as.factor(resultado_kmeans$cluster)

View(datos_kmeans)

# Para hacer las observaciones lo suyo es hacerlo con una reducci?n de dimensiones con PCA  
# de las columnas de tumores y visualizarlas con las agrupaciones 

datos_normalizados_limpios2 <- select(datos_normalizados_limpios, -muestra, -clase)

is.numeric(datos_normalizados_limpios2)

resultado_pca <- prcomp(datos_normalizados_limpios2, center = TRUE)

print(summary(resultado_pca)) # se cogen 2 dimensiones PC1 y PC2

pca_datos <- data.frame(resultado_pca$x[, 1:2], agrupaciones = as.factor(resultado_kmeans$cluster))


# Etiquetar los niveles de los grupos

etiquetas <- c("HPB", "CHC", "CFB", "AHG", "CGC")

levels(pca_datos$agrupaciones) <- etiquetas

# Gr?fica ggplot

cluster_plot <- ggplot(pca_datos, aes(x = PC1, y = PC2, color = agrupaciones)) + 
  geom_point() + 
  ggtitle(paste("Grupos Identificados con K-Means (K =", duda_k5, ")")) + 
  xlab("Componente Principal 1 (PC1)") + 
  ylab("Componente Principal 2 (PC2)")

print(cluster_plot)

# Se comprueba con duda_k3, duda_k4, duda_k5, duda_k6, duda_k7 tambi?n para ver cual se adapta mejor  
# Para K=5 3 agrupaciones claras (AHG, HPB, CGC) 2 con ruido (CHC, CFB)


  ######################################
 # DBSCAN (No supervisado- Clustering)#
######################################


# Preparaci?n de datos

datos_kmeans2 <- cbind(datos_normalizados_limpios2, clase_numerica)

# Configuraci?n DBSCAN

eps_value <- 0.5 # Distancia m?xima de 2 puntos para considerarse vecinos

minPts_value <- 5 # N?mero m?nimo de puntos para considerarse unn grupo

# Ejecuci?n DBSCAN

dbscan_result <- dbscan(datos_kmeans2[, -ncol(datos_kmeans2)], eps = eps_value, minPts = minPts_value)

# A?adir agrupaciones a datos iniciales

datos_kmeans$agrupaciones <- as.factor(dbscan_result$cluster)

# Gr?fica

cluster_plot2 <- ggplot(pca_datos, aes(x = PC1, y = PC2, color = agrupaciones)) + 
  geom_point() + 
  ggtitle(paste("Grupos Identificados con DBSCAN (eps =", eps_value, ", minpts =", minPts_value, ")")) + 
  xlab("Componente Principal 1 (PC1)") + 
  ylab("Componente Principal 2 (PC2)")


print(cluster_plot2)


  ########################
 # M?TODOS SUPERVISADOS #
########################


# Dividir el dataset en conjuntos de entrenamiento y prueba

set.seed(123)  # Para reproducibilidad

# Crear partici?n de datos

indice_entrenamiento <- createDataPartition(datos_pca$clase, p = 0.7, list = FALSE)
datos_entrenamiento <- datos_pca_reducidos[indice_entrenamiento, ]
datos_prueba <- datos_pca_reducidos[-indice_entrenamiento, ]

# Funci?n para entrenar y evaluar un modelo

evaluar_modelo <- function(funcion_modelo, nombre_modelo, datos_entrenamiento, datos_prueba) {
  modelo <- funcion_modelo(datos_entrenamiento)
  predicciones <- predict(modelo, datos_prueba)
  matriz_confusion <- confusionMatrix(as.factor(predicciones), as.factor(datos_prueba$clase))
  
  valor_kappa <- matriz_confusion$overall["Kappa"]
  precision <- matriz_confusion$overall["Accuracy"]
  
  cat("Resultados para", nombre_modelo, ":\n")
  print(matriz_confusion)
  cat("Kappa:", valor_kappa, "\n")
  cat("Precision:", precision, "\n\n")
  
  return(list(matriz_confusion = matriz_confusion, kappa = valor_kappa, precision = precision))
}

# Definir las funciones de los modelos

modelo_svm <- function(datos) {
  svm(clase ~ ., data = datos, kernel = "linear")
}

modelo_svm_gaussiano <- function(datos) {
  svm(clase ~ ., data = datos, kernel = "radial")
}

modelo_random_forest <- function(datos) {
  randomForest(clase ~ ., data = datos, importance = TRUE)
}

modelo_naive_bayes <- function(datos) {
  naiveBayes(clase ~ ., data = datos)
}

# Evaluar los modelos

resultados_svm <- evaluar_modelo(modelo_svm, "SVM", datos_entrenamiento, datos_prueba)
resultados_svm_gaussiano <- evaluar_modelo(modelo_svm_gaussiano, "SVM Gaussiano", datos_entrenamiento, datos_prueba)
resultados_random_forest <- evaluar_modelo(modelo_random_forest, "Random Forest", datos_entrenamiento, datos_prueba)
resultados_naive_bayes <- evaluar_modelo(modelo_naive_bayes, "Naive Bayes", datos_entrenamiento, datos_prueba)

# Comparar los resultados

resultados <- data.frame(
  Modelo = c("SVM", "SVM Gaussiano", "Random Forest", "Naive Bayes"),
  Precision = c(resultados_svm$precision, resultados_svm_gaussiano$precision, resultados_random_forest$precision, resultados_naive_bayes$precision),
  Kappa = c(resultados_svm$kappa, resultados_svm_gaussiano$kappa, resultados_random_forest$kappa, resultados_naive_bayes$kappa)
)

print(resultados)


