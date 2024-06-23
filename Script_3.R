# Cargar datos

rm(list=ls()) # Limpiar el entorno actual de trabajo 

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

# Realizar PCA y reducir a 3 dimensiones
pca <- prcomp(caracteristicas, center = TRUE, scale. = TRUE)
datos_pca_reducidos <- as.data.frame(pca$x[, 1:3])
datos_pca_reducidos$clase <- respuesta

# Dividir el dataset en conjuntos de entrenamiento y prueba

set.seed(123)  # Para reproducibilidad

# Crear partición de datos
indice_entrenamiento <- createDataPartition(datos_pca$clase, p = 0.7, list = FALSE)
datos_entrenamiento <- datos_pca_reducidos[indice_entrenamiento, ]
datos_prueba <- datos_pca_reducidos[-indice_entrenamiento, ]

# Función para entrenar y evaluar un modelo
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
