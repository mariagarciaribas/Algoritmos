#################
# Carga de Datos#
#################

rm(list=ls()) #Limpiar el entorno actual de trabajo 

# Cargar librerías necesarias
library(readr) # Para leer ficheros CSV y TXT
library(dplyr) # Librería para manipulación de datos

# Leer los nombres de columnas
column_names <- read_lines("column_names.txt")

# Leer datos desde el fichero gene_expression.csv
gene_data <- read_csv2("gene_expression.csv", col_names = FALSE)

# Asignar los nombres de columnas
colnames(gene_data) <- column_names

# Verificar datos
head(gene_data)

# Anadir columnas con el número de muestra (sample) y la clase (class)

# Leer el fichero con las columnas "sample" y "class"
sample_class_data <- read_csv2("classes.csv", col_names = FALSE)

# Asignar nombres de columna
colnames(sample_class_data) <- c("sample", "class")

# Combinar los dataframes gene_data y sample_class_data
gene_data$sample <- sample_class_data$sample
gene_data$class <- sample_class_data$class

# Comprobar datos
head(gene_data)

################################
# Preparar Datos para Limpieza #
################################

# Crear una dataframe para limpieza excluyendo las últimas dos columnas
data_to_clean = select(gene_data, -sample, -class)

# Verificar las dimensiones del conjunto de datos
print(paste("El conjunto de datos tiene", nrow(data_to_clean), "filas y", ncol(data_to_clean), "columnas."))

##########################################
# Identificación de Valores no Numéricos #
##########################################

# Identificar qué columnas no son numéricas
# Usamos sapply para aplicar una función a cada columna y luego identificamos las que no son numéricas
non_numeric_columns <- sapply(data_to_clean, function(x) !is.numeric(x))

# Imprimir los nombres de las columnas no numéricas
non_numeric_names <- names(data_to_clean)[non_numeric_columns]
print(non_numeric_names)

# Verificar los datos almacenados como texto
# data_to_clean$CFB

# Visualizar el tipo de cada una de las columnas no numéricas
non_numeric_types <- sapply(data_to_clean[non_numeric_columns], class)
print(non_numeric_types)

# Transformar números almacenados como texto a datos numéricos
# asumiendo que se necesita revisar todas las columnas y convertir donde sea aplicable
data_to_clean <- data_to_clean %>%
  mutate(across(where(is.character), ~parse_number(.))) # Convierte texto a numéricos

# Verificar los cambios
non_numeric_columns <- sapply(data_to_clean, function(x) !is.numeric(x))
non_numeric_names <- names(data_to_clean)[non_numeric_columns]
print(non_numeric_names)

#################################
# Identificación de Datos Nulos #
#################################

# Identificar valores nulos en el conjunto de datos
# Calcular el número de valores NA por columna
na_count_per_column <- sapply(data_to_clean, function(x) sum(is.na(x)))

# Para obtener una visión más general, podemos calcular el total de NAs en todo el conjunto de datos
total_nas <- sum(na_count_per_column)
print(paste("Total de valores NA en el conjunto de datos:", total_nas))

# Para ver las filas que contienen al menos un NA, puedes usar:
rows_with_na <- which(rowSums(is.na(data_to_clean)) > 0)
print(paste("Hay", length(rows_with_na), "filas con al menos un valor NA."))

#####################################
# Identificación de Datos Infinitos #
#####################################

# Identificar valores nulos en el conjunto de datos
# Calcular el número de valores NA por columna
inf_count_per_column <- sapply(data_to_clean, function(x) sum(is.infinite(x)))

# Para obtener una visión más general, podemos calcular el total de NAs en todo el conjunto de datos
total_inf <- sum(inf_count_per_column)
print(paste("Total de valores infinitos en el conjunto de datos:", total_inf))

##########################
# Normalización de Datos #
##########################

# Normalizar los datos
# Escalar cada columna para que tenga media ~0 y desviación estándar ~1
normalized_data <- as.data.frame(lapply(data_to_clean, scale))

# Añadir de nuevo las columnas "sample" y "class"
normalized_data$sample <- gene_data$sample
normalized_data$class <- gene_data$class

data_for_pca <- select(normalized_data, -MIER3, -ZCCHC12, -RPL22L1)




# Cargar las librerías necesarias
library(e1071)  # Para SVM, SVM Gaussiano y Naive Bayes
library(caret)  # Para la división de los datos y evaluación del modelo
library(randomForest)  # Para Random Forest
library(dplyr)  # Para la manipulación de datos


# Eliminar la columna "Sample"
data_pca <- data_for_pca[ , !(names(data_for_pca) %in% c("sample"))]


data_pca$class <- as.factor(ifelse(data_pca$class == "CGC", 1, 
                                   ifelse(data_pca$class == "CHC", 2, 
                                          ifelse(data_pca$class == "CFB", 3, 
                                                 ifelse(data_pca$class == "AGH", 4, 5)))))


# Separar las características y la variable objetivo
features <- data[ , !(names(data) %in% c("class"))]


response <- data$class

# Realizar PCA y reducir a 3 dimensiones
pca <- prcomp(features, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:3])
pca_data$class <- response

pca_data$class <- as.factor(pca_data$class)
data <- pca_data %>%
  mutate(class_num = as.numeric(factor(class, levels = c("CGC", "CHC", "CFB", "AGH", "HPB"))))

# Mostrar la varianza explicada por los componentes principales
explained_variance <- summary(pca)$importance[2, 1:3]
print(paste("Varianza explicada por los 3 primeros componentes:", round(explained_variance * 100, 2), "%"))

# Dividir el dataset en conjuntos de entrenamiento y prueba
set.seed(123)  # Para reproducibilidad
trainIndex <- createDataPartition(data$class, p = 0.7, list = FALSE)
data_train <- pca_data[trainIndex, ]
data_test <- pca_data[-trainIndex, ]

# Función para entrenar y evaluar un modelo
evaluate_model <- function(model_func, model_name, data_train, data_test) {
  model <- model_func(data_train)
  predictions <- predict(model, data_test)
  confusion_matrix <- confusionMatrix(as.factor(predictions), as.factor(data_test$class))
  
  kappa_value <- confusion_matrix$overall["Kappa"]
  accuracy <- confusion_matrix$overall["Accuracy"]
  
  cat("Resultados para", model_name, ":\n")
  print(confusion_matrix)
  cat("Kappa:", kappa_value, "\n")
  cat("Accuracy:", accuracy, "\n\n")
  
  return(list(confusion_matrix = confusion_matrix, kappa = kappa_value, accuracy = accuracy))
}


# Definir las funciones de los modelos
svm_model <- function(data) {
  svm(class ~ ., data = data, kernel = "linear")
}

svm_gaussian_model <- function(data) {
  svm(class ~ ., data = data, kernel = "radial")
}

random_forest_model <- function(data) {
  randomForest(class~ ., data = data, importance = TRUE)
}

naive_bayes_model <- function(data) {
  naiveBayes(class ~ ., data = data)
}

# Evaluar los modelos
results_svm <- evaluate_model(svm_model, "SVM", data_train, data_test)
results_svm_gaussian <- evaluate_model(svm_gaussian_model, "SVM Gaussiano", data_train, data_test)
results_random_forest <- evaluate_model(random_forest_model, "Random Forest", data_train, data_test)
results_naive_bayes <- evaluate_model(naive_bayes_model, "Naive Bayes", data_train, data_test)

# Comparar los resultados
results <- data.frame(
  Model = c("SVM", "SVM Gaussiano", "Random Forest", "Naive Bayes"),
  Accuracy = c(results_svm$accuracy, results_svm_gaussian$accuracy, results_random_forest$accuracy, results_naive_bayes$accuracy),
  Kappa = c(results_svm$kappa, results_svm_gaussian$kappa, results_random_forest$kappa, results_naive_bayes$kappa)
)

print(results)
