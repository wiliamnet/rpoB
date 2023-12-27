library(Biostrings)
library(DECIPHER)
library(ggplot2)


# Directorio donde se encuentran los archivos FNA
directorio <- "E:/fna"

# Lista para almacenar las secuencias
secuencias <- DNAStringSet()

# Leer los archivos FNA y almacenar las secuencias en la lista
fna_files <- list.files(path = directorio, pattern = "\\.fna$", full.names = TRUE)
for (archivo in fna_files) {
  secuencia <- readDNAStringSet(archivo)
  secuencias <- c(secuencias, secuencia)
}

# Objetivo 01 Alineación de secuencias, 
alineacion <- AlignSeqs(secuencias)

# Imprimir el resultado
print(alineacion)


# Cargar la biblioteca DECIPHER
library(DECIPHER)

# Visualizar la alineación

BrowseSeqs(alineacion)



##graficar un arbol filogeetico
## Carga la biblioteca 'ape'
library(ape)

# Supongamos que 'alineacion' es tu objeto de alineación obtenido previamente

# Convierte la alineación a un objeto DNAbin
alineacion_dna <- as.DNAbin(alineacion)

# Calcula la matriz de distancias
dist_matrix <- dist.dna(alineacion_dna)

# Continuar con la construcción del árbol filogenético...


# Supongamos que 'alineacion' es tu objeto de alineación obtenido previamente

# Calcula la matriz de distancias
dist_matrix <- dist.dna(alineacion_dna)


# Construye el árbol utilizando Neighbor-Joining
arbol <- nj(dist_matrix)

# Visualiza el árbol
plot(arbol, main = "Árbol Filogenético para patogenos", cex = 0.6)


##### MEJORA




# Verifica si el objeto es de clase DNAbin
if ("DNAbin" %in% class(secuencias_dna)) {
  # Calcula la matriz de distancias
  matriz_distancias <- dist.dna(secuencias_dna)
  
  # Construye el árbol utilizando Neighbor-Joining
  arbol <- nj(matriz_distancias)
  
  # Extrae solo el nombre del microorganismo (suponiendo que las etiquetas están en un formato específico)
  etiquetas <- gsub("^(\\w+)_.*", "\\1", dimnames(matriz_distancias)[[1]])
  
  # Visualiza el árbol con parámetros personalizados
  plot(
    arbol,
    main = "Árbol Filogenético para patógenos",
    cex = 0.6,  # Tamaño de los nodos
    col = "blue",  # Color de las ramas
    lwd = 2,  # Grosor de las ramas
    tip.label = etiquetas  # Etiquetas de las puntas del árbol
  )
} else {
  warning("La conversión a DNAbin no fue exitosa.")
}

# Convierte la alineación a un objeto DNAbin
secuencias_dna <- as.DNAbin(alineacion)

# Verifica si el objeto es de clase DNAbin
if ("DNAbin" %in% class(secuencias_dna)) {
  # Calcula la matriz de distancias
  matriz_distancias <- dist.dna(secuencias_dna)
  
  # Construye el árbol utilizando Neighbor-Joining
  arbol <- nj(matriz_distancias)
  
  # Visualiza el árbol con parámetros personalizados
  plot(
    arbol,
    main = "Árbol Filogenético para patógenos",
    cex = 0.6,  # Tamaño de los nodos
    col = "blue",  # Color de las ramas
    lwd = 2  # Grosor de las ramas
  )
} else {
  warning("La conversión a DNAbin no fue exitosa.")
}


##################encontrar mutaciones # Muestra los microorganismos con mutaciones en la posicion 526

if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 526:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
} else {
  print("No se encontraron mutaciones en la posición 526.")
}

# Identifica microorganismos con mutaciones comparando con la primera secuencia
microorganismos_con_mutaciones <- which(as.character(nucleotidos_en_posicion) != as.character(nucleotidos_en_posicion[1]))

# Muestra los microorganismos con mutaciones y los nucleótidos modificados
if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 526:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
  
  # Obtén los nucleótidos modificados
  nucleotidos_modificados <- nucleotidos_en_posicion[microorganismos_con_mutaciones]
  print(paste("Nucleótidos modificados:", paste(nucleotidos_modificados, collapse = ", ")))
} else {
  print("No se encontraron mutaciones en la posición 526.")
}




# Encuentra la base en la posición específica (526) para todas las secuencias
nucleotidos_en_posicion <- subseq(alineacion, start = 526, end = 526)

# Identifica microorganismos con mutaciones comparando con la primera secuencia
microorganismos_con_mutaciones <- which(as.character(nucleotidos_en_posicion) != as.character(nucleotidos_en_posicion[1]))

# Muestra los microorganismos con mutaciones y los nucleótidos modificados
if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 526:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
  
  # Itera sobre los microorganismos con mutaciones para obtener detalles de la sustitución
  for (microorganismo in microorganismos_con_mutaciones) {
    nucleotido_original <- as.character(nucleotidos_en_posicion[microorganismo])
    nucleotido_referencia <- as.character(nucleotidos_en_posicion[1])
    
    if (nucleotido_original != nucleotido_referencia) {
      print(paste("Microorganismo", microorganismo, ":",
                  "Nucleótido original:", nucleotido_original,
                  "Nucleótido de referencia:", nucleotido_referencia))
    }
  }
} else {
  print("No se encontraron mutaciones en la posición 526.")
}



# Encuentra la base en la posición específica (533) para todas las secuencias
nucleotidos_en_posicion <- subseq(alineacion, start = 533, end = 533)

# Identifica microorganismos con mutaciones comparando con la primera secuencia
microorganismos_con_mutaciones <- which(as.character(nucleotidos_en_posicion) != as.character(nucleotidos_en_posicion[1]))

# Muestra los microorganismos con mutaciones y los nucleótidos modificados
if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 533:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
  
  # Itera sobre los microorganismos con mutaciones para obtener detalles de la sustitución
  for (microorganismo in microorganismos_con_mutaciones) {
    nucleotido_original <- as.character(nucleotidos_en_posicion[microorganismo])
    nucleotido_referencia <- as.character(nucleotidos_en_posicion[1])
    
    if (nucleotido_original != nucleotido_referencia) {
      print(paste("Microorganismo", microorganismo, ":",
                  "Nucleótido original:", nucleotido_original,
                  "Nucleótido de referencia:", nucleotido_referencia))
    }
  }
} else {
  print("No se encontraron mutaciones en la posición 533.")
}



# Encuentra la base en la posición específica (550) para todas las secuencias
nucleotidos_en_posicion <- subseq(alineacion, start = 550, end = 550)

# Identifica microorganismos con mutaciones comparando con la primera secuencia
microorganismos_con_mutaciones <- which(as.character(nucleotidos_en_posicion) != as.character(nucleotidos_en_posicion[1]))

# Muestra los microorganismos con mutaciones y los nucleótidos modificados
if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 550:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
  
  # Itera sobre los microorganismos con mutaciones para obtener detalles de la sustitución
  for (microorganismo in microorganismos_con_mutaciones) {
    nucleotido_original <- as.character(nucleotidos_en_posicion[microorganismo])
    nucleotido_referencia <- as.character(nucleotidos_en_posicion[1])
    
    if (nucleotido_original != nucleotido_referencia) {
      print(paste("Microorganismo", microorganismo, ":",
                  "Nucleótido original:", nucleotido_original,
                  "Nucleótido de referencia:", nucleotido_referencia))
    }
  }
} else {
  print("No se encontraron mutaciones en la posición 550.")
}


# Encuentra la base en la posición específica (572) para todas las secuencias
nucleotidos_en_posicion <- subseq(alineacion, start = 572, end = 572)

# Identifica microorganismos con mutaciones comparando con la primera secuencia
microorganismos_con_mutaciones <- which(as.character(nucleotidos_en_posicion) != as.character(nucleotidos_en_posicion[1]))

# Muestra los microorganismos con mutaciones y los nucleótidos modificados
if (length(microorganismos_con_mutaciones) > 0) {
  print(paste("Microorganismos con mutaciones en la posición 572:", 
              paste(microorganismos_con_mutaciones, collapse = ", ")))
  
  # Itera sobre los microorganismos con mutaciones para obtener detalles de la sustitución
  for (microorganismo in microorganismos_con_mutaciones) {
    nucleotido_original <- as.character(nucleotidos_en_posicion[microorganismo])
    nucleotido_referencia <- as.character(nucleotidos_en_posicion[1])
    
    if (nucleotido_original != nucleotido_referencia) {
      print(paste("Microorganismo", microorganismo, ":",
                  "Nucleótido original:", nucleotido_original,
                  "Nucleótido de referencia:", nucleotido_referencia))
    }
  }
} else {
  print("No se encontraron mutaciones en la posición 572.")
}

