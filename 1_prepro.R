################ MALDI-TOF OTAMENDI ############################################
################ 1) PREPROCESAMIENTO  ##########################################
################ Primeras muestras  ############################################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library(binda)
library(here)
library(dplyr)
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)


### CARGA DE ESPECTROS #########################################################
################################################################################


# Creación de la ruta relativa de los archivos (Uno p/cada, corregir)
ruta_proyecto <- "C:/Users/urtea/OneDrive/Documents/Proyectos/MALDI_Otamendi/Datos/PBMC"
#ruta_proyecto <- "C:/Users/Facundo/Documents/Proyectos/MALDI_Vinchucas/Datos/PBMC"
ruta_datos <- file.path(ruta_proyecto)

# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)


### OBTENCIÓN DE METADATA DE ESPECTROS #########################################
################################################################################


# Creación de columnas vacías
col_grupo <- c() # Guardar grupo
col_ind <- c()
col_dia <- c() # Guardar dia de adquisición
col_rep_m <- c() # Guardar número de réplica bio
col_rep_t <- c() # Guardar número de réplica téc

# Patrones auxiliares para buscar el estado de la muestra
patron_ctl_sano <- "CTL_SANO"
patron_ctl_uti <- "CTL_UTI"
patron_sepsis <- "SEPSIS"
patron_shock <- "SHOCK"


# Ciclo que extrae dia, tipo, numero, well y réplica de cada muestra
for(i in 1:length(Spectra_list)) {
  nombre <- Spectra_list[[i]]@metaData$file
  nombre2 <- Spectra_list[[i]]@metaData[["name"]]
  
  if (grepl(patron_ctl_sano, nombre, ignore.case = TRUE)) {
    
    grupo <- patron_ctl_sano
    ind <- str_extract(nombre, "(?<=CTL_SANO\\\\)[^\\\\]+")
    dia <- "D1"
    
  }
  
  if (grepl(patron_ctl_uti, nombre, ignore.case = TRUE)) {
    
    grupo <- patron_ctl_uti
    
    # Extraer la parte después de "CTL_UTI\\"
    parte_relevante <- sub(".*CTL_UTI\\\\([^\\\\]+)\\\\.*", "\\1", nombre)
    
    # Dividir por "_"
    partes <- unlist(strsplit(parte_relevante, "_"))
    
    # Asignar valores a variables
    ind <- partes[1]
    dia <- partes[2]
    
  }
  
  if (grepl(patron_sepsis, nombre, ignore.case = TRUE)) {
    print(nombre)
    print(nombre2)
    grupo <- patron_sepsis
    
    # Extraer la parte después de "SEPSIS\\"
    parte_relevante <- sub(".*SEPSIS\\\\([^\\\\]+)\\\\.*", "\\1", nombre)
    
    # Dividir por "_"
    partes <- unlist(strsplit(parte_relevante, "_"))
    
    # Asignar valores a variables
    ind <- partes[1]
    dia <- partes[2]
  
  }
  
  if (grepl(patron_shock, nombre, ignore.case = TRUE)) {
    print(nombre)
    print(nombre2)
    grupo <- patron_shock
    
    # Extraer la parte después de "SEPSIS\\"
    parte_relevante <- sub(".*SHOCK\\\\([^\\\\]+)\\\\.*", "\\1", nombre)
    
    # Dividir por "_"
    partes <- unlist(strsplit(parte_relevante, "_"))
    
    # Asignar valores a variables
    ind <- partes[1]
    dia <- partes[2]
    
  }
  
  extraer_rep <- function(nombre) {
    pos_punto <- regexpr("\\.", nombre)  # Encuentra la posición del "."
    aux_t <- substr(nombre, pos_punto - 1, pos_punto - 1)  # Caracter anterior al "."
    aux_m <- substr(nombre, pos_punto - 3, pos_punto - 3)  # Tercer caracter antes del "."
    
    return(list(aux_m = as.numeric(aux_m), aux_t = as.numeric(aux_t)))
  }
  
  resultado <- extraer_rep(nombre2)
  
  # Asignar los valores a nuevas variables
  rep_m <- resultado$aux_m
  rep_t <- resultado$aux_t
    
  # Almacena los valores extraídos en sus respectivas columnas
  col_grupo <- c(col_grupo, grupo) # Guardar grupo
  col_ind <- c(col_ind, ind)
  col_dia <- c(col_dia, dia) # Guardar dia de adquisición
  col_rep_m <- c(col_rep_m, rep_m) # Guardar número de réplica bio
  col_rep_t <- c(col_rep_t, rep_t) # Guardar número de réplica téc

  # print(paste("Estado:", estado, "| Número:", numero_vinchuca, "| Sexo:", sexo,
  #             "| Réplica Muestra:", replica_muestra, "| Réplica Técnica:", replica_tecnica))

  # Guardo metadata en el espectro
  Spectra_list[[i]]@metaData$grupo <- grupo
  Spectra_list[[i]]@metaData$individuo <- ind
  Spectra_list[[i]]@metaData$dia <- dia
  Spectra_list[[i]]@metaData$rep_m <- rep_m
  Spectra_list[[i]]@metaData$rep_t <- rep_t

  # Data Frame con los datos limpios
  df_metadata <- data.frame(grupo = col_grupo, individuo = col_ind, dia = col_dia,
                            rep_m = col_rep_m, rep_t = col_rep_t)

  # Creación de factores de agrupamiento para su uso posterior
  df_metadata$factor_grupo_dia <- paste0(df_metadata$grupo, "_", df_metadata$dia)
  # Creación de factores de agrupamiento para su uso posterior
  df_metadata$paciente_rep_t <- paste0(df_metadata$individuo, "_", 
                                       df_metadata$dia, "_", df_metadata$rep_m, "_", df_metadata$rep_t)
  df_metadata$paciente_rep_m <- paste0(df_metadata$individuo, "_", 
                                       df_metadata$dia, "_", df_metadata$rep_m)
  df_metadata$factor_rep <- paste0(df_metadata$individuo, "_", 
                                       df_metadata$dia, "_", df_metadata$rep_m)
  df_metadata$factor_dia <- paste0(df_metadata$individuo, "_", 
                                   df_metadata$dia)

}


### CONTROL DE CALIDAD Y LIMPIEZA DE ESPECTROS #################################
################################################################################


# Screening inicial: Detección de espectros de baja calidad
sc.results <- screenSpectra(Spectra_list, meta = df_metadata)
summary(sc.results)
plot(sc.results, labels = TRUE)

# plot(Spectra_list[[253]]) # Ploteo de espectros ruidosos (no hay)

# Descartamos espectros defectuosos
Spectra_list_f1 <- sc.results$fspectra # Filtramos espectros
df_metadata_f1 <- sc.results$fmeta # Filtramos metadatos


### FILTRADO Y TRANSFORMACIÓN DE ESPECTROS #####################################
################################################################################

# Parámetros de procesamiento de espectros
thScale <- 10 # Smoothing
ite <- 105 # Baseline correction
SigNoi <- 2.5 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.03 # Peak binning

# Transformación/filtrado/corrección de espectros con parámetros definidos
# 1) Transformación de intensidad por medio de función sqrt
Spectra_list_f1 <- transfIntensity(Spectra_list_f1, fun = sqrt)
plot(Spectra_list_f1[[30]])
# 2) Suavizado del espectro mediante el método Wavelet
Spectra_list_f1 <- wavSmoothing(Spectra_list_f1, method = "Wavelet", n.levels = 4)
plot(Spectra_list_f1[[30]])
# Detección de la linea de base
baseline <- estimateBaseline(Spectra_list_f1[[30]], method = "SNIP",
                             iterations = ite)
plot(Spectra_list_f1[[30]])
lines(baseline, col="red", lwd=2)
# 3) Remoción de linea de base mediante método SNIP
Spectra_list_f2 <- removeBaseline(Spectra_list_f1, method = "SNIP",
                                  iterations = ite)
plot(Spectra_list_f2[[30]])
# 4) Calibración de intensidad mediante método PQN
Spectra_list_f2 <- calibrateIntensity(Spectra_list_f2, method = "PQN")
plot(Spectra_list_f2[[30]])
# 5) Alineación de espectros
Spectra_list_f3 <- alignSpectra(Spectra_list_f2,
                                halfWindowSize=20,
                                SNR=2,
                                tolerance=0.02, # Parámetro sensible
                                warpingMethod="lowess")
plot(Spectra_list_f3[[30]])


### PROMEDIO DE LECTURAS DE UNA MISMA RÉPLICA TÉCNICA #################
################################################################################


# Promedio de lecturas de una misma well
Spectra_list_prom_rep <- averageMassSpectra(Spectra_list_f3,
                                            labels = factor(df_metadata_f1$factor_rep),
                                            method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_rep <- df_metadata_f1 %>%
  distinct(df_metadata_f1$factor_rep, .keep_all = TRUE)

# Promedio de wells de una misma muestra
Spectra_list_prom_muestra <- averageMassSpectra(Spectra_list_prom_rep,
                                                labels = factor(df_metadata_prom_rep$factor_dia),
                                                method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_mue <- df_metadata_prom_rep %>%
  distinct(df_metadata_prom_rep$factor_dia, .keep_all = TRUE)


### EXTRACCIÓN DE PICOS Y ALINEACIÓN ###########################################
################################################################################

# A partir de acá probamos trabajar con Spectra_list_prom_rep

# Análisis de la SNR en espectros para chequear que utilizamos el valor correcto
noise <- estimateNoise(Spectra_list_prom_rep[[20]])
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue") # Se ve que es correcto el 2

# Detección de picos a partir del umbral definido de SNR
peaks <- detectPeaks(Spectra_list_prom_rep,
                     SNR = SigNoi,
                     halfWindowSize = 40)

# Ploteo de picos detectados en un espectro de ejemplo
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
points(peaks[[20]], col="red", pch=4)

# Alineado de picos
peaks <- alignPeaks(peaks,
                    minFreq = 0.8,
                    tolerance = tol)

#summaryPeaks(peaks[1:10])  # resumen estadistico de picos (primeros 10)

# Conteo de picos por perfil
cP <- countPeaks(peaks)

# Gráfico de picos
plot(cP, type = "n")
text(cP, label = 1:length(cP))


# Patrones de picos
peakPatterns(peaks)

# Filtrado de picos de baja frecuencia de aparición
picos_filtrados <- filterPeaks(peaks,
                               minFreq = 0.25,
                               labels = df_metadata_prom_rep$grupo ) #labels

# Patrones de picos
peakPatterns(picos_filtrados)

# Conteo de picos por perfil
cP2 <- countPeaks(picos_filtrados)

# Gráfico
plot(cP2, type = "n")
text(cP2, label = 1:length(cP2))

# Fusión de picos de la misma muestra
picos_fusion_muestra <- mergeMassPeaks(picos_filtrados,
                                       labels = df_metadata_prom_rep$factor_dia,
                                       method = "median")

# Patrones de picos
peakPatterns(picos_fusion_muestra)


### CREACIÓN DE MATRIZ DE INTENSIDADES Y DICOTÓMICA ############################
################################################################################


# Matriz de intensidades 19 individuos
matint_28_ind <- intensityMatrix(picos_fusion_muestra,
                                 Spectra_list_prom_muestra) # sin valores NA

# Matriz de intensidades de 80 muestras
# matint_na_51 <- intensityMatrix(picos_fusion_muestra) # con valores NA
matint_82_mue <- intensityMatrix(picos_filtrados,
                                 Spectra_list_prom_rep) # sin valores NA

# Definición de umbrales
thr1 <- optimizeThreshold(matint_28_ind,
                          df_metadata_prom_mue$grupo,
                          verbose = T)
thr2 <- optimizeThreshold(matint_82_mue,
                          df_metadata_prom_rep$grupo,
                          verbose = T)

# Dicotomización
matint_28_ind_dico <- dichotomize(matint_28_ind, thr1)
matint_82_mue_dico <- dichotomize(matint_82_mue, thr2)

# Agrego nombres a las filas de cada df
rownames(matint_28_ind_dico) <- df_metadata_prom_mue$factor_dia
rownames(matint_82_mue_dico) <- df_metadata_prom_rep$factor_rep
rownames(matint_28_ind) <- df_metadata_prom_mue$factor_dia
rownames(matint_82_mue) <- df_metadata_prom_rep$factor_rep


### GUARDAR DATOS ##############################################################
################################################################################

# Establecer el directorio de trabajo en la ubicación del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Se puede mejorar, subcarpeta

# Guardar matrices y metadata asociada como archivo .Rdata
save(matint_28_ind_dico,df_metadata_prom_mue, file = "matint_28_ind_dico.Rdata")
save(matint_82_mue_dico, df_metadata_prom_rep, file = "matint_82_mue_dico.Rdata")
save(matint_28_ind,df_metadata_prom_mue, file = "matint_28_ind.Rdata")
save(matint_82_mue, df_metadata_prom_rep, file = "matint_82_mue.Rdata")

# Guardar matrices y metadata asociada como archivo .csv
write.csv(matint_28_ind_dico, "matint_28_ind_dico.csv", row.names = TRUE)
write.csv(matint_28_ind, "matint_28_ind.csv", row.names = TRUE)
write.csv(matint_82_mue_dico, "matint_82_mue_dico.csv", row.names = TRUE)
write.csv(matint_82_mue, "matint_82_mue.csv", row.names = TRUE)
write.csv(df_metadata_prom_mue, "df_28.csv", row.names = TRUE)
write.csv(df_metadata_prom_rep, "df_82.csv", row.names = TRUE)

# #
# #
# #
# ### FIN ########################################################################
# ################################################################################
# 
# ##### ANEXO: GRAFICO DE PICOS PREPONDERANTES ###################################
# 
# # Cargar librería necesaria
# library(scales)  # Para agregar transparencia a los colores
# 
# # Definir los puntos de referencia y tolerancia
# highlight_positions <- c(2152, 3466, 5443, 8491, 6283)
# tolerance <- 20
# 
# # Definir los límites del eje x
# x_lim <- c(2000, 10000)
# 
# # Definir los colores con transparencia (50% de opacidad)
# colors <- c(alpha("blue", 0.3), alpha("blue", 0.3), alpha("red", 0.3), alpha("red", 0.3), alpha("red", 0.3))
# 
# # Ajustar la cantidad de espectros por gráfico
# spectra_per_plot <- 3
# total_spectra <- length(Spectra_list_prom_muestra)
# 
# # Calcular el número de gráficos necesarios
# num_plots <- ceiling(total_spectra / spectra_per_plot)
# 
# # Loop para generar los gráficos
# for (plot_idx in 1:num_plots) {
#   
#   # Definir los índices de los espectros que van en este gráfico
#   start_idx <- (plot_idx - 1) * spectra_per_plot + 1
#   end_idx <- min(plot_idx * spectra_per_plot, total_spectra)
#   
#   # Configurar un layout de 3 filas y 1 columna (3 subplots por gráfico)
#   par(mfrow=c(spectra_per_plot, 1), mar=c(4, 4, 2, 2))  # Márgenes ajustados
#   
#   # Graficar cada espectro en el rango de este gráfico
#   for (i in start_idx:end_idx) {
#     mass <- Spectra_list_prom_muestra[[i]]@mass
#     intensity <- Spectra_list_prom_muestra[[i]]@intensity
#     
#     # Graficar el espectro individual
#     plot(mass, intensity, type = "l", col = "black",
#          xlab = "Mass/Charge (m/z)", ylab = "Intensity",
#          main = paste("Espectro del individuo", i, "con regiones de interés"),
#          xlim = x_lim)
#     
#     # Dibujar las barras semitransparentes con tolerancia
#     for (j in 1:length(highlight_positions)) {
#       rect(highlight_positions[j] - tolerance, par("usr")[3],
#            highlight_positions[j] + tolerance, par("usr")[4],
#            col = colors[j], border = NA)
#     }
#     
#     # Dibujar el espectro encima para que la barra quede de fondo
#     lines(mass, intensity, col = "black")
#   }
#   
# }