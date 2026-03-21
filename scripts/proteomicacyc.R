#Instalación de proteoDA

install.packages("devtools")
devtools::install_github("ByrumLab/proteoDA", 
                         dependencies = TRUE, 
                         build_vignettes = TRUE)
library(proteoDA)
library(limma)
library(dplyr)
library(ggplot2)
library(here)

## ── Project root & analysis constants ────────────────────────────────────────
proj_dir     <- here::here()
FC_THRESH    <- 1     # log2 fold-change threshold for volcano plot coloring
LOG10_EPS    <- 1e-10 # epsilon to avoid -Inf in -log10(0)

#CARGAR ARCHIVOS NECESARIOS

#Matriz de expresión o abundancia
input_data <- read.csv2(file.path(proj_dir, "proteomicacyc.csv"), header = TRUE, row.names = 1)

#Anotación de las proteinas
annotation_data <- read.csv2(file.path(proj_dir, "annotation.csv"))
dim(annotation_data)

#Metadatos del experimento
sample_metadata <- read.csv2(file.path(proj_dir, "metadata.csv"), row.names = 1)
dim(input_data)
str(input_data)
input_data[] <- lapply(input_data, function(x) as.numeric(as.character(x)))

#Creación del objeto DA para llevar a cabo el análisis
raw <- DAList(data = input_data,
              annotation = annotation_data,
              metadata = sample_metadata,
              design = NULL,
              eBayes_fit = NULL,
              results = NULL,
              tags = NULL)

# Filter out unneeded samples and proteins with too much missing data
filtered <- raw |>
  filter_samples(group != "Pool") |>
  zero_to_missing() |>
  filter_proteins_by_proportion(min_prop = 0.33,
                                grouping_column = "group")

# Make the normalization report
write_norm_report(filtered,
                  grouping_column = "group")

# Normalize
normalized <- normalize_data(filtered, 
                             norm_method = "cycloess")

# Make the quality control report
write_qc_report(normalized,
                color_column = "group", filename = "group.pdf", overwrite = TRUE)

write_qc_report(normalized,
                color_column = "vph", filename = "vph.pdf", overwrite = TRUE)

write_qc_report(normalized,
                color_column = "condition", filename = "condition.pdf", overwrite = TRUE)

# Turn metadata column into a factor with desired levels
normalized$metadata$condition <- factor(normalized$metadata$condition, 
                                    levels = c("NORMAL", "TUMORAL"))
normalized$metadata$vph <- factor(normalized$metadata$vph, 
                                    levels = c("NEGATIVE", "POSITIVE"))
normalized$metadata$group <- factor(normalized$metadata$group, 
                                  levels = c("NORMAL_NEGATIVE", "TUMORAL_NEGATIVE", "NORMAL_POSITIVE", "TUMORAL_POSITIVE"))

# Add a statistical design, fit the model, and extract results
no_intercept <- add_design(normalized,
                        design_formula = ~0 + group)

no_intercept$design$design_matrix

# Otra manera de generar la matriz de contrastes (no correr)
#contrasts <- makeContrasts(
  #TVsS = (TUMORAL_NEGATIVE + TUMORAL_POSITIVE) - (NORMAL_NEGATIVE + NORMAL_POSITIVE),
  #TpVsTn = TUMORAL_POSITIVE - TUMORAL_NEGATIVE,
  #TpVsSp = TUMORAL_POSITIVE - NORMAL_POSITIVE,
  #levels = no_intercept$design$design_matrix)

#contrasts
#write.csv(contrasts, file = "contrasts.csv", row.names = FALSE)

# Adición de contrastes usando proteoDA

no_intercept <- add_contrasts(no_intercept, contrasts_vector = c("TVsS=(TUMORAL_NEGATIVE + TUMORAL_POSITIVE) - (NORMAL_NEGATIVE + NORMAL_POSITIVE)",
                                                                 "TpVsTn=TUMORAL_POSITIVE - TUMORAL_NEGATIVE",
                                                                 "TpVsSp = TUMORAL_POSITIVE - NORMAL_POSITIVE"))
no_intercept$design$contrast_matrix


final <- no_intercept |>
  fit_limma_model() |>
  extract_DA_results(pval_thresh = 0.05, 
                     lfc_thresh = 1,
                     adj_method = "BH",
                     extract_intercept = F)


# Export results
write_limma_tables(final, overwrite = TRUE)
write_limma_plots(final,
                  grouping_column = "group", overwrite = TRUE)


#___________________________________________________________________________________________

#Make volcano plot of DEGs result
#Determine point colors based on significance and sign of the logFC and the Padj.
res.df <- final$results$TVsS

library(org.Hs.eg.db)

res.df$gene_symbol <- mapIds(org.Hs.eg.db, key = row.names(res.df), keytype = "UNIPROT", column = "SYMBOL")

res.df <- res.df %>%
  mutate(point_color = case_when(
    adj.P.Val < 0.05 & logFC < -FC_THRESH ~ "down", # significantly down
    adj.P.Val < 0.05 & logFC > FC_THRESH  ~ "up",   # significantly up
    TRUE ~ "NS") # not significant
  )

#Create the Volcano plot
volcano_plot <- ggplot(data = res.df, aes(x = logFC, y = -log10(adj.P.Val), col = point_color)) +
  geom_hline(yintercept = -log10(0.05), col = "green4", linetype = 'dashed')+
  geom_point(color = ifelse(res.df$logFC < -1 & res.df$adj.P.Val<0.05, "blue", ifelse(res.df$logFC > 1 & res.df$adj.P.Val<0.05, "red", "grey")), size = 1.5) +
  geom_text(aes(label = ifelse(res.df$logFC < -2 & -log10(adj.P.Val) >2.45 | res.df$logFC > 2 & -log10(adj.P.Val) >2.45, res.df$gene_symbol, ""), size = 3),
            vjust = -0.5, hjust = 0, size = 3)+
  scale_color_manual(values = c("blue", "grey", "red"),
                     labels = c("Reprimido", "No significativo", "Sobre expresado"))+
  coord_cartesian(ylim = c(0, 3.8), xlim = c(-11, 11))+
  labs(color = '', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"Adjp-value"), size = 5)+
  ggtitle('Tumoral vs Sano') +
  guides(shape = guide_legend(override.aes = list(color = NULL, size = 6)))+
  scale_shape_manual(values = c(16, 17, 18))+
  theme(axis.line = element_line(),panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

print(volcano_plot)
ggsave("proteomicacyc_TvsS.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

#________________________________________________________________________________________

df <- read.csv2("david_enrichment.csv", header = TRUE)

df_filtered <- df[df$Benjamini < 0.05, ]
df_filtered$Benjamini <- as.numeric(df_filtered$Benjamini)
df_filtered$Count <- as.numeric(df_filtered$Count)

# Crea el gráfico
p <- ggplot(df_filtered, aes(x=Benjamini, y=reorder(Term, -Benjamini), size=Count, color=Benjamini)) +
  geom_point(stat="identity") +
  labs(x="Benjamini", y="Termino", title="Enriquecimiento en DAVID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),  # Aumenta el tamaño de la letra de los ejes
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=18),  # Aumenta el tamaño del título del eje x
        axis.title.y = element_text(size=18),  # Aumenta el tamaño del título del eje y
        legend.text = element_text(size=15),  # Aumenta el tamaño de la letra de las leyendas
        plot.title = element_text(hjust = 0.5, size=18)) +  # Centra el título y aumenta su tamaño
  scale_size_continuous(range = c(1, 10)) +  # Ajusta el tamaño de los puntos
  scale_color_gradient(low = "red", high = "blue")  # Ajusta el degradado de color
# Muestra el gráfico
print(p)

ggsave("Enrichment_todos.png", plot = p, width = 16, height = 14, dpi = 300)
