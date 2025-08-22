library(dplyr)
library(viridis)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(Rtsne)
library(hrbrthemes)
library(RColorBrewer)

data <- read.csv("finaltable.tsv", sep="\t", header = TRUE)

metadata <- read.csv("metadata_geoLocal.csv", header = TRUE)
metadata <- na.omit(metadata)

country <- metadata[,c(1,16,17)]

###########################################################################

family <- data %>% filter(Rank == "F")
#family <- left_join(family, country, by=c("sample_name"="Run"))
pcaFamily <- family[,c(1,2,7)]

class <- data %>% filter(Rank == "C")
#family <- left_join(family, country, by=c("sample_name"="Run"))
pcaClass <- class[,c(1,2,7)]

phylum <- data %>% filter(Rank == "P")
#phylum <- left_join(phylum, country, by=c("sample_name"="Run"))
pcaPhylum <- phylum[,c(1,2,7)]

specie <- data %>% filter(Rank == "S")
specie <- specie %>% filter(scientific_name != "Homosapiens")
#specie <- left_join(specie, country, by=c("sample_name"="Run"))
pcaSpecie <- specie[,c(1,2,7)]

order <- data %>% filter(Rank == "O")
#order <- left_join(order, country, by=c("sample_name"="Run"))
pcaOrder <- order[,c(1,2,7)]

genera <- data %>% filter(Rank == "G")
#genera <- left_join(genera, country, by=c("sample_name"="Run"))
pcaGenera <- genera[,c(1,2,7)]

kingdom <- data %>% filter(Rank == "K")
#kingdom <- left_join(kingdom, country, by=c("sample_name"="Run"))
pcaKingdom <- kingdom[,c(1,2,7)]



#### Diversity plot ###############

bioValues <- species %>%
  group_by(scientific_name, sample_name) #%>%
  #summarise(mean_percentage = mean(Percentage, na.rm = TRUE))

bioValues <- bioValues %>% filter(Percentage >= 0.05)
#bioValues <- bioValues %>% filter(scientific_name != "Chordata")
#bioValues <- bioValues %>% filter(scientific_name != "Proteobacteria")
bioValues <- bioValues %>% filter(scientific_name != "Chordata")
#bioValues <- bioValues %>% filter(scientific_name != "Escherechia")


annot <- read.csv("annot_metadata.tsv",sep="\t")
wild <- annot %>% filter(Local == "wild")
greenhouse <- annot %>% filter(Local == "greenhouse")
my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(29)
selvagem <- bioValues %>% filter(sample_name %in% wild$Run)


ggplot(bioValues %>% filter(sample_name %in% greenhouse$Run), aes(x = reorder(sample_name,scientific_name), y = Percentage, fill = scientific_name)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  coord_flip() +
  #labs(title = "Genus distribution by Library") +
  theme_ipsum() +
  scale_fill_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 2, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )


### selvagens e casa de vegetação

bioValues <- bioValues %>%
  mutate(Group = ifelse(sample_name %in% greenhouse$Run, "Greenhouse", "Wild"))

# Definir paleta de cores
#my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(bioValues$scientific_name)))


ggplot(bioValues, aes(x = reorder(sample_name, scientific_name), y = Percentage, fill = scientific_name)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +  # Divide os gráficos por grupo
  theme_ipsum() +
  scale_fill_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )




################################
pcaData <- pcaSpecie %>%
  pivot_wider(names_from = scientific_name, values_from = Percentage, values_fill = 0)

sample_names <- pcaData$sample_name
dados_numericos <- pcaData %>% select(-sample_name)
rownames(dados_numericos) <- sample_names

set.seed(123)  # Para reprodutibilidade
tsne_results <- Rtsne(dados_numericos, dims = 2, perplexity = 24, verbose = TRUE, max_iter = 1000)


# Criar um data frame com os resultados do t-SNE
tsne_df <- data.frame(
  X = tsne_results$Y[,1],
  Y = tsne_results$Y[,2],
  Category = as.factor(row.names(dados_numericos)) # Ajuste para uma coluna categórica, se aplicável
)

rownames(tsne_df) <- tsne_df$Category
df <- left_join(tsne_df, metadata, by = c("Category"="Run"))
rownames(df) <- df$Category
df <- na.omit(df)
#df <- df %>% filter(Country != "unknow")

# Plotar os resultados do t-SNE

my_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#17298a","#090a1321")

tsne_plot <- ggplot(df, aes(x = X, y = Y,color = Continent)) +
  geom_point(size = 4) +
  theme_ipsum() +
  scale_color_manual(values = my_palette) +
  labs(title = "Metagenome Species Scale of T. cacao", x = "Dim 1", y = "Dim 2") 

tsne_plot #+ BioProjecttsne_plot #+ theme(legend.Continenttsne_plot #+ BioProjecttsne_plot #+ theme(legend.position="none")



##### Complete map plot #####http://127.0.0.1:31861/graphics/plot_zoom_png?width=579&hehttp://127.0.0.1:31861/graphics/plot_zoom_png?width=463&height=349ight=436
library(ggiraph) # install.packages('ggiraph')http://127.0.0.1http://127.0.0.1:31861/graphics/plot_zoom_png?width=579&height=436:31861/graphics/plot_zoom_png?width=579&height=436
library(ggplot2) # install.packages('ggplot2')
library(dplyr) # install.packages('dplyr')
library(patchwork) # install.packages('patchwork')
library(tidyr) # install.packages('tidyr')
library(sf) # install.packages('sf')
set.seed(123)

# mapa mundi
world_df <- read_sf("https://raw.githubusercontent.com/holtzy/R-graph-gallery/master/DATA/world.geojson")
world_df <- world_df %>%
  filter(!name %in% c("Antarctica", "Greenland"))

# number of species per country
species_value <- as.data.frame(table(metadata$Country))
species_value <- species_value %>% filter(Var1 != "unknow")
species_value <- species_value[-1,]
df <- left_join(df, species_value, by=c("Country"="Var1"))

df <- left_join(world_df, df, by = c("name" = "Country"))

p1 <- ggplot(df, aes(
  x = X, y = Y,
  tooltip = name,
  data_id = name,
  color = name
)) +
  geom_point_interactive(size = 4) +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )



# Create the second chart (Bar plot)
p2 <- ggplot(species_value, aes(
  x = reorder(Var1,Freq),
  y = Freq,
  tooltip = Var1,
  data_id = Var1,
  fill = Var1
)) +
  geom_col_interactive() + #data = filter(df, !is.na(Freq))) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )
# Create the third chart (choropleth)
p3 <- ggplot() +
  geom_sf(data = df, fill = "lightgrey", color = "lightgrey") +
  geom_sf_interactive(
    data = filter(df, !is.na(Category)) ,
                  aes(geometry = geometry, fill = name, tooltip = name, data_id = name)) +
  coord_sf(crs = st_crs(3857)) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

# Combine the plots
combined_plot <- (p1 + p2) / p3 + plot_layout(heights = c(1, 2))


# Create the interactive plot
interactive_plot <- girafe(ggobj = combined_plot)

interactive_plot <- girafe_options(
  interactive_plot,
  opts_hover(css = "fill:red;stroke:black;")
)

# save as an html widget
htmltools::save_html(interactive_plot, "multiple-ggiraph-2.html")

pdf("static_plot2.pdf", width = 14, height = 8)
print(combined_plot)
dev.off()




#### Exclusives ####


# Filtrar apenas os gêneros (Rank == "G")
df_genera <- data %>% filter(Rank == "S")

# Contar quantas amostras possuem cada gênero
genera_counts <- df_genera %>%
  group_by(scientific_name) %>%
  summarise(n_samples = n_distinct(sample_name), .groups = "drop")

# Identificar os gêneros exclusivos (que aparecem em apenas uma amostra)
exclusive_genera <- genera_counts %>%
  filter(n_samples == 1) %>%
  pull(scientific_name)

# Filtrar o dataframe original para obter os gêneros exclusivos por amostra
df_exclusive_genera <- df_genera %>%
  filter(scientific_name %in% exclusive_genera)

# Exibir resultado
df_exclusive_genera





#### Perfis distintos ####

library(dplyr)
library(tidyr)
library(vegan)

# Filtrar apenas níveis taxonômicos desejados (por exemplo, Rank == "G" para gêneros)
df_genera <- data %>% filter(Rank == "C")
df_genera <- df_genera %>% filter(scientific_name != "Mammalia")

# Criar uma matriz de abundância com gêneros como colunas e amostras como linhas
abundance_matrix <- df_genera %>%
  select(sample_name, scientific_name, Percentage) %>%
  pivot_wider(names_from = scientific_name, values_from = Percentage, values_fill = 0)


## Para espécies - Itens com nomes repetidos
abundance_matrix <- df_genera %>%
  group_by(sample_name, scientific_name) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = "drop") %>%  # Somando os valores duplicados
  pivot_wider(names_from = scientific_name, values_from = Percentage, values_fill = list(Percentage = 0))  # Corrigindo `values_fill`



# Remover a coluna de sample_name para calcular distâncias
rownames(abundance_matrix) <- abundance_matrix$sample_name
annot <- abundance_matrix$sample_name
abundance_matrix <- abundance_matrix %>% select(-sample_name)
rownames(abundance_matrix) <- annot

# Calcular matriz de dissimilaridade de Bray-Curtis
bray_curtis_dist <- vegdist(abundance_matrix, method = "bray")

# Identificar as amostras mais distintas (com maior distância média das outras)
dist_summary <- apply(as.matrix(bray_curtis_dist), 1, mean)
distinct_samples <- names(sort(dist_summary, decreasing = TRUE)[1:5])  # Top 5 mais distintos

# Exibir amostras mais distintas
distinct_samples


#### PCA sem escala #####

# Aplicar PCA na matriz de abundância
abundance_matrix <- abundance_matrix[, apply(abundance_matrix, 2, var) > 0]
rownames(abundance_matrix) <- annot
pca_result <- prcomp(abundance_matrix, center = TRUE, scale. = TRUE)

# Converter resultados para data frame
pca_df <- data.frame(Sample = rownames(abundance_matrix),
                     PC1 = pca_result$x[,1],
                     PC2 = pca_result$x[,2])

# Plotar PCA
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(vjust = 1.5, size = 3) +
  theme_minimal() +
  ggtitle("PCA das amostras de microbioma")
##############################################








## PCA
outlier <- c("SRR23584467","SRR23584495")
newAnnot <- annot[!annot %in% outlier]
abundance_matrix <- subset(abundance_matrix,
                                    !rownames(abundance_matrix) %in% outlier)

# Remover colunas com variação zero
abundance_matrix_filtered <- abundance_matrix[, apply(abundance_matrix, 2, var) > 0]
rownames(abundance_matrix_filtered) <- newAnnot

# Verificar se ainda há colunas suficientes para o PCA
if (ncol(abundance_matrix_filtered) > 1) {
  # Aplicar PCA apenas se houver colunas suficientes
  pca_result <- prcomp(abundance_matrix_filtered, center = TRUE, scale. = TRUE)
  
  # Converter resultados para data frame
  pca_df <- data.frame(Sample = rownames(abundance_matrix_filtered),
                       PC1 = pca_result$x[,1],
                       PC2 = pca_result$x[,2])
  
  # Plotar PCA
  ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point() +
    geom_text(vjust = 1.5, size = 3) +
    theme_minimal() +
    ggtitle("PCA das amostras de microbioma")
  
} else {
  print("Erro: Poucas colunas com variação suficiente para o PCA.")
}





# Clusterização (K-means) para Detectar Perfis Diferentes
set.seed(42)  # Para reprodutibilidade
k_clusters <- kmeans(abundance_matrix, centers = 3)  # Ajuste `centers` conforme necessário

# Adicionar clusters ao data frame
pca_df$Cluster <- as.factor(k_clusters$cluster)

# Visualizar clusters no PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  theme_minimal() +
  ggtitle("Clusters das amostras de microbioma (K-means)")


########################################

pca_df$Label <- ifelse(pca_df$Sample %in% distinct_samples, pca_df$Sample, NA)

# Plot com rótulos apenas para amostras distintas
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Label)) +
  geom_point(size = 3) +
  #geom_text(vjust = 1.5, size = 2, na.rm = TRUE) +  # Apenas adiciona labels não-NA
  theme_minimal() +
  ggtitle("Clusters das amostras de microbioma (K-means)")


###################################################################

#### Most frequent species
specie <- data %>% filter(Rank == "S")
specie <- specie %>% filter(scientific_name != "Homosapiens")
specie <- specie %>% filter(Percentage >= 0.5)


genera <- data %>% filter(Rank == "G")
genera <- genera %>% filter(scientific_name != "Homo")
genera <- genera %>% filter(scientific_name != "Escherichia")
genera <- genera %>% filter(Percentage >= 0.99)


genera %>%
  ggplot(aes(x=scientific_name, y=Percentage, fill = scientific_name)) +
  geom_boxplot() +
  theme_ipsum() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(legend.position="none") +
  coord_flip()

###############################################################
# Phylum frequency
phylum <- data %>% filter(Rank == "P")
phylum <- phylum %>% filter(scientific_name != "Chordata")
#phylum <- phylum %>% filter(scientific_name != "Proteobacteria")
phylum <- phylum %>% filter(Percentage > 1)

phylum %>%
  ggplot(aes(x=scientific_name, y=Percentage, fill = scientific_name)) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(legend.position="none") +
  coord_flip()

# Genus freq
genera <- data %>% filter(Rank == "G")
genera <- genera %>% filter(scientific_name != "Homo")
genera <- genera %>% filter(scientific_name != "Escherichia")
genera <- genera %>% filter(Percentage > 0.7)

my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(30)

genera %>%
  ggplot(aes(x=scientific_name, y=Percentage, fill = scientific_name)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = my_colors) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(legend.position="none") +
  coord_flip()





#### t-SNE analysis #####

# Import data species
specie <- data %>% filter(Rank == "S")
specie <- specie %>% filter(scientific_name != "Homosapiens")
pcaSpecie <- specie[,c(1,2,7)]

# Genus
genera <- data %>% filter(Rank == "G")
genera <- genera %>% filter(scientific_name != "Homo")
pcaGenera <- genera[,c(1,2,7)]

####################################################################

# Annotation data
annot <- read.csv("annot_metadata.tsv",sep="\t")

wild <- annot %>% filter(Local == "wild")
greenhouse <- annot %>% filter(Local == "greenhouse")

# Color palette
my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(29)


pcaSpecie <- pcaSpecie %>%
  mutate(Percentage = ifelse(Percentage == 0, 0.0, as.numeric(Percentage)))

pcaData <- pcaSpecie %>%
  mutate(Percentage = as.double(Percentage))


pcaData <- pcaData %>%
  group_by(sample_name, scientific_name) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = scientific_name, values_from = Percentage, values_fill = list(Percentage = 0))


tsne_data <- pcaData %>%
  select(-sample_name) 



set.seed(42)
tsne_result <- Rtsne(tsne_data, dims = 2, perplexity = 15, verbose = TRUE, max_iter = 500)


tsne_df <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  sample_name = pcaData$sample_name  # Adicionar de volta os nomes das amostras
)


tsne_annot_df <- tsne_df %>%
  left_join(annot, by = c("sample_name" = "Run"))

ggplot(tsne_annot_df, aes(x = X, y = Y, color = Continent, shape = Local)) +
  geom_point(alpha = 0.85, size = 3.5) + 
  theme_ipsum() +
  labs(title = "Species relative abundance", x = "Dim 1", y = "Dim 2") +
  theme(legend.position = "right")  



























###########################################
#### Metagenome analysis of Phytophthora libraries #####


PRJNA785999 <- c("SRR17587545","SRR17587546","SRR17587547",
                 "SRR17587548","SRR17587549","SRR17587550",
                 "SRR17587551","SRR17587552","SRR17587553",
                 "SRR17587554","SRR17587555","SRR17587556",
                 "SRR17587557")

annot <- data.frame(
  libs = c("SRR17587545","SRR17587546","SRR17587547",
           "SRR17587548","SRR17587549","SRR17587550",
           "SRR17587551","SRR17587552","SRR17587553",
           "SRR17587554","SRR17587555","SRR17587556",
           "SRR17587557"),
  condition = c("P_palmivora","Control","Control","Control",
                "P_palmivora","P_palmivora","P_palmivora",
                "Control","Control","P_palmivora","P_palmivora",
                "Control","Control")
)





# Phylum data
phyto <- data %>% filter(sample_name %in% PRJNA785999)

phylum <- phyto %>% filter(Rank == "P")
phylum <- phylum %>% filter(scientific_name != "Chordata")
pcaPhylum <- phylum[,c(1,2,7)]

bioValues <- phylum %>%
  group_by(scientific_name, sample_name) 

bioValues <- bioValues %>% filter(Percentage > 0.05)

my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(8)

ggplot(bioValues, aes(x = reorder(sample_name,scientific_name), y = Percentage, fill = scientific_name)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  coord_flip() +
  theme_ipsum() +
  scale_fill_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 2, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )





# Genus data

genusData <- phyto %>% filter(Rank == "G")
genusData <- genusData %>% filter(scientific_name != "Homo")
pcagenusData <- genusData[,c(1,2,7)]

bioValues <- genusData %>%
  group_by(scientific_name, sample_name) 

bioValues <- bioValues %>% filter(Percentage > 0.01)

bioValues <- bioValues %>%
  mutate(Group = ifelse(sample_name %in% annot$condition,
                        "P_palmivora", "Control"))




my_colors <- colorRampPalette(brewer.pal(8, "Set1"))(300)

ggplot(bioValues%>%filter(Group==P_palmivora), aes(x = reorder(sample_name,scientific_name), y = Percentage, fill = scientific_name)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_ipsum() +
  scale_fill_manual(values = my_colors) +
  theme(
    strip.text = element_text(size = 2, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )



