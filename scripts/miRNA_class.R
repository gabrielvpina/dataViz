library(dplyr)
library(ggplot2)
library(viridis)

data <- read.table("result_bedtools_intersect.tsv", sep = "\t", header = TRUE)

table(data$Annot)

data$qPosition <- paste(data$qStart,"-",data$qEnd)

dados <- data.frame(
  class = c("Gene", "Exon", "Pseudogene", "Other"),
  freq = c(140, 22, 2, 29)
)


absolute <- unique(data$qPosition)

# ACHAR INTRONS
genes <- subset(data, Annot == "gene")
exon <- subset(data, Annot == 'exon')
cds <- subset(data, Annot == 'CDS')
region <- subset(data, Annot == 'region')
mRNA <- subset(data, Annot == 'mRNA')
pseudogene <- subset(data, Annot == 'pseudogene') 


region_positions <- unique(region$qPosition)
cds_positions <- unique(cds$qPosition)
exon_positions <- unique(exon$qPosition)
gene_positions <- unique(genes$qPosition)
mRNA_positions <- unique(mRNA$qPosition)
pGene_positions <- unique(pseudogene$qPosition)

intersect(gene_positions, region_positions)


miRNAs <- data.frame(
  class = c("exon","intron","pseudogene","other"),
  value = c(22,96,2,73)
)








introns <- subset(data, qPosition %in% gene_positions & 
                    Annot != "exon" & 
                    !qPosition %in% subset(data, Annot == "exon")$qPosition)
# Intron number
intron_number <- unique(introns$qPosition)

# OUTSIDE GENES
qPositions_region <- data %>%
  group_by(qPosition) %>%
  summarise(only_region = all(Annot == "region")) %>%
  filter(only_region) %>%
  pull(qPosition)

regions_only <- data %>%
  filter(qPosition %in% qPositions_region)










data <- miRNAs
# Compute percentages
data$fraction = data$value / sum(data$value)
# Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$class,": ", round(data$fraction*100),"%")
# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=class)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label, color=class), size=6) + # x here controls label position (inner / outer)
  scale_fill_viridis(discrete = TRUE, alpha=0.7, option = "D") +
  scale_color_viridis(discrete = TRUE, alpha=0.7, option = "D") +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")



#________________________________________
#Clusters


region <- subset(data, Annot == 'region')

clusters <- table(region$QueryScaffold)
df <- as.data.frame(clusters)
colnames(df) <- c("scaffold","EVEs")

write.table(region,"EVEs_by_scaffold.tsv", row.names = FALSE)
write.table(df,"EVEs_clusters.tsv", row.names = FALSE)

#_________________________________________
# Open table

clusters <- region[,c(1:3,8)]

check_neighbors <- function(data, threshold = 5000) {
  data <- data %>%
    arrange(qStart) %>% # Ordenar por qStart para facilitar a lógica
    mutate(
      has_neighbor = sapply(1:n(), function(i) {
        any(
          abs(qStart[i] - qEnd[-i]) < threshold |
            abs(qEnd[i] - qStart[-i]) < threshold
        )
      })
    )
  return(data)
}

# Aplicar função agrupando por QueryScaffold
result <- clusters %>%
  group_by(QueryScaffold) %>%
  group_modify(~ check_neighbors(.x)) %>%
  ungroup()


res <- as.data.frame(result)
res_clusters <- res %>% filter(has_neighbor == TRUE)

# write.table(res_clusters, "miRNA_clusters", row.names = FALSE)

res_5k <- as.data.frame(result)
res_clusters_5k <- res_5k %>% filter(has_neighbor == TRUE)



#________________________________________________

check_overlaps <- function(data) {
  data <- data %>%
    arrange(qStart) %>% # Ordenar por qStart para facilitar a lógica
    mutate(
      overlaps = sapply(1:n(), function(i) {
        any(
          (qStart[i] >= qStart[-i] & qStart[i] <= qEnd[-i]) | # Início sobrepõe outro intervalo
            (qEnd[i] >= qStart[-i] & qEnd[i] <= qEnd[-i])       # Final sobrepõe outro intervalo
        )
      })
    )
  return(data)
}

# Aplicar função agrupando por QueryScaffold
result <- res_clusters %>%
  group_by(QueryScaffold) %>%
  group_modify(~ check_overlaps(.x)) %>%
  ungroup()

overlap <- as.data.frame(result)
overlap <- overlap %>% filter(overlaps == TRUE)

write.table(result, "result_cluster_miRNA.tsv", sep = "\t", row.names = FALSE)
