library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
library(ggExtra)

# Create a color palette - RColorBrewer
my_palette <- colorRampPalette(brewer.pal(8, "Set1"))(19)


data <- read.csv("BLAST-allContigs.csv")
data <- unique(data)

viralMaterial <- data[grepl("virus|phage", data$BLASTx_Organism_Name, ignore.case = TRUE), ]
HostMaterial <- data[grepl("Theobroma cacao", data$BLASTx_Organism_Name, ignore.case = TRUE), ]
metaMaterial <- data[!grepl("Theobroma cacao", data$BLASTx_Organism_Name, ignore.case = TRUE), ]

speciesList <- unique(metaMaterial$BLASTx_Organism_Name)

#### PLOT 01) Viral Alignments ####

hits <- data.frame(
  Type <- c("Meta", "Viral", "Host"),
  Value <- c(110338, 253, 37185)
)

hits$log10 <- round(log10(hits$Value),1)

ggplot(hits, aes(x=Type,y=log10, fill=Type)) +
  geom_bar(stat = "identity", ) +
  theme_classic() +
  #geom_text(aes(label=paste(Percentage)), 
  #          position=position_dodge(width=0.9), 
  #          hjust=-0.25, size=3.5) +
  scale_fill_manual(values = my_palette) +
  ylim(0,5) +
  #coord_flip() +
  xlab("")+
  #ggtitle("Contig BLASTx Alignment")+
  theme(legend.position="none",text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))


#### viral-like alignments ####

polyprotein <- viralMaterial[grepl("polyprotein", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
coat <- viralMaterial[grepl("coat", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
capsid <- viralMaterial[grepl("capsid", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
rdrp <- viralMaterial[grepl("RdRp|RNA-dependent RNA polymerase|Polymerase|Replicase", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
hypo <- viralMaterial[grepl("hypothetical protein", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
RT <- viralMaterial[grepl("RVT|Transcriptase", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
MP <- viralMaterial[grepl("MP|movement protein", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]
NS_poly <- viralMaterial[grepl("non-structural polyprotein", viralMaterial$BLASTx_Subject_Title, ignore.case = TRUE), ]


hit_freq <- data.frame(
  type = c("Capsid Protein","Coat Protein","Hypothetical Protein",
           "Movement Protein","Non-structural Polyprotein",
           "Polyprotein","RdRp","Reverse Transcriptase"),
  value = c(25,7,29,11,5,79,21,22)
) 




ggplot(hit_freq, aes(x=reorder(type,value) ,y=value, fill="darkred")) +
  geom_bar(stat = "identity", ) +
  theme_classic() +
  scale_fill_manual(values = my_palette) +
  coord_flip() +
  xlab("")+
  ggtitle("Viral-like BLASTx Alignment")+
  theme(legend.position="none")







#### BLAST scatter plot ####

# Filogenia
filo <- read.csv("filogenia.csv")

viralMaterial$BLASTx_GenomeComposition[is.na(viralMaterial$BLASTx_GenomeComposition)] <- NA

###
filo$domains <- ifelse(filo$Contig.Domains == filo$BestHit.Domains, "Complete", "Fragment")


p <- ggplot(viralMaterial, aes(x=BLASTx_Cover ,y=BLASTx_Ident, )) +
  geom_point(alpha=0.5, size=4, color="darkred") +
  theme_classic() +
  scale_fill_manual(values = my_palette) +
  coord_flip()
#theme(legend.position="none")
p

p2 <- ggMarginal(p, type = "boxplot", color="darkred")
p2


#### Viral analysis ####

# viral contigs plot

ggplot(filo, aes(x = domains, y = ContigLength, fill = domains)) +
  geom_boxplot() +
  coord_flip() +  # Para facilitar a leitura
  labs(x = "Conserved domains", y = "Contig lenght (nt)")+
  theme_minimal()


# Contagem e cálculo da porcentagem
domain_counts <- filo %>%
  count(domains) %>%
  mutate(perc = (n / sum(n)) * 100)  # Calcula a porcentagem

# Criando o gráfico de pizza
ggplot(domain_counts, aes(x = "", y = n, fill = domains)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  # Adiciona a porcentagem
            position = position_stack(vjust = 0.5)) +  # Centraliza os textos
  labs(fill = "Conserved domains structure") +
  theme_void()





#### boxplot ####

ggplot(filo, aes(x = reorder(family, BLASTn.Ident), y = BLASTn.Ident)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +  # Melhor para visualizar muitas categorias
  labs(x = "", y = "Nucleotide Identity (%)")





### Solanum tuberosum contamination ###
species <- data.frame(data$Sample_name, data$BLASTx_Organism_Name)
colnames(species) <- c("sample","species")
potato <- species %>% filter(species == "Solanum tuberosum")

sample_counts <- potato %>%
  count(sample) %>%
  arrange(desc(n))

ggplot(sample_counts, aes(x = reorder(sample, -n), y = n)) +
  geom_bar(stat = "identity", fill = "#619CFF", color = "black") +
  labs(
    title = "Frequence of S. tuberosum contigs in samples",
    x = "Sample",
    y = "Number of assembled contigs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
