library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pheatmap)

# Create a color palette - RColorBrewer
my_palette <- colorRampPalette(brewer.pal(8, "Set1"))(5)  # 5 colors

# Create a color palette - Viridis (A,B,C,D,E)
my_palette <- viridis(5, option = "B") 

data <- read.csv("tables/total_quantData.csv")
data <- data[-1,]
data <- data %>% select(-c(SRR23584496))
index <- read.csv("tables/quant.sf",sep="\t")


#write table
#write.csv(data,"tables/final_count.csv")


rownames(data) <- index$Name

viralquant <- data[c(1:1353),]
viralquant[] <- sapply(viralquant, as.numeric)

# Annot info
metadata <- read.csv("tables/annot_pheatmap.tsv",sep="\t")

# Create matrix
data_total <- as.matrix(viralquant[rowSums(viralquant)>0,])

# Annotation - Local
colinfo <- metadata[,c(1,6)]
rownames(colinfo) <- colinfo$Run
colinfo <- as.data.frame(colinfo[,-1])

my_palette <- viridis(100, option = "B")

pheatmap(log10(data_total + 1), 
         scale="none", 
         fontsize = 6, 
         #cutree_cols = 6,
         cluster_rows = T,
         show_colnames = FALSE,
         cluster_cols = T,
         border_color = 'NA',
         color = rev(my_palette))
         #annotation_col = colinfo)









#### Bioproject Quant ####


PRJEB35419 <- c("ERR3665546","ERR3665547","ERR3665548",
                "ERR3665549","ERR3665550","ERR3665551",
                "ERR3665552","ERR3665553","ERR3665554",
                "ERR3665555","ERR3665556","ERR3665557")

PRJNA604260 <- c("SRR11060111","SRR11060112","SRR11060113")

PRJNA613342 <- c("SRR11389076","SRR11389077","SRR11389078",
                 "SRR11389079","SRR11389080","SRR11389081",
                 "SRR11389082","SRR11389083","SRR11389084")

PRJNA742476 <- c("SRR15010977","SRR15010978","SRR15010979",
                 "SRR15010980","SRR15010981","SRR15010982",
                 "SRR15010983","SRR15010984")

PRJNA785999 <- c("SRR17587545","SRR17587546","SRR17587547",
                 "SRR17587548","SRR17587549","SRR17587550",
                 "SRR17587551","SRR17587552","SRR17587553",
                 "SRR17587554","SRR17587555","SRR17587556",
                 "SRR17587557")

PRJNA558793 <- c("SRR23455105","SRR23455106","SRR23455107",
                 "SRR23455108","SRR23455109","SRR23455110",
                 "SRR23455111","SRR23455112","SRR23455113",
                 "SRR23455114","SRR23455115","SRR23455116",
                 "SRR23455117","SRR23455118","SRR23566373",
                 "SRR23566374","SRR23566375","SRR23566376",
                 "SRR23566377","SRR23566381","SRR23566382",
                 "SRR23566393","SRR23566404","SRR23566405",
                 "SRR23584464","SRR23584465","SRR23584466",
                 "SRR23584467","SRR23584468","SRR23584473",
                 "SRR23584484","SRR23584495")

PRJNA314774 <- c("SRR3217276","SRR3217277","SRR3217278",
                 "SRR3217279","SRR3217280","SRR3217292",
                 "SRR3217294","SRR3217297","SRR3217298",
                 "SRR3217299","SRR3217301","SRR3217304",
                 "SRR3217315","SRR3217317","SRR3217318",
                 "SRR3217319","SRR3217280","SRR3217281",
                 "SRR3217282","SRR3217283","SRR3217284",
                 "SRR3217292","SRR3217294")

PRJNA421343 <- c("SRR7235732","SRR7235734")

PRJNA326055 <- c("SRR9203098","SRR9203099","SRR9203100",
                 "SRR9203101","SRR9203102","SRR9203103",
                 "SRR9203104","SRR9203105","SRR9203106",
                 "SRR9203107")



# Select BioProject
df <- viralquant %>% select(PRJEB35419)

# Create matrix
data_subset <- as.matrix(df[rowSums(df)>0,])

# Annotation - Local
#colinfo <- metadata[,c(1,5)]
#rownames(colinfo) <- colinfo$Run
#colinfo <- colinfo[PRJNA314774,]
#colinfo <- colinfo %>% select(2)

my_palette <- viridis(100, option = "B")
pheatmap(data_subset, scale = "row", fontsize = 12, 
         cluster_rows = T, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "YlGn"))(100),
         color = rev(my_palette),
         cluster_cols = T,
         #annotation_col = colinfo,
         fontsize_col = 10,
         angle_col = "45")











#### Known and Novel viruses ####

my_palette <- viridis(20, option = "B")


myData <- read.csv("tables/final_count.csv")
rownames(myData) <- myData$X
myData <- myData[,-1]
rownames(myData)

viralquant <- myData[c(1:1343),]
tail(viralquant)
viralquant[] <- sapply(viralquant, as.numeric)

rownames(viralquant)
df <- viralquant[rowSums(viralquant)>0,]
df$name <- rownames(df)

# import name table
ids <- read.csv("contig_names.tsv",sep="\t")

# Merge com preenchimento de valores ausentes
df <- left_join(df, ids, by = "name") %>%
  mutate(virus = coalesce(virus, name)) 
df <- df %>% relocate(virus)


# virus circulantes de cacau
cacaoVirus <- df[c(2:20),]
rownames(cacaoVirus) <- cacaoVirus$virus
cacaoVirus <- cacaoVirus[,-1]
cacaoVirus <- cacaoVirus[,-110]
cacaoVirus <- as.matrix(cacaoVirus)

pheatmap(cacaoVirus, 
         scale="row", 
         fontsize = 5, 
         #cutree_cols = 6,
         cluster_rows = T,
         show_colnames = T,
         cluster_cols = T,
         border_color = 'NA',
         color = rev(my_palette))



# virus novos encontrados
newVirus <- df[c(21:49,53:56),]
rownames(newVirus) <- newVirus$virus
newVirus <- newVirus[,-1]
newVirus <- newVirus[,-110]
exclude <- c("Cacao leafroll virus – contig19",
             "Badnavirus cocoa (BVC) – contig3324",
             "Badnavirus cocoa (BVC) – contig3706",
             "Citrus tatter leaf virus – contig13",
             "Red mite associated cystovirus – contig4",
             "Citrus tatter leaf virus – contig13",
             "Capillovirus mali –contig6875",
             "Capillovirus mali – contig3404",
             "apillovirus mali – contig6328")
newVirus <- newVirus[!(rownames(newVirus) %in% exclude), ]
newVirus <- as.matrix(newVirus)

my_palette <- colorRampPalette(c("yellow", "black"))(30)

pheatmap(newVirus, 
         scale="row", 
         fontsize = 6, 
         cutree_rows = 4,
         cluster_rows = T,
         show_colnames = T,
         cluster_cols = T,
         border_color = 'NA',
         color = my_palette)



# reference genes
ref <- df[c(53,55:58),]
rownames(ref) <- ref$virus
ref <- ref[,-1]
ref <- ref[,-110]
ref <- as.matrix(ref)



pheatmap(ref, 
         scale="row", 
         fontsize = 6, 
         #cutree_cols = 6,
         cluster_rows = T,
         show_colnames = FALSE,
         cluster_cols = T,
         border_color = 'NA',
         color = rev(my_palette))




# Virus quant
bvc <- df %>% filter(virus == "Badnavirus cocoa (BVC)")
bvc <- as.data.frame(t(bvc))
bvc$name <- rownames(bvc)
bvc <- bvc[-1,]
bvc <- bvc[-110,]
bvc$V1 <- as.numeric(bvc$V1)
rownames(bvc) <- NULL

bvc %>%
  filter(V1 > 0) %>%
    ggplot(aes(x=name,y=V1)) +
    geom_bar(stat = "identity") +
    theme_classic()


  

########## track viral origins ###########
library(tidyverse)


TcDV_1 <- myData %>% filter(X == "Contig28")
TcDV1_long <- TcDV_1 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)

ggplot(TcDV1_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()


TcDV_2 <- myData %>% filter(X == "Contig17")
TcDV2_long <- TcDV_2 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcDV2_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcDV2_long$sample)


TcDV_3 <- myData %>% filter(X == "SRR17587555-contig_4059-USA")
TcDV3_long <- TcDV_3 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcDV3_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcDV3_long$sample)


TcDV_4 <- myData %>% filter(X == "SRR17587555-contig_4318-USA")
TcDV4_long <- TcDV_4 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcDV4_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcDV4_long$sample)


TcFV_1 <- myData %>% filter(X == "Contig8")
TcFV1_long <- TcFV_1 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcFV1_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcFV1_long$sample)


TcFV_2 <- myData %>% filter(X == "SRR11060111-contig_16197-India")
TcFV2_long <- TcFV_2 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcFV2_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcFV2_long$sample)


TcKV_1 <- myData %>% filter(X == "SRR11060111-contig_422-India")
TcKV1_long <- TcKV_1 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcKV1_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcKV1_long$sample)


TcKV_2 <- myData %>% filter(X == "SRR23566404-contig_1591-USA")
TcKV2_long <- TcKV_2 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcKV2_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcKV2_long$sample)


TcPV_1 <- myData %>% filter(X == "Contig26")
TcPV1_long <- TcPV_1 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcPV1_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcPV1_long$sample)


TcPV_2 <- myData %>% filter(X == "Contig20")
TcPV2_long <- TcPV_2 %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(TcPV2_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(TcPV2_long$sample)



BVC <- myData %>% filter(X == "Contig21+23")
BVC_long <- BVC %>%
  pivot_longer(-X, names_to = "sample", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  filter(expression > 0)
ggplot(BVC_long, aes(y = expression, x = sample)) +
  geom_bar(stat = "identity", fill = "lightblue", color="grey") +
  theme_minimal()
print(BVC_long$sample)




