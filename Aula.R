################################################################
## Script - Visualizacao de dados no R
## Autor: Anderson Freitas
## Dicas: - nao usar espacos nos nomes de arquivos e objetos
##        - nao usar cedilhas nem acentos
##        - atencao ao uso de maiusculas e minusculas
################################################################

#Instalando Pacotes

pacotes <- c("tidyverse","ggplot2", "devtools", 
             "microeco", "dplyr", "magrittr",
             "ggpubr", "vegan", "randomForest")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ALDEx2")
install.packages("file2meco", repos = BiocManager::repositories())


####

getwd()

#setwd("C:/Users/Dinos/OneDrive/Documentos")


#Chamando as bibliotecas
library(tidyverse)
library(phyloseq)
library(microeco)
library(ggplot2)
library(magrittr)
library(file2meco)
library(ggpubr)
library(ALDEx2)
library(randomForest)

####

#Importando arquivos
otu2 <- read_delim("otu2.csv", delim = "\t", 
                  escape_double = FALSE, trim_ws = TRUE)
otu_mat <- otu2 %>%
  tibble::column_to_rownames("OTU_ID") 

tax2 <- read_delim("tax2.csv", delim = "\t", 
                escape_double = FALSE, trim_ws = TRUE)
tax_mat <- tax2 %>% 
  tibble::column_to_rownames("OTU_ID")

mapp <- read_csv("mapp.csv")
samples_df <- mapp %>% 
  tibble::column_to_rownames("Sample_ID") 

####

#Criando um objeto R6 do microeco
dataset <- microtable$new(sample_table = samples_df,
                          otu_table    = otu_mat,
                          tax_table    = tax_mat)

dataset$tidy_dataset()
dataset
dataset$tax_table %<>% tidy_taxonomy

#Criando um objeto phyloseq
physeq <- meco2phyloseq(dataset)
physeq
microbiome::summarize_phyloseq(physeq)

dataset.r <- clone(dataset)
dataset.r$rarefy_samples(sample.size = 82270)
dataset.r$tidy_dataset()
dataset.r$tax_table %<>% tidy_taxonomy

####

#Abundancia
t1 <- trans_abund$new(dataset = dataset.r, taxrank = "Phylum", ntaxa = 8) #separando por filo
t1$plot_bar(others_color = "grey70", facet = "Substrate", xtext_keep = FALSE, legend_text_italic = FALSE)

t1 <- trans_abund$new(dataset = dataset.r, taxrank = "Family", ntaxa = 10) #separando por família
t1$plot_bar(others_color = "pink", facet = "Sand", xtext_keep = FALSE, legend_text_italic = TRUE)

# O parametro groupmean pode ser usado para fazer uma barra por grupo
t1 <- trans_abund$new(dataset = dataset.r, taxrank = "Phylum", ntaxa = 10, groupmean = "Substrate")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1
g1 + theme_classic() + 
  theme(axis.title.y = element_text(size = 10)) 


####

#Heatmap
t1 <- trans_abund$new(dataset = dataset.r, taxrank = "Genus", ntaxa = 6)
t1$plot_heatmap(facet = "Substrate", xtext_keep = FALSE, withmargin = FALSE)
lara.teste <- t1$plot_heatmap(facet = "Substrate", xtext_keep = FALSE, withmargin = FALSE)
lara.teste + theme(axis.title.y = element_text(size = 10))


#Donut plot
t1 <- trans_abund$new(dataset = dataset.r, taxrank = "Phylum", ntaxa = 8, groupmean = "Substrate")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)

teste.aula <- t1$plot_donut(label = TRUE)
teste.aula

ggplot2::ggsave(filename = "Donut_Plot.svg", 
                plot = teste.aula,
                device = "svg",
                dpi = 1200,
                width = 20,
                height = 12,
                units = "cm")

#Diagrama de Venn
dataset.venn <- dataset$merge_samples("Substrate")
t1 <- trans_venn$new(dataset.venn, ratio = "seqratio")
t1$plot_venn()

####

#Alfa diversidade
physeq.r <- meco2phyloseq(dataset.r)
diversity <- microbiome::alpha(physeq.r, index = "all")
meta=microbiome::meta(physeq.r)
alpha= cbind(diversity,meta)

#Grafico
alpha$Substrate <- factor(alpha$Substrate, levels = c("Control", "ADE"))

alpha.plot = 
  ggplot(data = alpha, aes(x=Substrate, y=observed, fill = Substrate)) +
  geom_boxplot() +
  labs(x = "Substrate", y= "Number of Taxa") +
  labs(color='Substrate') + guides(color = "none")+
  scale_fill_manual(values = c("#cb6751", "#7aa457")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

alpha.plot

#Grafico com valor de p
comparison = list(c("Control", "ADE"))

alpha.p = 
  ggplot(data = alpha, aes(x=Substrate, y=observed, fill = Substrate)) +
  geom_boxplot() +
  labs(x = "Substrate", y= "Number of Taxa") +
  labs(color='Substrate') + guides(color = "none")+
  scale_fill_manual(values = c("#cb6751", "#7aa457")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test")

alpha.p


alpha.bar = 
  ggplot(data = alpha, aes(x=Substrate, y=observed, fill = Substrate)) +
  geom_bar(stat="identity", color="black") +
  labs(x = "Substrate", y= "Number of Taxa") +
  labs(color='Substrate') + guides(color = "none")+
  scale_fill_manual(values = c("#cb6751", "#7aa457")) +
  theme_bw() 

alpha.bar

###

#Beta diversidade
physeq.clr = microbiome::transform(physeq, "clr")
input_ord = ordinate(physeq.clr, method = "PCoA" , "euclidean") 

p4 = plot_ordination(physeq.clr, input_ord, color = "Substrate")
p4
p1 = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right")+
  theme_bw()
p1

df        = as(sample_data(physeq.clr), "data.frame")
ds        = phyloseq::distance(physeq.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate, data = df, permutations = 999)
permanova

p2 = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -60, y = 100, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.27  |  p = 0.03'), size = 3)+
  theme_bw()
p2

p3 = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 70, y = -200, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.27  |  p = 0.03'), size = 3)+
  stat_ellipse(type = "t") +
  theme_bw()
p3
       
####

#Abundância diferencial
#ALDEx2

genera = microbiome::aggregate_rare(physeq, level = "Genus",
                                     detection = 1/100, prevalence = 1/100)
mi    = as.data.frame((otu_table(genera)))
var   = sample_data(genera)
treat = var$Substrate
x <- aldex(mi, treat, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=TRUE, denom="iqlr", verbose=TRUE)
aldex.plot(x, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")
#Separando os mais importantes
names  <- rownames(x)
raball <- cbind(x, names)
newdata1 <- raball[order(raball$wi.ep),] 
res1=(newdata1[(newdata1$wi.ep<="0.05"),])
res1

#Random Forest
#Vamos usar o mesmo arquivo de entrada da RandomForest
genera.clr = microbiome::transform(genera, "clr")
predictors <- t(otu_table(genera))
dim(predictors)
response <- as.factor(sample_data(genera.clr)$Substrate)
rf.data <- data.frame(response, predictors)
set.seed(27)
erie.classify <- randomForest(response~., data = rf.data, ntree = 10000)
print(erie.classify)
# Make a data frame with predictor names and their importance
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top 10 predictors
imp.20 <- imp.sort[1:20, ]
# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  ggtitle("Preditores") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

####

#Networks
#ADE
group_ADE <- clone(dataset)
group_ADE$sample_table <- subset(group_ADE$sample_table, Substrate == "ADE")
group_ADE$tidy_dataset()
group_ADE
group_ADE$cal_abund()
set.seed(27)
t1 <- trans_network$new(dataset = group_ADE, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.005)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
t1$res_network_attr
t1$save_network(filepath = "./ADE_network.gexf")

#Control
group_Control <- clone(dataset)
group_Control$sample_table <- subset(group_Control$sample_table, Substrate != "ADE")
group_Control$tidy_dataset()
group_Control
group_Control$cal_abund()
set.seed(27)
t1 <- trans_network$new(dataset = group_Control, 
                        cor_method = "sparcc",
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.005)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
t1$res_network_attr
t1$save_network(filepath = "Control_network.gexf")


t1 <- trans_network$new(dataset = group_ADE, 
                        cor_method = "spearman", 
                        use_WGCNA_pearson_spearman = TRUE, 
                        filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.7)
t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
t1$res_network_attr
t1$save_network(filepath = "ADE_WGCNA_network.gexf")

####

#Analisando dados que não são de NGS

shapiro.test(samples_df$Root_size)
#p>0.05 = normal; p<0.05= não normal

kruskal.test(x = samples_df$Root_size, g = samples_df$Substrate)

dunn.test::dunn.test(x = samples_df$Root_size, g = samples_df$Substrate)

#ANOVA

one.way <- aov(Root_size ~ Substrate, data = samples_df)
summary(one.way)

samples_df[is.na(samples_df)] <- 0
colnames(samples_df)
my.variables <- samples_df[c(3:15)]

for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = samples_df$Substrate)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)
  } else {
    print(colnames(my.variables[i]))
    print("No difference")
  }
}

#PCA

colnames(samples_df)

chem   <- prcomp(samples_df[,-c(1:2)], scale=TRUE)
summary(chem)
print(chem)   
biplot <- data.frame(chem$x, Substrate = samples_df$Substrate)
#biplot
ggplot(data = biplot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Substrate, shape = Substrate), size = 4) +
  labs(
    x = "PC1 (73.38%)",
    y = "PC2 (12.04%)",
    color = "Substrate") +
  scale_color_manual(values=c("red", "green"),
                     labels = c("Control", "ADE")) +
  theme(legend.position="bottom") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) + 
  theme_bw()

dev.print(tiff, "./biplotPCA.tiff", compression = "lzw", res=600, height=6, width=9, units="in")



factoextra::fviz_pca_biplot(chem, ggtheme = theme_bw(), repel = T)
#dev.print(tiff, "./Figures/biplotPCA_Color.tiff", compression = "lzw", res=600, height=6, width=9, units="in")




### Extra: grafico de barras com desvio
sum.mb = Rmisc::summarySE(data = Mapfile, measurevar = "Conc", groupvars = c("Treatment", "Source", "Depth"))
sum.mb = sum.mb %>% mutate(Source=factor(Source, levels=c("ADE Forest", "ADE Cassava", "Oxisol Cassava")))
sum.mb$Depth <- as.character(sum.mb$Depth)
sum.mb$Conc <- as.numeric(sum.mb$Conc)
mb = 
  ggplot(data = sum.mb, aes(x=Depth, y=Conc, fill = Source))  +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Conc, ymax=Conc+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  labs(x = "", y= "DNA Concentration \n (ng/uL)", title = "Microbial Biomass") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  guides(color = "none")+
  scale_fill_manual(values = c("darkgreen" ,"#7aa457", "#cb6751"))
mb
