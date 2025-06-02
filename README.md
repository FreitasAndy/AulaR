# AulaR
Scripts e dados para a aula de "Análise e visualização de dados no R", da disciplina "CEN5807 - Bioinformática aplicada ao estudo de comunidades microbianas" do CENA/USP.

## Primeiros passos

Instale o R e o RStudio.

Na sequência, instale os pacotes que utilizaremos com os seguintes comandos:

pacotes <- c("tidyverse","ggplot2", "devtools", "microeco", "dplyr", "magrittr", "ggpubr", 
             "vegan", "randomForest")

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

## Seguindo a aula

Baixe todo o conteúdo aqui do GitHub e salve-o em uma única pasta separada.
Abra o arquivo Aula Dados R.Rproj no RStudio e siga os passos em aula.
