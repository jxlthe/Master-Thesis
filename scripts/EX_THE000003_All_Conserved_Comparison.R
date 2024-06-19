#load data ####################
library("DESeq2")
library(dplyr)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

all_HOGs <- HOG_DE.a2tea[HOG_DE.a2tea$HOG != "singleton",]

#all
Gene_Ta_all <- all_HOGs$gene[grepl("Traes", all_HOGs$gene)]
Gene_Os_all <- all_HOGs$gene[grepl("Os", all_HOGs$gene)]
Gene_Ms_all <- all_HOGs$gene[grepl("Misin", all_HOGs$gene)]
Gene_Sb_all <- all_HOGs$gene[grepl("SORBI", all_HOGs$gene)]
Gene_Zm_all <- all_HOGs$gene[grepl("Zm", all_HOGs$gene)]

conserved_HOGs <- read.csv(".\\conserved_HOGs.csv")

conserved_data <- all_HOGs[all_HOGs$HOG %in% conserved_HOGs$x, ]

Gene_Ta_conserved <- conserved_data$gene[grepl("Traes", conserved_data$gene)]
Gene_Os_conserved <- conserved_data$gene[grepl("Os", conserved_data$gene)]
Gene_Ms_conserved <- conserved_data$gene[grepl("Misin", conserved_data$gene)]
Gene_Sb_conserved <- conserved_data$gene[grepl("SORBI", conserved_data$gene)]
Gene_Zm_conserved <- conserved_data$gene[grepl("Zm", conserved_data$gene)]

species <- c("Triticum_aestivum", "Oryza_sativa", "Miscanthus_sinensis", "Sorghum_bicolor", "Zea_mays")
all_ProtCount <- c(length(Gene_Ta_all), length(Gene_Os_all), length(Gene_Ms_all), length(Gene_Sb_all), length(Gene_Zm_all))
all_ProtCount_ratio <- all_ProtCount / length(all_HOGs$gene)
conserved_ProtCount <- c(length(Gene_Ta_conserved), length(Gene_Os_conserved), length(Gene_Ms_conserved), length(Gene_Sb_conserved), length(Gene_Zm_conserved))
conserved_ProtCount_ratio <- conserved_ProtCount / length(conserved_data$gene)

df <- data.frame(Species = species, ProtInAll = all_ProtCount, ProtInAllRatio = all_ProtCount_ratio,
                 ProtInConserved = conserved_ProtCount, ProtInConservedRatio = conserved_ProtCount_ratio)


#this is the number of HOGs in all conserved
length(unique(HOG_DE.a2tea$HOG))
length(conserved_HOGs$x)
length(conserved_HOGs$x) / length(unique(HOG_DE.a2tea$HOG))

#this is the number of Proteins in all and conserved
length(all_HOGs$gene)
length(conserved_data$gene)
length(conserved_data$gene) / length(all_HOGs$gene)
