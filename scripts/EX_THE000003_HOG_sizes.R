library("DESeq2")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

conserved_HOGs <- read.csv(".\\conserved_hogs.csv")
conserved_HOGs <- conserved_HOGs$x
conserved_HOGs <- conserved_HOGs[conserved_HOGs != "singleton"]
conserved_HOGs <- HOG_DE.a2tea[HOG_DE.a2tea$HOG %in% conserved_HOGs,]$HOG

all_HOGs <- HOG_DE.a2tea$HOG
all_HOGs <- all_HOGs[all_HOGs != "singleton"]

conserved_sizes <- table(conserved_HOGs)
all_sizes <- table(all_HOGs)

summary(as.numeric(conserved_sizes))
summary(as.numeric(all_sizes))
