library("DESeq2")
library(tidyr)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

which_hypotheses <- function(hog)
{
  
  print(paste(hog, "is in:"))
  
  if(HOG_level_list$hypothesis_1[which(HOG_level_list$hypothesis_1$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 1")
  
  if(HOG_level_list$hypothesis_2[which(HOG_level_list$hypothesis_2$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 2")
  
  if(HOG_level_list$hypothesis_3[which(HOG_level_list$hypothesis_3$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 3")
  
  if(HOG_level_list$hypothesis_4[which(HOG_level_list$hypothesis_4$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 4")
  
  if(HOG_level_list$hypothesis_5[which(HOG_level_list$hypothesis_5$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 5")
  
  if(HOG_level_list$hypothesis_6[which(HOG_level_list$hypothesis_6$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 6")
  
  if(HOG_level_list$hypothesis_7[which(HOG_level_list$hypothesis_7$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 7")
  
  if(HOG_level_list$hypothesis_8[which(HOG_level_list$hypothesis_8$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 8")
  
  if(HOG_level_list$hypothesis_9[which(HOG_level_list$hypothesis_9$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 9")
  
  if(HOG_level_list$hypothesis_10[which(HOG_level_list$hypothesis_10$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 10")
  
  if(HOG_level_list$hypothesis_11[which(HOG_level_list$hypothesis_11$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 11")
  
  if(HOG_level_list$hypothesis_12[which(HOG_level_list$hypothesis_12$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 12")
  
}

which_hypotheses("N0.HOG0011860")

