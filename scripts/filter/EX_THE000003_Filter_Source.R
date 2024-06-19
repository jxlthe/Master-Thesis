#load data ####################
library("DESeq2")
library(dplyr)


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

all_hogs <- unique(HOG_DE.a2tea$HOG)

#get GO Data
library(GO.db)
GO_TERMS <- Term(GOTERM)
GO_ONTOLOGY <- Ontology(GOTERM)
GO_DEFINITION <- Definition(GOTERM)

#function to extract GO Terms and related information from HOG list

extract_GO <- function(HOGs_list)
{
  
  #init vectors
  Ta_GO <- c()
  Os_GO <- c()
  Mi_GO <- c()
  SO_GO <- c()
  Zm_GO <- c()
  
  #get every unique GO Term list for every HOG
  for(hog in HOGs_list)
  {
    Ta_GO <- append(Ta_GO, SFA$Triticum_aestivum[SFA$Triticum_aestivum$HOG == hog,]$'Gene-Ontology-Term')
    Os_GO <- append(Os_GO, SFA$Oryza_sativa[SFA$Oryza_sativa$HOG == hog,]$'Gene-Ontology-Term')
    Mi_GO <- append(Mi_GO, SFA$Miscanthus_sinensis[SFA$Miscanthus_sinensis$HOG == hog,]$'Gene-Ontology-Term')
    SO_GO <- append(SO_GO, SFA$Sorghum_bicolor[SFA$Sorghum_bicolor$HOG == hog,]$'Gene-Ontology-Term')
    Zm_GO <- append(Zm_GO, SFA$Zea_mays[SFA$Zea_mays$HOG == hog,]$'Gene-Ontology-Term')
  }
  
  
  # split GO-Term list and combine all unique GO-Term into one vector
  Ta_GO <- unlist(strsplit(Ta_GO, ", "))
  Os_GO <- unlist(strsplit(Os_GO, ", "))
  Mi_GO <- unlist(strsplit(Mi_GO, ", "))
  SO_GO <- unlist(strsplit(SO_GO, ", "))
  Zm_GO <- unlist(strsplit(Zm_GO, ", "))
  
  Ta_GO_ID_count <- table(Ta_GO[!is.na(Ta_GO)])
  Os_GO_ID_count <- table(Os_GO[!is.na(Os_GO)])
  Mi_GO_ID_count <- table(Mi_GO[!is.na(Mi_GO)])
  SO_GO_ID_count <- table(SO_GO[!is.na(SO_GO)])
  Zm_GO_ID_count <- table(Zm_GO[!is.na(Zm_GO)])
  
  
  #combine and kick out NA
  all_GO_ID <- c(Ta_GO, Os_GO, Mi_GO, SO_GO, Zm_GO)
  all_GO_ID <- all_GO_ID[!is.na(all_GO_ID)] # remove NA
  all_GO_ID_total_count <- table(all_GO_ID) # unique and calculate the count of each GO ID this is the unsorted one
  all_GO_ID <- unique(all_GO_ID) # unique
  
  # get corresponding data for each GO_ID
  all_GO_Terms <- GO_TERMS[all_GO_ID]
  
  all_GO_total_count <- all_GO_ID_total_count[all_GO_ID] # count of how many times GO-Term was found. This one is the sorted version
  
  Ta_GO_count <- all_GO_total_count
  Ta_GO_count[] <- 0
  Ta_GO_count[rownames(Ta_GO_ID_count)] <- Ta_GO_ID_count
  Ta_GO_count <- as.array(Ta_GO_count)
  
  Os_GO_count <- all_GO_total_count
  Os_GO_count[] <- 0
  Os_GO_count[rownames(Os_GO_ID_count)] <- Os_GO_ID_count
  Os_GO_count <- as.array(Os_GO_count)
  
  Mi_GO_count <- all_GO_total_count
  Mi_GO_count[] <- 0
  Mi_GO_count[rownames(Mi_GO_ID_count)] <- Mi_GO_ID_count
  Mi_GO_count <- as.array(Mi_GO_count)
  
  SO_GO_count <- all_GO_total_count
  SO_GO_count[] <- 0
  SO_GO_count[rownames(SO_GO_ID_count)] <- SO_GO_ID_count
  SO_GO_count <- as.array(SO_GO_count)
  
  Zm_GO_count <- all_GO_total_count
  Zm_GO_count[] <- 0
  Zm_GO_count[rownames(Zm_GO_ID_count)] <- Zm_GO_ID_count
  Zm_GO_count <- as.array(Zm_GO_count)
  
  all_GO_Ontology <- GO_ONTOLOGY[all_GO_ID] # Ontology
  
  all_GO_Definition <- GO_DEFINITION[all_GO_ID] # definition of the GO-Term
  
  
  #count HOGs which have a particular GO Term
  all_HOGs_with_GO_count <- rep(0, length(all_GO_ID)) # definition of the GO-Term
  all_HOGs_with_GO <- rep("NA", length(all_GO_ID))
  
  SFA_bound <- bind_rows(SFA$Oryza_sativa, SFA$Triticum_aestivum, SFA$Sorghum_bicolor, SFA$Miscanthus_sinensis, SFA$Zea_mays)
  #filter for only my HOGs_list
  SFA_bound <- SFA_bound[SFA_bound$HOG %in% HOGs_list,]
  
  n <- 1
  for(id in all_GO_ID)
  {
    h <- SFA_bound[grepl(id, SFA_bound$`Gene-Ontology-Term`),]
    all_HOGs_with_GO_count[n] <- length(unique(h$HOG))
    all_HOGs_with_GO[n] <- paste(unique(h$HOG), collapse = " ")
    n <- n + 1
  }
  
  #combine to dataframe and sort along Ontology
  df <- data.frame(GOID = all_GO_ID, GOTERM = all_GO_Terms,
                   ONTOLOGY = all_GO_Ontology,
                   HOGs_with_GO_COUNT = all_HOGs_with_GO_count,
                   HOGs_with_GO = all_HOGs_with_GO,
                   TOTAL_COUNT = all_GO_total_count,
                   Ta_COUNT = Ta_GO_count, Os_COUNT = Os_GO_count, Ms_COUNT = Mi_GO_count,
                   Sb_COUNT = SO_GO_count, Zm_COUNT = Zm_GO_count,
                   DEFINITION = all_GO_Definition)
  
  df <- df[order(df$ONTOLOGY),]
  rownames(df) <- NULL
  return(df)
}
