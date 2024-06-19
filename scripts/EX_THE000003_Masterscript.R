This it the Masterscript with every script used in EX_THE000003. This is made to represent the workflow order, and is not meant to be a working R script
These are the filtering scripts used to filter the A2TEA RData returned by the A2TEA Workflow in EX_THE000003
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
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

# Zm expanded, All Zm DEG, no DEG in C3

ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs = c()
ZmExp_AtLeastOneZmDEG_NoDEGC3_genes = c()
AtLeastOneZmDEG_NoDEGC3_hogs = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for everything that has no C3 plant significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        # #check for NA in Zm Log2FoldChanges
        # if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
        # {
          #check for Zm that have at least one significant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
          {
            AtLeastOneZmDEG_NoDEGC3_hogs <- append(AtLeastOneZmDEG_NoDEGC3_hogs, hog)
            #filter for HOGs, that have Zm expanded
            if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs <- append(ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs, hog)
              ZmExp_AtLeastOneZmDEG_NoDEGC3_genes <- append(ZmExp_AtLeastOneZmDEG_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }

        # }
      }
    }
    
  }
}

write.csv(ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
write.csv(ZmExp_AtLeastOneZmDEG_NoDEGC3_genes, ".\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_genes")
write.csv(extract_GO(ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs), ".\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_unique_GO.csv")
write.csv(AtLeastOneZmDEG_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\AtLeastOneZmDEG_NoDEGC3_hogs.csv")



# Zm expanded, All Zm Down, no DEG in C3

ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs = c()
ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_genes = c()
AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for everything that has no C3 plant significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        # #check for NA in Zm Log2FoldChanges
        # if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
        # {
          #check for Zm that have at least one significant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
          {
            #no significant downregulated
            if(sum(genes[genes$gene == "Zm" & genes$significant == "yes",]$log2FoldChange > 0) == 0)
            {
              AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- append(AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs, hog)
              #filter for HOGs, that have Zm expanded
              if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- append(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs, hog)
                ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_genes <- append(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          # }
          
        }
      }
    }
    
  }
}

write.csv(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
write.csv(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_genes, ".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_genes")
write.csv(extract_GO(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs), ".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_unique_GO.csv")
write.csv(AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")




# Zm expanded, All Zm Up, no DEG in C3

ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs = c()
ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_genes = c()
AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for everything that has no C3 plant significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        # #check for NA in Zm Log2FoldChanges
        # if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
        # {
          #check for Zm that have at least one significant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
          {
            #no significant upregulated
            if(sum(genes[genes$gene == "Zm" & genes$significant == "yes",]$log2FoldChange < 0) == 0)
            {
              AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- append(AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, hog)
              #filter for HOGs, that have Zm expanded
              if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- append(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, hog)
                ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_genes <- append(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
            
          # }
          
        }
      }
    }
    
  }
}

write.csv(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
write.csv(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_genes, ".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_genes")
write.csv(extract_GO(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs), ".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_unique_GO.csv")
write.csv(AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")


#Are there HOGs where no C4 is regulated, but expanded and all C3 are regulated?

# AtLeastOneC4Exp_NoDEGC4_AllC3DEG C3 plants are DEGregulated and at least one C4 is expanded

AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs = c()
AtLeastOneC4Exp_NoDEGC4_AllC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for HOGs, that have more than one Zm Mi or SO/expanded
    if(count["Zm"] > 1 || count["Mi"] > 2 || count["SO"] > 1)
    {
      #filter for HOGs which have at least one Os and three Ta
      if(count["Tr"] >= 3 && count["Os"] >= 1)
      {
        #filter for everything that has no C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              #check for NA in Tr and Os Log2FoldChanges
              na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
              na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
              if(na_sum_Tr == 0 && na_sum_Os == 0)
              {
                #check for Tr, Os that are all not nonsignificant
                notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
                notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
                if(notsig_Tr == 0 && notsig_Os == 0)
                {
                  print("hit")
                  AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs <- append(AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs, hog)
                  AtLeastOneC4Exp_NoDEGC4_AllC3DEG_genes <- append(AtLeastOneC4Exp_NoDEGC4_AllC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                }
              }
            }
          }
        }
      }
    }
  }
}

write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3DEG_genes, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_genes.csv")
write.csv(extract_GO(AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs), ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Down\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_unique_GO.csv")


# AtLeastOneC4Exp_NoDEGC4_AllC3Up C3 plants are upregulated and at least one C4 is expanded

AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs = c()
AtLeastOneC4Exp_NoDEGC4_AllC3Up_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for HOGs, that have more than one Zm Mi or SO/expanded
    if(count["Zm"] > 1 || count["Mi"] > 2 || count["SO"] > 1)
    {
      #filter for HOGs which have at least one Os and three Ta
      if(count["Tr"] >= 3 && count["Os"] >= 1)
      {
        #filter for everything that has no C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              #check for NA in Tr and Os Log2FoldChanges
              na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
              na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
              if(na_sum_Tr == 0 && na_sum_Os == 0)
              {
                #filter for all Zm, Mi and SO, that are not downregulated
                down_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange < 0)
                down_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange < 0)
                if(down_Tr == 0 && down_Os == 0)
                {
                  #check for Tr, Os that are all not nonsignificant
                  notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
                  notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
                  if(notsig_Tr == 0 && notsig_Os == 0)
                  {
                    print("hit")
                    AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs <- append(AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs, hog)
                    AtLeastOneC4Exp_NoDEGC4_AllC3Up_genes <- append(AtLeastOneC4Exp_NoDEGC4_AllC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Up\\AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs.csv")
write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3Up_genes, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Up\\AtLeastOneC4Exp_NoDEGC4_AllC3Up_genes.csv")
write.csv(extract_GO(AtLeastOneC4Exp_NoDEGC4_AllC3Up_hogs), ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Down\\AtLeastOneC4Exp_NoDEGC4_AllC3Up_unique_GO.csv")



# AtLeastOneC4Exp_NoDEGC4_AllC3Down C3 plants are Downregulated and at least one C4 is expanded

AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs = c()
AtLeastOneC4Exp_NoDEGC4_AllC3Down_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for HOGs, that have more than one Zm Mi or SO/expanded
    if(count["Zm"] > 1 || count["Mi"] > 2 || count["SO"] > 1)
    {
      #filter for HOGs which have at least one Os and three Ta
      if(count["Tr"] >= 3 && count["Os"] >= 1)
      {
        #filter for everything that has no C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              #check for NA in Tr and Os Log2FoldChanges
              na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
              na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
              if(na_sum_Tr == 0 && na_sum_Os == 0)
              {
                #filter for all Zm, Mi and SO, that are not downregulated
                up_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange > 0)
                up_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange > 0)
                if(up_Tr == 0 && up_Os == 0)
                {
                  #check for Tr, Os that are all not nonsignificant
                  notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
                  notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
                  if(notsig_Tr == 0 && notsig_Os == 0)
                  {
                    print("hit")
                    AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs <- append(AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs, hog)
                    AtLeastOneC4Exp_NoDEGC4_AllC3Down_genes <- append(AtLeastOneC4Exp_NoDEGC4_AllC3Down_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Down\\AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs.csv")
write.csv(AtLeastOneC4Exp_NoDEGC4_AllC3Down_genes, ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Down\\AtLeastOneC4Exp_NoDEGC4_AllC3Down_genes.csv")
write.csv(extract_GO(AtLeastOneC4Exp_NoDEGC4_AllC3Down_hogs), ".\\AtLeastOneC4Exp_NoDEGC4_AllC3Down\\AtLeastOneC4Exp_NoDEGC4_AllC3Down_unique_GO.csv")




# AllC4Exp_NoDEGC4_AllC3DEG

AllC4Exp_NoDEGC4_AllC3DEG_hogs = c()
AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs = c()

AllC4Exp_NoDEGC4_AllC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
      #filter for everything that has no C4 plant significantly regulated
      if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
      {
        if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {

              #check for Zm, Mi and SO that are all not nonsignificant
              notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
              notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
              if(notsig_Tr == 0 && notsig_Os == 0)
              {
                AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs <- append(AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs, hog)
                
                #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                {
                  print("hit")
                  AllC4Exp_NoDEGC4_AllC3DEG_hogs <- append(AllC4Exp_NoDEGC4_AllC3DEG_hogs, hog)
                  AllC4Exp_NoDEGC4_AllC3DEG_genes <- append(AllC4Exp_NoDEGC4_AllC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                }
              }
            
            
          }
        }
      }
    

  }
}

write.csv(AllC4Exp_NoDEGC4_AllC3DEG_hogs, ".\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_hogs, ".\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_genes, ".\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_genes.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AllC3DEG_hogs), ".\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_unique_GO.csv")




# AllC4Exp_NoDEGC4_AllC3Up

AllC4Exp_NoDEGC4_AllC3Up_hogs = c()

AllC4Exp_NoDEGC4_AllC3Up_genes = c()

NoDEGC4_AllC3Up_hogs = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

        #filter for everything that has no C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              
              lfc_Tr <- genes[genes$gene == "Tr",]$log2FoldChange
              lfc_Tr[is.na(lfc_Tr)] <- 0
              lfc_Os <- genes[genes$gene == "Os",]$log2FoldChange
              lfc_Os[is.na(lfc_Os)] <- 0
              #filter for all Tr, Os, that are not downregulated
              down_Tr <- sum(lfc_Tr < 0)
              down_Os <- sum(lfc_Os < 0)
              #print(paste(down_Tr, down_Os))
              if(down_Tr == 0 && down_Os == 0)
              {
                #check for Zm, Mi and SO that are all not nonsignificant
                notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
                notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
                if(notsig_Tr == 0 && notsig_Os == 0)
                {
                  
                  #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                  if(hog %in% HOG_level_list$hypothesis_12$HOG)
                  {
                    
                    NoDEGC4_AllC3Up_hogs <- append(NoDEGC4_AllC3Up_hogs, hog)
                    
                    if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                    {
                      print("hit")
                      AllC4Exp_NoDEGC4_AllC3Up_hogs <- append(AllC4Exp_NoDEGC4_AllC3Up_hogs, hog)
                      AllC4Exp_NoDEGC4_AllC3Up_genes <- append(AllC4Exp_NoDEGC4_AllC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                
                    }
                  }
                }
              }
              
            }
          }
        }
        


  }
}

write.csv(AllC4Exp_NoDEGC4_AllC3Up_hogs, ".\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AllC3Up_genes, ".\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_genes.csv")
write.csv(NoDEGC4_AllC3Up_hogs, ".\\AllC4Exp_NoDEGC4_AllC3Up\\NoDEGC4_AllC3Up_hogs.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AllC3Up_hogs), ".\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_unique_GO.csv")




# AllC4Exp_NoDEGC4_AllC3Down

AllC4Exp_NoDEGC4_AllC3Down_hogs = c()
NoDEGC4_AllC3Down_hogs = c()

AllC4Exp_NoDEGC4_AllC3Down_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    
        #filter for everything that has no C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              
              lfc_Tr <- genes[genes$gene == "Tr",]$log2FoldChange
              lfc_Tr[is.na(lfc_Tr)] <- 0
              lfc_Os <- genes[genes$gene == "Os",]$log2FoldChange
              lfc_Os[is.na(lfc_Os)] <- 0
              #filter for all Tr, Os, that are not downregulated
              up_Tr <- sum(lfc_Tr > 0)
              up_Os <- sum(lfc_Os > 0)
              
              if(up_Tr == 0 && up_Os == 0)
              {
                #check for Zm, Mi and SO that are all not nonsignificant
                notsig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "no")
                notsig_Os <- sum(genes[genes$gene == "Os",]$significant == "no")
                if(notsig_Tr == 0 && notsig_Os == 0)
                {
                  #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                  if(hog %in% HOG_level_list$hypothesis_12$HOG)
                  {
                    NoDEGC4_AllC3Down_hogs <- append(NoDEGC4_AllC3Down_hogs, hog)
                    #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                    if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                    {
                      print("hit")
                      AllC4Exp_NoDEGC4_AllC3Down_hogs <- append(AllC4Exp_NoDEGC4_AllC3Down_hogs, hog)
                      AllC4Exp_NoDEGC4_AllC3Down_genes <- append(AllC4Exp_NoDEGC4_AllC3Down_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                
                    }
                  }
                }
              }
              
            }
          }
        }
        

  }
}

write.csv(AllC4Exp_NoDEGC4_AllC3Down_hogs, ".\\AllC4Exp_NoDEGC4_AllC3Down\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
write.csv(NoDEGC4_AllC3Down_hogs, ".\\AllC4Exp_NoDEGC4_AllC3Down\\NoDEGC4_AllC3Down_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AllC3Down_genes, ".\\AllC4Exp_NoDEGC4_AllC3Down\\AllC4Exp_NoDEGC4_AllC3Down_genes.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AllC3Down_hogs), ".\\AllC4Exp_NoDEGC4_AllC3Down\\AllC4Exp_NoDEGC4_AllC3Down_unique_GO.csv")











# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG

AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs = c()
NoDEGC4_AtLeastOneOfEachC3DEG_hogs = c()
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

      
      #filter for everything that has no C4 plant significantly regulated
      if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
      {
        if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            # #check for NA in Tr and Os Log2FoldChanges
            # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
            # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
            # if(na_sum_Tr == 0 && na_sum_Os == 0)
            # {
              #check for Tr and Os which have at least one significant DEG for each of the species
              sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
              sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
              if(sig_Tr > 0 && sig_Os > 0)
              {
                NoDEGC4_AtLeastOneOfEachC3DEG_hogs <- append(NoDEGC4_AtLeastOneOfEachC3DEG_hogs, hog)
                
                #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                {
                  print("hit")
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs, hog)
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_genes <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                }
              }
              
            # }
          }
        }
      }
    
    
  }
}

write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(NoDEGC4_AtLeastOneOfEachC3DEG_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG\\NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_genes, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_genes.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs), ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_unique_GO.csv")




# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp

AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

      
      #filter for everything that has no C4 plant significantly regulated
      if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
      {
        if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            # #check for NA in Tr and Os Log2FoldChanges
            # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
            # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
            # if(na_sum_Tr == 0 && na_sum_Os == 0)
            # {
            #check for Tr and Os which have at least one significant DEG for each of the species
            sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
            sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
            if(sig_Tr > 0 && sig_Os > 0)
            {

              #filter for all Tr, Os, that are not downregulated
              down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange < 0)
              down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange < 0)
              if(down_Tr == 0 && down_Os == 0)
              {
                NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
                
                #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                {
                  print("hit")
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_genes <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                }
                
              }
              
            }
            
            # }
          }
        }
      }
    
    
  }
}

write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_genes, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_genes.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs), ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_unique_GO.csv")




# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown

AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

      #filter for everything that has no C4 plant significantly regulated
      if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
      {
        if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            # #check for NA in Tr and Os Log2FoldChanges
            # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
            # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
            # if(na_sum_Tr == 0 && na_sum_Os == 0)
            # {
            #check for Tr and Os which have at least one significant DEG for each of the species
            sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
            sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
            if(sig_Tr > 0 && sig_Os > 0)
            {
              
              #filter for all Tr, Os, that are not downregulated
              up_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange > 0)
              up_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange > 0)
              if(up_Tr == 0 && up_Os == 0)
              {
                NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
                
                #filter for HOGs, that have more than one Zm Mi or SO/all expanded
                if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
                {
                  print("hit")
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
                  AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_genes <- append(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                }
                
              }
              
            }
            
            # }
          }
        }
      
    }
    
  }
}

write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_genes, ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_genes.csv")
write.csv(extract_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs), ".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_unique_GO.csv")




setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")


# is something regulated in C4 which is not regulated in C3?


#No DEG C3 , All C4 expanded, At least one DEG per C4 plant

AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs = c()

AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    #filter for everything that has no C3 plants significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        #filter for HOGs, that have more than one Zm Mi and SO/all C4 expanded
        if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
        {
          
          #filter for everything that has at least one regulated gene per species of C4 plants
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
          {
            if(sum(genes[genes$gene == "SO",]$significant == "yes") > 0)
            {
              if(sum(genes[genes$gene == "Mi",]$significant == "yes") > 0)
              {
                print("hit")
                AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs <- append(AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs, hog)
                AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_genes <- append(AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
                
              }
            }
          }
        }
      }
    }
  }
}

write.csv(AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs, ".\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs.csv")
write.csv(AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_genes, ".\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_genes.csv")
#write.csv(extract_GO(AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_hogs), ".\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3\\AllC4Exp_AtLeastOneDEGPerC4_NoDEGC3_unique_GO.csv")
# the list is empty


#No DEG C3 , All C4 expanded, DEG in Any C4 Plant

AllC4Exp_DEGAnyC4_NoDEGC3_hogs = c()
AllC4Exp_DEGAnyC4_NoDEGC3_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    #filter for everything that has no C3 plants significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        #filter for HOGs, that have more than one Zm Mi and SO /all C4 expanded
        if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
        {
          
          #filter for everything that has at least one regulated gene in any C4 plant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0
             || sum(genes[genes$gene == "Mi",]$significant == "yes") > 0
             || sum(genes[genes$gene == "SO",]$significant == "yes") > 0)
          {
            print("hit")
            AllC4Exp_DEGAnyC4_NoDEGC3_hogs <- append(AllC4Exp_DEGAnyC4_NoDEGC3_hogs, hog)
            AllC4Exp_DEGAnyC4_NoDEGC3_genes <- append(AllC4Exp_DEGAnyC4_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
          }
        }
      }
    }
  }
}

write.csv(AllC4Exp_DEGAnyC4_NoDEGC3_hogs, ".\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
write.csv(AllC4Exp_DEGAnyC4_NoDEGC3_genes, ".\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_genes.csv")
write.csv(extract_GO(AllC4Exp_DEGAnyC4_NoDEGC3_hogs), ".\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_unique_GO.csv")


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
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

conserved_hogs = c()
AllC4Exp_hogs = c()
ZmExp_hogs = c()
MsExp_hogs = c()
SbExp_hogs = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    conserved_hogs <- append(conserved_hogs, hog)
    
    #Ms Expanded
    if(hog %in% HOG_level_list$hypothesis_3$HOG)
    {
      if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
      {
        MsExp_hogs <- append(MsExp_hogs, hog)
      }
    }
    
    #Sb Expanded
    if(hog %in% HOG_level_list$hypothesis_6$HOG)
    {
      if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
      {
        SbExp_hogs <- append(SbExp_hogs, hog)
      }
    }
    
    #Zm Expanded
    if(hog %in% HOG_level_list$hypothesis_9$HOG)
    {
      if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
      {
        ZmExp_hogs <- append(ZmExp_hogs, hog)
      }
    }
    
    if(hog %in% HOG_level_list$hypothesis_9$HOG)
    {
      #AllC4 Expanded
      if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog,]$expansion == "yes")
      {
        AllC4Exp_hogs <- append(AllC4Exp_hogs, hog)
      }
    }
    
    
    
    
  }   
}

write.csv(conserved_hogs, ".\\conserved_hogs.csv")
write.csv(AllC4Exp_hogs, ".\\AllC4Exp_hogs.csv")
write.csv(ZmExp_hogs, ".\\ZmExp_hogs.csv")
write.csv(MsExp_hogs, ".\\MsExp_hogs.csv")
write.csv(SbExp_hogs, ".\\SbExp_hogs.csv")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# Mi expanded, No Mi DEG, At least One Of Each C3 DEG

MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs = c()
NoDEGMi_AtLeastOneOfEachC3DEG_hogs = c()
MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      #filter for everything where At least one of each C3 species is sig
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
          
        #check for Mi that are all nonsignificant
        if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
        {
          NoDEGMi_AtLeastOneOfEachC3DEG_hogs <- append(NoDEGMi_AtLeastOneOfEachC3DEG_hogs, hog)
          #filter for HOGs, that have Mi expanded
          if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
          {
            print("hit")
            MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs, hog)
            MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_genes <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
          }
        }
        
      }
      
    # }
    
  }
}

write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(NoDEGMi_AtLeastOneOfEachC3DEG_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG\\NoDEGMi_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_genes, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_hogs), ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_unique_GO.csv")



# Mi expanded, No Mi DEG, At least One Of Each C3 Down

MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange > 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange > 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for Mi that are all nonsignificant
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
            #filter for HOGs, that have Mi expanded
            if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
              MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_genes <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown\\NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_genes, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_hogs), ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyDown_unique_GO.csv")




# Mi expanded, No Mi DEG, At least One Of Each C3 Up

MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
# 
#     #check for NA in Tr and Os Log2FoldChanges
#     na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
#     na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
#     if(na_sum_Tr == 0 && na_sum_Os == 0)
#     {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange < 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange < 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for Mi that are all nonsignificant
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
            #filter for HOGs, that have Mi expanded
            if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
              MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_genes <- append(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp\\NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_genes, ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_hogs), ".\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp\\MiExp_NoDEGMi_AtLeastOneOfEachC3DEG_OnlyUp_unique_GO.csv")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# Mi expanded, No Mi DEG, All C3 DEG

MiExp_NoDEGMi_AllC3DEG_hogs = c()
NoDEGMi_AllC3DEG_hogs = c()
MiExp_NoDEGMi_AllC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
      #filter for everything that has no C3 plant not significantly regulated
      if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
      {
        if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
        {
          #check for Mi that are all nonsignificant
          if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
          {
            NoDEGMi_AllC3DEG_hogs <- append(NoDEGMi_AllC3DEG_hogs, hog)
            #filter for HOGs, that have Mi expanded
            if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              MiExp_NoDEGMi_AllC3DEG_hogs <- append(MiExp_NoDEGMi_AllC3DEG_hogs, hog)
              MiExp_NoDEGMi_AllC3DEG_genes <- append(MiExp_NoDEGMi_AllC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            
            }
          }
        }
      }
    
  }
}

write.csv(MiExp_NoDEGMi_AllC3DEG_hogs, ".\\MiExp_NoDEGMi_AllC3DEG\\MiExp_NoDEGMi_AllC3DEG_hogs.csv")
write.csv(NoDEGMi_AllC3DEG_hogs, ".\\MiExp_NoDEGMi_AllC3DEG\\NoDEGMi_AllC3DEG_hogs.csv")
write.csv(MiExp_NoDEGMi_AllC3DEG_genes, ".\\MiExp_NoDEGMi_AllC3DEG\\MiExp_NoDEGMi_AllC3DEG_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AllC3DEG_hogs), ".\\MiExp_NoDEGMi_AllC3DEG\\MiExp_NoDEGMi_AllC3DEG_unique_GO.csv")



# Mi expanded, No Mi DEG, All C3 Up

MiExp_NoDEGMi_AllC3Up_hogs = c()
NoDEGMi_AllC3Up_hogs = c()
MiExp_NoDEGMi_AllC3Up_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          down_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange < 0)
          down_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange < 0)
          if(down_Tr == 0 && down_Os == 0)
          {
        
            #check for Mi that are all nonsignificant
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              NoDEGMi_AllC3Up_hogs <- append(NoDEGMi_AllC3Up_hogs, hog)
              #filter for HOGs, that have Mi expanded
              if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                MiExp_NoDEGMi_AllC3Up_hogs <- append(MiExp_NoDEGMi_AllC3Up_hogs, hog)
                MiExp_NoDEGMi_AllC3Up_genes <- append(MiExp_NoDEGMi_AllC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
    
  }
}

write.csv(MiExp_NoDEGMi_AllC3Up_hogs, ".\\MiExp_NoDEGMi_AllC3Up\\MiExp_NoDEGMi_AllC3Up_hogs.csv")
write.csv(NoDEGMi_AllC3Up_hogs, ".\\MiExp_NoDEGMi_AllC3Up\\NoDEGMi_AllC3Up_hogs.csv")
write.csv(MiExp_NoDEGMi_AllC3Up_genes, ".\\MiExp_NoDEGMi_AllC3Up\\MiExp_NoDEGMi_AllC3Up_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AllC3Up_hogs), ".\\MiExp_NoDEGMi_AllC3Up\\MiExp_NoDEGMi_AllC3Up_unique_GO.csv")



# Mi expanded, No Mi Down, All C3 Down

MiExp_NoDEGMi_AllC3Down_hogs = c()
NoDEGMi_AllC3Down_hogs = c()
MiExp_NoDEGMi_AllC3Down_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          up_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange > 0)
          up_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange > 0)
          if(up_Tr == 0 && up_Os == 0)
          {
 
            #check for Mi that are all nonsignificant
            if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
            {
              NoDEGMi_AllC3Down_hogs <- append(NoDEGMi_AllC3Down_hogs, hog)
              #filter for HOGs, that have Mi expanded
              if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                MiExp_NoDEGMi_AllC3Down_hogs <- append(MiExp_NoDEGMi_AllC3Down_hogs, hog)
                MiExp_NoDEGMi_AllC3Down_genes <- append(MiExp_NoDEGMi_AllC3Down_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
  
  }
}

write.csv(MiExp_NoDEGMi_AllC3Down_hogs, ".\\MiExp_NoDEGMi_AllC3Down\\MiExp_NoDEGMi_AllC3Down_hogs.csv")
write.csv(NoDEGMi_AllC3Down_hogs, ".\\MiExp_NoDEGMi_AllC3Down\\NoDEGMi_AllC3Down_hogs.csv")
write.csv(MiExp_NoDEGMi_AllC3Down_genes, ".\\MiExp_NoDEGMi_AllC3Down\\MiExp_NoDEGMi_AllC3Down_genes.csv")
write.csv(extract_GO(MiExp_NoDEGMi_AllC3Down_hogs), ".\\MiExp_NoDEGMi_AllC3Down\\MiExp_NoDEGMi_AllC3Down_unique_GO.csv")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# SO expanded, No SO DEG, At least One Of Each C3 DEG

SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs = c()
NoDEGSO_AtLeastOneOfEachC3DEG_hogs = c()
SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      #filter for everything where At least one of each C3 species is sig
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
          
        #check for SO that are all nonsignificant
        if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
        {
          NoDEGSO_AtLeastOneOfEachC3DEG_hogs <- append(NoDEGSO_AtLeastOneOfEachC3DEG_hogs, hog)
          #filter for HOGs, that have SO expanded
          if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
          {
            print("hit")
            SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs, hog)
            SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_genes <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
          }
        }
        
      }
      
    # }
    
  }
}

write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(NoDEGSO_AtLeastOneOfEachC3DEG_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG\\NoDEGSO_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_genes, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_hogs), ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_unique_GO.csv")



# SO expanded, No SO DEG, At least One Of Each C3 Down

SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange > 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange > 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for SO that are all nonsignificant
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
            #filter for HOGs, that have SO expanded
            if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
              SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_genes <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown\\NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_genes, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_hogs), ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyDown_unique_GO.csv")




# SO expanded, No SO DEG, At least One Of Each C3 Up

SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
# 
#     #check for NA in Tr and Os Log2FoldChanges
#     na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
#     na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
#     if(na_sum_Tr == 0 && na_sum_Os == 0)
#     {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange < 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange < 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for SO that are all nonsignificant
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
            #filter for HOGs, that have SO expanded
            if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
              SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_genes <- append(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp\\NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_genes, ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_hogs), ".\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp\\SOExp_NoDEGSO_AtLeastOneOfEachC3DEG_OnlyUp_unique_GO.csv")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# SO expanded, No SO DEG, All C3 DEG

SOExp_NoDEGSO_AllC3DEG_hogs = c()
NoDEGSO_AllC3DEG_hogs = c()
SOExp_NoDEGSO_AllC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
      #filter for everything that has no C3 plant not significantly regulated
      if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
      {
        if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
        {
          #check for SO that are all nonsignificant
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            NoDEGSO_AllC3DEG_hogs <- append(NoDEGSO_AllC3DEG_hogs, hog)
            #filter for HOGs, that have SO expanded
            if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              SOExp_NoDEGSO_AllC3DEG_hogs <- append(SOExp_NoDEGSO_AllC3DEG_hogs, hog)
              SOExp_NoDEGSO_AllC3DEG_genes <- append(SOExp_NoDEGSO_AllC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            
            }
          }
        }
      }
    
  }
}

write.csv(SOExp_NoDEGSO_AllC3DEG_hogs, ".\\SOExp_NoDEGSO_AllC3DEG\\SOExp_NoDEGSO_AllC3DEG_hogs.csv")
write.csv(NoDEGSO_AllC3DEG_hogs, ".\\SOExp_NoDEGSO_AllC3DEG\\NoDEGSO_AllC3DEG_hogs.csv")
write.csv(SOExp_NoDEGSO_AllC3DEG_genes, ".\\SOExp_NoDEGSO_AllC3DEG\\SOExp_NoDEGSO_AllC3DEG_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AllC3DEG_hogs), ".\\SOExp_NoDEGSO_AllC3DEG\\SOExp_NoDEGSO_AllC3DEG_unique_GO.csv")



# SO expanded, No SO DEG, All C3 Up

SOExp_NoDEGSO_AllC3Up_hogs = c()
NoDEGSO_AllC3Up_hogs = c()
SOExp_NoDEGSO_AllC3Up_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          down_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange < 0)
          down_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange < 0)
          if(down_Tr == 0 && down_Os == 0)
          {
        
            #check for SO that are all nonsignificant
            if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
            {
              NoDEGSO_AllC3Up_hogs <- append(NoDEGSO_AllC3Up_hogs, hog)
              #filter for HOGs, that have SO expanded
              if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                SOExp_NoDEGSO_AllC3Up_hogs <- append(SOExp_NoDEGSO_AllC3Up_hogs, hog)
                SOExp_NoDEGSO_AllC3Up_genes <- append(SOExp_NoDEGSO_AllC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
    
  }
}

write.csv(SOExp_NoDEGSO_AllC3Up_hogs, ".\\SOExp_NoDEGSO_AllC3Up\\SOExp_NoDEGSO_AllC3Up_hogs.csv")
write.csv(NoDEGSO_AllC3Up_hogs, ".\\SOExp_NoDEGSO_AllC3Up\\NoDEGSO_AllC3Up_hogs.csv")
write.csv(SOExp_NoDEGSO_AllC3Up_genes, ".\\SOExp_NoDEGSO_AllC3Up\\SOExp_NoDEGSO_AllC3Up_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AllC3Up_hogs), ".\\SOExp_NoDEGSO_AllC3Up\\SOExp_NoDEGSO_AllC3Up_unique_GO.csv")



# SO expanded, No SO Down, All C3 Down

SOExp_NoDEGSO_AllC3Down_hogs = c()
NoDEGSO_AllC3Down_hogs = c()
SOExp_NoDEGSO_AllC3Down_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          up_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange > 0)
          up_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange > 0)
          if(up_Tr == 0 && up_Os == 0)
          {
 
            #check for SO that are all nonsignificant
            if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
            {
              NoDEGSO_AllC3Down_hogs <- append(NoDEGSO_AllC3Down_hogs, hog)
              #filter for HOGs, that have SO expanded
              if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                SOExp_NoDEGSO_AllC3Down_hogs <- append(SOExp_NoDEGSO_AllC3Down_hogs, hog)
                SOExp_NoDEGSO_AllC3Down_genes <- append(SOExp_NoDEGSO_AllC3Down_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
  
  }
}

write.csv(SOExp_NoDEGSO_AllC3Down_hogs, ".\\SOExp_NoDEGSO_AllC3Down\\SOExp_NoDEGSO_AllC3Down_hogs.csv")
write.csv(NoDEGSO_AllC3Down_hogs, ".\\SOExp_NoDEGSO_AllC3Down\\NoDEGSO_AllC3Down_hogs.csv")
write.csv(SOExp_NoDEGSO_AllC3Down_genes, ".\\SOExp_NoDEGSO_AllC3Down\\SOExp_NoDEGSO_AllC3Down_genes.csv")
write.csv(extract_GO(SOExp_NoDEGSO_AllC3Down_hogs), ".\\SOExp_NoDEGSO_AllC3Down\\SOExp_NoDEGSO_AllC3Down_unique_GO.csv")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# Zm expanded, No Zm DEG, At least One Of Each C3 DEG

ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs = c()
NoDEGZm_AtLeastOneOfEachC3DEG_hogs = c()
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      #filter for everything where At least one of each C3 species is sig
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
          
        #check for Zm that are all nonsignificant
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
        {
          NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- append(NoDEGZm_AtLeastOneOfEachC3DEG_hogs, hog)
          #filter for HOGs, that have Zm expanded
          if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
          {
            print("hit")
            ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs, hog)
            ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_genes <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
          }
        }
        
      }
      
    # }
    
  }
}

write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(NoDEGZm_AtLeastOneOfEachC3DEG_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_genes, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs), ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_unique_GO.csv")



# Zm expanded, No Zm DEG, At least One Of Each C3 Down

ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs = c()
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    
    # #check for NA in Tr and Os Log2FoldChanges
    # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
    # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
    # if(na_sum_Tr == 0 && na_sum_Os == 0)
    # {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange > 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange > 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for Zm that are all nonsignificant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
          {
            NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
            #filter for HOGs, that have Zm expanded
            if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs, hog)
              ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_genes <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_genes, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs), ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_unique_GO.csv")




# Zm expanded, No Zm DEG, At least One Of Each C3 Up

ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs = c()
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
# 
#     #check for NA in Tr and Os Log2FoldChanges
#     na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
#     na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
#     if(na_sum_Tr == 0 && na_sum_Os == 0)
#     {
      sig_Tr <- sum(genes[genes$gene == "Tr",]$significant == "yes")
      sig_Os <- sum(genes[genes$gene == "Os",]$significant == "yes")
      if(sig_Tr > 0 && sig_Os > 0)
      {
        #filter for all Tr, Os, that are not downregulated
        down_Tr <- sum(genes[genes$gene == "Tr" & genes$significant == "yes",]$log2FoldChange < 0)
        down_Os <- sum(genes[genes$gene == "Os"  & genes$significant == "yes",]$log2FoldChange < 0)
        if(down_Tr == 0 && down_Os == 0)
        {
          
          #check for Zm that are all nonsignificant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
          {
            NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
            #filter for HOGs, that have Zm expanded
            if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs, hog)
              ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_genes <- append(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
      
    # }
    
  }
}

write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_genes, ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs), ".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_unique_GO.csv")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")



# Zm expanded, No Zm DEG, All C3 DEG

ZmExp_NoDEGZm_AllC3DEG_hogs = c()
NoDEGZm_AllC3DEG_hogs = c()
ZmExp_NoDEGZm_AllC3DEG_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
      #filter for everything that has no C3 plant not significantly regulated
      if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
      {
        if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
        {
          #check for Zm that are all nonsignificant
          if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
          {
            NoDEGZm_AllC3DEG_hogs <- append(NoDEGZm_AllC3DEG_hogs, hog)
            #filter for HOGs, that have Zm expanded
            if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
            {
              print("hit")
              ZmExp_NoDEGZm_AllC3DEG_hogs <- append(ZmExp_NoDEGZm_AllC3DEG_hogs, hog)
              ZmExp_NoDEGZm_AllC3DEG_genes <- append(ZmExp_NoDEGZm_AllC3DEG_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            
            }
          }
        }
      }
    
  }
}

write.csv(ZmExp_NoDEGZm_AllC3DEG_hogs, ".\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
write.csv(NoDEGZm_AllC3DEG_hogs, ".\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs.csv")
write.csv(ZmExp_NoDEGZm_AllC3DEG_genes, ".\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AllC3DEG_hogs), ".\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_unique_GO.csv")



# Zm expanded, No Zm DEG, All C3 Up

ZmExp_NoDEGZm_AllC3Up_hogs = c()
NoDEGZm_AllC3Up_hogs = c()
ZmExp_NoDEGZm_AllC3Up_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          down_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange < 0)
          down_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange < 0)
          if(down_Tr == 0 && down_Os == 0)
          {
        
            #check for Zm that are all nonsignificant
            if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
            {
              NoDEGZm_AllC3Up_hogs <- append(NoDEGZm_AllC3Up_hogs, hog)
              #filter for HOGs, that have Zm expanded
              if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                ZmExp_NoDEGZm_AllC3Up_hogs <- append(ZmExp_NoDEGZm_AllC3Up_hogs, hog)
                ZmExp_NoDEGZm_AllC3Up_genes <- append(ZmExp_NoDEGZm_AllC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
    
  }
}

write.csv(ZmExp_NoDEGZm_AllC3Up_hogs, ".\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
write.csv(NoDEGZm_AllC3Up_hogs, ".\\ZmExp_NoDEGZm_AllC3Up\\NoDEGZm_AllC3Up_hogs.csv")
write.csv(ZmExp_NoDEGZm_AllC3Up_genes, ".\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AllC3Up_hogs), ".\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_unique_GO.csv")



# Zm expanded, No Zm Down, All C3 Down

ZmExp_NoDEGZm_AllC3Down_hogs = c()
NoDEGZm_AllC3Down_hogs = c()
ZmExp_NoDEGZm_AllC3Down_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant not significantly regulated
    if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
      {
        # #check for NA in Tr and Os Log2FoldChanges
        # na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
        # na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
        # if(na_sum_Tr == 0 && na_sum_Os == 0)
        # {
          #filter for all Tr, Os, that are not downregulated
          up_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange > 0)
          up_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange > 0)
          if(up_Tr == 0 && up_Os == 0)
          {
 
            #check for Zm that are all nonsignificant
            if(sum(genes[genes$gene == "Zm",]$significant == "yes") == 0)
            {
              NoDEGZm_AllC3Down_hogs <- append(NoDEGZm_AllC3Down_hogs, hog)
              #filter for HOGs, that have Zm expanded
              if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
              {
                print("hit")
                ZmExp_NoDEGZm_AllC3Down_hogs <- append(ZmExp_NoDEGZm_AllC3Down_hogs, hog)
                ZmExp_NoDEGZm_AllC3Down_genes <- append(ZmExp_NoDEGZm_AllC3Down_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              }
            }
          }
        # }
      }
    }
  
  }
}

write.csv(ZmExp_NoDEGZm_AllC3Down_hogs, ".\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
write.csv(NoDEGZm_AllC3Down_hogs, ".\\ZmExp_NoDEGZm_AllC3Down\\NoDEGZm_AllC3Down_hogs.csv")
write.csv(ZmExp_NoDEGZm_AllC3Down_genes, ".\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_genes.csv")
write.csv(extract_GO(ZmExp_NoDEGZm_AllC3Down_hogs), ".\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_unique_GO.csv")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

# Zm expanded, All Zm DEG, no DEG in C3

ZmExp_AllZmDEG_NoDEGC3_hogs = c()
ZmExp_AllZmDEG_NoDEGC3_genes = c()

x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    if(hog %in% HOG_level_list$hypothesis_9$HOG)
    {
      #filter for HOGs, that have Zm expanded
      if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
      {
        #filter for everything that has no C3 plant significantly regulated
        if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
          {
            # #check for NA in Zm Log2FoldChanges
            # if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
            #   
            #check for Zm that are all not nonsignificant
            if(sum(genes[genes$gene == "Zm",]$significant == "no") == 0)
            {
              print("hit")
              ZmExp_AllZmDEG_NoDEGC3_hogs <- append(ZmExp_AllZmDEG_NoDEGC3_hogs, hog)
              ZmExp_AllZmDEG_NoDEGC3_genes <- append(ZmExp_AllZmDEG_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
            }
          }
        }
      }
    }
    
  }
}

write.csv(ZmExp_AllZmDEG_NoDEGC3_hogs, ".\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
write.csv(ZmExp_AllZmDEG_NoDEGC3_genes, ".\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_genes.csv")
write.csv(extract_GO(ZmExp_AllZmDEG_NoDEGC3_hogs), ".\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_unique_GO.csv")


# Zm expanded, only any Zm DEG, no DEG in C3

ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs = c()
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs = c()
ZmExp_OnlyAnyZmDEG_NoDEGC3_genes = c()
x = 0

for(hog in all_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for conserved HOGs
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {

    #filter for everything that has no C3 plant significantly regulated and no Mi and SO
    if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
    {
      if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
      {
        if(sum(genes[genes$gene == "Mi",]$significant == "yes") == 0)
        {
          if(sum(genes[genes$gene == "SO",]$significant == "yes") == 0)
          {
            # #check for NA in Zm Log2FoldChanges
            # if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
            # {
              #check for Zm where there is at least one significant
              if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
              {
                ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs <- append(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs, hog)
                #filter for HOGs, that have Zm expanded
                if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
                {
                  print("hit")
                  ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs <- append(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs, hog)
                  ZmExp_OnlyAnyZmDEG_NoDEGC3_genes <- append(ZmExp_OnlyAnyZmDEG_NoDEGC3_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
              
                }
              }
            # }  
          }    
        }
      }
    }
    
  }
}

write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs, ".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs, ".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs.csv")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_genes, ".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_genes.csv")
write.csv(extract_GO(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs), ".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_unique_GO.csv")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
source("EX_THE000003_Filter_Source.R")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

#filter hypotheses 12


# MSZ_expanded_deg_any_exp_species_hogs ################################################

MSZ_expanded_deg_any_exp_species_hogs = c()

MSZ_expanded_deg_any_exp_species_genes = c()


x = 0
c = 0

hypo12_hogs = HOG_level_list$hypothesis_12$HOG

for(hog in hypo12_hogs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(all_hogs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  # filter for every HOGs that has at least one of each plant/ is conserved
  if(sum(c("Mi", "SO", "Zm", "Tr", "Os") %in% unique(genes$gene)) == 5)
  {
    #filter for expanded HOGs (expansion criteria is all three C4 need to be expanded)
    if(HOG_level_list$hypothesis_12[HOG_level_list$hypothesis_12$HOG == hog, ]$expansion == "yes")
    {
      #filter for HOGs which have at least one Os and one Ta
      if(count["Tr"] >= 1 && count["Os"] >= 1)
      {
        #filter for everything that has at least one C4 plant significantly regulated
        if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0 || sum(genes[genes$gene == "SO",]$significant == "yes") > 0 || sum(genes[genes$gene == "Mi",]$significant == "yes") > 0)
        {
            c = c + 1
            cat("hit ", c, "\n")
            MSZ_expanded_deg_any_exp_species_hogs <- append(MSZ_expanded_deg_any_exp_species_hogs, hog)
            MSZ_expanded_deg_any_exp_species_genes <- append(MSZ_expanded_deg_any_exp_species_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
        }
      }
    }
  }
}

write.csv(MSZ_expanded_deg_any_exp_species_hogs, ".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_hogs.csv")
write.csv(MSZ_expanded_deg_any_exp_species_genes, ".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_genes.csv")
write.csv(extract_GO(MSZ_expanded_deg_any_exp_species_hogs), ".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_unique_GO.csv")

#not necessary

# #filter for stuff that is only upregulated in C3
# 
# 
# MSZ_expanded_deg_any_exp_species_onlyallC3Up_hogs = c()
# 
# MSZ_expanded_deg_any_exp_species_onlyallC3Up_genes = c()
# 
# x = 0
# 
# for(hog in MSZ_expanded_deg_any_exp_species_hogs)
# {
#   
#   x = x + 1
#   if(x %% 1000 == 1)
#   {
#     cat(x, "of", length(all_hogs), "\n")
#   }
#   genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
#   genes$gene <- substr(genes$gene, 1, 2)
#   count <- table(genes$gene)
#   
#   #filter for everything that has all C3 plants significantly regulated
#   if(sum(genes[genes$gene == "Tr",]$significant == "no") == 0)
#   {
#     if(sum(genes[genes$gene == "Os",]$significant == "no") == 0)
#     {
#       #check for NA in Tr and Os Log2FoldChanges
#       na_sum_Tr <- sum(is.na(genes[genes$gene == "Tr",]$log2FoldChange))
#       na_sum_Os <- sum(is.na(genes[genes$gene == "Os",]$log2FoldChange))
#       if(na_sum_Tr == 0 && na_sum_Os == 0)
#       {
#         #filter for all Ta Os, that are not downregulated
#         down_Tr <- sum(genes[genes$gene == "Tr",]$log2FoldChange < 0)
#         down_Os <- sum(genes[genes$gene == "Os",]$log2FoldChange < 0)
#         if(down_Tr == 0 && down_Os == 0)
#         {
# 
#           print("hit")
#           MSZ_expanded_deg_any_exp_species_onlyallC3Up_hogs <- append(MSZ_expanded_deg_any_exp_species_onlyallC3Up_hogs, hog)
#           MSZ_expanded_deg_any_exp_species_onlyallC3Up_genes <- append(MSZ_expanded_deg_any_exp_species_onlyallC3Up_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
#           
#         }
#       }
#     }
#   }
# }
# 
# write.csv(MSZ_expanded_deg_any_exp_species_onlyallC3Up_hogs, ".\\MSZ_expanded_deg_any_exp_species_onlyallC3Up_hogs.csv")
# write.csv(MSZ_expanded_deg_any_exp_species_onlyallC3Up_genes, ".\\MSZ_expanded_deg_any_exp_species_onlyallC3Up_genes.csv")
# 
# 
# 
# #filter for stuff that is not regulated in C3 but at least once regulated in all C4
# 
# MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_hogs = c()
# 
# MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_genes = c()
# 
# x = 0
# 
# for(hog in MSZ_expanded_deg_any_exp_species_hogs)
# {
#   
#   x = x + 1
#   if(x %% 1000 == 1)
#   {
#     cat(x, "of", length(all_hogs), "\n")
#   }
#   genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
#   genes$gene <- substr(genes$gene, 1, 2)
#   count <- table(genes$gene)
#   
#   #filter for everything that has no C3 plants significantly regulated
#   if(sum(genes[genes$gene == "Tr",]$significant == "yes") == 0)
#   {
#     if(sum(genes[genes$gene == "Os",]$significant == "yes") == 0)
#     {
#       
#       #filter for everything that has at least one regulated gene per species of C4 plants
#       if(sum(genes[genes$gene == "Zm",]$significant == "yes") > 0)
#       {
#         if(sum(genes[genes$gene == "SO",]$significant == "yes") > 0)
#         {
#           if(sum(genes[genes$gene == "Tr",]$significant == "yes") > 0)
#           {
#             print("hit")
#             MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_hogs <- append(MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_hogs, hog)
#             MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_genes <- append(MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_genes, HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
#             
#           }
#         }
#       }
#     }
#   }
# }
# 
# write.csv(MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_hogs, ".\\MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_hogs.csv")
# write.csv(MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_genes, ".\\MSZ_expanded_deg_any_exp_species_noDEGC3__at_least_one_DEG_in_each_C4_genes.csv")



The script to start the A2TEA WebApp to parse throught EX_THE000003_A2TEA_finished.RData
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")
library(A2TEA.WebApp)
A2TEA_App()script to calculate go term enrichment with topGO
library(topGO)

#The analysis was done for sets with more than 20 HOGs

#load data ####################
library("DESeq2")
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(forcats)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

#make a list of conserved HOGs
conserved_HOGs <- c()

x <- 0
for(hog in unique(HOG_DE.a2tea$HOG))
{
  x <- x + 1
  if(x %% 1000 == 0)
    cat(x, " of ", length(unique(HOG_DE.a2tea$HOG)), "\n")
  # you need to look at the gene ids and not the species column, because sometime
  # there is no entry for species for genes
  species <- unique(HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
  species <- species[!is.na(species)]
  species <- unique(substr(species, 1, 2))
  if(length(species) == 5)
  {
    conserved_HOGs <- append(conserved_HOGs, hog)
  }
  
}

#get all HOGs/GO list from conserved set
SFA_bound <- bind_rows(SFA$Oryza_sativa, SFA$Triticum_aestivum, SFA$Sorghum_bicolor, SFA$Miscanthus_sinensis, SFA$Zea_mays)
SFA_conserved <- SFA_bound[SFA_bound$HOG %in% conserved_HOGs,]


#make a HOG-GO annotation list

#split DF by HOG
HOG2GO <- rep("NA", length(unique(SFA_conserved$HOG[!is.na(SFA_conserved$HOG)])))
names(HOG2GO) <- unique(SFA_conserved$HOG[!is.na(SFA_conserved$HOG)])


#get HOG and GO pair
annotation <- SFA_conserved[,c("HOG", "Gene-Ontology-Term")]
# split the strings of the GO Terms
annotation$`Gene-Ontology-Term` <- str_split(SFA_conserved$`Gene-Ontology-Term`, ", ")
#split the dataframe by HOG number
annotation <- split(annotation, annotation$HOG)
#unlist each GOlist of each HOG independently
for(i in 1:length(annotation))
{
  HOG2GO[unique(annotation[[i]]$HOG)] <- list(unique(unlist(annotation[[i]]$`Gene-Ontology-Term`)))
}


#default background is the conserved HOGs
test_GO_enrich <- function(interest, background = conserved_HOGs, nTop = 200)
{
  #my HOGs are my genes
  
  #list of HOG Of Interest
  HOGs_interest <- interest
  
  #filter for background so you can have different backgrounds too
  SFA_background <- SFA_conserved[SFA_conserved$HOG %in% background,]
  
  #make the gene universe (HOGs_all) and mark everything that is of interest
  #init HOGs all with not of interest
  HOGs_all <- rep(0, length(unique(SFA_background$HOG[!is.na(SFA_background$HOG)])))
  #set HOGs as names
  names(HOGs_all) <- unique(SFA_background$HOG[!is.na(SFA_background$HOG)])
  #set values (0 = not interesting, 1 = interesting)
  HOGs_all[names(HOGs_all) %in% HOGs_interest] <- 1
  #convert to factor
  HOGs_all <- as.factor(HOGs_all)
  
  #object which includes enrichment test result of all three ontologies
  allRes <- list(BP = c(), MF = c(), CC = c())
  
  #go through all three ontologies
  for(ontology in c("BP", "MF", "CC"))
  {
    #make a topGodata object to analyse
    GOdata <- new("topGOdata", ontology = ontology, allGenes = HOGs_all,
                  annot = annFUN.gene2GO, gene2GO = HOG2GO)
    print(GOdata)
    #use the fishers weight01 function
    result_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    #make a results table
    res <- GenTable(GOdata, weight01 = result_weight01, orderBy = "weight01", ranksOf = "classic", topNodes = nTop)
    res <- res[res$Significant > 1,]
    res$Ontology <- rep(ontology, length(res$GO.ID))
    allRes[ontology] <- list(res)
  }
  
  allRes <- bind_rows(allRes$BP, allRes$MF, allRes$CC)
  
  allRes$weight01 <- as.numeric(allRes$weight01)
  
  return(allRes)
  
}


plot_GO_graph <- function(Terms, top)
{
  #combine ID with Term
  Terms$Term <- paste(Terms$GO.ID, Terms$Term)
  
  #only plot top n terms and filter out NA
  Terms <- Terms[1:top,]
  Terms <- Terms[!is.na(Terms$GO.ID),]
  
  if(dim(Terms)[1] == 0)
  {
    ggplot() + ggtitle("\tNo Data")
  }
  else
  {
    maxn <- max(Terms$Significant)
    Terms_order <- Terms[order(Terms$weight01, decreasing=T),]
    ggplot(Terms_order) + 
      geom_point(aes(x = Terms_order$weight01, y = fct_inorder(Terms_order$Term), size = Terms_order$Significant)) +
      geom_vline(xintercept=0.05) +
      geom_vline(xintercept=0.01) +
      geom_vline(xintercept=0.001) +
      scale_size_binned_area(name = "Significant\nAnnotated", limits = c(0, maxn),
                             breaks = seq(0, maxn, maxn/5)) +
      scale_x_reverse() +
      xlab("p-Value") +
      ylab("GO Term")
  }
  
}

plot_GO <- function(data, title, pdf, top = 10, width=10, height=10)
{
  bp <- plot_GO_graph(data[data$Ontology == "BP",], top = top)
  mf <- plot_GO_graph(data[data$Ontology == "MF",], top = top)
  cc <- plot_GO_graph(data[data$Ontology == "CC",], top = top)
  
  grid <- plot_grid(bp, mf, cc, labels = c("BP", "MF", "CC"), ncol=1)
  t <- ggdraw() + draw_label(title)
  plot <- plot_grid(t, grid, ncol=1, rel_heights=c(0.1, 1))
  ggsave(filename = pdf, plot = plot, width=width, height=height)
}

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

#51 HOGs
# MSZ_expanded_deg_any_exp_species_hogs ################################################
MSZ_expanded_deg_any_exp_species_hogs <- read.csv(".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_HOGs.csv")
MSZ_expanded_deg_any_exp_species_hogs <- MSZ_expanded_deg_any_exp_species_hogs$x
MSZ_expanded_deg_any_exp_species_enrich <- test_GO_enrich(MSZ_expanded_deg_any_exp_species_hogs)
plot_GO(MSZ_expanded_deg_any_exp_species_enrich, "Ms, Sb, Zm vs Ta, Os, Conserved, Expanded, DEG in Any Expanded Species\ntopGO Analysis",
        ".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.pdf")
write.csv(MSZ_expanded_deg_any_exp_species_enrich, file = ".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")

#20 HOGs
# AllC4Exp_NoDEGC4_AllC3DEG_hogs ################################################
AllC4Exp_NoDEGC4_AllC3DEG_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
AllC4Exp_NoDEGC4_AllC3DEG_hogs <- AllC4Exp_NoDEGC4_AllC3DEG_hogs$x
AllC4Exp_NoDEGC4_AllC3DEG_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3DEG_hogs)
plot_GO(AllC4Exp_NoDEGC4_AllC3DEG_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
#enrichemnt test with background conserved, AllC3DEG, NoDEGC4
AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs <- read.csv(".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs.csv")
AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs <- AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs$x
AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3DEG_hogs, background = AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs)
plot_GO(AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG vs Conserved, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.csv")


# # AllC4Exp_NoDEGC4_AllC3Up_hogs ################################################
# AllC4Exp_NoDEGC4_AllC3Up_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
# AllC4Exp_NoDEGC4_AllC3Up_hogs <- AllC4Exp_NoDEGC4_AllC3Up_hogs$x
# AllC4Exp_NoDEGC4_AllC3Up_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3Up_hogs)
# plot_GO(AllC4Exp_NoDEGC4_AllC3Up_enrich, "All C4 Expanded, No DEG in C4, All C3 Upregulated\ntopGO Analysis",
#         ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs_moreThan1Sig_topGO.pdf")
# write.csv(AllC4Exp_NoDEGC4_AllC3Up_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs_moreThan1Sig_topGO.csv")

# # AllC4Exp_DEGAnyC4_NoDEGC3 ################################################
# AllC4Exp_DEGAnyC4_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
# AllC4Exp_DEGAnyC4_NoDEGC3_hogs <- AllC4Exp_DEGAnyC4_NoDEGC3_hogs$x
# AllC4Exp_DEGAnyC4_NoDEGC3_enrich <- test_GO_enrich(AllC4Exp_DEGAnyC4_NoDEGC3_hogs)
# plot_GO(AllC4Exp_DEGAnyC4_NoDEGC3_enrich, "All C4 Expanded, DEG in Any C4, No DEG in C3\ntopGO Analysis",
#         ".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_moreThan1Sig_topGO.pdf")
# write.csv(AllC4Exp_DEGAnyC4_NoDEGC3_enrich, file = ".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_moreThan1Sig_topGO.csv")

# # ZmExp_AllZmDEG_NoDEGC3_hogs ################################################
# ZmExp_AllZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
# ZmExp_AllZmDEG_NoDEGC3_hogs <- ZmExp_AllZmDEG_NoDEGC3_hogs$x
# ZmExp_AllZmDEG_NoDEGC3_enrich <- test_GO_enrich(ZmExp_AllZmDEG_NoDEGC3_hogs)
# plot_GO(ZmExp_AllZmDEG_NoDEGC3_enrich, "Zm Expanded, All Zm DEG, No DEG in C3\ntopGO Analysis",
#         ".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
# write.csv(ZmExp_AllZmDEG_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")

#48 HOGs
# ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs ################################################
ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs <- ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs$x
ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich <- test_GO_enrich(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs)
plot_GO(ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich, "Zm Expanded, Only Any Zm DEG, No DEG in C3\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
#enrichment test with background: Conserved Only Any Zm DEG, No DEG C3
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs <- read.csv(".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs.csv")
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs <- ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs$x
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich <- test_GO_enrich(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs, background=ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs)
plot_GO(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich, "Zm Expanded, Only Any Zm DEG, No DEG in C3 vs Conserved, Only Any Zm DEG, No DEG in C3\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_bg_topGO.pdf")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich, file = ".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_bg_topGO.csv")


#513 HOGs
# ZmExp_NoDEGZm_AllC3DEG_hogs ################################################
ZmExp_NoDEGZm_AllC3DEG_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
ZmExp_NoDEGZm_AllC3DEG_hogs <- ZmExp_NoDEGZm_AllC3DEG_hogs$x
ZmExp_NoDEGZm_AllC3DEG_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3DEG_enrich, "Zm Expanded, No Zm DEG, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3DEG_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG DEG
ZmExp_NoDEGZm_AllC3DEG_bg_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs.csv")
ZmExp_NoDEGZm_AllC3DEG_bg_hogs <- ZmExp_NoDEGZm_AllC3DEG_bg_hogs$x
ZmExp_NoDEGZm_AllC3DEG_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_hogs, background=ZmExp_NoDEGZm_AllC3DEG_bg_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3DEG_bg_enrich, "Zm Expanded, No Zm DEG, All C3 DEGregulated vs Conserved, No Zm DEG, All C3 DEG\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.csv")

NoDEGZm_AllC3DEG_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_bg_hogs)
plot_GO(NoDEGZm_AllC3DEG_bg_enrich, "No Zm DEG, All C3 DEGregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(NoDEGZm_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs_moreThan1Sig_bg_topGO.csv")



#204 HOGs
# ZmExp_NoDEGZm_AllC3Up_hogs ################################################
ZmExp_NoDEGZm_AllC3Up_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
ZmExp_NoDEGZm_AllC3Up_hogs <- ZmExp_NoDEGZm_AllC3Up_hogs$x
ZmExp_NoDEGZm_AllC3Up_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Up_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Up_enrich, "Zm Expanded, No Zm DEG, All C3 Upregulated\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Up_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG Up
NoDEGZm_AllC3Up_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\NoDEGZm_AllC3Up_hogs.csv")
NoDEGZm_AllC3Up_hogs <- NoDEGZm_AllC3Up_hogs$x
ZmExp_NoDEGZm_AllC3Up_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Up_hogs, background = NoDEGZm_AllC3Up_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Up_bg_enrich, "Zm Expanded, No Zm DEG, All C3 Upregulated vs Conserved, No Zm DEG, All C3 Upregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Up_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_bg_hogs_moreThan1Sig_bg_topGO.csv")

#105 HOGs
# ZmExp_NoDEGZm_AllC3Down_hogs ################################################
ZmExp_NoDEGZm_AllC3Down_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
ZmExp_NoDEGZm_AllC3Down_hogs <- ZmExp_NoDEGZm_AllC3Down_hogs$x
ZmExp_NoDEGZm_AllC3Down_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Down_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Down_enrich, "Zm Expanded, No Zm DEG, All C3 Downregulated\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Down_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG Down
NoDEGZm_AllC3Down_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\NoDEGZm_AllC3Down_hogs.csv")
NoDEGZm_AllC3Down_hogs <- NoDEGZm_AllC3Down_hogs$x
ZmExp_NoDEGZm_AllC3Down_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Down_hogs, background = NoDEGZm_AllC3Down_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Down_bg_enrich, "Zm Expanded, No Zm DEG, All C3 Downregulated vs Conserved, No Zm DEG, All C3 Downregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Down_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_bg_hogs_moreThan1Sig_bg_topGO.csv")

After filtering and calculating the topThe GO Terms and HOG counts were put into ReviGO for clustering and the significant terms were higlighted in a TreeMap. The following are all the Treemap scripts in EX_THE000003
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9337110928954263,0.04668641,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.9074954800738403,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9260545876803968,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.9562618187433756,0.06960742,"protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9728461942621843,0.04505158,"calmodulin binding"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9791953685095711,-0,"phosphatidylinositol transfer activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,6,0.9303677761411089,0.0830445,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9369434002734686,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.931141723353588,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9456232781319497,0.05997119,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9481756669577021,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9463330863851006,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9414959979979491,0.04472286,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9603176088388498,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,2,0.8295433660701551,0,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8716740809455901,0.35450436,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.7962202662229924,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9205756783951948,0.03322883,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.905059530511471,-0,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9136241846650163,0.31015716,"iron-sulfur cluster binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,4,0.8634289110966092,0.26746377,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.9129423235538945,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.8989871378015214,0.37942352,"iron-sulfur cluster binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9323615679242189,0.10462209,"iron-sulfur cluster binding"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,1,0.8447059239259296,0.03812803,"protein serine kinase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.866922162628014,0.42844795,"protein serine kinase activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8530913035056045,0.25738438,"protein serine kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9308067313434633,0.18603615,"protein serine kinase activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,1,0.7838849590751893,0.37283533,"protein serine kinase activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9157122115697741,0.05781294,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.935748566292739,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9712716402102628,0.01879944,"pterocarpan synthase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9593242515245328,0.23977693,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, " Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}






# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}



dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000502","proteasome complex",0.37652382533874423,1,0.8271175970868048,-0,"proteasome complex"),
                     c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,1,0.9999075555427188,4.999E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9878249124521837,2.721E-05,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,8,0.876325831792117,0.09354521,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,6,0.8517450403470601,0.17160779,"cytoplasm"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7954947398274242,0,"vacuole"),
                     c("GO:0005634","nucleus",16.5161752456724,5,0.722821401188212,0.35145602,"vacuole"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7846459032454033,0.21968674,"vacuole"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7604285311621104,0.25147478,"vacuole"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7266424122076163,0.18309297,"vacuole"),
                     c("GO:0005777","peroxisome",0.7476308639759881,1,0.8084648576917414,0.18435889,"vacuole"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7583222273090043,0.22703414,"vacuole"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.8098307472795678,0.18302188,"vacuole"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7808471979237221,0.21139748,"vacuole"),
                     c("GO:0005886","plasma membrane",17.177321395000487,5,0.9808601938612675,0.06338235,"plasma membrane"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999273299581207,4.05E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,7,0.9998444774268712,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999407554378136,3.357E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, " Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)





# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,2,0.859657561830051,-0,"monoatomic ion transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8452765335777909,0.4188028,"monoatomic ion transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.846647839368461,0.30464429,"monoatomic ion transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.840916109634737,0.37247286,"monoatomic ion transport"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,2,0.8772005106481118,-0,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.9292980500454672,0.12277613,"lipid catabolic process"),
                     c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8925541854984664,0.46381088,"lipid catabolic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.9254223629521837,0.07741388,"phosphorylation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9474388241629303,0.11172681,"phosphorylation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,2,0.7760089427909355,0,"intracellular signal transduction"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.8131845409871619,0.37522238,"intracellular signal transduction"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8258431389896448,0.47274124,"intracellular signal transduction"),
                     c("GO:0009733","response to auxin",0.1256056236300358,1,0.8526513918430121,0.49934907,"intracellular signal transduction"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.9428240885966082,0.39926133,"intracellular signal transduction"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8276119997819956,0.46171648,"intracellular signal transduction"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9605487851922724,0.21912257,"intracellular signal transduction"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8305467899052755,0.42391909,"intracellular signal transduction"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9803322186234948,0.11943324,"intracellular signal transduction"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8115566104700873,0.42151038,"intracellular signal transduction"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8876779238953187,0.07105297,"DNA demethylation"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8589276766539276,0.36848795,"DNA demethylation"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8610960350815848,0.27735821,"DNA demethylation"),
                     c("GO:0006508","proteolysis",5.2622572267907,1,0.8629236559306371,0.44475214,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9180290003370447,0.15925058,"DNA demethylation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9344204935933396,0.23039539,"DNA demethylation"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.9016549919916592,0.15970956,"DNA demethylation"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9067097943832549,0.1131676,"DNA demethylation"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9622253596069233,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.8988350893712908,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9569296911885713,0.43353892,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9555047392627222,0.35287398,"acquisition of seed longevity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)






# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}



dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,9,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,30,0.9999298415024575,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9766709014171154,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,204,0.9006467084488148,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,86,0.8711955549649424,0.17160779,"cytoplasm"),
                     c("GO:0009504","cell plate",0.0007255512462560189,2,0.9395571191954046,0.05708197,"cell plate"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,4,0.9802540780053186,2.997E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,128,0.9365244522321454,0.34881718,"plant-type cell wall"),
                     c("GO:0009898","cytoplasmic side of plasma membrane",0.2657902935257323,1,0.9550365138301414,0.23284101,"plant-type cell wall"),
                     c("GO:0009925","basal plasma membrane",0.13578741268972233,1,0.9463298922396219,0.18839421,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,25,0.9764118381717231,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,106,0.7353647898965974,0,"chloroplast"),
                     c("GO:0000138","Golgi trans cisterna",0.007541260384887046,2,0.7868433828079426,0.37607379,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,11,0.8161229772690657,0.13280054,"chloroplast"),
                     c("GO:0000785","chromatin",1.087181392930179,4,0.6154547100122836,0.45756389,"chloroplast"),
                     c("GO:0001673","male germ cell nucleus",0.006684016617906476,1,0.8400958319844106,0.21436526,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,279,0.7434541920863925,0.35145602,"chloroplast"),
                     c("GO:0005730","nucleolus",1.2114693153196516,25,0.6574244508381973,0.20773614,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,49,0.7694728451028692,0.2503098,"chloroplast"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,12,0.6508642725908665,0.20731432,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,12,0.7958642379536632,0.21057364,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,9,0.8012526236944911,0.17354881,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,59,0.707932123912852,0.19658734,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,4,0.7748705589208366,0.38610723,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,10,0.7214407684248566,0.44106058,"chloroplast"),
                     c("GO:0005905","clathrin-coated pit",0.09808111076528711,3,0.8806413166374408,0.46642691,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,100,0.8057497688689464,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,2,0.8190122800115598,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,18,0.7710315941262023,0.39454155,"chloroplast"),
                     c("GO:0012511","monolayer-surrounded lipid storage body",0.0046142077544432435,1,0.8242011923160282,0.28396619,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,2,0.6825635969117441,0.39861008,"chloroplast"),
                     c("GO:0031977","thylakoid lumen",0.009951481990600534,4,0.7865019059678459,0.29998046,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,6,0.8028487298659475,0.22189912,"chloroplast"),
                     c("GO:0032044","DSIF complex",0.014580598332295611,1,0.6977513011508626,0.49997569,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.788446534024034,0.28709921,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7810150470710108,0.31004802,"chloroplast"),
                     c("GO:0042470","melanosome",0.049131249288425556,1,0.7757444321211127,0.3330034,"chloroplast"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7274865822451015,0.486126,"chloroplast"),
                     c("GO:0043661","peribacteroid membrane",1.987811633578134E-05,1,0.7980277845638879,0.30368923,"chloroplast"),
                     c("GO:0070971","endoplasmic reticulum exit site",0.06951625759076932,1,0.77292085047214,0.45185985,"chloroplast"),
                     c("GO:0009524","phragmoplast",0.002879842104146322,3,0.9348019928088872,0.06218567,"phragmoplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999411529088681,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8857840483760228,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,9,0.9999250488252842,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,378,0.9999011439048419,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,6,1,-0,"protein-containing complex"),
                     c("GO:0034045","phagophore assembly site membrane",0.08896450966078939,1,0.9108104635559916,0.07998094,"phagophore assembly site membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,3,0.999932835127444,5.465E-05,"cell projection"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999462429024586,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,5,0.9965175146401201,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,3,0.9120657479855075,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0070469","respirasome",0.5984977809313304,3,0.968203139496543,0.0715415,"respirasome"),
                     c("GO:0090406","pollen tube",0.001702063711251277,4,0.9947034317355141,2.369E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.995454037968868,0.27222822,"pollen tube"),
                     c("GO:0098552","side of membrane",0.7241349304670945,3,0.9677007152814887,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,17,0.7938877775225776,-0,"ribonucleoprotein complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,2,0.7836928687770224,0.25406433,"ribonucleoprotein complex"),
                     c("GO:0000148","1,3-beta-D-glucan synthase complex",0.012428792238947283,2,0.7603198759962925,0.41849798,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7658118738451125,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000178","exosome (RNase complex)",0.10685978389207654,1,0.7658931000593635,0.34132399,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,4,0.7616170306174925,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0000974","Prp19 complex",0.0635900941581645,1,0.8509571458457826,0.24779784,"ribonucleoprotein complex"),
                     c("GO:0005667","transcription regulator complex",0.7880231963702957,1,0.8214255995273086,0.32589738,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7261137000342258,0.4447901,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,12,0.6207102252827627,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.813311779148287,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7984871088787376,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005950","anthranilate synthase complex",0.00014660110797638738,1,0.8267219358756437,0.15701804,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.771298490982265,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.801646627727481,0.43971777,"ribonucleoprotein complex"),
                     c("GO:0010330","cellulose synthase complex",0.0019629639881584074,1,0.8288151557394036,0.37568323,"ribonucleoprotein complex"),
                     c("GO:0016272","prefoldin complex",0.03795726314317447,1,0.8558414471368166,0.23619435,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,2,0.7844583658805595,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0016602","CCAAT-binding factor complex",0.03082101937862897,2,0.7042485204988549,0.48430429,"ribonucleoprotein complex"),
                     c("GO:0030289","protein phosphatase 4 complex",0.015278817168589934,1,0.818783800666579,0.29979865,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7595578125779482,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.8016531148164839,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,3,0.7618393545642829,0.47090253,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7817374230316501,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0045252","oxoglutarate dehydrogenase complex",0.07266196950090673,1,0.7324805062638388,0.46971574,"ribonucleoprotein complex"),
                     c("GO:0045273","respiratory chain complex II",0.05444864540824706,1,0.7586647868595982,0.39997785,"ribonucleoprotein complex"),
                     c("GO:0046930","pore complex",0.06221601936645362,3,0.7967552314432623,0.24728353,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7158452399269217,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071439","clathrin complex",0.011926869801468804,1,0.7703665716373392,0.36298129,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.728818219916503,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,9,0.7813093676340951,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.7137965779612403,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7966516257371146,0.40438888,"ribonucleoprotein complex"),
                     c("GO:0160064","multi-pass translocon complex",0.0006336149582030303,1,0.8478591733073518,0.17222646,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000003","reproduction",1.1154428097959028,1,1,-0,"reproduction"),
                     c("GO:0002376","immune system process",0.9427113541805001,2,1,-0,"immune system process"),
                     c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,3,0.9686358216302681,-0,"intracellular iron ion homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.981020579886502,0.46430576,"intracellular iron ion homeostasis"),
                     c("GO:0006952","defense response",1.1604144588756624,44,0.9134724519311507,-0,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.879827777764831,0.44951817,"defense response"),
                     c("GO:0009414","response to water deprivation",0.019128848432972426,16,0.8785314214004789,0.38192823,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,8,0.9254516172618769,0.45786068,"defense response"),
                     c("GO:0009625","response to insect",0.000539778096485113,2,0.9378941685822639,0.49043501,"defense response"),
                     c("GO:0009646","response to absence of light",0.0015108857221249963,3,0.9253577229366186,0.4554337,"defense response"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,11,0.8531535145513461,0.23262547,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8982731987736398,0.25765651,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,2,0.933346065833084,0.4902216,"defense response"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9350745269183776,0.4570295,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9329884633724373,0.29372182,"defense response"),
                     c("GO:0010447","response to acidic pH",0.005513623752681268,1,0.9263956476869073,0.48753241,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8634333996745125,0.38492887,"defense response"),
                     c("GO:0019236","response to pheromone",0.03625386311050012,1,0.9148117746363615,0.4993889,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,3,0.8757630073551622,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,4,0.9108509718162283,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,6,0.8413285625884777,0.40205681,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9203174924393598,0.46171648,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,13,0.9025456246148672,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,2,0.9820782649233186,0.48259105,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9298786822575796,0.49002119,"defense response"),
                     c("GO:0071000","response to magnetism",7.887168533115806E-05,1,0.9402222825420283,0.39597089,"defense response"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9390096720600496,0.14771546,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9284492141524632,0.44641368,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9248163774630318,0.39541577,"defense response"),
                     c("GO:1901562","response to paraquat",0.00022182661499388205,1,0.9351793434848085,0.37754219,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.934434494507892,0.31063873,"defense response"),
                     c("GO:0007017","microtubule-based process",1.4522840599239462,1,0.990143487369369,0.01367267,"microtubule-based process"),
                     c("GO:0007049","cell cycle",2.8073119376140743,9,0.989492087424824,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,3,0.9757722735243382,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,9,0.9959005435852354,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,12,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,19,1,-0,"metabolic process"),
                     c("GO:0008219","cell death",0.4651211168388386,1,0.9910898639992425,0.01191292,"cell death"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,2,0.9920904124444596,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9901838028127955,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,8,0.9563835628896364,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,7,0.9540238308470911,0.21768792,"biosynthetic process"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,3,0.9530524036319646,0.27529491,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9280437098504655,0.16846091,"biosynthetic process"),
                     c("GO:0071704","organic substance metabolic process",50.871943734517124,1,0.9521326283355345,0.31960258,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9431382619255397,0.17227516,"biosynthetic process"),
                     c("GO:0009695","jasmonic acid biosynthetic process",0.005493905831348478,4,0.8820553239475333,0.06061134,"jasmonic acid biosynthetic process"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,4,0.9107687165294583,0.26380719,"jasmonic acid biosynthetic process"),
                     c("GO:0006694","steroid biosynthetic process",0.2949431320360321,3,0.8643322450747801,0.40616444,"jasmonic acid biosynthetic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,1,0.8889571051236331,0.46657264,"jasmonic acid biosynthetic process"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9322091152525966,0.38978709,"jasmonic acid biosynthetic process"),
                     c("GO:0019432","triglyceride biosynthetic process",0.09997725537774264,2,0.8930989767583173,0.49870367,"jasmonic acid biosynthetic process"),
                     c("GO:0019760","glucosinolate metabolic process",0.005523482713347663,2,0.8984957731210338,0.22547604,"jasmonic acid biosynthetic process"),
                     c("GO:0031407","oxylipin metabolic process",0.012508556345488349,1,0.9241968951629272,0.27726473,"jasmonic acid biosynthetic process"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,2,0.9032366017649239,0.31813021,"jasmonic acid biosynthetic process"),
                     c("GO:0033387","putrescine biosynthetic process from ornithine",0.01868765994315126,1,0.8776589989161953,0.2843421,"jasmonic acid biosynthetic process"),
                     c("GO:0033491","coniferin metabolic process",1.7253181166190823E-05,1,0.9325722453908004,0.3924192,"jasmonic acid biosynthetic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,2,0.9100973444117998,0.3067044,"jasmonic acid biosynthetic process"),
                     c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1,0.9270117887483821,0.26504175,"jasmonic acid biosynthetic process"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9080846791615016,0.48869939,"jasmonic acid biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,2,0.962069788115057,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,19,0.8822204902848627,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,7,0.9155596328357238,0.4203896,"flower development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,4,0.9119539320899355,0.38141807,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,2,0.9356805668650227,0.48720779,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.8583572508885404,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9190179371314678,0.29044182,"flower development"),
                     c("GO:0019827","stem cell population maintenance",0.017275363827690213,2,0.9245142188070515,0.36342025,"flower development"),
                     c("GO:0021987","cerebral cortex development",0.017519373104183483,1,0.9190617509063173,0.44335853,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.9261666268625263,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.920758101491622,0.42767441,"flower development"),
                     c("GO:0055047","generative cell mitosis",7.147746483136198E-05,1,0.8937433537239626,0.40521407,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9379660468836518,0.28605225,"flower development"),
                     c("GO:0090708","specification of plant organ axis polarity",0.002785156388256519,2,0.9048766078118994,0.42388435,"flower development"),
                     c("GO:1905177","tracheary element differentiation",0.0002661919379926584,1,0.9383884191347898,0.2910313,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,2,0.9295369247237566,0.4298524,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.9938703007851393,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9593105824517872,0.04613723,"chlororespiration"),
                     c("GO:0015031","protein transport",3.093438694074183,39,0.8856881222322535,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,28,0.9299296200807329,0.42145531,"protein transport"),
                     c("GO:0006833","water transport",0.07339010320064257,2,0.9512306966604464,0.25687781,"protein transport"),
                     c("GO:0006897","endocytosis",0.6399968963991822,5,0.9251490553342954,0.32211657,"protein transport"),
                     c("GO:0009852","auxin catabolic process",3.4506362332381646E-05,1,0.918609329134268,0.44733804,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,2,0.9523856179247027,0.21356797,"protein transport"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9433319652038283,0.48373949,"protein transport"),
                     c("GO:0010966","regulation of phosphate transport",0.0010204024289718573,1,0.9540785580000554,0.48139262,"protein transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.951568450969917,0.18579698,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9498793467403328,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,14,0.9257846983818245,0.38566904,"protein transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,8,0.9031255984668667,0.44154323,"protein transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9561906223935729,0.30124551,"protein transport"),
                     c("GO:0045037","protein import into chloroplast stroma",0.006780500198312994,1,0.9230660938123257,0.49429905,"protein transport"),
                     c("GO:0050821","protein stabilization",0.10143638155636905,1,0.9390658629772833,0.43402974,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9224341182431762,0.44310395,"protein transport"),
                     c("GO:0051646","mitochondrion localization",0.04394385243028803,1,0.9470244054068786,0.2322362,"protein transport"),
                     c("GO:0060918","auxin transport",0.015500750907739157,3,0.9005987910456148,0.22426661,"protein transport"),
                     c("GO:0061635","regulation of protein complex stability",0.0017622892191180627,1,0.9519447151888122,0.33196934,"protein transport"),
                     c("GO:0072699","protein localization to cortical microtubule cytoskeleton",0.00015774337066231612,1,0.9467319327152923,0.36351639,"protein transport"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,2,0.9278062822324343,0.3859558,"protein transport"),
                     c("GO:0120010","intermembrane phospholipid transfer",0.007554428610624982,1,0.8902307770751422,0.42470759,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9314889060544345,0.49534489,"protein transport"),
                     c("GO:1904823","purine nucleobase transmembrane transport",0.0731485586643159,1,0.9290749747825838,0.49367033,"protein transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1,0.9297548608646681,0.30569128,"protein transport"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9337466975912142,0.46789196,"protein transport"),
                     c("GO:2000694","regulation of phragmoplast microtubule organization",0.0019791863537787476,1,0.9443902205137444,0.42590434,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,9,0.9550108711717666,0.04477571,"photosynthesis"),
                     c("GO:0016477","cell migration",0.49093927008395993,1,0.991049257897926,0.01198612,"cell migration"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,47,0.8875503224800294,0,"protein ubiquitination"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,3,0.889910754013842,0.4354711,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8756768221539741,0.45078305,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,20,0.9424227992072893,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9378176150799266,0.43492456,"protein ubiquitination"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,1,0.9207936575917155,0.34119038,"protein ubiquitination"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,2,0.8961847229904666,0.46608863,"protein ubiquitination"),
                     c("GO:0006076","(1->3)-beta-D-glucan catabolic process",2.4647401665986893E-06,1,0.9369468412709423,0.49567686,"protein ubiquitination"),
                     c("GO:0006182","cGMP biosynthetic process",0.05174475505757287,1,0.8495420820947562,0.45290106,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,3,0.8740544592960273,0.4022316,"protein ubiquitination"),
                     c("GO:0006283","transcription-coupled nucleotide-excision repair",0.06059070751549558,1,0.8366645313379071,0.46757848,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,24,0.8417882101772437,0.19691418,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,18,0.8296522788341744,0.48035482,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,20,0.8674236250172993,0.40866092,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,10,0.8693454815537701,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,4,0.9423874205568638,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,32,0.9411424466613101,0.12094833,"protein ubiquitination"),
                     c("GO:0006747","FAD biosynthetic process",0.030434611577160615,1,0.8585945486528581,0.43629805,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,3,0.9201546126255483,0.45718665,"protein ubiquitination"),
                     c("GO:0008652","amino acid biosynthetic process",2.679426429329935,12,0.8381579891861963,0.49906638,"protein ubiquitination"),
                     c("GO:0009820","alkaloid metabolic process",0.01601588160255828,2,0.9485491751247405,0.14390223,"protein ubiquitination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8965124392057823,0.367811,"protein ubiquitination"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9437965846692097,0.32178408,"protein ubiquitination"),
                     c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8647900867621974,0.36610043,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,4,0.8811583262165019,0.20835326,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,77,0.9098459003605229,0.44883397,"protein ubiquitination"),
                     c("GO:0016554","cytidine to uridine editing",0.012560315888986921,1,0.8967884109527346,0.3977518,"protein ubiquitination"),
                     c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,2,0.8745940089841773,0.46779521,"protein ubiquitination"),
                     c("GO:0018345","protein palmitoylation",0.0003918936864891916,1,0.9136811946983828,0.37471442,"protein ubiquitination"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9399102664551204,0.33268781,"protein ubiquitination"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9241297699418284,0.27916289,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,8,0.8875206936688738,0.25012422,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,2,0.8710247194655414,0.28833637,"protein ubiquitination"),
                     c("GO:0033358","UDP-L-arabinose biosynthetic process",2.21826614993882E-05,1,0.9080094382781345,0.46513452,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9496007741312238,0.45767193,"protein ubiquitination"),
                     c("GO:0042245","RNA repair",0.022010129687726296,1,0.9033176951860613,0.35967549,"protein ubiquitination"),
                     c("GO:0042726","flavin-containing compound metabolic process",0.22289631222618586,1,0.9078809166818013,0.23046044,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,2,0.9061961605967513,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,3,0.9234150366642541,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.888352038284189,0.31626431,"protein ubiquitination"),
                     c("GO:0046166","glyceraldehyde-3-phosphate biosynthetic process",0.03550458209985412,1,0.9017590795849828,0.33432601,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,9,0.9142523755930573,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,2,0.9025197296743325,0.41108531,"protein ubiquitination"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.8229889350754062,0.42789517,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8736344757012352,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.895328418162861,0.45115494,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,3,0.9364364843856474,0.45105203,"protein ubiquitination"),
                     c("GO:0140040","mitochondrial polycistronic RNA processing",3.4506362332381646E-05,1,0.9158728863742371,0.41070687,"protein ubiquitination"),
                     c("GO:1900871","chloroplast mRNA modification",0.00043625900948796796,1,0.9104502567020882,0.47298773,"protein ubiquitination"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.8391362254428162,0.46167258,"protein ubiquitination"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,2,0.8812813452757842,0.47967538,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,4,0.8947764310896186,0.31086297,"protein ubiquitination"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8768921732476757,0.49989178,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9920083214141256,0.01031083,"establishment of cell polarity"),
                     c("GO:0032259","methylation",2.6278542060840238,15,0.9678323230608921,0.05843136,"methylation"),
                     c("GO:0032502","developmental process",4.156005433076044,1,1,-0,"developmental process"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9772237518509119,0.03897815,"collagen metabolic process"),
                     c("GO:0040011","locomotion",0.5148398554794674,1,1,-0,"locomotion"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,16,0.8915756535066078,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,5,0.9097910614586415,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009786","regulation of asymmetric cell division",0.004426673339211246,1,0.955108677528184,0.14531513,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,5,0.924008625745951,0.38244242,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,2,0.9428171668433051,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010075","regulation of meristem growth",0.0026200187970944065,2,0.8605703832674918,0.39444345,"positive regulation of DNA-templated transcription"),
                     c("GO:0010089","xylem development",0.003164726373912717,1,0.9289878202466676,0.47349401,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,3,0.9515244205675838,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9498995339125141,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010366","negative regulation of ethylene biosynthetic process",0.0002045734338276912,1,0.9350142978377944,0.45396999,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9178782014905182,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0010928","regulation of auxin mediated signaling pathway",0.004086539196220627,1,0.9331502998788527,0.40226117,"positive regulation of DNA-templated transcription"),
                     c("GO:0032784","regulation of DNA-templated transcription elongation",0.17948484367188314,2,0.9295724228994686,0.37447493,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9399437687794638,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042659","regulation of cell fate specification",0.0026569898995933866,2,0.9372355320335743,0.39470747,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,4,0.948889587906444,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9380556035294847,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9355785026747696,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,10,0.8940374526132623,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,3,0.9293109081356303,0.42989004,"positive regulation of DNA-templated transcription"),
                     c("GO:0048509","regulation of meristem development",0.004372449055546074,1,0.9374968356560581,0.4381267,"positive regulation of DNA-templated transcription"),
                     c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,2,0.9385062593712483,0.39214155,"positive regulation of DNA-templated transcription"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,2,0.9105462894204667,0.45733164,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,3,0.9437530455211556,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,5,0.9361829156562242,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0060147","regulation of post-transcriptional gene silencing",0.021610841780737307,1,0.9400609881569504,0.26975857,"positive regulation of DNA-templated transcription"),
                     c("GO:0060195","negative regulation of antisense RNA transcription",0.00012077226816333578,1,0.9398665110039087,0.26135971,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9324036814431228,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,6,0.9497056054756448,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:0141005","retrotransposon silencing by heterochromatin formation",0.001355607091629279,1,0.9428929429040672,0.29342052,"positive regulation of DNA-templated transcription"),
                     c("GO:1900057","positive regulation of leaf senescence",0.0003204162216578296,1,0.9301175786325859,0.35419608,"positive regulation of DNA-templated transcription"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,5,0.9255845936231454,0.138781,"positive regulation of DNA-templated transcription"),
                     c("GO:1900457","regulation of brassinosteroid mediated signaling pathway",0.001054908791304239,1,0.9375717953283268,0.34087448,"positive regulation of DNA-templated transcription"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,3,0.9277197793630935,0.32535433,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9387344175913066,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1902326","positive regulation of chlorophyll biosynthetic process",0.00019717921332789513,1,0.9378805987248041,0.43671847,"positive regulation of DNA-templated transcription"),
                     c("GO:1902448","positive regulation of shade avoidance",4.190058283217772E-05,1,0.9348739791521216,0.30806662,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9232814504088713,0.48792077,"positive regulation of DNA-templated transcription"),
                     c("GO:1903329","regulation of iron-sulfur cluster assembly",0.0002464740166598689,1,0.949749293047066,0.40066627,"positive regulation of DNA-templated transcription"),
                     c("GO:1903527","positive regulation of membrane tubulation",0.0019298915504467734,1,0.9308269458228253,0.46246533,"positive regulation of DNA-templated transcription"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9341913688746851,0.40819627,"positive regulation of DNA-templated transcription"),
                     c("GO:1904966","positive regulation of vitamin E biosynthetic process",1.9717921332789515E-05,1,0.9401339564910491,0.38536455,"positive regulation of DNA-templated transcription"),
                     c("GO:1905038","regulation of membrane lipid metabolic process",0.009506502822571143,1,0.9404618961423398,0.25241702,"positive regulation of DNA-templated transcription"),
                     c("GO:1905157","positive regulation of photosynthesis",0.00010598382716374363,1,0.9382733536899409,0.39735289,"positive regulation of DNA-templated transcription"),
                     c("GO:1905183","negative regulation of protein serine/threonine phosphatase activity",1.7253181166190823E-05,1,0.9393315556699844,0.24023597,"positive regulation of DNA-templated transcription"),
                     c("GO:2000022","regulation of jasmonic acid mediated signaling pathway",0.011707515791343773,1,0.9292449834339215,0.42518681,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,5,0.9395119566225884,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,4,0.9278259249908258,0.37396531,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,3,0.9351126910533618,0.36714483,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,2,0.9306424609003232,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000049","positive regulation of cell-cell adhesion mediated by cadherin",0.0027210731439249528,1,0.9389922460317134,0.41411697,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,3,0.9369698607582899,0.39204232,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,2,0.9337478031918454,0.49708545,"positive regulation of DNA-templated transcription"),
                     c("GO:2000134","negative regulation of G1/S transition of mitotic cell cycle",0.02027248787027422,1,0.9295991590062869,0.16280015,"positive regulation of DNA-templated transcription"),
                     c("GO:2000185","regulation of phosphate transmembrane transport",0.00017253181166190822,1,0.9556565290595118,0.11823538,"positive regulation of DNA-templated transcription"),
                     c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,1,0.9425165946111713,0.2795365,"positive regulation of DNA-templated transcription"),
                     c("GO:2000280","regulation of root development",0.008382581306602141,1,0.93530346561912,0.45378622,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,2,0.9400621691337849,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.9320504846101171,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,2,0.9564176021664041,0.12578285,"positive regulation of DNA-templated transcription"),
                     c("GO:2000762","regulation of phenylpropanoid metabolic process",0.005555524335513445,1,0.9435717501109856,0.23331889,"positive regulation of DNA-templated transcription"),
                     c("GO:2001006","regulation of cellulose biosynthetic process",0.002356291599268347,2,0.9450457632488612,0.22862385,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,3,1,-0,"rhythmic process"),
                     c("GO:0048544","recognition of pollen",0.04932191547380636,1,0.9924948346594075,0.00950258,"recognition of pollen"),
                     c("GO:0050896","response to stimulus",17.567785530535815,12,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,3,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,6,0.9900713261441032,0.01381155,"cell division"),
                     c("GO:0060361","flight",9.612486649734888E-05,1,0.9992041163488156,-0,"flight"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,24,0.9115489367113405,0.01286368,"cell wall organization"),
                     c("GO:0000741","karyogamy",0.004604134631206351,1,0.9409217300776628,0.24960876,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,18,0.909556730491422,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,2,0.9283815069843394,0.39503669,"cell wall organization"),
                     c("GO:0007005","mitochondrion organization",0.8276597479438399,1,0.9175331228084536,0.4723023,"cell wall organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,3,0.9205428040419101,0.38546222,"cell wall organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,2,0.9238270729386129,0.42383478,"cell wall organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9263348218637661,0.40792275,"cell wall organization"),
                     c("GO:0009658","chloroplast organization",0.0906284959258338,7,0.92573827825937,0.31305597,"cell wall organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9565340798418022,0.2057032,"cell wall organization"),
                     c("GO:0010275","NAD(P)H dehydrogenase complex assembly",0.0022404488114382086,1,0.9342787929768167,0.41693818,"cell wall organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.929799954548471,0.45272646,"cell wall organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,5,0.9160784415903175,0.35075221,"cell wall organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9274077633240769,0.39384196,"cell wall organization"),
                     c("GO:0034462","small-subunit processome assembly",0.005688620304509775,1,0.9171810923498835,0.44451454,"cell wall organization"),
                     c("GO:0043622","cortical microtubule organization",0.008850881938255893,3,0.9340381291657605,0.26122046,"cell wall organization"),
                     c("GO:0048564","photosystem I assembly",0.005678761343843379,1,0.8947685122755664,0.44445982,"cell wall organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9258000800482549,0.48565697,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9293688456269509,0.4524177,"cell wall organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,2,0.9505782803397702,0.20466634,"cell wall organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9461820487839615,0.22063843,"cell wall organization"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.957105431386808,0.08070223,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9352773509874961,0.08958079,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
#pdf( file="test_BP.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")

    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")

    }
  )
}

# 
# 
# topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP",]
# 
# terms <- intersect(t$tm$description, topGO_data$Term)
# 
# 
# 
# p <- topGO_data[topGO_data$Term %in% terms,]$weight01
# l <- log10(p)
# 
# 
# colfunc<-colorRampPalette(c("blue", "yellow"))
# 
# colors <- colfunc(100)
# col_index <- round(l / log10(0.01) * 100)
# col_plot <- colors[col_index]
# 
# 
# with(
#   # t$tm is a data frame with one row per rectangle.  Filter to the group we
#   # want to highlight.
#   t$tm %>%
#     filter(description %in% terms),
#   {
#     # Use grid.rect to add a rectangle on top of the treemap.
#     grid.rect(x = x0 + (w / 2),
#               y = y0 + (h / 2),
#               width = w,
#               height = h,
#               gp = gpar(col = col_plot, fill = NA, lwd = 3),
#               vp = "data")
#     
#   }
# )




dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000995","RNA polymerase III general transcription initiation factor activity",0.007712441507517546,1,0.9841871246288971,0.00722361,"RNA polymerase III general transcription initiation factor activity"),
                     c("GO:0003674","molecular_function",100,32,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,10,0.9689397995337248,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,42,0.9521522430374141,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,2,0.9641766955594843,0.42104634,"mRNA binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,6,0.9791567156099864,0.36500624,"mRNA binding"),
                     c("GO:0001006","RNA polymerase III type 3 promoter sequence-specific DNA binding",0.0190615719374988,1,0.9608605950710434,0.36290179,"mRNA binding"),
                     c("GO:0001217","DNA-binding transcription repressor activity",0.09621481442486451,1,0.9835529779383934,0.49577888,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,96,0.9454417456386672,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,3,0.966190227513994,0.19519762,"mRNA binding"),
                     c("GO:0003684","damaged DNA binding",0.33181460177613226,1,0.9582563633962112,0.4175852,"mRNA binding"),
                     c("GO:0003697","single-stranded DNA binding",0.40260363864918647,2,0.9575327891197729,0.42694223,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.959081297510557,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9591098450616944,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,11,0.9324790327408994,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,2,0.9636243827418178,0.428381,"mRNA binding"),
                     c("GO:0008327","methyl-CpG binding",0.009309036644208929,1,0.9540229178639519,0.34771157,"mRNA binding"),
                     c("GO:0019237","centromeric DNA binding",0.004625690909914204,1,0.9647495718932401,0.20059314,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,9,0.9680097931257933,0.13060477,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.964953719942075,0.41077037,"mRNA binding"),
                     c("GO:0035198","miRNA binding",0.029599099839661934,1,0.9644971613561167,0.41680219,"mRNA binding"),
                     c("GO:0042134","rRNA primary transcript binding",0.01629192287398833,1,0.9659286451165663,0.39794793,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9587398118533373,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,21,0.9458079172850538,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,6,0.9959828945185376,-0,"structural constituent of ribosome"),
                     c("GO:0003774","cytoskeletal motor activity",0.3926027441124113,1,1,-0,"cytoskeletal motor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,34,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,62,0.7991312682332853,0,"protein kinase activity"),
                     c("GO:0000234","phosphoethanolamine N-methyltransferase activity",0.0010222643861315665,1,0.9005656057222542,0.45971198,"protein kinase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,4,0.8861130940300226,0.4205441,"protein kinase activity"),
                     c("GO:0003975","UDP-N-acetylglucosamine-dolichyl-phosphate N-acetylglucosaminephosphotransferase activity",0.005718915079898721,1,0.9032388864135855,0.35348317,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8893158783180034,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8898915598540006,0.48555173,"protein kinase activity"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8899083254406018,0.46089842,"protein kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,2,0.8913250950727402,0.25688843,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8465994648490534,0.49507895,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,10,0.8545818422016799,0.41376722,"protein kinase activity"),
                     c("GO:0004801","transaldolase activity",0.029270910839342038,1,0.9113911394817553,0.21824965,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,22,0.8235425141772773,0.47753726,"protein kinase activity"),
                     c("GO:0008146","sulfotransferase activity",0.16979789278712867,1,0.8972810614263171,0.25770873,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,15,0.8467149998284565,0.35677507,"protein kinase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,7,0.8646788483075069,0.30425857,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,7,0.8647495345145602,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,3,0.871128081007062,0.35123485,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,6,0.8791448240759565,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8934390285663193,0.40329157,"protein kinase activity"),
                     c("GO:0008798","beta-aspartyl-peptidase activity",0.023041972313000234,1,0.8906679336776527,0.44886191,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.8506214055277963,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,11,0.8716753777217396,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,37,0.876628926563111,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,2,0.8859296476106466,0.31489566,"protein kinase activity"),
                     c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.04143767537877,2,0.8620848457185354,0.43304972,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.8988388139109963,0.49223587,"protein kinase activity"),
                     c("GO:0033855","nicotianamine aminotransferase activity",6.6524797362141E-06,1,0.9356878940872976,0.45711636,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8743065572322031,0.49253483,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9150530150150219,0.30679413,"protein kinase activity"),
                     c("GO:0047326","inositol tetrakisphosphate 5-kinase activity",0.0016542499610719059,1,0.9024036089977014,0.40502966,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,2,0.8849374424809161,0.48733868,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9222526746562715,0.29376839,"protein kinase activity"),
                     c("GO:0050734","hydroxycinnamoyltransferase activity",0.0034969868480032116,1,0.90309679849696,0.35899295,"protein kinase activity"),
                     c("GO:0051753","mannan synthase activity",0.005627997856837128,2,0.8976585684333968,0.49140545,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9039044082734732,0.47557181,"protein kinase activity"),
                     c("GO:0052667","phosphomethylethanolamine N-methyltransferase activity",0.00011530964876104439,1,0.9103636574993704,0.40189746,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9213773782071678,0.48656541,"protein kinase activity"),
                     c("GO:0070012","oligopeptidase activity",0.014826159838775825,1,0.8933720677284,0.31386666,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.9013348578283663,0.18579171,"protein kinase activity"),
                     c("GO:0102732","myo-inositol-1,2,3,4,6-heptakisphosphate 5-kinase activity",0.0006741179466030287,1,0.9086267095387811,0.38119148,"protein kinase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.8843734136147638,0.45431535,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8832947747905774,0.39355689,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,6,0.9746143375880698,-0,"GTPase activator activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,173,0.9765859638010934,0.07766851,"protein binding"),
                     c("GO:0097159","organic cyclic compound binding",39.22969739688024,1,0.9707619643887931,0.13134756,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,30,0.9548331501228443,0.05738816,"zinc ion binding"),
                     c("GO:0000035","acyl binding",0.03278342014006308,1,0.9779977426252987,0.11604884,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,148,0.9381184188920123,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,6,0.9622465842144501,0.3830934,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,2,0.9718546830497732,0.21260669,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,52,0.9580118674808056,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,12,0.9632929394898434,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,20,0.9471629003822944,0.19635374,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,5,0.9562520289102728,0.24645234,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9551793223880397,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.95672267465334,0.29697106,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,2,0.9651036010951621,0.26177282,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,3,0.9505238324807589,0.31800282,"zinc ion binding"),
                     c("GO:0030955","potassium ion binding",0.05969270067304912,1,0.9720836061172309,0.26193125,"zinc ion binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9821582854632487,0.2693746,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,2,0.956819596612351,0.2885161,"zinc ion binding"),
                     c("GO:0043325","phosphatidylinositol-3,4-bisphosphate binding",0.013358179310317913,1,0.9737742438971143,0.211085,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9775729510461595,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,7,0.9506049035463411,0.33990845,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,1,0.9564209130429503,0.34531063,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,5,0.955636599785464,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,12,0.9623105778002776,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9605050879749866,0.28523818,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,3,0.9590380608318773,0.29755636,"zinc ion binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,4,0.9552335844763332,0.30272555,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,6,0.9811618864945685,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,3,1,-0,"electron transfer activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,2,0.9870190666322757,0.02769872,"auxin binding"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,6,0.9143108305903176,-0,"antiporter activity"),
                     c("GO:0005338","nucleotide-sugar transmembrane transporter activity",0.0598301852542642,2,0.9275650340655242,0.3876045,"antiporter activity"),
                     c("GO:0015079","potassium ion transmembrane transporter activity",0.4915428577358782,3,0.9006833655355473,0.39988553,"antiporter activity"),
                     c("GO:0015267","channel activity",1.6825096949913436,3,0.9041584076811204,0.45334576,"antiporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9440732240270157,0.46922631,"antiporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9296392077666529,0.35182404,"antiporter activity"),
                     c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,2,0.9245893969014204,0.41895885,"antiporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,3,0.9555101850473171,0.20645071,"antiporter activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,2,0.9235065769550767,0.36857202,"antiporter activity"),
                     c("GO:0042910","xenobiotic transmembrane transporter activity",0.2653474592383718,1,0.9291784951893058,0.37757719,"antiporter activity"),
                     c("GO:0051119","sugar transmembrane transporter activity",0.13162596406073218,2,0.9178488318907265,0.35505294,"antiporter activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9343558497731906,0.4373877,"antiporter activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,7,0.8810902927317588,0.0531541,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004322","ferroxidase activity",0.06266857660838222,1,0.9179378657781734,0.34042691,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004324","ferredoxin-NADP+ reductase activity",0.018482806200448173,1,0.9239712537810244,0.30749203,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,5,0.8975107343242501,0.46480825,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,4,0.8877923428856705,0.41437051,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0010242","oxygen evolving activity",0.0031710153409287207,1,0.9312468395343216,0.26980743,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0015930","glutamate synthase activity",0.041380641452497105,2,0.9158182898465297,0.32846804,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9200321084537753,0.46927032,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,3,0.8964682192013551,0.41891471,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,2,0.8998981679879802,0.44913707,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016630","protochlorophyllide reductase activity",0.0022374506846133423,1,0.92786850794828,0.26342111,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,2,0.9090472749922742,0.39197077,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016651","oxidoreductase activity, acting on NAD(P)H",0.7914588191768638,1,0.90165282556963,0.4378256,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,6,0.8961932035919564,0.47359913,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9079736615097322,0.49506653,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.924129806024186,0.41198066,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0019139","cytokinin dehydrogenase activity",0.0057233500663895305,1,0.928972415421427,0.2813582,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0045543","gibberellin 2-beta-dioxygenase activity",0.006071496505918068,1,0.918140273601625,0.48062226,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0047037","salutaridine reductase (NADPH) activity",1.77399459632376E-05,1,0.9382493475303578,0.38914098,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0047501","(+)-neomenthol dehydrogenase activity",2.66099189448564E-05,1,0.9371855713015793,0.3973742,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0048529","magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase activity",0.002614424536332141,1,0.9256782569963189,0.33893217,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,6,0.9021566538274115,0.43460917,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102111","gibberellin A20,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102652","gibberellin A9,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102721","ubiquinol:oxygen oxidoreductase activity",0.0020112663735820627,1,0.9300274199055578,0.26152924,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102772","sphingolipid C4-monooxygenase activity",4.4349864908093996E-05,1,0.9387728340281034,0.20804942,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102924","gibberellin A44,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.35384313,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0106292","superoxide-generating NADPH oxidase activity",0.0027186467188661623,1,0.9270682711279755,0.26695083,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:1990136","linoleate 9S-lipoxygenase activity",0.00026609918944856393,1,0.9338910724870338,0.23015728,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,199,0.94795157278954,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,59,0.9514904785786901,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,96,0.9474747332767991,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,19,0.9568951737037287,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,6,0.9299383755313524,0.04651897,"carboxy-lyase activity"),
                     c("GO:0009978","allene oxide synthase activity",9.09172230615927E-05,1,0.9582979753288584,0.33685723,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.949980592318686,0.45277577,"carboxy-lyase activity"),
                     c("GO:0047769","arogenate dehydratase activity",0.02556769711951619,2,0.9441430404368263,0.47917335,"carboxy-lyase activity"),
                     c("GO:0106099","2-keto-3-deoxy-L-rhamnonate aldolase activity",0.00019513940559561356,1,0.9563425525954127,0.48677249,"carboxy-lyase activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9642486139462519,0.28175041,"carboxy-lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,18,0.9585373428327074,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,15,0.9573490421877254,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,28,0.8699916892191324,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9214819569250295,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004301","epoxide hydrolase activity",0.06371079843372243,1,0.9311110808384813,0.22109293,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9392140747595851,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,11,0.8900502416316907,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008568","microtubule severing ATPase activity",0.013493446398287598,1,0.9039060921890452,0.41923201,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.9014398579567025,0.38812549,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,11,0.9114266549752635,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9116529749120368,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9136811150660402,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,2,0.9316529135975979,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0033919","glucan 1,3-alpha-glucosidase activity",0.000915824710352141,1,0.9276848533170794,0.4867618,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,2,0.9322888130169236,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047668","amygdalin beta-glucosidase activity",2.4392425699451695E-05,1,0.9382288266081679,0.44127145,"ATP hydrolysis activity"),
                     c("GO:0047701","beta-L-arabinosidase activity",4.878485139890339E-05,1,0.9401401842035029,0.41027052,"ATP hydrolysis activity"),
                     c("GO:0047782","coniferin beta-glucosidase activity",2.66099189448564E-05,1,0.9380013667085084,0.44258599,"ATP hydrolysis activity"),
                     c("GO:0047840","dCTP diphosphatase activity",0.00782775115627859,1,0.9364963382573893,0.42129781,"ATP hydrolysis activity"),
                     c("GO:0047884","FAD diphosphatase activity",0.00013970207446049608,1,0.9456770136428667,0.31970647,"ATP hydrolysis activity"),
                     c("GO:0050224","prunasin beta-glucosidase activity",1.55224527178329E-05,1,0.9393837076440363,0.43456671,"ATP hydrolysis activity"),
                     c("GO:0080079","cellobiose glucosidase activity",7.982975683456919E-05,1,0.9349762316979057,0.45988357,"ATP hydrolysis activity"),
                     c("GO:0080083","beta-gentiobiose beta-glucosidase activity",7.76122635891645E-05,1,0.9350575361915731,0.42073911,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9402006720937442,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9387706512195393,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140326","ATPase-coupled intramembrane lipid transporter activity",0.056601515088954966,1,0.9035448513096077,0.47545686,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,2,0.9359460693424815,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,8,0.981632546991765,0.04901613,"carbohydrate binding"),
                     c("GO:0030674","protein-macromolecule adaptor activity",1.2206258094127533,3,0.9818006831662439,-0,"protein-macromolecule adaptor activity"),
                     c("GO:0032556","pyrimidine deoxyribonucleotide binding",9.09172230615927E-05,1,0.9763354793621403,0.07629645,"pyrimidine deoxyribonucleotide binding"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9477041521979613,0.03310522,"DNA demethylase activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,4,0.98155668749802,-0,"signaling receptor activity"),
                     c("GO:0009885","transmembrane histidine kinase cytokinin receptor activity",3.769738517187989E-05,1,0.8767415148035875,0.38178873,"signaling receptor activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,15,0.9437887856978878,0.04095538,"protein homodimerization activity"),
                     c("GO:0000149","SNARE binding",0.26943873427614345,2,0.942792002641888,0.40354822,"protein homodimerization activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,11,0.942883750636101,0.40281405,"protein homodimerization activity"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,13,0.9354328706465737,0.38709915,"protein homodimerization activity"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,2,0.9493128327061413,0.3218275,"protein homodimerization activity"),
                     c("GO:0019901","protein kinase binding",0.4024661540679714,5,0.9328854808754286,0.41867396,"protein homodimerization activity"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.9441555178342247,0.39266904,"protein homodimerization activity"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,3,0.9463011211239948,0.34139999,"protein homodimerization activity"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9481717059123352,0.32923903,"protein homodimerization activity"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,3,0.9449938446832636,0.34990764,"protein homodimerization activity"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9484686286287121,0.32731003,"protein homodimerization activity"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,2,0.9624225318942571,0.2370668,"protein homodimerization activity"),
                     c("GO:0032050","clathrin heavy chain binding",0.02621520514717436,2,0.9500809474712852,0.30626295,"protein homodimerization activity"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,2,0.9434110870812179,0.39860027,"protein homodimerization activity"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9521192316858261,0.48770765,"protein homodimerization activity"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,2,0.9455109213449056,0.33872441,"protein homodimerization activity"),
                     c("GO:0042393","histone binding",0.41579772345934446,6,0.9407521339487543,0.41995229,"protein homodimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,9,0.9398351958719575,0.42737733,"protein homodimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,4,0.9438266869026394,0.3952865,"protein homodimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,6,0.9490761271825151,0.32336445,"protein homodimerization activity"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,9,0.935226970665565,0.4679269,"protein homodimerization activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,11,0.9390050009843673,0.43412799,"protein homodimerization activity"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,2,0.9422377809164217,0.40798975,"protein homodimerization activity"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,3,0.9529726326090413,0.28789489,"protein homodimerization activity"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,2,0.9490154651867558,0.32375837,"protein homodimerization activity"),
                     c("GO:1990935","splicing factor binding",0.0007095978385295039,1,0.961086045252152,0.24567361,"protein homodimerization activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.980052909297772,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9969861288752431,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9808115239171371,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9862914642874989,0.49249984,"ER retention sequence binding"),
                     c("GO:1904408","melatonin binding",0.00021066185831344647,1,0.9793104873168712,0.43339648,"ER retention sequence binding"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,2,0.9623917017344634,0.02253254,"oxalate-CoA ligase activity"),
                     c("GO:0003972","RNA ligase (ATP) activity",0.012360307349885798,1,0.9369948658263713,0.4132865,"oxalate-CoA ligase activity"),
                     c("GO:0003989","acetyl-CoA carboxylase activity",0.08698117255099935,1,0.9486761697043934,0.40308719,"oxalate-CoA ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.955587030443094,0.28284902,"oxalate-CoA ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.963824697656016,0.2351635,"oxalate-CoA ligase activity"),
                     c("GO:0106286","(E)-caffeate-CoA ligase activity",2.4392425699451695E-05,1,0.9642654643218527,0.39917988,"oxalate-CoA ligase activity"),
                     c("GO:0050373","UDP-arabinose 4-epimerase activity",8.20472500799739E-05,1,0.9612338253353215,0.02256811,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0016871","cycloartenol synthase activity",0.0006563780006397912,1,0.9581354922059703,0.28013495,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:2001070","starch binding",0.024815966909324,1,0.9847875682274669,0.03482451,"starch binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,5,0.8693142036127602,0,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,2,0.8166078163982737,0.44154323,"monoatomic ion transport"),
                     c("GO:0035672","oligopeptide transmembrane transport",0.23299928216907387,1,0.8627103740851154,0.43869347,"monoatomic ion transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8501483265362694,0.4188028,"monoatomic ion transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.8648015434858708,0.30464429,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,3,0.8381661356544009,-0,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.8545317476844199,0.4259419,"defense response"),
                     c("GO:0010039","response to iron ion",0.05924742412469929,1,0.8678249553799107,0.46879008,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8452169740379526,0.46097507,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,2,0.7888309661612708,0.37522238,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8500399690180803,0.46171648,"defense response"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8262373553023209,0.42151038,"defense response"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0009820","alkaloid metabolic process",0.01601588160255828,2,0.9436353225576943,0.05107695,"alkaloid metabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8654486740322442,0.36848795,"alkaloid metabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8677155725862653,0.27735821,"alkaloid metabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,2,0.8680221978172313,0.44475214,"alkaloid metabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.91015480010716,0.11683693,"alkaloid metabolic process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,2,0.8973857222610786,-0,"lipid catabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9353800814722798,0.15163831,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.9338811362930068,0.12277613,"lipid catabolic process"),
                     c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.9085864140814529,0.46381088,"lipid catabolic process"),
                     c("GO:0016137","glycoside metabolic process",0.12090536413233209,1,0.9517454394802891,0.06048027,"glycoside metabolic process"),
                     c("GO:0033491","coniferin metabolic process",1.7253181166190823E-05,1,0.9442566754145938,0.44674727,"glycoside metabolic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.9338745778886478,0.07741388,"phosphorylation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9507663781277531,0.11172681,"phosphorylation"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9798385633655547,-0,"chromosome condensation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,0.9829765143020287,-0,"carbohydrate homeostasis"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9455995749748067,-0,"root development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9004928829494658,0.29362419,"root development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9484621531515232,0.43353892,"root development"),
                     c("GO:0048316","seed development",0.03280076213709535,1,0.9264830084500535,0.46716202,"root development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8845783069732632,0.07105297,"DNA demethylation"),
                     c("GO:0006076","(1->3)-beta-D-glucan catabolic process",2.4647401665986893E-06,1,0.9449009301204617,0.22339313,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9011708594395732,0.15925058,"DNA demethylation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9307160410633457,0.23039539,"DNA demethylation"),
                     c("GO:0016926","protein desumoylation",0.02830261133305275,1,0.8950802362388776,0.33184604,"DNA demethylation"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.893709612099068,0.15970956,"DNA demethylation"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9670039793939942,0.08184898,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0045739","positive regulation of DNA repair",0.0330817425160876,1,0.9461004114183711,0.41456906,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:1902456","regulation of stomatal opening",0.0015108857221249963,1,0.9714778756693939,-0,"regulation of stomatal opening"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,1,0.9208134133543636,0.38658972,"regulation of stomatal opening"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.94918833865427,0.11892868,"regulation of stomatal opening"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9277321826146957,0.04070747,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,3,0.9637918578215985,0.07993224,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,3,0.934805058457749,0,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9300274339252134,0.31015716,"zinc ion binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,4,0.8963589951387093,0.26746377,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.9237444654140653,0.37942352,"zinc ion binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9504269344941726,0.12083997,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.9298237490314022,0.18168089,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,1,0.9726340124201915,0.05834827,"lipid binding"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.9292742341183012,-0,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9578984371331922,0.29459182,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.951676000791839,0.35458321,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9497163330768981,0.35290783,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,9,0.9359613624514527,0.0830445,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,4,0.9417497750842001,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.936645291325742,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9493991577398884,0.05997119,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9516616283484537,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9500276011425831,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9490108836223107,0.04472286,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9662945343000091,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,2,0.8210042296808621,-0,"ATP hydrolysis activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.8792213438440544,0.42844795,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8730328155525962,0.35450436,"ATP hydrolysis activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8368571495135634,0.46981862,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,2,0.788574523776987,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.8094718727251017,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,2,0.7802928649580687,0.3656957,"ATP hydrolysis activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8054879923335817,0.45595338,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.854669789283242,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,1,0.866099721971034,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0047782","coniferin beta-glucosidase activity",2.66099189448564E-05,1,0.8784272638229184,0.3973634,"ATP hydrolysis activity"),
                     c("GO:0070139","SUMO-specific endopeptidase activity",0.006703482080858407,1,0.8474368226275719,0.18227358,"ATP hydrolysis activity"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,1,0.8634077729955099,0.37283533,"ATP hydrolysis activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,2,0.983830955048966,-0,"signaling receptor activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9578211223943313,0.05189795,"identical protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9599078328645053,0.39898091,"identical protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9612831540639439,0.38355794,"identical protein binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.925570731080364,0.06961566,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,2,0.9133923051843015,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9377979962929239,0.33073906,"sequence-specific DNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9435955158068591,0.28841983,"sequence-specific DNA binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9220990254326026,0.03322883,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9127134136818348,0.30799911,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0047037","salutaridine reductase (NADPH) activity",1.77399459632376E-05,1,0.9517261712274524,0.38914098,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0047501","(+)-neomenthol dehydrogenase activity",2.66099189448564E-05,1,0.9509287428710922,0.3973742,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9292986410001567,0.05781294,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.947165212862748,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9746685872764549,0.01879944,"pterocarpan synthase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9635594141554411,0.23977693,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)


if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,2,0.9999102062499213,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9907586081680615,2.598E-05,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,10,0.8713615860984057,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,6,0.8466445742055034,0.17160779,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,8,0.985458991941732,0.06338235,"plasma membrane"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.9150687312533304,-0,"riboflavin synthase complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,1,0.8480183102844111,0.34749127,"riboflavin synthase complex"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7560655280977256,0,"chloroplast"),
                     c("GO:0000793","condensed chromosome",0.4223503377863461,1,0.8081396817689885,0.44699462,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,10,0.6984134357909436,0.35145602,"chloroplast"),
                     c("GO:0005694","chromosome",2.5245804089932373,1,0.7711269266077386,0.15404313,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7401638138584955,0.25147478,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.776927919424568,0.18302188,"chloroplast"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7293294274784751,0.18309297,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,1,0.7906266379318838,0.18435889,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7213112123835315,0.22703414,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,1,0.7901992768240182,0.18475413,"chloroplast"),
                     c("GO:0016607","nuclear speck",0.16023998054727537,1,0.7668007904886469,0.14877708,"chloroplast"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7615264882153783,0.21139748,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.999928779234615,3.782E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,13,0.9998529851438804,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999416045492258,3.171E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.07909597668631853,3,0.9217149565546587,0.05769476,"malate metabolic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,1,0.9229464063040455,0.3067044,"malate metabolic process"),
                     c("GO:0045487","gibberellin catabolic process",0.00595974172283563,2,0.8762033265144406,0.26507989,"malate metabolic process"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,10,0.9122748819072783,-0,"monoatomic ion transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9432813994604033,0.22103265,"monoatomic ion transport"),
                     c("GO:0015031","protein transport",3.093438694074183,4,0.8694349327399881,0.42145531,"monoatomic ion transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9349354269439446,0.19142099,"monoatomic ion transport"),
                     c("GO:0015783","GDP-fucose transmembrane transport",0.013410651246463469,2,0.916147426399219,0.42609188,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,4,0.8762118526965018,0.44154323,"monoatomic ion transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9449366276391375,0.16068432,"monoatomic ion transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9062455303621825,0.44310395,"monoatomic ion transport"),
                     c("GO:0051646","mitochondrion localization",0.04394385243028803,1,0.9298546266470985,0.24160098,"monoatomic ion transport"),
                     c("GO:0072583","clathrin-dependent endocytosis",0.06990989008540521,2,0.9308555064685979,0.26649209,"monoatomic ion transport"),
                     c("GO:0072699","protein localization to cortical microtubule cytoskeleton",0.00015774337066231612,1,0.9433731536634987,0.15867185,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,20,0.8876272230726193,0,"defense response"),
                     c("GO:0000165","MAPK cascade",0.14246691110973742,2,0.8588402128071748,0.45634702,"defense response"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8956770006647753,0.4896406,"defense response"),
                     c("GO:0006955","immune response",0.6433982378290884,2,0.9111344210371278,0.30172547,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8637303629479638,0.44951817,"defense response"),
                     c("GO:0007300","ovarian nurse cell to oocyte transport",0.00031548674132463224,1,0.8738949051222475,0.46428331,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,2,0.9035040983547749,0.45786068,"defense response"),
                     c("GO:0009641","shade avoidance",0.0003204162216578296,1,0.9195167266622747,0.49975625,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,3,0.8965158682779218,0.22702886,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,1,0.9222402904495552,0.48975475,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9207574345423947,0.20406574,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8872073311603659,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,3,0.8183183765622872,0.37522238,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,6,0.882390566826153,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9813132397214319,0.48259105,"defense response"),
                     c("GO:0048573","photoperiodism, flowering",0.0014763793597926149,2,0.8198788496978948,0.46247226,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9183645836820575,0.26100449,"defense response"),
                     c("GO:0071456","cellular response to hypoxia",0.026404761404771757,4,0.859551573207622,0.39162542,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9085482240172137,0.40112537,"defense response"),
                     c("GO:0007049","cell cycle",2.8073119376140743,2,0.9904146740988637,0.01433045,"cell cycle"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,1,0.9955431861403087,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,4,1,-0,"metabolic process"),
                     c("GO:0008219","cell death",0.4651211168388386,1,0.9918986929793846,0.01044746,"cell death"),
                     c("GO:0015979","photosynthesis",0.228607115192195,2,0.9607037950593874,0.06311855,"photosynthesis"),
                     c("GO:0016477","cell migration",0.49093927008395993,1,0.9918610395838999,0.01050372,"cell migration"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,12,0.9027159603986793,-0,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,7,0.9465354534227138,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.928854251766739,0.44424972,"protein ubiquitination"),
                     c("GO:0006004","fucose metabolic process",0.06850745293061056,2,0.9071417508253294,0.43756433,"protein ubiquitination"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.8923138496616365,0.4768145,"protein ubiquitination"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,3,0.8025323564716771,0.41940643,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,2,0.9033179992280823,0.21165279,"protein ubiquitination"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8896195004172376,0.37832162,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.8915054911061395,0.37292242,"protein ubiquitination"),
                     c("GO:0006508","proteolysis",5.2622572267907,5,0.9046114465127817,0.44941082,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9465016164309483,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.945311378466218,0.12094833,"protein ubiquitination"),
                     c("GO:0006694","steroid biosynthetic process",0.2949431320360321,2,0.8659733844576685,0.33431465,"protein ubiquitination"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,3,0.9575436664180312,0.12823919,"protein ubiquitination"),
                     c("GO:0008380","RNA splicing",0.9742871404547958,1,0.8911240835883005,0.39771984,"protein ubiquitination"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9596258578399225,0.21768792,"protein ubiquitination"),
                     c("GO:0009073","aromatic amino acid family biosynthetic process",0.5154831526629496,1,0.8632746428917868,0.41955763,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.898647303381967,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,24,0.9201586061779947,0.46324041,"protein ubiquitination"),
                     c("GO:0018345","protein palmitoylation",0.0003918936864891916,1,0.927482938228154,0.37471442,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,4,0.8903472165832508,0.14129671,"protein ubiquitination"),
                     c("GO:0031146","SCF-dependent proteasomal ubiquitin-dependent protein catabolic process",0.12335531585793119,2,0.9054990986254112,0.46794754,"protein ubiquitination"),
                     c("GO:0033314","mitotic DNA replication checkpoint signaling",0.042962885843981745,1,0.8389238091593698,0.4495745,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9474209773171544,0.45767193,"protein ubiquitination"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8653016314022828,0.4941006,"protein ubiquitination"),
                     c("GO:0042853","L-alanine catabolic process",0.0402368832197236,1,0.8984918709082395,0.42825993,"protein ubiquitination"),
                     c("GO:0046488","phosphatidylinositol metabolic process",0.5598188987797269,2,0.8680597043051811,0.47629275,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9247723766968365,0.40604291,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9264348474927936,0.20437656,"protein ubiquitination"),
                     c("GO:0055047","generative cell mitosis",7.147746483136198E-05,1,0.8982388243261324,0.45916025,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9161100139273048,0.45115494,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9444518331173675,0.45105203,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,1,0.9094524352519981,0.21493739,"protein ubiquitination"),
                     c("GO:1990918","double-strand break repair involved in meiotic recombination",0.009479390680738558,1,0.7765808641048081,0.40319107,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9927490803456864,0.00919456,"establishment of cell polarity"),
                     c("GO:0030154","cell differentiation",2.279586420543629,3,0.8990410185370589,0.01240155,"cell differentiation"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9143084839834023,0.48841036,"cell differentiation"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.9219472873147521,0.48986076,"cell differentiation"),
                     c("GO:0055046","microgametogenesis",0.0014049018949612527,1,0.9243587561470956,0.40770388,"cell differentiation"),
                     c("GO:0032259","methylation",2.6278542060840238,1,0.97034765211817,0.05843136,"methylation"),
                     c("GO:0040011","locomotion",0.5148398554794674,1,1,-0,"locomotion"),
                     c("GO:0048544","recognition of pollen",0.04932191547380636,1,0.9931984117235034,0.00854634,"recognition of pollen"),
                     c("GO:0050896","response to stimulus",17.567785530535815,4,1,-0,"response to stimulus"),
                     c("GO:0060361","flight",9.612486649734888E-05,1,1,-0,"flight"),
                     c("GO:0070814","hydrogen sulfide biosynthetic process",0.03377679924306844,1,0.9402148354754484,0.03784896,"hydrogen sulfide biosynthetic process"),
                     c("GO:0000103","sulfate assimilation",0.10156701278519879,1,0.9563161258332863,0.48363934,"hydrogen sulfide biosynthetic process"),
                     c("GO:0042350","GDP-L-fucose biosynthetic process",0.02421607213683212,1,0.8835380103040614,0.34490714,"hydrogen sulfide biosynthetic process"),
                     c("GO:0046166","glyceraldehyde-3-phosphate biosynthetic process",0.03550458209985412,1,0.9116980138604488,0.11775102,"hydrogen sulfide biosynthetic process"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,8,0.9197627632454867,-0,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,3,0.9293852555636286,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.942052264292397,0.37383592,"cell wall organization"),
                     c("GO:0009662","etioplast organization",0.0006778035458146395,1,0.9517053599085655,0.22083297,"cell wall organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9401213080094034,0.32062999,"cell wall organization"),
                     c("GO:0051013","microtubule severing",0.012493767904488756,1,0.9408521597191466,0.26779102,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9495513901672242,0.4524177,"cell wall organization"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9798816262649623,-0,"intracellular auxin homeostasis"),
                     c("GO:0098771","inorganic ion homeostasis",1.0570555799893464,1,0.9677054967424179,0.45565197,"intracellular auxin homeostasis"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,3,0.915753539696631,-0,"regulation of defense response to fungus"),
                     c("GO:0010105","negative regulation of ethylene-activated signaling pathway",0.005397780964851129,1,0.9119840385339306,0.43825323,"regulation of defense response to fungus"),
                     c("GO:1900457","regulation of brassinosteroid mediated signaling pathway",0.001054908791304239,1,0.9314577900279192,0.34087448,"regulation of defense response to fungus"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9228448013437632,0.32535433,"regulation of defense response to fungus"),
                     c("GO:1901001","negative regulation of response to salt stress",0.00037710524548959944,1,0.9192906002379682,0.46916415,"regulation of defense response to fungus"),
                     c("GO:2000022","regulation of jasmonic acid mediated signaling pathway",0.011707515791343773,1,0.9222577951337201,0.42518681,"regulation of defense response to fungus"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,1,0.9233401420289508,0.37952976,"regulation of defense response to fungus"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,2,0.9238018124888107,0.37732875,"regulation of defense response to fungus"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,1,0.9572395386714166,0.08761323,"regulation of secondary cell wall biogenesis"),
                     c("GO:1903527","positive regulation of membrane tubulation",0.0019298915504467734,1,0.9338975824675887,0.46246533,"regulation of secondary cell wall biogenesis"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9364922075877135,0.40819627,"regulation of secondary cell wall biogenesis"),
                     c("GO:2001006","regulation of cellulose biosynthetic process",0.002356291599268347,2,0.947090576042956,0.0943689,"regulation of cellulose biosynthetic process"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9214809266030887,0.38244242,"regulation of cellulose biosynthetic process"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9488957450726474,0.37175537,"regulation of cellulose biosynthetic process"),
                     c("GO:0031047","regulatory ncRNA-mediated gene silencing",0.26322685557224024,2,0.9070631389466772,0.33664393,"regulation of cellulose biosynthetic process"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,2,0.899438868614809,0.23293868,"regulation of cellulose biosynthetic process"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9414309705668448,0.18590747,"regulation of cellulose biosynthetic process"),
                     c("GO:0051510","regulation of unidimensional cell growth",0.00424921204721614,1,0.9265414734291049,0.45200802,"regulation of cellulose biosynthetic process"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9331342057273261,0.13410528,"regulation of cellulose biosynthetic process"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,1,0.9506863636529148,0.15074542,"regulation of cellulose biosynthetic process"),
                     c("GO:1901141","regulation of lignin biosynthetic process",0.00025140349699306625,1,0.9510187473965689,0.14312882,"regulation of cellulose biosynthetic process"),
                     c("GO:1905038","regulation of membrane lipid metabolic process",0.009506502822571143,1,0.9380563272407281,0.25241702,"regulation of cellulose biosynthetic process"),
                     c("GO:1905157","positive regulation of photosynthesis",0.00010598382716374363,1,0.9456428144310399,0.39735289,"regulation of cellulose biosynthetic process"),
                     c("GO:1905183","negative regulation of protein serine/threonine phosphatase activity",1.7253181166190823E-05,1,0.9463200482012567,0.1289891,"regulation of cellulose biosynthetic process"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9433671437873972,0.49983935,"regulation of cellulose biosynthetic process"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9442925564043755,0.347328,"regulation of cellulose biosynthetic process"),
                     c("GO:2000306","positive regulation of photomorphogenesis",0.00015281389032911874,1,0.9191016589735107,0.33080078,"regulation of cellulose biosynthetic process"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9296953924094261,0.14897448,"regulation of cellulose biosynthetic process"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9396520077476884,0.26299316,"regulation of cellulose biosynthetic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}





dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,9,0.9999164107002521,4.598E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,37,0.888698863151965,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,18,0.8576699668760611,0.17160779,"cytoplasm"),
                     c("GO:0009504","cell plate",0.0007255512462560189,1,0.9373247065460408,0.05708197,"cell plate"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,1,0.9709662236773764,2.997E-05,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,9,0.9501447151073252,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,22,0.7106448688607407,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,2,0.804497965675289,0.13280054,"chloroplast"),
                     c("GO:0001673","male germ cell nucleus",0.006684016617906476,1,0.8360702504197768,0.1163509,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,59,0.7074440060782706,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,3,0.7432982904889285,0.2503098,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,2,0.7873356189967988,0.17354881,"chloroplast"),
                     c("GO:0005788","endoplasmic reticulum lumen",0.11812570632538061,4,0.7056582812642246,0.47469365,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,17,0.649438687571789,0.19658734,"chloroplast"),
                     c("GO:0005797","Golgi medial cisterna",0.015012947362598856,2,0.738297847695811,0.39669825,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,5,0.700122434902463,0.15783714,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,3,0.6941870352481858,0.49095344,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,18,0.7869536915557038,0.17236349,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,3,0.7599599827038843,0.43399423,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.7759755773508121,0.22189912,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7519763330461235,0.30090248,"chloroplast"),
                     c("GO:0012505","endomembrane system",6.893425119210306,5,0.9999099064894741,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,116,0.9998758904921088,9.537E-05,"membrane"),
                     c("GO:0019005","SCF ubiquitin ligase complex",0.14169866753507532,3,0.9091563726356336,-0,"SCF ubiquitin ligase complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,2,0.8257573386651481,0.19052614,"SCF ubiquitin ligase complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7852224787846288,0.38761376,"SCF ubiquitin ligase complex"),
                     c("GO:0005886","plasma membrane",17.177321395000487,44,0.9129233654725646,0.34881718,"SCF ubiquitin ligase complex"),
                     c("GO:0009898","cytoplasmic side of plasma membrane",0.2657902935257323,1,0.9378338399477214,0.20588623,"SCF ubiquitin ligase complex"),
                     c("GO:0010330","cellulose synthase complex",0.0019629639881584074,1,0.9333903733494184,0.38954564,"SCF ubiquitin ligase complex"),
                     c("GO:0016602","CCAAT-binding factor complex",0.03082101937862897,1,0.7691215119887675,0.41728323,"SCF ubiquitin ligase complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.7508330451999053,0.19254961,"SCF ubiquitin ligase complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,2,0.8719318836064335,0.21998964,"SCF ubiquitin ligase complex"),
                     c("GO:0045273","respiratory chain complex II",0.05444864540824706,1,0.8811186657595379,0.46599171,"SCF ubiquitin ligase complex"),
                     c("GO:0046930","pore complex",0.06221601936645362,1,0.8913857435781342,0.47090253,"SCF ubiquitin ligase complex"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.99992041171445,5.465E-05,"cell projection"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9906499036338272,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9064298792259897,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9913127437403347,2.369E-05,"pollen tube"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,7,0.9419628176739984,0.066121,"transcription cis-regulatory region binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,5,0.9499660975376847,0.22593798,"transcription cis-regulatory region binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,24,0.94105616324916,0.45223656,"transcription cis-regulatory region binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9665228957965674,0.31120183,"transcription cis-regulatory region binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,16,0.9864138522601702,0.40628023,"transcription cis-regulatory region binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,3,0.9562940102124945,0.31664803,"transcription cis-regulatory region binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9430915716595631,0.27764537,"transcription cis-regulatory region binding"),
                     c("GO:0019237","centromeric DNA binding",0.004625690909914204,1,0.9667754581017608,0.49059442,"transcription cis-regulatory region binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9651443006653624,0.14801235,"transcription cis-regulatory region binding"),
                     c("GO:0003674","molecular_function",100,10,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,1,0.9738454051056479,0.0533165,"chromatin binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,11,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,23,0.7744996333649337,0,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,2,0.8845724263170435,0.34856872,"protein kinase activity"),
                     c("GO:0003975","UDP-N-acetylglucosamine-dolichyl-phosphate N-acetylglucosaminephosphotransferase activity",0.005718915079898721,1,0.8976980242425416,0.35348317,"protein kinase activity"),
                     c("GO:0004106","chorismate mutase activity",0.044957458057334886,1,0.9475028729130611,0.47467998,"protein kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.8991077512603162,0.25688843,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8337664879965059,0.49507895,"protein kinase activity"),
                     c("GO:0004781","sulfate adenylyltransferase (ATP) activity",0.04230755362907627,1,0.8822900625187642,0.41713542,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,1,0.8725292313409576,0.35677507,"protein kinase activity"),
                     c("GO:0008236","serine-type peptidase activity",1.2663216927208079,2,0.8316623479874524,0.49481639,"protein kinase activity"),
                     c("GO:0008373","sialyltransferase activity",0.03638684666384572,1,0.8914279426635238,0.43994578,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,2,0.8750557656849572,0.30161168,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,1,0.8443049032880112,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,6,0.8672270504187616,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,12,0.8727657548687641,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,3,0.8851166610554807,0.21201958,"protein kinase activity"),
                     c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.04143767537877,1,0.8564528272494051,0.43304972,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.8910492855322953,0.37780808,"protein kinase activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9559696927130285,0.38075906,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,1,0.8730438925023081,0.48733868,"protein kinase activity"),
                     c("GO:0050734","hydroxycinnamoyltransferase activity",0.0034969868480032116,1,0.8981091164002102,0.18416899,"protein kinase activity"),
                     c("GO:0051753","mannan synthase activity",0.005627997856837128,1,0.8996393290260886,0.47286456,"protein kinase activity"),
                     c("GO:0061630","ubiquitin protein ligase activity",0.6665452071699717,5,0.8056777694326984,0.45681379,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,1,0.9952182751637113,-0,"GTPase activator activity"),
                     c("GO:0005515","protein binding",8.610051728351934,46,0.9716427288876132,0.07766851,"protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,6,0.9166999090790998,0.05255315,"microtubule binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,5,0.9309748412630602,0.40281405,"microtubule binding"),
                     c("GO:0019901","protein kinase binding",0.4024661540679714,2,0.92603792544991,0.41867396,"microtubule binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9335958346372112,0.3860137,"microtubule binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9443981783814701,0.31827679,"microtubule binding"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,1,0.9357296086277153,0.37244819,"microtubule binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9283267326886423,0.41995229,"microtubule binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,3,0.9271877990120283,0.42737733,"microtubule binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,5,0.9307049517133249,0.38709915,"microtubule binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9321461934406067,0.3952865,"microtubule binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9386618294523559,0.35396093,"microtubule binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,5,0.9214704781073605,0.4651985,"microtubule binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9261568333118312,0.4679269,"microtubule binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.930172341394325,0.40798975,"microtubule binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,11,0.9405145586980309,-0,"zinc ion binding"),
                     c("GO:0000035","acyl binding",0.03278342014006308,1,0.9712133733070926,0.11604884,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,1,0.9497634817826637,0.3830934,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,3,0.9512389854106326,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,2,0.9321986694039569,0.35257211,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,3,0.9451137395441458,0.13873458,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9527322215068251,0.20283352,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,1,0.9381654281197742,0.31639809,"zinc ion binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9635077978033577,0.17458707,"zinc ion binding"),
                     c("GO:0043325","phosphatidylinositol-3,4-bisphosphate binding",0.013358179310317913,1,0.964338608907905,0.17101403,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,1,0.9341467725870477,0.31135007,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,1,0.944874429694901,0.3462071,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,2,0.9438052293524416,0.16293278,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9466715608305454,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9479750821545041,0.28982977,"zinc ion binding"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,1,0.9821431196277406,0.02993406,"auxin binding"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,4,0.9090083182364842,-0,"antiporter activity"),
                     c("GO:0005338","nucleotide-sugar transmembrane transporter activity",0.0598301852542642,2,0.9193689970476834,0.33273016,"antiporter activity"),
                     c("GO:0008553","P-type proton-exporting transporter activity",0.01860920331543624,1,0.8966669279706063,0.47571857,"antiporter activity"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.8847836857049981,0.39050958,"antiporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9238202754911319,0.35182404,"antiporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.9525787792498923,0.20645071,"antiporter activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.9109513286038611,0.3876045,"antiporter activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9055820684265667,0.43240123,"antiporter activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9263030778312501,0.28474685,"antiporter activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,3,0.8796902445130788,0.0531541,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004324","ferredoxin-NADP+ reductase activity",0.018482806200448173,1,0.9281374100439902,0.30749203,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9024394673204892,0.46480825,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9268905757653264,0.31440091,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,1,0.9058668394137597,0.41891471,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9011585081526731,0.47359913,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.906941954769149,0.49506653,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0045543","gibberellin 2-beta-dioxygenase activity",0.006071496505918068,1,0.9167721931620054,0.48062226,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9069600712008533,0.43460917,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102111","gibberellin A20,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102652","gibberellin A9,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102772","sphingolipid C4-monooxygenase activity",4.4349864908093996E-05,1,0.9380050977762479,0.28012766,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102924","gibberellin A44,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.20747889,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,55,0.9437836522476242,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,18,0.9479842952457235,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,25,0.9432159328868454,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,7,0.9543080049624433,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,7,0.9562022310468056,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,5,0.9548328537234818,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,5,0.871929271156829,0.05997119,"ATP hydrolysis activity"),
                     c("GO:0004301","epoxide hydrolase activity",0.06371079843372243,1,0.924843077902784,0.22109293,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,3,0.8864327502617124,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8607554481442179,0.454689,"ATP hydrolysis activity"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,2,0.86706935078231,0.32020834,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,4,0.9048402323081485,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9050971095270066,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9058211895932066,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0140326","ATPase-coupled intramembrane lipid transporter activity",0.056601515088954966,1,0.9049973218123009,0.47545686,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,4,0.9781376477450502,0.05647972,"carbohydrate binding"),
                     c("GO:0030527","structural constituent of chromatin",0.19884040182219404,1,1,-0,"structural constituent of chromatin"),
                     c("GO:0030674","protein-macromolecule adaptor activity",1.2206258094127533,1,1,-0,"protein-macromolecule adaptor activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,4,0.9740485483191771,-0,"signaling receptor activity"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9966774087712642,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8780311780345175,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,1,0.9719635630397723,0.02253254,"oxalate-CoA ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9645846351874537,0.44358297,"oxalate-CoA ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9539150239545147,0.33857975,"oxalate-CoA ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9712579065613889,0.4447539,"oxalate-CoA ligase activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9331203890479748,0.0567159,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9455587216231844,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9688129977521718,0.01868194,"pterocarpan synthase activity"),
                     c("GO:0004383","guanylate cyclase activity",0.041591303310810554,1,0.9425222402969877,0.49730148,"pterocarpan synthase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,1,0.9387243910696051,0.48814337,"pterocarpan synthase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9507451098580505,0.32867735,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
c("GO:0005576","extracellular region",3.0008202208427486,4,0.9998992004886821,5.282E-05,"extracellular region"),
c("GO:0009505","plant-type cell wall",0.013780776369110043,1,0.9872638457957056,3.275E-05,"plant-type cell wall"),
c("GO:0005886","plasma membrane",14.184672303876209,6,0.9722381493674755,0.29250802,"plant-type cell wall"),
c("GO:0009506","plasmodesma",0.02618690220789438,1,0.9999359590919521,3.43E-05,"plasmodesma"),
c("GO:0009536","plastid",0.5104294469183077,7,0.6823513310670152,0,"plastid"),
c("GO:0005634","nucleus",12.339605699963672,13,0.578123451304178,0.44231542,"plastid"),
c("GO:0005737","cytoplasm",25.080076249313628,5,0.7543622687518001,0.22725924,"plastid"),
c("GO:0005739","mitochondrion",2.7614025547937224,3,0.6289338728049949,0.33771423,"plastid"),
c("GO:0005764","lysosome",0.2904358594155184,1,0.6782503799265944,0.22601794,"plastid"),
c("GO:0005773","vacuole",0.6250852016776067,1,0.6767722261747079,0.28151423,"plastid"),
c("GO:0005783","endoplasmic reticulum",2.123313387573007,4,0.5768419972738681,0.27526608,"plastid"),
c("GO:0005811","lipid droplet",0.07544965542347233,1,0.7019048550948056,0.4118361,"plastid"),
c("GO:0005829","cytosol",1.696303522837096,3,0.7330563829423168,0.22154616,"plastid"),
c("GO:0005840","ribosome",3.532158825827589,2,0.591537393537252,0.22161604,"plastid"),
c("GO:0009507","chloroplast",0.4142191413803319,5,0.6059529774655197,0.2334722,"plastid"),
c("GO:0010369","chromocenter",0.0029358879747399397,1,0.7574317221125776,0.31480507,"plastid"),
c("GO:0012511","monolayer-surrounded lipid storage body",0.00680851841612842,1,0.7451383986254116,0.33527212,"plastid"),
c("GO:0022625","cytosolic large ribosomal subunit",0.11811713241045162,2,0.6404092045880214,0.43014036,"plastid"),
c("GO:0030018","Z disc",0.03448430804052515,1,0.6884421737859988,0.3833401,"plastid"),
c("GO:0031966","mitochondrial membrane",1.347511654066338,1,0.55442047181033,0.30801698,"plastid"),
c("GO:0043661","peribacteroid membrane",3.807896205888378E-05,1,0.7718186759491326,0.12528944,"plastid"),
c("GO:0016020","membrane",61.57163680895251,19,0.9998372235275705,0.00011624,"membrane"),
c("GO:1990904","ribonucleoprotein complex",4.368586074828206,2,0.9163188984692716,-0,"ribonucleoprotein complex"),
c("GO:0005854","nascent polypeptide-associated complex",0.01439765555446396,1,0.7808724333896949,0.26990927,"ribonucleoprotein complex"),
c("GO:0009317","acetyl-CoA carboxylase complex",0.08877729214408166,1,0.7457109324693386,0.3836939,"ribonucleoprotein complex"),
c("GO:0019005","SCF ubiquitin ligase complex",0.043566140491568935,1,0.9321001607798933,0.29740373,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])
hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)
dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0001671","ATPase activator activity",0.0341266115541096,1,1,-0,"ATPase activator activity"),
c("GO:0003700","DNA-binding transcription factor activity",4.532603540907692,4,0.974517641048897,-0,"DNA-binding transcription factor activity"),
c("GO:0003735","structural constituent of ribosome",2.2675249562595297,2,0.9876859144497072,-0,"structural constituent of ribosome"),
c("GO:0003989","acetyl-CoA carboxylase activity",0.09778675834076506,1,0.9680756504817798,0.03241463,"acetyl-CoA carboxylase activity"),
c("GO:0005515","protein binding",5.073796515930947,5,0.9671970182257469,0.07857259,"protein binding"),
c("GO:0008270","zinc ion binding",3.776138036817838,2,0.9118799298452052,0.06655862,"zinc ion binding"),
c("GO:0000287","magnesium ion binding",1.7504660637586555,1,0.9330218898914651,0.37813207,"zinc ion binding"),
c("GO:0005509","calcium ion binding",1.4167951420399318,1,0.9346426100333597,0.36751657,"zinc ion binding"),
c("GO:0008198","ferrous iron binding",0.04409462238483419,1,0.9422771150598274,0.47821211,"zinc ion binding"),
c("GO:0043169","cation binding",18.51293079024014,1,0.9318775708616098,0.27606386,"zinc ion binding"),
c("GO:0016491","oxidoreductase activity",11.721830439531795,7,0.948754590514374,0.09355129,"oxidoreductase activity"),
c("GO:0016740","transferase activity",21.216634542374862,13,0.9444826589956185,0.07040944,"transferase activity"),
c("GO:0016787","hydrolase activity",21.340557702986732,6,0.9444375928461459,0.11166207,"transferase activity"),
c("GO:0016757","glycosyltransferase activity",2.3960575386602994,5,0.8658684166356257,0,"glycosyltransferase activity"),
c("GO:0004185","serine-type carboxypeptidase activity",0.2130250273078274,1,0.8738354742290061,0.37435024,"glycosyltransferase activity"),
c("GO:0004478","methionine adenosyltransferase activity",0.0400000828170345,1,0.9046567569670109,0.21140702,"glycosyltransferase activity"),
c("GO:0004674","protein serine/threonine kinase activity",1.268552650275868,3,0.8008849369801089,0.30164209,"glycosyltransferase activity"),
c("GO:0004721","phosphoprotein phosphatase activity",0.5987290200478247,1,0.8497554998303689,0.41513428,"glycosyltransferase activity"),
c("GO:0008168","methyltransferase activity",2.684244473270584,1,0.8475540584513593,0.34937253,"glycosyltransferase activity"),
c("GO:0016746","acyltransferase activity",3.4236888973740025,1,0.8609030877987566,0.34378837,"glycosyltransferase activity"),
c("GO:0016763","pentosyltransferase activity",0.5429555612510428,1,0.8178482066765651,0.4643706,"glycosyltransferase activity"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.40741643979158,1,0.8467183847227898,0.41755188,"glycosyltransferase activity"),
c("GO:0016776","phosphotransferase activity, phosphate group as acceptor",0.3442981997320215,1,0.8730271763372275,0.45182766,"glycosyltransferase activity"),
c("GO:0050505","hydroquinone glucosyltransferase activity",0.00010352129313132951,1,0.8705521572792748,0.38255741,"glycosyltransferase activity"),
c("GO:0050734","hydroxycinnamoyltransferase activity",0.004307030643174526,1,0.9138157539858599,0.3128147,"glycosyltransferase activity"),
c("GO:0051753","mannan synthase activity",0.002838662827443299,1,0.859925513385278,0.34619078,"glycosyltransferase activity"),
c("GO:0052636","arabinosyltransferase activity",0.0025662383718345365,1,0.8681219548068662,0.17080964,"glycosyltransferase activity"),
c("GO:0102406","omega-hydroxypalmitate O-sinapoyl transferase activity",9.534855946306665E-05,1,0.9291050517417418,0.13884408,"glycosyltransferase activity"),
c("GO:0016829","lyase activity",3.4809062057855114,1,0.9557116593592981,0.05103459,"lyase activity"),
c("GO:0016838","carbon-oxygen lyase activity, acting on phosphates",0.12872327951969606,1,0.9306632613945891,0.03327247,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0008948","oxaloacetate decarboxylase activity",0.024681655678153825,1,0.9436743477893873,0.43322627,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0050552","(4S)-limonene synthase activity",1.3621222780438093E-05,1,0.9436599512926328,0.40007135,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0102903","gamma-terpinene synthase activity",2.7242445560876187E-06,1,0.9484835211799921,0.35927317,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0016887","ATP hydrolysis activity",3.291701972876113,3,0.8903472934655432,0.04836706,"ATP hydrolysis activity"),
c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.7350856328136782,1,0.9135470321751943,0.29181845,"ATP hydrolysis activity"),
c("GO:0052793","pectin acetylesterase activity",0.0008881037252845636,1,0.9327902018153947,0.16150403,"ATP hydrolysis activity"),
c("GO:0032453","histone H3K4 demethylase activity",0.010000701765397649,1,0.9133527628479083,0.02670362,"histone H3K4 demethylase activity"),
c("GO:0042910","xenobiotic transmembrane transporter activity",0.1732101931206069,1,0.967766881976088,0.01034163,"xenobiotic transmembrane transporter activity"),
c("GO:0015171","amino acid transmembrane transporter activity",0.2472497116659562,1,0.9670156354346987,0.35414066,"xenobiotic transmembrane transporter activity"),
c("GO:0043565","sequence-specific DNA binding",1.6398235953577127,3,0.9421291030608252,-0,"sequence-specific DNA binding"),
c("GO:0000166","nucleotide binding",19.572068037245,9,0.9207263260943546,0.47112053,"sequence-specific DNA binding"),
c("GO:0000976","transcription cis-regulatory region binding",0.37429485653910227,2,0.9502676073579898,0.38815494,"sequence-specific DNA binding"),
c("GO:0003677","DNA binding",11.827531128307996,6,0.9284325863898248,0.41192319,"sequence-specific DNA binding"),
c("GO:0003729","mRNA binding",0.29463794571910035,2,0.9535506502136918,0.25887051,"sequence-specific DNA binding"),
c("GO:0005525","GTP binding",1.9488618978002925,1,0.923165453996413,0.34132961,"sequence-specific DNA binding"),
c("GO:0020037","heme binding",1.466676059861896,3,0.9579581837897007,0.13335963,"sequence-specific DNA binding"),
c("GO:0031418","L-ascorbic acid binding",0.07127168607636428,1,0.947832325379137,0.20191161,"sequence-specific DNA binding"),
c("GO:0051287","NAD binding",0.9436456232940781,1,0.9445296330886084,0.12665645,"sequence-specific DNA binding"),
c("GO:0051087","protein-folding chaperone binding",0.14845770708399478,2,0.9298292504190665,0.04628281,"protein-folding chaperone binding"),
c("GO:0030544","Hsp70 protein binding",0.04413276180861943,1,0.9354711628598626,0.39749369,"protein-folding chaperone binding"),
c("GO:0031072","heat shock protein binding",0.10834320599560461,1,0.9313805060735177,0.42425405,"protein-folding chaperone binding"),
c("GO:0043621","protein self-association",0.01613842475026305,1,0.9395344681131841,0.37126244,"protein-folding chaperone binding"),
c("GO:0051082","unfolded protein binding",0.4145292243879603,2,0.9243039485814422,0.47170041,"protein-folding chaperone binding"),
c("GO:0140297","DNA-binding transcription factor binding",0.06685568565094625,1,0.9336380105890657,0.40943701,"protein-folding chaperone binding"),
c("GO:0051213","dioxygenase activity",0.6317060003992653,2,0.881211831865022,0.03928991,"dioxygenase activity"),
c("GO:0004497","monooxygenase activity",1.2382263598775005,1,0.8749229607046218,0.43585558,"dioxygenase activity"),
c("GO:0004601","peroxidase activity",0.4509741680593005,2,0.8613915519395966,0.36163411,"dioxygenase activity"),
c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.3673746215479463,2,0.8542871256794812,0.40450551,"dioxygenase activity"),
c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.4015829604387384,1,0.8736960818213074,0.44216553,"dioxygenase activity"),
c("GO:0016712","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",0.05945663743661228,1,0.8691652305490726,0.37156298,"dioxygenase activity"),
c("GO:0033771","licodione synthase activity",1.3621222780438093E-05,1,0.9075974940010506,0.38603396,"dioxygenase activity"),
c("GO:0047890","flavanone 4-reductase activity",1.906971189261333E-05,1,0.9195822764721437,0.41114057,"dioxygenase activity"),
c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.451820100478857E-05,1,0.9055300158092532,0.25783003,"dioxygenase activity"),
c("GO:0102469","naringenin 2-hydroxylase activity",2.7242445560876187E-06,1,0.912862992625522,0.36819953,"dioxygenase activity"),
c("GO:0102472","eriodictyol 2-hydroxylase activity",2.7242445560876187E-06,1,0.912862992625522,0.16834365,"dioxygenase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])
hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)




library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006123","mitochondrial electron transport, cytochrome c to oxygen",0.047058804746273615,1,0.9708155666584158,0.05401416,"mitochondrial electron transport, cytochrome c to oxygen"),
c("GO:0006457","protein folding",1.0590625918025878,2,0.9915181497788998,0.01051422,"protein folding"),
c("GO:0006730","one-carbon metabolic process",0.3470171077448388,1,0.9520586567451994,0.02817366,"one-carbon metabolic process"),
c("GO:0006108","malate metabolic process",0.054951831581367426,1,0.9516425963768287,0.2431141,"one-carbon metabolic process"),
c("GO:0006865","amino acid transport",0.6097621009115332,2,0.9183573676658028,-0,"amino acid transport"),
c("GO:0006811","monoatomic ion transport",4.401907471506606,1,0.9333500838758619,0.33548775,"amino acid transport"),
c("GO:0042908","xenobiotic transport",0.24723501249579571,1,0.9503034467272828,0.24847794,"amino acid transport"),
c("GO:0098655","monoatomic cation transmembrane transport",2.923892461650467,1,0.9104411394803963,0.41450816,"amino acid transport"),
c("GO:0006952","defense response",1.0410213876080876,4,0.8625541555957023,-0,"defense response"),
c("GO:0009739","response to gibberellin",0.015047641551607423,1,0.8875402733962205,0.47319047,"defense response"),
c("GO:0010044","response to aluminum ion",0.00414442117004926,1,0.8964838355643986,0.4322357,"defense response"),
c("GO:0098869","cellular oxidant detoxification",0.8170663139318542,2,0.8466371269544272,0.34943531,"defense response"),
c("GO:0007017","microtubule-based process",0.8354933550892162,1,0.9917066152132382,0.01025273,"microtubule-based process"),
c("GO:0009653","anatomical structure morphogenesis",0.7824540774185857,2,0.969406932589886,-0,"anatomical structure morphogenesis"),
c("GO:0009877","nodulation",0.0006951717692939769,1,0.9758601487640157,0.45533059,"anatomical structure morphogenesis"),
c("GO:0055046","microgametogenesis",0.0015100860443036628,1,0.9658151639493533,0.47747485,"anatomical structure morphogenesis"),
c("GO:0009813","flavonoid biosynthetic process",0.015656332622328944,2,0.9411979066705443,0.03204325,"flavonoid biosynthetic process"),
c("GO:0006556","S-adenosylmethionine biosynthetic process",0.05360805457278003,1,0.9428016116093211,0.12361764,"flavonoid biosynthetic process"),
c("GO:0030244","cellulose biosynthetic process",0.03770891075830534,1,0.9329711041803669,0.11447671,"flavonoid biosynthetic process"),
c("GO:0010345","suberin biosynthetic process",0.0015433478514469152,1,0.9237423665884029,0.09640071,"suberin biosynthetic process"),
c("GO:0016310","phosphorylation",7.474579996508841,4,0.9484773228342563,0.08716535,"phosphorylation"),
c("GO:0016567","protein ubiquitination",0.8261966799926771,2,0.9297890717699825,-0,"protein ubiquitination"),
c("GO:0006260","DNA replication",1.444034747678592,1,0.9254773538525394,0.13081923,"protein ubiquitination"),
c("GO:0006412","translation",5.085673767131161,2,0.8731225150522767,0.35814401,"protein ubiquitination"),
c("GO:0006486","protein glycosylation",0.628352124923897,1,0.8976389749369134,0.45949342,"protein ubiquitination"),
c("GO:0006508","proteolysis",5.350747086797883,1,0.9245589314391058,0.46208387,"protein ubiquitination"),
c("GO:0042744","hydrogen peroxide catabolic process",0.1277353179722325,1,0.9346249079263118,0.47893235,"protein ubiquitination"),
c("GO:0046777","protein autophosphorylation",0.0549352006777958,1,0.926755760387361,0.36847714,"protein ubiquitination"),
c("GO:0051603","proteolysis involved in protein catabolic process",0.9636876860000256,1,0.8995835722780419,0.28918282,"protein ubiquitination"),
c("GO:0019375","galactolipid biosynthetic process",0.0015799358393044929,1,0.9301884269657771,0.09615063,"galactolipid biosynthetic process"),
c("GO:0016114","terpenoid biosynthetic process",0.30528351832219996,1,0.905325448217654,0.41553037,"galactolipid biosynthetic process"),
c("GO:0019953","sexual reproduction",0.35416507009992376,2,1,-0,"sexual reproduction"),
c("GO:0032259","methylation",3.1286355155207577,1,0.9807346094481942,0.03567345,"methylation"),
c("GO:0045893","positive regulation of DNA-templated transcription",0.723806859063603,3,0.8384389936430329,-0,"positive regulation of DNA-templated transcription"),
c("GO:0006357","regulation of transcription by RNA polymerase II",1.917576443615649,1,0.880174015900248,0.4674003,"positive regulation of DNA-templated transcription"),
c("GO:0010888","negative regulation of lipid storage",0.0023216741385990235,1,0.9219438450941123,0.41955325,"positive regulation of DNA-templated transcription"),
c("GO:0031047","RNA-mediated gene silencing",0.1393037744966557,1,0.8771826635315785,0.32456012,"positive regulation of DNA-templated transcription"),
c("GO:0040010","positive regulation of growth rate",0.00047897002286283576,1,0.9122786726752258,0.44241796,"positive regulation of DNA-templated transcription"),
c("GO:0042691","positive regulation of crystal cell differentiation",0.00015633049357328667,1,0.914839848767494,0.44286692,"positive regulation of DNA-templated transcription"),
c("GO:0050821","protein stabilization",0.036474897713290676,1,0.9366171809747011,0.17581023,"positive regulation of DNA-templated transcription"),
c("GO:0051726","regulation of cell cycle",0.5059719579017281,1,0.9184944712765444,0.24768465,"positive regulation of DNA-templated transcription"),
c("GO:0060966","regulation of gene silencing by RNA",0.009459657951541007,1,0.9222941851133702,0.26306049,"positive regulation of DNA-templated transcription"),
c("GO:1902290","positive regulation of defense response to oomycetes",0.00011974250571570894,1,0.8890214439396393,0.40937928,"positive regulation of DNA-templated transcription"),
c("GO:2000652","regulation of secondary cell wall biogenesis",0.0007616953835804819,1,0.9461555405419374,0.15137793,"positive regulation of DNA-templated transcription"),
c("GO:0050896","response to stimulus",14.674014998147983,1,1,-0,"response to stimulus"),
c("GO:0051085","chaperone cofactor-dependent protein refolding",0.016697427185912744,1,0.9817854419792842,0.00726942,"chaperone cofactor-dependent protein refolding"),
c("GO:0071555","cell wall organization",0.8657349901438613,7,0.8699903461137262,0,"cell wall organization"),
c("GO:0007005","mitochondrion organization",0.3670107800186479,1,0.9149178070077862,0.45467033,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006325","chromatin organization",1.3384303174082528,6,0.9037388951885091,0.01352937,"chromatin organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9333411730790718,0.39503669,"chromatin organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9279374065611924,0.35552134,"chromatin organization"),
                     c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.8820687665564247,0.33525093,"chromatin organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9579440959878853,0.21045842,"chromatin organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9256312676607006,0.45272646,"chromatin organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,2,0.9115965269804759,0.36480707,"chromatin organization"),
                     c("GO:0051017","actin filament bundle assembly",0.08024701034412013,1,0.9174737417503959,0.46988587,"chromatin organization"),
                     c("GO:0051289","protein homotetramerization",0.01763275115184702,1,0.9215957100689234,0.48332609,"chromatin organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.921292474100556,0.48565697,"chromatin organization"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,2,0.9264926397767697,0.4063955,"chromatin organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9543856736003223,0.20937319,"chromatin organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9499810451975432,0.22611843,"chromatin organization"),
                     c("GO:0006412","translation",4.38869169324396,8,0.7926785473312071,0.01596295,"translation"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,1,0.8677237896585783,0.33478559,"translation"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8520901314386792,0.39091196,"translation"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.934003880398652,0.15163831,"translation"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.9172294311001228,0.49989178,"translation"),
                     c("GO:0006368","transcription elongation by RNA polymerase II",0.1500706345236944,1,0.8393995899047425,0.41096565,"translation"),
                     c("GO:0006397","mRNA processing",1.242261085587905,8,0.8084440392980892,0.48035482,"translation"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.835936321514754,0.47678922,"translation"),
                     c("GO:0006572","tyrosine catabolic process",0.043362173750970734,1,0.8371953264870066,0.47697306,"translation"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.9324316490173363,0.14656365,"translation"),
                     c("GO:0006744","ubiquinone biosynthetic process",0.2844852395091539,1,0.8859151989131273,0.49014462,"translation"),
                     c("GO:0006747","FAD biosynthetic process",0.030434611577160615,1,0.8516574668367464,0.49383913,"translation"),
                     c("GO:0006914","autophagy",0.44134623319182764,1,0.8875365731855603,0.45718665,"translation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9486550285633618,0.23095606,"translation"),
                     c("GO:0009072","aromatic amino acid metabolic process",0.7377977862098182,1,0.8480798733104356,0.32628769,"translation"),
                     c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.8559943480922595,0.44469555,"translation"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8758672995845684,0.36563143,"translation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9357301655135775,0.23899188,"translation"),
                     c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8582844465109717,0.49160318,"translation"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,4,0.8540507583586597,0.11630623,"translation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,7,0.9024381441247201,0.41915608,"translation"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,7,0.8716226622696978,0.4385246,"translation"),
                     c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,1,0.8549385071091905,0.39483568,"translation"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9307485795768821,0.24495457,"translation"),
                     c("GO:0022900","electron transport chain",0.8207017864535318,2,0.9333472366556608,0.105264,"translation"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9078036745139673,0.31208096,"translation"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8912983901142879,0.27892474,"translation"),
                     c("GO:0032544","plastid translation",0.006610433126817685,1,0.8419312670230298,0.46103483,"translation"),
                     c("GO:0033611","oxalate catabolic process",0.006792823899145988,1,0.8951568888908242,0.3986442,"translation"),
                     c("GO:0042364","water-soluble vitamin biosynthetic process",1.0506004254930246,1,0.8683508000478061,0.2442778,"translation"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,1,0.8974335419432055,0.39728749,"translation"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8558057009916047,0.435373,"translation"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,2,0.9445158923116669,0.11952856,"translation"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,1,0.8494579766589408,0.43903701,"translation"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8945424620674443,0.25604839,"translation"),
                     c("GO:0046473","phosphatidic acid metabolic process",0.11459562930583944,1,0.8854244492520411,0.46395125,"translation"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9169646294953643,0.20101185,"translation"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9040388038942404,0.23016052,"translation"),
                     c("GO:0046940","nucleoside monophosphate phosphorylation",0.18958288413443797,1,0.8374780159780952,0.37868979,"translation"),
                     c("GO:0071040","nuclear polyadenylation-dependent antisense transcript catabolic process",0.008039982423444923,1,0.8356113046566447,0.40606321,"translation"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8447728017865375,0.49920408,"translation"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.814945114628871,0.31056631,"translation"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9340960408288801,0.20646638,"translation"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,1,0.8616523601909517,0.41679514,"translation"),
                     c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,1,0.859860881669243,0.37088281,"translation"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.951378202843326,0.11478523,"translation"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8638776036655083,0.33283067,"translation"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9295141548961214,0.11660613,"translation"),
                     c("GO:0006915","apoptotic process",0.3972840732335429,1,0.9854642682152334,0.01170435,"apoptotic process"),
                     c("GO:0007049","cell cycle",2.8073119376140743,2,0.9886804497897039,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9829888947796039,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,2,0.9951870710351908,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,4,1,-0,"metabolic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9533893830047817,0.08207669,"toxin metabolic process"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,4,0.8994341575880751,-0,"response to salt stress"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,3,0.9046004999747203,0.44531989,"response to salt stress"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.9156091405660307,0.38955838,"response to salt stress"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,3,0.8908940499207098,0.49159695,"response to salt stress"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,3,0.8393855848938463,0.31063873,"response to salt stress"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8921612280569129,0.17275947,"response to salt stress"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9304287573707649,0.23963034,"response to salt stress"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8561489198440102,0.38492887,"response to salt stress"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8683695017144614,0.42625684,"response to salt stress"),
                     c("GO:0034059","response to anoxia",0.0003401341429906191,1,0.9239237437324571,0.45207873,"response to salt stress"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,1,0.9003240931273356,0.46099612,"response to salt stress"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8286460161569248,0.46171648,"response to salt stress"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9116752148000451,0.27995591,"response to salt stress"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9784320132590988,0.48813486,"response to salt stress"),
                     c("GO:0051607","defense response to virus",0.17287441054506547,3,0.8905814624060662,0.36620957,"response to salt stress"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9346344614831905,0.43018112,"response to salt stress"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.9262375341932602,0.16037728,"response to salt stress"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9591183995889442,0.06965752,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,5,0.9065471644188652,-0,"flower development"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.9430738021597689,0.4229884,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,4,0.9356400235231461,0.4203896,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9464499525346355,0.48835672,"flower development"),
                     c("GO:0009944","polarity specification of adaxial/abaxial axis",0.0024376280247661035,2,0.9434678225980172,0.38300054,"flower development"),
                     c("GO:0010073","meristem maintenance",0.03536902139069119,1,0.9377404342222467,0.39934449,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9293437414221889,0.29044182,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.941796478243767,0.39539149,"flower development"),
                     c("GO:0048367","shoot system development",0.08373708242002387,1,0.9303744015693652,0.48899493,"flower development"),
                     c("GO:0061137","bud dilation",9.612486649734888E-05,1,0.9538573069389013,0.29194579,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9539044179192605,0.28605225,"flower development"),
                     c("GO:0015031","protein transport",3.093438694074183,12,0.8911444677877459,0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,3,0.9371998749958568,0.42145531,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9565182844571332,0.21356797,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9416038378199687,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,6,0.9317059101071705,0.38566904,"protein transport"),
                     c("GO:0048250","iron import into the mitochondrion",0.008732574410259155,1,0.9473693396780738,0.4861705,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9386774339596206,0.49534489,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,3,0.9277070548179531,0.35175462,"protein transport"),
                     c("GO:0032259","methylation",2.6278542060840238,4,0.9630381097319232,0.06915605,"methylation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,0.9890168330900395,-0,"carbohydrate homeostasis"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,4,0.8817156359903194,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.00031548674132463224,1,0.9386595886188672,0.44967237,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9327569894298242,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,4,0.89929908295257,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9283701550229393,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0009926","auxin polar transport",0.004261535748049133,1,0.9151122521821295,0.4567727,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9409409043106762,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,1,0.9499354952025029,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9441846315281429,0.19548819,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9090817532027623,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0032876","negative regulation of DNA endoreduplication",0.0004264000488215732,1,0.9325263834384296,0.48325207,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.937124668130686,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,1,0.9467761303273207,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9349975528089349,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9340602359715098,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,3,0.8886162915959869,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,2,0.924959055124064,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,1,0.932876854775429,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.935688272889682,0.33352541,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.9507973047224163,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,1,0.9278657196743837,0.47223923,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9370307576793555,0.38691897,"positive regulation of DNA-templated transcription"),
                     c("GO:1901000","regulation of response to salt stress",0.0019915100546117406,1,0.9391213556476719,0.47172255,"positive regulation of DNA-templated transcription"),
                     c("GO:1902183","regulation of shoot apical meristem development",0.0016217990296219374,1,0.9357692310781289,0.38065657,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9278599682057307,0.454995,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,2,0.9374785388067236,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,1,0.9395637552639283,0.3403934,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9403832387750986,0.12485439,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0051301","cell division",1.5693197819947182,1,0.9893374613953471,0.01381155,"cell division"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,4,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,5,0.9999224840330628,6.017E-05,"extracellular region"),
                     c("GO:0005739","mitochondrion",4.856981674016684,20,0.7429688158951828,0,"mitochondrion"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,4,0.7849669552020982,0.15656526,"mitochondrion"),
                     c("GO:0000785","chromatin",1.087181392930179,2,0.6267578627158689,0.45756389,"mitochondrion"),
                     c("GO:0005634","nucleus",16.5161752456724,75,0.7108148502548343,0.35145602,"mitochondrion"),
                     c("GO:0005730","nucleolus",1.2114693153196516,8,0.6525206369991348,0.22801288,"mitochondrion"),
                     c("GO:0005737","cytoplasm",43.076660800468886,64,0.8844138649186254,0.12447513,"mitochondrion"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,5,0.625087928442393,0.22750481,"mitochondrion"),
                     c("GO:0005773","vacuole",1.35235298008496,4,0.7746528732709345,0.23143591,"mitochondrion"),
                     c("GO:0005777","peroxisome",0.7476308639759881,3,0.7769793003190435,0.21411812,"mitochondrion"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,13,0.692503713783265,0.26090848,"mitochondrion"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,2,0.7819153382655257,0.38610723,"mitochondrion"),
                     c("GO:0005829","cytosol",14.048080989011824,27,0.8532490738086973,0.19224382,"mitochondrion"),
                     c("GO:0009507","chloroplast",0.6990388085931706,14,0.7294336041196315,0.21231674,"mitochondrion"),
                     c("GO:0009536","plastid",0.7624897559369845,10,0.7863656568074683,0.21465144,"mitochondrion"),
                     c("GO:0009579","thylakoid",0.26089530737804617,1,0.7779106746861505,0.39454155,"mitochondrion"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6909463303728424,0.40758865,"mitochondrion"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,4,0.6894561784958844,0.25147478,"mitochondrion"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.7825608346576599,0.2461316,"mitochondrion"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,1,0.7832725491827927,0.4052793,"mitochondrion"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7358403967214099,0.486126,"mitochondrion"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9030105125468458,0.10673752,"mitochondrion"),
                     c("GO:0005886","plasma membrane",17.177321395000487,23,0.959322598108123,0.00015294,"plasma membrane"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,5,0.9795750356957036,3.812E-05,"plasmodesma"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999360816469345,4.693E-05,"cell surface"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.999916561578383,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,69,0.9998855005769441,0.00010118,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999420397195233,4.148E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9999461291247308,3.787E-05,"apoplast"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9957105861594396,2.697E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.9957188065510704,0.27222822,"pollen tube"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,10,0.7917625287215776,-0,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7433502552256501,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,3,0.7515619719814391,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,5,0.5914476605766589,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8088404021989811,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7931079870659599,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005942","phosphatidylinositol 3-kinase complex",0.07251039886384637,1,0.7616223441935326,0.46964777,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.772242956869616,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0030131","clathrin adaptor complex",0.06838817448871372,1,0.7570390531220923,0.40377927,"ribonucleoprotein complex"),
                     c("GO:0033179","proton-transporting V-type ATPase, V0 domain",0.05699304429922707,1,0.7650237145721609,0.32669992,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7278739987734285,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7991344382389276,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7445202547292793,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.778691591901441,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7012162813786493,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.7152465682937388,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,5,0.7773723415349952,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.6989975611723888,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7945554103286121,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,4,1,-0,"molecular_function"),
                     c("GO:0003729","mRNA binding",1.136762432364793,22,0.919189989329702,0,"mRNA binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,32,0.9055575952249222,0.48160605,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9396185885442078,0.42104634,"mRNA binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,3,0.9162007849864184,0.31664803,"mRNA binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,27,0.9332895993562562,0.20275684,"mRNA binding"),
                     c("GO:0003713","transcription coactivator activity",0.24395086691346182,2,0.9753831979375028,0.29461398,"mRNA binding"),
                     c("GO:0003723","RNA binding",6.099813894661886,31,0.9194448330560689,0.40176997,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9309513931976233,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9302157238216323,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,6,0.9123267158409292,0.24940265,"mRNA binding"),
                     c("GO:0003919","FMN adenylyltransferase activity",0.02766322823642363,1,0.9008570455896044,0.2955927,"mRNA binding"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8946958237802933,0.3505363,"mRNA binding"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8978500615159202,0.33108637,"mRNA binding"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8971382099233871,0.3394079,"mRNA binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,2,0.9460972935657948,0.3693843,"mRNA binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,7,0.9200294102727018,0.13371594,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9386791003004495,0.428381,"mRNA binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,5,0.9345942704456097,0.19635374,"mRNA binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9360695411920246,0.29697106,"mRNA binding"),
                     c("GO:0016780","phosphotransferase activity, for other substituted phosphate groups",0.32920461222629094,1,0.8811699811123633,0.35102513,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9526208862981007,0.13769902,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9409402074530071,0.41077037,"mRNA binding"),
                     c("GO:0032934","sterol binding",0.11902838493358808,1,0.9520409820555689,0.10144903,"mRNA binding"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8792496643222873,0.27195941,"mRNA binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9713245237085395,0.2693746,"mRNA binding"),
                     c("GO:0044183","protein folding chaperone",0.39473153762799984,1,0.9817725532787365,0.32443613,"mRNA binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.927582774792719,0.33990845,"mRNA binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,3,0.9412561574070847,0.16455763,"mRNA binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,1,0.9407046839109325,0.29755636,"mRNA binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,1,0.9346852271099148,0.30272555,"mRNA binding"),
                     c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9550440586924542,0.25647876,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,3,0.9940744567948091,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,10,1,-0,"catalytic activity"),
                     c("GO:0005515","protein binding",8.610051728351934,38,0.9644651205302437,0.0711246,"protein binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,3,0.9716392748624806,0.05041733,"lipid binding"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9621416322696915,0.0291488,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0004165","delta(3)-delta(2)-enoyl-CoA isomerase activity",0.021527424426388823,1,0.958026799223716,0.44447646,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0004619","phosphoglycerate mutase activity",0.0367083831844294,1,0.9536718610095938,0.40080388,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,1,0.9377674971012107,-0,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9003538444090412,0.39044024,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.930714653561947,0.47124298,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9532029424934207,0.04318065,"CoA-ligase activity"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9179776234311781,0.04634246,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9395229810745828,0.2770188,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9191237651479918,0.44913707,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9258026679233834,0.39379271,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9328564423684245,0.28092429,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,28,0.9473893227859949,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,12,0.9514633183740351,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,27,0.9468379296218662,0.12712323,"transferase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,4,0.8829813457657951,0.05667901,"glycosyltransferase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9212959877218342,0.20205602,"glycosyltransferase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9074792207022762,0.24793639,"glycosyltransferase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,4,0.8474485826623853,0.33973871,"glycosyltransferase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,1,0.8922989155723076,0.2917808,"glycosyltransferase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8828842450274552,0.28934569,"glycosyltransferase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8778123282706567,0.36085078,"glycosyltransferase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.9005106240285965,0.49223587,"glycosyltransferase activity"),
                     c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.915281639990702,0.17480439,"glycosyltransferase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.879379968175537,0.45431535,"glycosyltransferase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8812924459963372,0.39355689,"glycosyltransferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,2,0.9575628370270132,0.05997119,"lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9643519416613686,0.03194325,"strictosidine synthase activity"),
                     c("GO:0004300","enoyl-CoA hydratase activity",0.09504176049804544,1,0.9562995134925896,0.39952982,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,3,0.9593802156632314,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,3,0.9580668665713258,0.0589874,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,13,0.859843567619985,-0,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.899100324036126,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,3,0.8471534642533737,0.41870234,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9334241023351099,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004430","1-phosphatidylinositol 4-kinase activity",0.022294677089298852,1,0.8963787319871025,0.37585892,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8988180044165047,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0004714","transmembrane receptor protein tyrosine kinase activity",0.12685170110337585,1,0.8401190602818799,0.47011379,"ATP hydrolysis activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,2,0.8206292145514233,0.48775459,"ATP hydrolysis activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,9,0.8243494089336002,0.3656957,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8882722935976386,0.39048386,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.9004105219525277,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,2,0.8774789740202524,0.20698894,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9233918951177618,0.23079392,"ATP hydrolysis activity"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,3,0.8445907312449777,0.37283533,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9221043593761462,0.3204136,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9372040574477077,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030295","protein kinase activator activity",0.14413927844455088,1,0.9873128641687623,-0,"protein kinase activator activity"),
                     c("GO:0043022","ribosome binding",0.5329368041478477,2,0.9513019234838447,0.04582217,"ribosome binding"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9699060444846261,0.054576,"protein-containing complex binding"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9666702280639078,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9757673014373125,0.49249984,"ER retention sequence binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.96753560753899,0.03635569,"ubiquinone binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9183948627609377,0.0463004,"unfolded protein binding"),
                     c("GO:0003779","actin binding",0.7421262469463454,2,0.9105620046594465,0.44654025,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9320630371182815,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,2,0.9155094327542793,0.42444044,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,1,0.9280562478047915,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,2,0.930543892534813,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9263198924398789,0.38789043,"unfolded protein binding"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9309390668773161,0.36031388,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9495734550393145,0.25391209,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9207010176079609,0.42217449,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.926138112437574,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9247713362959409,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9317478271669182,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9134297138568106,0.47948026,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9478734934103231,0.26331899,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.9390444789707422,0.31312228,"unfolded protein binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006952","defense response",1.1604144588756624,21,0.9014415233402141,-0,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8752042541957586,0.44951817,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,2,0.9150452051149945,0.45786068,"defense response"),
                     c("GO:0009625","response to insect",0.000539778096485113,1,0.9291554169291139,0.49043501,"defense response"),
                     c("GO:0009646","response to absence of light",0.0015108857221249963,2,0.9120105516585164,0.48888107,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,7,0.8891454207552074,0.4259419,"defense response"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,7,0.8436645479161755,0.24614308,"defense response"),
                     c("GO:0009751","response to salicylic acid",0.010802956150202055,2,0.9066274414130375,0.48305816,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8945114961461071,0.26818052,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,1,0.9216148131432202,0.4902216,"defense response"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9286548257423337,0.41832444,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9225360330506323,0.29372182,"defense response"),
                     c("GO:0010343","singlet oxygen-mediated programmed cell death",0.0019717921332789512,3,0.8891696093474535,0.38975058,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8541037605927543,0.39975389,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,2,0.8704782250522954,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,4,0.8995467814880478,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,4,0.8347303631233572,0.42828273,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9093418202069201,0.46171648,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,5,0.8888638884058782,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9800865143688429,0.48259105,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9232817742053951,0.44100881,"defense response"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.8100738023559206,0.384161,"defense response"),
                     c("GO:0071000","response to magnetism",7.887168533115806E-05,1,0.9303898364080667,0.31824892,"defense response"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9304188450813394,0.14771546,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9219237996016767,0.40537115,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9133752749105614,0.26816179,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.9277672192478896,0.32241615,"defense response"),
                     c("GO:0007017","microtubule-based process",1.4522840599239462,1,0.9899491638152037,0.01367267,"microtubule-based process"),
                     c("GO:0007018","microtubule-based movement",0.6202173565622278,2,0.9888196900259673,0.01231347,"microtubule-based movement"),
                     c("GO:0007049","cell cycle",2.8073119376140743,4,0.9892705821896848,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9790446456918097,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,6,0.9959533166598166,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,7,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,9,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9919661296944029,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9910094743951549,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,6,0.9528457169537842,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,4,0.9502312409889643,0.21768792,"biosynthetic process"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,3,0.9491604339430597,0.27529491,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9229285056035131,0.16846091,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9387905680641739,0.17227516,"biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9595194773116903,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,10,0.8825629182428915,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,4,0.9126127560427858,0.4203896,"flower development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9137796479654485,0.3897097,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9372901018299922,0.48835672,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,2,0.9342080187589882,0.3161767,"flower development"),
                     c("GO:0010077","maintenance of inflorescence meristem identity",0.0001355607091629279,2,0.9279314546127306,0.29658706,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.8517435466529647,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9240734370776816,0.29044182,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.9235019774120778,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.917819379933899,0.42767441,"flower development"),
                     c("GO:0048262","determination of dorsal/ventral asymmetry",0.0029848003417510126,2,0.9201281255410524,0.38718574,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.937038912692902,0.28605225,"flower development"),
                     c("GO:0090626","plant epidermis morphogenesis",0.016804598455869863,1,0.9153740668740785,0.46768707,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,1,0.9278642103639095,0.48720779,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.9937958592974828,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9585695496235732,0.04613723,"chlororespiration"),
                     c("GO:0010960","magnesium ion homeostasis",0.024344238625495253,2,0.9731856934415024,-0,"magnesium ion homeostasis"),
                     c("GO:0010268","brassinosteroid homeostasis",0.011392029050019141,1,0.9758982367300477,0.46517239,"magnesium ion homeostasis"),
                     c("GO:0060586","multicellular organismal-level iron ion homeostasis",0.00353443739890252,1,0.9233250804375102,0.48479779,"magnesium ion homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9798449387920019,0.37586433,"magnesium ion homeostasis"),
                     c("GO:0015031","protein transport",3.093438694074183,20,0.9020310464464982,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,10,0.9408560440643576,0.42145531,"protein transport"),
                     c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.889803447393566,0.42491424,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9587691104465781,0.21356797,"protein transport"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9479669069238335,0.48373949,"protein transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9606534356467149,0.39622712,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9524081710689201,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,8,0.9355613033604673,0.38566904,"protein transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9617817562793647,0.30124551,"protein transport"),
                     c("GO:0048250","iron import into the mitochondrion",0.008732574410259155,1,0.9460649327887085,0.4861705,"protein transport"),
                     c("GO:0060918","auxin transport",0.015500750907739157,3,0.9119859198185138,0.22426661,"protein transport"),
                     c("GO:0120010","intermembrane phospholipid transfer",0.007554428610624982,1,0.9010444626280496,0.42470759,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.942342927030281,0.49534489,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,4,0.9248036519119001,0.35175462,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,7,0.9523107492733556,0.04477571,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,21,0.8867418503322062,0,"protein ubiquitination"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,2,0.8857069838859444,0.17872097,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8719000511074585,0.40302693,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,8,0.938475339077507,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9378623204198039,0.43492456,"protein ubiquitination"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,1,0.9162879551515956,0.41634329,"protein ubiquitination"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,3,0.9076858446580361,0.40337462,"protein ubiquitination"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",18.82728500884891,2,0.8598648721274008,0.31102513,"protein ubiquitination"),
                     c("GO:0006221","pyrimidine nucleotide biosynthetic process",0.5101395959817636,2,0.8177561276725704,0.39188311,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,9,0.8362351537008291,0.48035482,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,15,0.8179568353632322,0.4385246,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,14,0.8564446603318003,0.47678922,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,5,0.8552032043842711,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,3,0.9384366369838889,0.15178462,"protein ubiquitination"),
                     c("GO:0006541","glutamine metabolic process",0.5478057552077248,3,0.8657742027539862,0.48867994,"protein ubiquitination"),
                     c("GO:0006572","tyrosine catabolic process",0.043362173750970734,2,0.8529629947027482,0.46616537,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,19,0.9370739768470074,0.12862028,"protein ubiquitination"),
                     c("GO:0006631","fatty acid metabolic process",1.936945636803579,6,0.8375113762255482,0.1033813,"protein ubiquitination"),
                     c("GO:0006744","ubiquinone biosynthetic process",0.2844852395091539,2,0.8806614315335874,0.31493359,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,2,0.9105803245525947,0.45718665,"protein ubiquitination"),
                     c("GO:0009695","jasmonic acid biosynthetic process",0.005493905831348478,1,0.8734302267163684,0.46296398,"protein ubiquitination"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9271848026810958,0.37903497,"protein ubiquitination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8926605170285802,0.36563143,"protein ubiquitination"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9433089891320287,0.32178408,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,4,0.8758501588520232,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,37,0.9022626348237903,0.46958224,"protein ubiquitination"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9400872737254753,0.33268781,"protein ubiquitination"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9190334711922511,0.27916289,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,4,0.887881288036158,0.14129671,"protein ubiquitination"),
                     c("GO:0031407","oxylipin metabolic process",0.012508556345488349,1,0.92142645830669,0.34584249,"protein ubiquitination"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,2,0.8988273425569252,0.34563628,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,2,0.8658033640727484,0.28833637,"protein ubiquitination"),
                     c("GO:0033358","UDP-L-arabinose biosynthetic process",2.21826614993882E-05,1,0.9044380961315857,0.46513452,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9461418414882614,0.45767193,"protein ubiquitination"),
                     c("GO:0042245","RNA repair",0.022010129687726296,1,0.9008423249166922,0.26040271,"protein ubiquitination"),
                     c("GO:0042726","flavin-containing compound metabolic process",0.22289631222618586,1,0.9038603215380578,0.18149234,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,2,0.906545601370787,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,2,0.915799458724189,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.8857711662541896,0.26614879,"protein ubiquitination"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.9013344167394665,0.16428511,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,6,0.9116976020514135,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,2,0.8999673429304007,0.41108531,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9041969061294882,0.44155167,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.86697030996335,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8960007461728071,0.43527038,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9380390776392794,0.45105203,"protein ubiquitination"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.8377417472723272,0.33478559,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,3,0.8855880125377597,0.32141303,"protein ubiquitination"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8702808049653483,0.49716341,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9918814181577328,0.01031083,"establishment of cell polarity"),
                     c("GO:0032259","methylation",2.6278542060840238,7,0.9655103582987943,0.05843136,"methylation"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9757454380382784,0.03897815,"collagen metabolic process"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,9,0.8991566082954726,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9428039972514977,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,4,0.9135546939424173,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9311564383206755,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,2,0.9474527886383333,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,3,0.9547254021333066,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9529935368317014,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9219415586561889,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0010928","regulation of auxin mediated signaling pathway",0.004086539196220627,1,0.9378189029717181,0.40226117,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9436478537127534,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9326499349548192,0.31416992,"positive regulation of DNA-templated transcription"),
                     c("GO:0042548","regulation of photosynthesis, light reaction",0.0070861279789712316,1,0.9443667897618462,0.23746411,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,3,0.9521274585792716,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9418312780824202,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9412598393586022,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,6,0.9022931592291114,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,2,0.9345363630100864,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.947303756171555,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9400267042149186,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9380809659054875,0.33352541,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,4,0.953751044933043,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9373416751386906,0.46789196,"positive regulation of DNA-templated transcription"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,3,0.9237448337854424,0.34078171,"positive regulation of DNA-templated transcription"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9343407088175688,0.39696279,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9435583442134363,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1901529","positive regulation of anion channel activity",1.2323700832993446E-05,1,0.9516654506605777,0.45462725,"positive regulation of DNA-templated transcription"),
                     c("GO:1902183","regulation of shoot apical meristem development",0.0016217990296219374,1,0.9439574409215425,0.38065657,"positive regulation of DNA-templated transcription"),
                     c("GO:1902448","positive regulation of shade avoidance",4.190058283217772E-05,1,0.9403999884252786,0.30806662,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.929191830567516,0.45547566,"positive regulation of DNA-templated transcription"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9434116934817975,0.40819627,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,3,0.9443940753676973,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,3,0.9324739861219885,0.36922831,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,3,0.9394989877255865,0.13634361,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.93548146543556,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000032","regulation of secondary shoot formation",0.007307954593965113,1,0.9380226816663005,0.49488849,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9414798963672286,0.42989004,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9379112523451213,0.46345207,"positive regulation of DNA-templated transcription"),
                     c("GO:2000185","regulation of phosphate transmembrane transport",0.00017253181166190822,1,0.959784708305648,0.44318481,"positive regulation of DNA-templated transcription"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9431411244988532,0.15519254,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,2,0.9435362881050269,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.9371495805560641,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,1,0.9618554096420213,0.12578285,"positive regulation of DNA-templated transcription"),
                     c("GO:2000762","regulation of phenylpropanoid metabolic process",0.005555524335513445,1,0.9473593492243972,0.23331889,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,3,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,7,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,2,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,2,0.9898740850883301,0.01381155,"cell division"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,10,0.9223546881051858,0.01286368,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,8,0.9160040513190786,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9355201785624562,0.4007549,"cell wall organization"),
                     c("GO:0007005","mitochondrion organization",0.8276597479438399,1,0.9257385419234894,0.4723023,"cell wall organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9287212961507355,0.42383478,"cell wall organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,2,0.9314090845663998,0.39090469,"cell wall organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9336721204720532,0.41402298,"cell wall organization"),
                     c("GO:0009658","chloroplast organization",0.0906284959258338,3,0.9309591048054646,0.31305597,"cell wall organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9604528369759061,0.2057032,"cell wall organization"),
                     c("GO:0009920","cell plate formation involved in plant-type cell wall biogenesis",0.00015527863049571742,1,0.9296977845087128,0.45904668,"cell wall organization"),
                     c("GO:0010207","photosystem II assembly",0.0070861279789712316,1,0.8968887601852811,0.45155549,"cell wall organization"),
                     c("GO:0010275","NAD(P)H dehydrogenase complex assembly",0.0022404488114382086,1,0.9392468742119799,0.41693818,"cell wall organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9350898198040947,0.45272646,"cell wall organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,2,0.9225334499289222,0.35075221,"cell wall organization"),
                     c("GO:0034462","small-subunit processome assembly",0.005688620304509775,1,0.9233327415361822,0.44451454,"cell wall organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9313780013526536,0.48565697,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9344573376323081,0.4524177,"cell wall organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,2,0.9555833813546034,0.20466634,"cell wall organization"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9541501166999655,0.03891167,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9346695357018706,0.08958079,"intrachromosomal DNA recombination"),
                     c("GO:0006270","DNA replication initiation",0.18677061560434888,1,0.8915880117660419,0.25665205,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,7,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,20,0.9999272053216985,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9770013299882434,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,120,0.897241834188223,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,54,0.8676759780409261,0.17160779,"cytoplasm"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,2,0.9871015694677739,2.997E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,65,0.9402851631532837,0.3208057,"plant-type cell wall"),
                     c("GO:0009925","basal plasma membrane",0.13578741268972233,1,0.9499567674091494,0.18839421,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,14,0.9670987319565381,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,66,0.7267748257471457,0,"chloroplast"),
                     c("GO:0000138","Golgi trans cisterna",0.007541260384887046,1,0.7847126229203457,0.37607379,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,8,0.8097345480148721,0.13280054,"chloroplast"),
                     c("GO:0000785","chromatin",1.087181392930179,2,0.6310218292716503,0.45756389,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,147,0.7351725457812942,0.35145602,"chloroplast"),
                     c("GO:0005730","nucleolus",1.2114693153196516,11,0.6593840512562242,0.22801288,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,32,0.7631637114037019,0.21465144,"chloroplast"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,5,0.6487140970015214,0.22750481,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,8,0.7910318848950383,0.23143591,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,7,0.7950367744269272,0.17354881,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,30,0.7058158314534835,0.2503098,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,2,0.7805719258777527,0.38610723,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,5,0.7346995442253187,0.44106058,"chloroplast"),
                     c("GO:0005905","clathrin-coated pit",0.09808111076528711,1,0.8843034674173224,0.46642691,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,62,0.801418427855384,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,1,0.8262819239280932,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,9,0.7767429599132744,0.39454155,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6864100932783908,0.39861008,"chloroplast"),
                     c("GO:0031977","thylakoid lumen",0.009951481990600534,2,0.7872051218152879,0.29998046,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,3,0.7976194431299694,0.2461316,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.7927351528465045,0.28709921,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7823312993438182,0.31004802,"chloroplast"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.737053290975465,0.486126,"chloroplast"),
                     c("GO:0043661","peribacteroid membrane",1.987811633578134E-05,1,0.8037970927187932,0.30368923,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999393415696792,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8842947842593032,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,5,0.9999220057920939,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,189,0.9998955392513816,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0034045","phagophore assembly site membrane",0.08896450966078939,1,0.9106158180456937,0.07998094,"phagophore assembly site membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.9999304352273388,5.465E-05,"cell projection"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999447452210687,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,3,0.9951486090250801,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,2,0.9107301932146584,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9926089325091942,2.369E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.993655243626757,0.27222822,"pollen tube"),
                     c("GO:0098552","side of membrane",0.7241349304670945,1,0.9680390950268467,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,12,0.7977674519751797,-0,"ribonucleoprotein complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,1,0.7854921537208012,0.25406433,"ribonucleoprotein complex"),
                     c("GO:0000159","protein phosphatase type 2A complex",0.06136622989309897,1,0.812960587135565,0.32835456,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7621299398938972,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,3,0.7627026089730917,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7270950971474648,0.4447901,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,5,0.6140553343058295,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8160882818678357,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.8011306332882022,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.7764207224685161,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.7822423676723407,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0030131","clathrin adaptor complex",0.06838817448871372,1,0.7599273560528329,0.47444895,"ribonucleoprotein complex"),
                     c("GO:0033179","proton-transporting V-type ATPase, V0 domain",0.05699304429922707,1,0.7570345167166008,0.46766199,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7521800807563347,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.8054714902444956,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,2,0.7711540233433036,0.30930488,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7700242866731067,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.787512812863161,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7166713683204003,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.7298382771511723,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,6,0.7869422711005678,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.7145901294300672,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.8005408034760669,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,17,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,4,0.9649126240562929,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,28,0.9429737859643986,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9573796320794878,0.42104634,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,46,0.9368707036252987,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9625789804256175,0.19519762,"mRNA binding"),
                     c("GO:0003714","transcription corepressor activity",0.25663049329068593,2,0.9768100542757125,0.29590577,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9512767475415395,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9511270322037229,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,8,0.9313035085598319,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9567183725746465,0.428381,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,6,0.9647944226121699,0.13060477,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9583097999647859,0.41077037,"mRNA binding"),
                     c("GO:0035198","miRNA binding",0.029599099839661934,1,0.9577632770366261,0.41680219,"mRNA binding"),
                     c("GO:0042134","rRNA primary transcript binding",0.01629192287398833,1,0.9594766474822041,0.39794793,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9556931252544846,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,12,0.9395673469359026,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,3,0.996953857194905,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,22,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,29,0.7966955598053842,0,"protein kinase activity"),
                     c("GO:0000234","phosphoethanolamine N-methyltransferase activity",0.0010222643861315665,1,0.9101961182379926,0.45971198,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,3,0.9041483348680777,0.34856872,"protein kinase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.917911386475061,0.20796201,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8943318146483252,0.45866285,"protein kinase activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,4,0.8580116616013308,0.41115081,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8945077899471767,0.48555173,"protein kinase activity"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.89424154920464,0.46089842,"protein kinase activity"),
                     c("GO:0004619","phosphoglycerate mutase activity",0.0367083831844294,1,0.9534616077201333,0.4684783,"protein kinase activity"),
                     c("GO:0004659","prenyltransferase activity",0.25118876486646274,1,0.8953678773240584,0.26852283,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,5,0.8356403557433328,0.41376722,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,7,0.8225227726558936,0.47753726,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,7,0.860152910741148,0.35677507,"protein kinase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,2,0.86787225776421,0.30425857,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8606159888939169,0.35123485,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,4,0.8785653070669813,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8974818261195343,0.40329157,"protein kinase activity"),
                     c("GO:0008798","beta-aspartyl-peptidase activity",0.023041972313000234,1,0.880325975233364,0.42550609,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.8540627345701525,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,7,0.8760345470755091,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,19,0.8809090870526365,0.35580262,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.908880760441058,0.49223587,"protein kinase activity"),
                     c("GO:0033855","nicotianamine aminotransferase activity",6.6524797362141E-06,1,0.9354078780944137,0.45711636,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8810299774872884,0.49253483,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.8991384315701102,0.49478365,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9171381390504977,0.42478125,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,2,0.8896139113429542,0.48733868,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9251136312087593,0.26985858,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9069933912781002,0.47557181,"protein kinase activity"),
                     c("GO:0052667","phosphomethylethanolamine N-methyltransferase activity",0.00011530964876104439,1,0.9190552469379513,0.40189746,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9277004896234339,0.49777323,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.9042916666933052,0.35248487,"protein kinase activity"),
                     c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.9072941853896729,0.31861303,"protein kinase activity"),
                     c("GO:0106262","1-acylglycerophosphoethanolamine O-acyltransferase activity",0.0003703213719825849,1,0.9122456680807022,0.44573481,"protein kinase activity"),
                     c("GO:0106263","1-acylglycerophosphoserine O-acyltransferase activity",0.0004324111828539164,1,0.9115852939314706,0.15964368,"protein kinase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.892179336596866,0.45431535,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.891768907063418,0.39355689,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,3,0.9777415698833525,-0,"GTPase activator activity"),
                     c("GO:0005515","protein binding",8.610051728351934,88,0.9721879680375627,0.07766851,"protein binding"),
                     c("GO:0097159","organic cyclic compound binding",39.22969739688024,1,0.9650038045512186,0.13134756,"protein binding"),
                     c("GO:0005543","phospholipid binding",0.6999096105403309,3,0.9638905596726693,0.04714527,"phospholipid binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,15,0.9479648068525763,0.05738816,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,88,0.9273531578387632,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,4,0.95568184563194,0.3830934,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,2,0.9688147268833677,0.21260669,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,34,0.9506227862818133,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,6,0.9569447823796046,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,13,0.9382587711667798,0.19635374,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9480095235583085,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9499164392995864,0.29697106,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,3,0.9439997884366288,0.31800282,"zinc ion binding"),
                     c("GO:0030955","potassium ion binding",0.05969270067304912,1,0.9674813988533001,0.26193125,"zinc ion binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9785687783374003,0.2693746,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,2,0.9484806002806669,0.2885161,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9738119879027932,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,5,0.9428110545651036,0.33990845,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,3,0.9486566541871637,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,6,0.9538370946043542,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9546900860383924,0.28523818,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,3,0.9529625669172973,0.29755636,"zinc ion binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,3,0.9483293052956454,0.30272555,"zinc ion binding"),
                     c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9653524532195201,0.25647876,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,4,0.9777629687338549,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,3,0.9931998046540986,-0,"photoreceptor activity"),
                     c("GO:0009882","blue light photoreceptor activity",0.009566265860675875,1,0.9807881196363898,0.45612612,"photoreceptor activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,2,0.9872072771190255,0.02769872,"auxin binding"),
                     c("GO:0015293","symporter activity",0.6379683717164414,3,0.923821166373216,-0,"symporter activity"),
                     c("GO:0005254","chloride channel activity",0.1887840699542837,2,0.9242613427716713,0.3683835,"symporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9161173207841773,0.46407837,"symporter activity"),
                     c("GO:0015165","pyrimidine nucleotide-sugar transmembrane transporter activity",0.046425438585792796,1,0.9438418196992429,0.32778933,"symporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9445611913696884,0.47124298,"symporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9386913249698924,0.3537519,"symporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.96521610515402,0.20748322,"symporter activity"),
                     c("GO:0042910","xenobiotic transmembrane transporter activity",0.2653474592383718,1,0.9394495999257795,0.3797985,"symporter activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9246894291754755,0.43531693,"symporter activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,5,0.9032465887001682,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,4,0.9045234927405373,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,1,0.9010181384485876,0.40806599,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0010242","oxygen evolving activity",0.0031710153409287207,1,0.9368825710449135,0.26712027,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9251361833992277,0.46927032,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,4,0.8900399704432891,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9055060857322698,0.41247223,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.906835546147857,0.44173969,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016630","protochlorophyllide reductase activity",0.0022374506846133423,1,0.9324827698857283,0.26085905,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,1,0.9156698087779124,0.3863248,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016639","oxidoreductase activity, acting on the CH-NH2 group of donors, NAD or NADP as acceptor",0.09697097962154752,1,0.9173384851662528,0.34937139,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9178622094095388,0.48602252,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0019139","cytokinin dehydrogenase activity",0.0057233500663895305,1,0.9347256405517566,0.27843728,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0019797","procollagen-proline 3-dioxygenase activity",0.004459378916508851,1,0.8891486089271369,0.3991339,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0048529","magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase activity",0.002614424536332141,1,0.9342795792957843,0.26361814,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9090203772696658,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,95,0.9463340990000323,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,33,0.9501098369352733,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,61,0.9458248194786081,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,12,0.9558474527626463,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,3,0.9312451570699782,0.04651897,"carboxy-lyase activity"),
                     c("GO:0009978","allene oxide synthase activity",9.09172230615927E-05,1,0.9586286303151453,0.33685723,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9507275672268762,0.45277577,"carboxy-lyase activity"),
                     c("GO:0018812","3-hydroxyacyl-CoA dehydratase activity",0.012180690397008016,1,0.9470672300218994,0.4947816,"carboxy-lyase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9445868285842907,0.48814337,"carboxy-lyase activity"),
                     c("GO:0106099","2-keto-3-deoxy-L-rhamnonate aldolase activity",0.00019513940559561356,1,0.9569514944109843,0.48677249,"carboxy-lyase activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9644978415543609,0.28175041,"carboxy-lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,10,0.9575820482613279,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,11,0.9563272943543439,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,20,0.86154668418708,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9066338370600061,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9303409957872354,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,4,0.872994834358962,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8863418508493948,0.38812549,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,4,0.8980043197460424,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.898268516476489,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.8995965306598724,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,1,0.9213380368262046,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,2,0.9236516297613782,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047668","amygdalin beta-glucosidase activity",2.4392425699451695E-05,1,0.9297547974380621,0.44127145,"ATP hydrolysis activity"),
                     c("GO:0047701","beta-L-arabinosidase activity",4.878485139890339E-05,1,0.9311589858280174,0.41027052,"ATP hydrolysis activity"),
                     c("GO:0047884","FAD diphosphatase activity",0.00013970207446049608,1,0.9378284698963498,0.31970647,"ATP hydrolysis activity"),
                     c("GO:0050224","prunasin beta-glucosidase activity",1.55224527178329E-05,1,0.9310740794270776,0.43456671,"ATP hydrolysis activity"),
                     c("GO:0080079","cellobiose glucosidase activity",7.982975683456919E-05,1,0.9260382380344921,0.45988357,"ATP hydrolysis activity"),
                     c("GO:0080083","beta-gentiobiose beta-glucosidase activity",7.76122635891645E-05,1,0.9261311566013791,0.42073911,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9289246491016637,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9258252545029961,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9432055679184503,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,4,0.9783323314588592,0.04901613,"carbohydrate binding"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9477681450997801,0.03310522,"DNA demethylase activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9764184523015997,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9955908316513356,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8901095611423163,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9775829143049821,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9839743343245027,0.49249984,"ER retention sequence binding"),
                     c("GO:0050373","UDP-arabinose 4-epimerase activity",8.20472500799739E-05,1,0.9657380793298711,0.02256811,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9591450397867443,0.44049754,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,9,0.9291812586835051,0.0463004,"unfolded protein binding"),
                     c("GO:0000149","SNARE binding",0.26943873427614345,1,0.9335923777889343,0.40559978,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,4,0.93369926377063,0.40485813,"unfolded protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,5,0.9263343769046012,0.43412799,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9411896537168161,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,3,0.9240850272295623,0.42444044,"unfolded protein binding"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.9351809171549877,0.39461121,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,2,0.9376807483820682,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9398601683828287,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,2,0.9361576356371167,0.38789043,"unfolded protein binding"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9402061067010586,0.36031388,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9564547862130692,0.25391209,"unfolded protein binding"),
                     c("GO:0032050","clathrin heavy chain binding",0.02621520514717436,1,0.9439845079074578,0.33497261,"unfolded protein binding"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,1,0.9343136186262746,0.4006017,"unfolded protein binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9439935311809833,0.31955157,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,3,0.931216122123806,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,6,0.9301481221369102,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,8,0.9359980831746438,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,3,0.934797809663528,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,5,0.9409138808360964,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,3,0.9247830160760373,0.4679269,"unfolded protein binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.9329467250005615,0.41008684,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.953210853019769,0.26331899,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,2,0.9449917264125887,0.48770765,"unfolded protein binding"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,1,0.9408432064229358,0.35601456,"unfolded protein binding"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.963076983842748,0.02238092,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0003972","RNA ligase (ATP) activity",0.012360307349885798,1,0.9363054715674739,0.36034906,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9546792169747752,0.28284902,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,1,0.9633250098724788,0.2351635,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:2001070","starch binding",0.024815966909324,1,0.9816736548671297,0.03482451,"starch binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000911","cytokinesis by cell plate formation",0.008712856488926366,1,0.9943739120321408,0.00821801,"cytokinesis by cell plate formation"),
                     c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,1,1,-0,"intracellular iron ion homeostasis"),
                     c("GO:0007155","cell adhesion",1.1091577223710762,1,0.9917207501957832,0.01321071,"cell adhesion"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,1,1,-0,"metabolic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,2,0.8123060199978509,-0,"steroid metabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,1,0.9366026159428562,0.15163831,"steroid metabolic process"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,1,0.8058495451850985,0.48590569,"steroid metabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9352545421856103,0.10901633,"steroid metabolic process"),
                     c("GO:0008203","cholesterol metabolic process",0.2200421431132646,2,0.790056955917325,0.45443823,"steroid metabolic process"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.7666090970461414,0.49264795,"steroid metabolic process"),
                     c("GO:0042759","long-chain fatty acid biosynthetic process",0.04469313344093403,1,0.7993405984293194,0.3981042,"steroid metabolic process"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.7771081847394969,0.42785472,"steroid metabolic process"),
                     c("GO:0009793","embryo development ending in seed dormancy",0.0206816347379296,2,0.8483925358219894,-0,"embryo development ending in seed dormancy"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.8854406117318956,0.40853039,"embryo development ending in seed dormancy"),
                     c("GO:0010015","root morphogenesis",0.013467340270295237,1,0.8560617638199274,0.43918043,"embryo development ending in seed dormancy"),
                     c("GO:0010090","trichome morphogenesis",0.006709022733481632,1,0.861821252260259,0.35446469,"embryo development ending in seed dormancy"),
                     c("GO:0015031","protein transport",3.093438694074183,3,0.927458842668181,0,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,1,0.9427287900065358,0.38566904,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9405307129667297,0.41259033,"protein transport"),
                     c("GO:0016135","saponin biosynthetic process",9.858960666394757E-06,1,0.9554118934769257,0.03120971,"saponin biosynthetic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,3,0.8780480233264657,0.05779371,"phosphorylation"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0070417","cellular response to cold",0.007012185773973271,1,0.9078628220050563,0.00808117,"cellular response to cold"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8945346792326945,0.4896406,"cellular response to cold"),
                     c("GO:0006952","defense response",1.1604144588756624,1,0.8982880972287655,0.45786068,"cellular response to cold"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.9123593013600194,0.31137091,"cellular response to cold"),
                     c("GO:0009624","response to nematode",0.0011042035946362127,1,0.9274399796732652,0.34321484,"cellular response to cold"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9110428880140717,0.14965038,"cellular response to cold"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8330004593099042,0.30465566,"cellular response to cold"),
                     c("GO:0043434","response to peptide hormone",0.19767216136121488,1,0.8792714644655645,0.18156262,"cellular response to cold"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9918913945769161,0.01286368,"cell wall organization"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9562316419110493,0.05291306,"macromolecule depalmitoylation"),
                     c("GO:0006013","mannose metabolic process",0.050640551462936674,1,0.9316211055970944,0.33438746,"macromolecule depalmitoylation"),
                     c("GO:0006412","translation",4.38869169324396,1,0.8436745944312076,0.25211911,"macromolecule depalmitoylation"),
                     c("GO:0045489","pectin biosynthetic process",0.012338489273993038,1,0.9108732274597405,0.10510175,"macromolecule depalmitoylation"),
                     c("GO:1905421","regulation of plant organ morphogenesis",0.005429822587016912,1,0.9545010570790557,0.09345563,"regulation of plant organ morphogenesis"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1,0.9470816734848149,0.1515582,"regulation of plant organ morphogenesis"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,3,0.9233290044708715,0.38614556,"regulation of plant organ morphogenesis"),
                     c("GO:0010030","positive regulation of seed germination",0.0013630013121290752,1,0.9522524048671861,0.39302598,"regulation of plant organ morphogenesis"),
                     c("GO:0032012","regulation of ARF protein signal transduction",0.05576721100946194,1,0.9542546721464553,0.14116446,"regulation of plant organ morphogenesis"),
                     c("GO:0040008","regulation of growth",0.2209540969749061,1,0.9512874129782445,0.12053119,"regulation of plant organ morphogenesis"),
                     c("GO:0045927","positive regulation of growth",0.06802189911779062,1,0.9501629219113467,0.14325584,"regulation of plant organ morphogenesis"),
                     c("GO:1902066","regulation of cell wall pectin metabolic process",0.00012816648866313183,1,0.9645507253381704,0.16810102,"regulation of plant organ morphogenesis"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9424827065100794,-0,"regulation of auxin polar transport"),
                     c("GO:0010540","basipetal auxin transport",0.0002908393396586453,1,0.9285531757265586,0.45061093,"regulation of auxin polar transport"),
                     c("GO:0048209","regulation of vesicle targeting, to, from or within Golgi",2.4647401665986893E-06,1,0.9661947727929658,0.33066454,"regulation of auxin polar transport"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9449895798698879,0.3459764,"regulation of auxin polar transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}







dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000015","phosphopyruvate hydratase complex",0.041644653723461905,1,0.9086144088931674,0.07573438,"phosphopyruvate hydratase complex"),
                     c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,2,0.9998986033747166,4.647E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,6,0.8824264242600797,0.08417172,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,2,0.8410762486926184,0.17160779,"cytoplasm"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,1,0.9561687738809086,3.018E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,3,0.9120875054971473,0.28538862,"plant-type cell wall"),
                     c("GO:0009536","plastid",0.7624897559369845,5,0.801612826414364,0,"plastid"),
                     c("GO:0000326","protein storage vacuole",0.0009591191132014498,1,0.8819440264886808,0.10321831,"plastid"),
                     c("GO:0005634","nucleus",16.5161752456724,5,0.7191363823434868,0.35145602,"plastid"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.754500473358835,0.2503098,"plastid"),
                     c("GO:0005773","vacuole",1.35235298008496,1,0.7891237398913662,0.21057364,"plastid"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,2,0.6295475531624943,0.19658734,"plastid"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7635287132890735,0.17236349,"plastid"),
                     c("GO:0017119","Golgi transport complex",0.0552089833580907,1,0.714689753644367,0.44260647,"plastid"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.9998889604884231,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,9,0.9998340429243633,9.537E-05,"membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.9999043808194399,5.465E-05,"cell projection"),
                     c("GO:0045202","synapse",1.0801172073373506,1,0.8209852136732425,4.855E-05,"synapse"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}





dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,3,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9675380968275328,0.03247589,"phosphopyruvate hydratase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.775754863567482,0.03898172,"protein serine/threonine phosphatase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8871755471748216,0.45866285,"protein serine/threonine phosphatase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8899807862151725,0.48555173,"protein serine/threonine phosphatase activity"),
                     c("GO:0004559","alpha-mannosidase activity",0.07698914798720577,1,0.908677723345562,0.17913273,"protein serine/threonine phosphatase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,2,0.7403563373645067,0.41376722,"protein serine/threonine phosphatase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8171674534539154,0.44701273,"protein serine/threonine phosphatase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8254607650217357,0.3938196,"protein serine/threonine phosphatase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.8781281494576612,0.32827501,"protein serine/threonine phosphatase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8690471915223179,0.25703449,"protein serine/threonine phosphatase activity"),
                     c("GO:0005085","guanyl-nucleotide exchange factor activity",0.4464279576581196,1,1,-0,"guanyl-nucleotide exchange factor activity"),
                     c("GO:0005515","protein binding",8.610051728351934,4,0.9467882816804609,0.06960742,"protein binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9327176990849406,-0,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9335678771959217,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.9562486155148083,0.27075548,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,6,0.9350402801182242,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9412670072024919,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.9357744210752793,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9494660574333623,0.05562821,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.951874095225543,0.04962823,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9510245110061127,0.03229967,"intramolecular lyase activity"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9326437823545031,0.05445618,"heme binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,3,0.8792267248252142,0.48153711,"heme binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,2,0.986595269489872,0.3967737,"heme binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9119890358065028,0.33073906,"heme binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8995508020656475,0.28189022,"heme binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.8897098336288308,0.15481958,"heme binding"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,1,0.9604459476809007,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9423192183734204,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9212899096591148,0.05306855,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9268058925071231,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9330046203362361,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.931722119464749,0.42311414,"protein dimerization activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9253715940546006,0.03095884,"sphingosine N-acyltransferase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9054205927040189,0.1977699,"sphingosine N-acyltransferase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8987044811614465,0.19133236,"sphingosine N-acyltransferase activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9029540873696541,0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,4,0.8524816754298336,0.26407481,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9047778837205973,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.8992681387269431,0.34784993,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8847467633007552,0.18168089,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,4,0.9999045028833746,4.647E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,8,0.8807663144139185,0.08417172,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,2,0.8308974626973006,0.17160779,"cytoplasm"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,1,0.9653046351644721,3.018E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,8,0.9028962103285528,0.28538862,"plant-type cell wall"),
                     c("GO:0009536","plastid",0.7624897559369845,5,0.7919367026122143,0,"plastid"),
                     c("GO:0000326","protein storage vacuole",0.0009591191132014498,1,0.8744842036391205,0.10321831,"plastid"),
                     c("GO:0005634","nucleus",16.5161752456724,13,0.7099293946394818,0.35145602,"plastid"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.744532265502735,0.2503098,"plastid"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,1,0.6245014656964339,0.20731432,"plastid"),
                     c("GO:0005773","vacuole",1.35235298008496,1,0.7793048162302162,0.21057364,"plastid"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,4,0.6374782178787609,0.19658734,"plastid"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7623733963062687,0.17236349,"plastid"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.9998957949175861,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,16,0.9998470848141064,9.537E-05,"membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.9999097497667101,5.465E-05,"cell projection"),
                     c("GO:0045202","synapse",1.0801172073373506,1,0.8571940349807317,4.855E-05,"synapse"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,1,0.9195873050105273,-0,"ribonucleoprotein complex"),
                     c("GO:0000015","phosphopyruvate hydratase complex",0.041644653723461905,1,0.8682150166326739,0.23819845,"ribonucleoprotein complex"),
                     c("GO:0030117","membrane coat",0.29445950681101296,1,0.7782435040005375,0.29013538,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}



dev.off()



# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,1,1,-0,"intracellular iron ion homeostasis"),
                     c("GO:0007155","cell adhesion",1.1091577223710762,1,0.9918552211007083,0.01321071,"cell adhesion"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,2,0.8534247524356472,-0,"steroid metabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9480000801318849,0.15163831,"steroid metabolic process"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,1,0.8377826242161245,0.48590569,"steroid metabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9468844395095692,0.10901633,"steroid metabolic process"),
                     c("GO:0008203","cholesterol metabolic process",0.2200421431132646,2,0.837604392634995,0.45443823,"steroid metabolic process"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.8177359333941377,0.49264795,"steroid metabolic process"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9321240528154487,0.15725843,"steroid metabolic process"),
                     c("GO:0042759","long-chain fatty acid biosynthetic process",0.04469313344093403,1,0.8440047120986782,0.3981042,"steroid metabolic process"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.8282662942985284,0.42785472,"steroid metabolic process"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,3,0.7958593720914472,0.00987901,"auxin-activated signaling pathway"),
                     c("GO:0002238","response to molecule of fungal origin",0.007682595099288114,2,0.8561750466668394,0.47317835,"auxin-activated signaling pathway"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8798364605318134,0.4896406,"auxin-activated signaling pathway"),
                     c("GO:0006952","defense response",1.1604144588756624,3,0.876429535144067,0.24614308,"auxin-activated signaling pathway"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8443530521025387,0.44951817,"auxin-activated signaling pathway"),
                     c("GO:0009624","response to nematode",0.0011042035946362127,1,0.9061332936456427,0.34321484,"auxin-activated signaling pathway"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.8802376082010779,0.4259419,"auxin-activated signaling pathway"),
                     c("GO:0010030","positive regulation of seed germination",0.0013630013121290752,1,0.9330367249858501,0.39302598,"auxin-activated signaling pathway"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.8960262652098409,0.17080919,"auxin-activated signaling pathway"),
                     c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,1,0.9239894723365195,0.33236762,"auxin-activated signaling pathway"),
                     c("GO:0032012","regulation of ARF protein signal transduction",0.05576721100946194,1,0.9122333943999571,0.46282204,"auxin-activated signaling pathway"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.7926893926416296,0.42828273,"auxin-activated signaling pathway"),
                     c("GO:0040008","regulation of growth",0.2209540969749061,1,0.9296589724366165,0.15374745,"auxin-activated signaling pathway"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8906418163906689,0.46472215,"auxin-activated signaling pathway"),
                     c("GO:0045927","positive regulation of growth",0.06802189911779062,1,0.9098907835228487,0.13260447,"auxin-activated signaling pathway"),
                     c("GO:0048209","regulation of vesicle targeting, to, from or within Golgi",2.4647401665986893E-06,1,0.9556732001053965,0.33066454,"auxin-activated signaling pathway"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,1,0.9131070277074601,0.42863919,"auxin-activated signaling pathway"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9279811128506251,0.1517125,"auxin-activated signaling pathway"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9165286189212768,0.3459764,"auxin-activated signaling pathway"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.9093436558702916,0.32241615,"auxin-activated signaling pathway"),
                     c("GO:1902882","regulation of response to oxidative stress",0.01609475328788944,1,0.9185275199482019,0.12060378,"auxin-activated signaling pathway"),
                     c("GO:1905421","regulation of plant organ morphogenesis",0.005429822587016912,1,0.9397297178437349,0.11290105,"auxin-activated signaling pathway"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9247599954019007,0.10716362,"auxin-activated signaling pathway"),
                     c("GO:0009793","embryo development ending in seed dormancy",0.0206816347379296,2,0.8812401528239566,-0,"embryo development ending in seed dormancy"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.9164516963298328,0.40853039,"embryo development ending in seed dormancy"),
                     c("GO:0010015","root morphogenesis",0.013467340270295237,1,0.8895735484931658,0.43918043,"embryo development ending in seed dormancy"),
                     c("GO:0010047","fruit dehiscence",0.002082705440775892,1,0.9484688509873693,0.46834396,"embryo development ending in seed dormancy"),
                     c("GO:0010090","trichome morphogenesis",0.006709022733481632,1,0.897637974036138,0.35446469,"embryo development ending in seed dormancy"),
                     c("GO:0015031","protein transport",3.093438694074183,4,0.9022621150389559,0,"protein transport"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.1438718130046987,2,0.9110832468027599,0.2741348,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,2,0.9319397887570026,0.38566904,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9292971562616306,0.41259033,"protein transport"),
                     c("GO:0016135","saponin biosynthetic process",9.858960666394757E-06,1,0.9597459419647616,0.03120971,"saponin biosynthetic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,6,0.9001037157078688,0.05779371,"phosphorylation"),
                     c("GO:0044843","cell cycle G1/S phase transition",0.05461124787132716,1,0.9857270870712774,0.00959067,"cell cycle G1/S phase transition"),
                     c("GO:0000911","cytokinesis by cell plate formation",0.008712856488926366,1,0.9864782001597643,0.45372294,"cell cycle G1/S phase transition"),
                     c("GO:0050896","response to stimulus",17.567785530535815,2,1,-0,"response to stimulus"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,2,0.9746267827770626,0.01286368,"cell wall organization"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9526808770000171,0.05291306,"macromolecule depalmitoylation"),
                     c("GO:0006013","mannose metabolic process",0.050640551462936674,1,0.9464052548391848,0.33438746,"macromolecule depalmitoylation"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1,0.8772859143564911,0.4946872,"macromolecule depalmitoylation"),
                     c("GO:0006412","translation",4.38869169324396,1,0.8626879215497344,0.25211911,"macromolecule depalmitoylation"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,1,0.8766733169127588,0.43192082,"macromolecule depalmitoylation"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",0.0004411884898211653,1,0.9205940057980978,0.26905503,"macromolecule depalmitoylation"),
                     c("GO:0045489","pectin biosynthetic process",0.012338489273993038,1,0.9231630072637175,0.10510175,"macromolecule depalmitoylation"),
                     c("GO:1902066","regulation of cell wall pectin metabolic process",0.00012816648866313183,1,0.9466438616787388,0.09252635,"regulation of cell wall pectin metabolic process"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1,0.9193648063327486,0.17609066,"regulation of cell wall pectin metabolic process"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9397706342454366,0.12668641,"regulation of cell wall pectin metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,3,1,-0,"molecular_function"),
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.91463056486156,0.0528646,"mRNA binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,2,0.8983401523761327,0.31664803,"mRNA binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9647688759912599,0.35731485,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,8,0.88239745029626,0.45223656,"mRNA binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8932722122445766,0.24450519,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,1,0.9360712416699458,0.13060477,"mRNA binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9677288198492044,0.03485749,"phosphopyruvate hydratase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,5,0.7343532747717424,0,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8731727486406581,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8737872867029048,0.48555173,"protein kinase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9017348796980814,0.2197454,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.7781290636399422,0.41376722,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8142778824417033,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8227360115156692,0.35123485,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8790819384399416,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8892380938671747,0.31489566,"protein kinase activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9184520922391256,0.21182638,"protein kinase activity"),
                     c("GO:0005085","guanyl-nucleotide exchange factor activity",0.4464279576581196,1,1,-0,"guanyl-nucleotide exchange factor activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,9,0.9504356501827331,0.06960742,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,8,0.9360752783255815,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9421263523713012,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,8,0.9367913538542104,0.12712323,"hydrolase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8827212565347011,0.05402108,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8732006896593303,0.29096067,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8739323593148979,0.32827501,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.8957649409643302,0.17478095,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9500792730745891,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9524195269800306,0.05541399,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9543983189341552,0.03465455,"intramolecular lyase activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9630106320612827,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9487311608463957,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9222392509874762,0.05306855,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9277179683923644,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.933832199442618,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9325696250840736,0.42311414,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9390304801000545,0.37611117,"protein dimerization activity"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.9246415584092138,0.0474854,"dioxygenase activity"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,1,0.9207737630804183,0.4568902,"dioxygenase activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,1,0.9196885870423279,0.42767891,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.9231628623785453,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9532195581978973,0.45929075,"dioxygenase activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9112094691370106,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,7,0.8592894703414051,0.28451694,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9078469500025347,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,2,0.9047787738628846,0.34784993,"iron-sulfur cluster binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.9107662553807204,0.3693843,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8902244543675969,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9107682460071076,0.11631471,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

if(length(topGO_data$GO.ID) > 0)
{
  library(grid)
  library(dplyr)
  
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()



# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.9154451308330057,-0,"monoatomic ion transport"),
                     c("GO:0070588","calcium ion transmembrane transport",0.4224811119566813,1,0.8666153222276416,0.3228373,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.7015775537782857,0,"defense response"),
                     c("GO:0009416","response to light stimulus",0.2589924319660237,1,0.7854965966676558,0.27535725,"defense response"),
                     c("GO:0009733","response to auxin",0.1256056236300358,1,0.7978122908701059,0.25746396,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.675481479067557,0.37522238,"defense response"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,1,0.9411691623050732,0.07105297,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9332400755998402,0.12277613,"lipid catabolic process"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,1,-0,"gametophyte development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8441020893355486,-0,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9062375099417653,0.15925058,"DNA demethylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998480528425548,0,"extracellular region"),
c("GO:0005737","cytoplasm",43.076660800468886,2,0.8388123989769932,0.00024363,"cytoplasm"),
c("GO:0005829","cytosol",14.048080989011824,1,0.8751259862390451,0.17160779,"cytoplasm"),
c("GO:0043231","intracellular membrane-bounded organelle",29.853178741885195,1,0.7517918558885712,0.23467685,"cytoplasm"),
c("GO:0005886","plasma membrane",17.177321395000487,1,0.9997853522349387,0.00015294,"plasma membrane"),
c("GO:0016020","membrane",49.2542153160787,2,0.9997056601254689,9.537E-05,"membrane"),
c("GO:0048046","apoplast",0.0960361495472436,1,0.9999114077929783,3.703E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003677","DNA binding",12.252016067246458,1,0.9175114432932019,0.07367259,"DNA binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9702347536865078,0.04505158,"calmodulin binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,2,0.9426878477320969,0.0830445,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9533385771708953,0.05997119,"lyase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.953905728294811,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.923513779183175,0.04301538,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9394411123873799,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8548448400924897,0.05781294,"ATP hydrolysis activity"),
                     c("GO:0005388","P-type calcium transporter activity",0.05516679695917812,1,0.9553370293207659,0.4743178,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.8780630243161062,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.8435989880930058,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,2,0.8257511538368354,0.38577556,"iron-sulfur cluster binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,2,0.872864697242199,0.23204379,"iron-sulfur cluster binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.8758674032257879,0.01849091,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.8942184621241172,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9821612378176074,0,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998480528425548,0,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8388123989769932,0.00024363,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,1,0.8751259862390451,0.17160779,"cytoplasm"),
                     c("GO:0043231","intracellular membrane-bounded organelle",29.853178741885195,1,0.7517918558885712,0.23467685,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,1,0.9997853522349387,0.00015294,"plasma membrane"),
                     c("GO:0016020","membrane",49.2542153160787,2,0.9997056601254689,9.537E-05,"membrane"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999114077929783,3.703E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.9154451308330057,-0,"monoatomic ion transport"),
c("GO:0070588","calcium ion transmembrane transport",0.4224811119566813,1,0.8666153222276416,0.3228373,"monoatomic ion transport"),
c("GO:0006952","defense response",1.1604144588756624,2,0.7015775537782857,0,"defense response"),
c("GO:0009416","response to light stimulus",0.2589924319660237,1,0.7854965966676558,0.27535725,"defense response"),
c("GO:0009733","response to auxin",0.1256056236300358,1,0.7978122908701059,0.25746396,"defense response"),
c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.675481479067557,0.37522238,"defense response"),
c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
c("GO:0016042","lipid catabolic process",1.4093137798594644,1,0.9411691623050732,0.07105297,"lipid catabolic process"),
c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9332400755998402,0.12277613,"lipid catabolic process"),
c("GO:0048229","gametophyte development",0.017253181166190824,1,1,-0,"gametophyte development"),
c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8441020893355486,-0,"DNA demethylation"),
c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9062375099417653,0.15925058,"DNA demethylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
c("GO:0003677","DNA binding",12.252016067246458,1,0.9175114432932019,0.07367259,"DNA binding"),
c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9702347536865078,0.04505158,"calmodulin binding"),
c("GO:0016787","hydrolase activity",22.243200201100027,2,0.9426878477320969,0.0830445,"hydrolase activity"),
c("GO:0016829","lyase activity",3.6221510367468346,1,0.9533385771708953,0.05997119,"lyase activity"),
c("GO:0016874","ligase activity",3.248266153118885,1,0.953905728294811,0.0589874,"ligase activity"),
c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.923513779183175,0.04301538,"acid-amino acid ligase activity"),
c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9394411123873799,0.4447539,"acid-amino acid ligase activity"),
c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8548448400924897,0.05781294,"ATP hydrolysis activity"),
c("GO:0005388","P-type calcium transporter activity",0.05516679695917812,1,0.9553370293207659,0.4743178,"ATP hydrolysis activity"),
c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.8780630243161062,0.26945252,"ATP hydrolysis activity"),
c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.8435989880930058,-0,"iron-sulfur cluster binding"),
c("GO:0000166","nucleotide binding",18.476222463002568,2,0.8257511538368354,0.38577556,"iron-sulfur cluster binding"),
c("GO:0046872","metal ion binding",18.074696526070646,2,0.872864697242199,0.23204379,"iron-sulfur cluster binding"),
c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.8758674032257879,0.01849091,"catalytic activity, acting on DNA"),
c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.8942184621241172,0.41904136,"catalytic activity, acting on DNA"),
c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9821612378176074,0,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",0.9427113541805001,2,1,-0,"immune system process"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,7,0.9186356768617372,0.01352937,"chromatin organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,2,0.9390813746700919,0.36382358,"chromatin organization"),
                     c("GO:0007032","endosome organization",0.07057290519022028,1,0.9401718046072549,0.34948919,"chromatin organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,1,0.9350750787398051,0.41402298,"chromatin organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9372819581737211,0.38649952,"chromatin organization"),
                     c("GO:0010027","thylakoid membrane organization",0.08263041408522105,3,0.929376668550639,0.32157946,"chromatin organization"),
                     c("GO:0010405","arabinogalactan protein metabolic process",0.00041161160782198105,1,0.9116607404941047,0.48519003,"chromatin organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9384416533666503,0.37383592,"chromatin organization"),
                     c("GO:0042254","ribosome biogenesis",2.136224808753416,5,0.9175878132704369,0.40498561,"chromatin organization"),
                     c("GO:0051017","actin filament bundle assembly",0.08024701034412013,1,0.9292945375911447,0.45180791,"chromatin organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9401720574956033,0.4524177,"chromatin organization"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,5,0.9262728972101378,0.4063955,"chromatin organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9584379981320287,0.20937319,"chromatin organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9543707006702298,0.22611843,"chromatin organization"),
                     c("GO:0006952","defense response",1.1604144588756624,20,0.8957636459317844,-0,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,3,0.8903484054146903,0.47274124,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,4,0.8943590811635584,0.49159695,"defense response"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,6,0.835229802251021,0.23262547,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8867081186019515,0.25765651,"defense response"),
                     c("GO:0010188","response to microbial phytotoxin",0.00018485551249490168,1,0.9336715752723834,0.27319058,"defense response"),
                     c("GO:0010196","nonphotochemical quenching",0.0005693549784842973,1,0.9231741544045722,0.37622087,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9226374812104545,0.20406574,"defense response"),
                     c("GO:0010343","singlet oxygen-mediated programmed cell death",0.0019717921332789512,1,0.8852394962381888,0.44410187,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8480469496818611,0.38492887,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,2,0.8595221172259188,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8927239928546409,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,5,0.8212463874064555,0.40205681,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9073299609359704,0.46171648,"defense response"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,3,0.8633771504728337,0.35855228,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,7,0.8864965580958846,0.45104908,"defense response"),
                     c("GO:0051103","DNA ligation involved in DNA repair",0.04254388001565997,1,0.8144245704915166,0.43114998,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9236295310126758,0.4243543,"defense response"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.7957322799946543,0.41244237,"defense response"),
                     c("GO:0071461","cellular response to redox state",0.0002957688199918427,1,0.9516903714411864,0.16681679,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.9260116244231166,0.31063873,"defense response"),
                     c("GO:1990414","replication-born double-strand break repair via sister chromatid exchange",0.027676567330736677,1,0.8186751991451701,0.39308191,"defense response"),
                     c("GO:0007049","cell cycle",2.8073119376140743,3,0.9887995062983475,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9821067489482855,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,3,0.9968868658947184,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,7,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,8,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9916705180602594,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9901453047902254,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,3,0.9485961331730296,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,2,0.9455384672030147,0.21768792,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.916406665515398,0.16846091,"biosynthetic process"),
                     c("GO:0071704","organic substance metabolic process",50.871943734517124,1,0.9430877798011215,0.29474266,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9338959248200909,0.17227516,"biosynthetic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9530544271030384,0.06359719,"toxin metabolic process"),
                     c("GO:0009813","flavonoid biosynthetic process",0.06199314467029023,3,0.9280455846162399,0.05649883,"flavonoid biosynthetic process"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,2,0.9003447070443624,0.20397686,"flavonoid biosynthetic process"),
                     c("GO:0009693","ethylene biosynthetic process",0.004658358914871523,1,0.9362773595837025,0.10463112,"flavonoid biosynthetic process"),
                     c("GO:0009697","salicylic acid biosynthetic process",0.021366832504244038,2,0.8527037207100914,0.31042729,"flavonoid biosynthetic process"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,1,0.8928899941250885,0.34289816,"flavonoid biosynthetic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,1,0.9006858251553228,0.3067044,"flavonoid biosynthetic process"),
                     c("GO:0042372","phylloquinone biosynthetic process",0.033601802691239926,2,0.8853206771987169,0.11806969,"flavonoid biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9568121087752541,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,9,0.8908392724438149,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,3,0.9230354927540166,0.4203896,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9421605704579668,0.48835672,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,1,0.9381000734135579,0.48720779,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.852418032883764,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9294022137406913,0.29044182,"flower development"),
                     c("GO:0022038","corpus callosum development",0.002092564401442287,1,0.9302073632677357,0.39347824,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.928925869834166,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9236609156308953,0.42767441,"flower development"),
                     c("GO:0048262","determination of dorsal/ventral asymmetry",0.0029848003417510126,2,0.9271829185931671,0.38718574,"flower development"),
                     c("GO:0090626","plant epidermis morphogenesis",0.016804598455869863,1,0.9183486446093306,0.47435352,"flower development"),
                     c("GO:1905392","plant organ morphogenesis",0.021231271795081108,1,0.9034616319762211,0.44850026,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,1,0.9321525162187911,0.34979162,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.993597558808241,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9549422791824785,0.04613723,"chlororespiration"),
                     c("GO:0015031","protein transport",3.093438694074183,14,0.9103422405483595,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,7,0.9460843053691552,0.42145531,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9507515579714438,0.35562027,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,7,0.9398786887052012,0.38566904,"protein transport"),
                     c("GO:0097339","glycolate transmembrane transport",0.0003376694028240204,1,0.951599441018805,0.41757087,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9485023667597943,0.49534489,"protein transport"),
                     c("GO:1901975","glycerate transmembrane transport",0.00010351908699714494,1,0.952832221721286,0.26680356,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,4,0.9283663106447313,0.35175462,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,3,0.9492701579522043,0.04477571,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,20,0.8787817727983146,0,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8608063849588752,0.40302693,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,5,0.9350254592557123,0.15163831,"protein ubiquitination"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",18.82728500884891,2,0.8436208381558299,0.27206623,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,8,0.8195192732168812,0.48035482,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,12,0.8074755633777055,0.4385246,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,7,0.8493320361736643,0.47678922,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,7,0.8528425168304801,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9349826835902781,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,9,0.9334748490371114,0.12094833,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,2,0.907561429820102,0.45718665,"protein ubiquitination"),
                     c("GO:0008652","amino acid biosynthetic process",2.679426429329935,4,0.8187448752720609,0.49906638,"protein ubiquitination"),
                     c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.8641141306992849,0.46840296,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.8636457055124153,0.16277389,"protein ubiquitination"),
                     c("GO:0018258","protein O-linked glycosylation via hydroxyproline",0.0002957688199918427,1,0.8965164020303701,0.36854741,"protein ubiquitination"),
                     c("GO:0018364","peptidyl-glutamine methylation",1.9717921332789515E-05,1,0.9407597825211361,0.31815889,"protein ubiquitination"),
                     c("GO:0030212","hyaluronan metabolic process",0.023962203899672456,1,0.916362091979296,0.14846587,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,3,0.8955529307937089,0.14129671,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,1,0.8509575240184966,0.28833637,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,1,0.9126631931540495,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,2,0.9141606037501719,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.8740122109151391,0.26614879,"protein ubiquitination"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8937726619685422,0.16428511,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,4,0.9073240796496012,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,1,0.9074103442869037,0.41108531,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8553631478800025,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.88211423643748,0.43527038,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,1,0.8748603118957694,0.32141303,"protein ubiquitination"),
                     c("GO:0032259","methylation",2.6278542060840238,4,0.9630452431500669,0.05843136,"methylation"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9743453680080014,0.03897815,"collagen metabolic process"),
                     c("GO:0045454","cell redox homeostasis",0.337647220162521,2,0.9782047884873202,-0,"cell redox homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9847712188746586,0.42883749,"cell redox homeostasis"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,9,0.8795963347957488,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9328147254894259,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,5,0.8997569399837453,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0006450","regulation of translational fidelity",0.29561600610151356,1,0.9347178755935446,0.18254504,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9451145884433402,0.38949219,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,2,0.9481091027917375,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9453653482900994,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9080353016145205,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0030512","negative regulation of transforming growth factor beta receptor signaling pathway",0.013122276646971421,1,0.9070792414433485,0.42782958,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.935007086825738,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9208022567925225,0.31416992,"positive regulation of DNA-templated transcription"),
                     c("GO:0042548","regulation of photosynthesis, light reaction",0.0070861279789712316,1,0.9351907319277292,0.23746411,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,2,0.9449907287239626,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0042981","regulation of apoptotic process",0.6083841390223874,1,0.9291369098141573,0.22270387,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9328435266119133,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0043484","regulation of RNA splicing",0.23003419974865566,1,0.9199613319251218,0.37266406,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9322503081007453,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,4,0.8854071590941985,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,1,0.9418206371631411,0.41192877,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9393485456220058,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9306901176114808,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080038","positive regulation of cytokinin-activated signaling pathway",0.0001996439534944938,1,0.9225478867658441,0.35146486,"positive regulation of DNA-templated transcription"),
                     c("GO:0090070","positive regulation of ribosome biogenesis",0.00022675609532707942,1,0.9433122010325001,0.35407617,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.946765199196527,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,3,0.9118699682618242,0.37396531,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9412071089323403,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,2,0.912096614000361,0.48792077,"positive regulation of DNA-templated transcription"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9439336984968786,0.49983935,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,3,0.9143430387136154,0.13990144,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,2,0.9270459558020184,0.36922831,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.9224841775726509,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000032","regulation of secondary shoot formation",0.007307954593965113,1,0.93652840146289,0.46431004,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9403446562038746,0.41762823,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9225784256688021,0.49708545,"positive regulation of DNA-templated transcription"),
                     c("GO:2000280","regulation of root development",0.008382581306602141,1,0.93764660082577,0.1419752,"positive regulation of DNA-templated transcription"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9397099562114293,0.15519254,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9340937151495133,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.926396799635331,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,2,0.8901036103505938,0.07681302,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,2,0.9113751604809034,0.3693131,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006221","pyrimidine nucleotide biosynthetic process",0.5101395959817636,2,0.8020872555179069,0.41640622,"phosphatidylinositol dephosphorylation"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9355979539909652,0.31764343,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,35,0.896548239281793,0.47540491,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016311","dephosphorylation",0.5642283189377719,1,0.9168936756356868,0.32207069,"phosphatidylinositol dephosphorylation"),
                     c("GO:1901137","carbohydrate derivative biosynthetic process",5.186311188037295,1,0.8715682360898135,0.47912132,"phosphatidylinositol dephosphorylation"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,5,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,1,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,2,0.9894456293343421,0.01381155,"cell division"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9518651679962883,0.09043985,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9263950558989176,0.08958079,"intrachromosomal DNA recombination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8815458283936144,0.15020607,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,7,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,10,0.9999257530914285,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9777373502571378,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,93,0.8953425503676568,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,41,0.8659805935478573,0.17160779,"cytoplasm"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,11,0.9749079730975561,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,46,0.7318287361495825,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,3,0.8115344680290907,0.13280054,"chloroplast"),
                     c("GO:0000786","nucleosome",0.17397327417075828,3,0.6654943895171002,0.41609548,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,126,0.7380939282232801,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,17,0.7647495322700621,0.26090848,"chloroplast"),
                     c("GO:0005742","mitochondrial outer membrane translocase complex",0.03996992242217233,1,0.6251152373967991,0.45499251,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,6,0.792921125079359,0.21802429,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,9,0.7956693369993979,0.17354881,"chloroplast"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,19,0.7078727052934302,0.2030659,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,1,0.7855235907068862,0.42381054,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,7,0.7293296072016058,0.15783714,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,3,0.7418478568689448,0.49095344,"chloroplast"),
                     c("GO:0008278","cohesin complex",0.03346232408674592,1,0.7116417682771674,0.4995934,"chloroplast"),
                     c("GO:0009522","photosystem I",0.033362933505067006,2,0.6778045196054654,0.35622153,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,47,0.8033821939658199,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,1,0.8308007132050742,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,5,0.78174159519074,0.43399423,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6933246414373261,0.40758865,"chloroplast"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,5,0.7125904709330525,0.23571918,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.8024356179131102,0.23050241,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.7995452269663076,0.4052793,"chloroplast"),
                     c("GO:0032586","protein storage vacuole membrane",0.00013169252072455135,1,0.8005647844849978,0.39770966,"chloroplast"),
                     c("GO:0033097","amyloplast membrane",7.4542936259180024E-06,1,0.7956395472568772,0.45278409,"chloroplast"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,1,0.7625147011790306,0.44880082,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999383493815162,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8871293438012026,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,4,0.9999203225276694,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,135,0.9998923458366435,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999439252066005,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9999477722348837,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,3,0.9101336775332723,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0062023","collagen-containing extracellular matrix",0.3197519250837527,1,0.9736314571567926,3.527E-05,"collagen-containing extracellular matrix"),
                     c("GO:0005886","plasma membrane",17.177321395000487,50,0.9340226189443323,0.35740582,"collagen-containing extracellular matrix"),
                     c("GO:0019897","extrinsic component of plasma membrane",0.23434808301161014,1,0.9444410808444494,0.22734667,"collagen-containing extracellular matrix"),
                     c("GO:0090404","pollen tube tip",0.00040998614942549014,1,0.9999633644803905,2.175E-05,"pollen tube tip"),
                     c("GO:0098552","side of membrane",0.7241349304670945,1,0.9689450210341369,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,11,0.785911280813043,-0,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7508173070633546,0.42313164,"ribonucleoprotein complex"),
                     c("GO:0005667","transcription regulator complex",0.7880231963702957,1,0.8151175856854667,0.32589738,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7250098958054277,0.38761376,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8106290470898265,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7951516237973681,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005853","eukaryotic translation elongation factor 1 complex",0.01002354016231774,1,0.8216842407429149,0.21073066,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.7697752543865839,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.7715584015735832,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.6829025228810525,0.25767523,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7937233990166496,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7507656679905836,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7743188138755523,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,5,0.7722807812800847,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7888424149240479,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,13,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,4,0.9611351788227029,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,23,0.9397627743394634,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9552139747970011,0.42104634,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,46,0.9304571034993059,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9582993692262955,0.19519762,"mRNA binding"),
                     c("GO:0003713","transcription coactivator activity",0.24395086691346182,2,0.9646157548972838,0.29461398,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9486811596841636,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.948376383859577,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,5,0.9273718332035783,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9545072611422076,0.428381,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,5,0.961702143109721,0.13060477,"mRNA binding"),
                     c("GO:0032934","sterol binding",0.11902838493358808,1,0.9602753423560685,0.10144903,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9500987454666109,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,11,0.9323975125272463,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,4,0.9963047271514212,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,20,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,28,0.8012889512105662,0,"protein kinase activity"),
                     c("GO:0003878","ATP citrate synthase activity",0.010080724293609766,1,0.90111202967131,0.32478201,"protein kinase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9196246748807231,0.20796201,"protein kinase activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,4,0.8692663224412006,0.41115081,"protein kinase activity"),
                     c("GO:0004659","prenyltransferase activity",0.25118876486646274,1,0.8956877844769857,0.26852283,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8507841213030327,0.49507895,"protein kinase activity"),
                     c("GO:0004721","phosphoprotein phosphatase activity",0.6445676316147657,3,0.827746251708972,0.4549882,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,8,0.8195039001097916,0.47753726,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,4,0.8652158530088412,0.35677507,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8946693927575179,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8985318456786945,0.40329157,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.854758657386247,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,5,0.877182289485356,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,9,0.8821684863712335,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,2,0.9012515924548455,0.21201958,"protein kinase activity"),
                     c("GO:0019786","protein-phosphatidylethanolamide deconjugating activity",0.010599617713034465,1,0.903575015657121,0.30543886,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.9146798450832966,0.49223587,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8797085259364104,0.49253483,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9147899061265138,0.30679413,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9229099827044737,0.29376839,"protein kinase activity"),
                     c("GO:0051741","2-methyl-6-phytyl-1,4-benzoquinone methyltransferase activity",0.0016742074002805483,1,0.9169995905573177,0.47516484,"protein kinase activity"),
                     c("GO:0051742","2-methyl-6-solanyl-1,4-benzoquinone methyltransferase activity",1.1087466227023499E-05,1,0.9344914229851815,0.12943404,"protein kinase activity"),
                     c("GO:0052630","UDP-N-acetylgalactosamine diphosphorylase activity",4.4349864908093996E-05,1,0.9201347359742593,0.25791543,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9138386017595478,0.3846796,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9282637979198199,0.49777323,"protein kinase activity"),
                     c("GO:0070204","2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylic-acid synthase activity",0.007958583257757468,1,0.9201019614889142,0.19601592,"protein kinase activity"),
                     c("GO:0070569","uridylyltransferase activity",0.1030402586342202,1,0.8813259139456104,0.45345646,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.8989056130651917,0.18579171,"protein kinase activity"),
                     c("GO:0102550","2-methyl-6-geranylgeranyl-1,4-benzoquinol methyltransferase activity",7.76122635891645E-05,1,0.9297380499560842,0.39293218,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.9000121762031924,0.39355689,"protein kinase activity"),
                     c("GO:1990714","hydroxyproline O-galactosyltransferase activity",0.0006475080276581724,1,0.9193389090433758,0.42279362,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,73,0.9709642464979593,0.07766851,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,13,0.9431261818443263,0.05738816,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,65,0.9198304274296102,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,3,0.951980980576073,0.3830934,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,26,0.9461399564400542,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,6,0.9533692929286738,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,7,0.9363475014305243,0.19635374,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9451857612930475,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9468571600713641,0.29697106,"zinc ion binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,2,0.9474421564492345,0.24057124,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.9450352440119275,0.2885161,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9722946711923405,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.9385909232649782,0.33990845,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,3,0.9443192341157636,0.34531063,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,2,0.9432731516283231,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,5,0.9495443006076708,0.18168089,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,2,0.9493950186232514,0.29755636,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,5,0.9768257194973491,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,1,0.9942543233776147,-0,"photoreceptor activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,6,0.8968261985078174,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,6,0.8982143348026331,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,2,0.8907659990398528,0.40806599,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9189858187898816,0.31075808,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,2,0.890425659736503,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.8961949187017461,0.41247223,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9007265043525489,0.44173969,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,1,0.9103083375357702,0.3863248,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016639","oxidoreductase activity, acting on the CH-NH2 group of donors, NAD or NADP as acceptor",0.09697097962154752,1,0.9115284848589319,0.34937139,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016711","flavonoid 3'-monooxygenase activity",8.20472500799739E-05,1,0.937270470200862,0.21343278,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9078332520701534,0.4215727,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0045486","naringenin 3-dioxygenase activity",0.00033705897330151436,1,0.9297661353959318,0.38276908,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0050664","oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor",0.059595130970251306,1,0.915726359435088,0.3347022,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,3,0.9030988232273004,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:1990136","linoleate 9S-lipoxygenase activity",0.00026609918944856393,1,0.9349402817373502,0.22819902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,74,0.9460377115566114,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,26,0.9499681698786101,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,38,0.9455067431342304,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,8,0.9559105526427042,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,2,0.9419042058348917,0.04651897,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9569626007807175,0.45277577,"carboxy-lyase activity"),
                     c("GO:0043748","O-succinylbenzoate synthase activity",0.006195676127660732,1,0.9555538612351807,0.43317182,"carboxy-lyase activity"),
                     c("GO:0070205","2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase activity",0.0028051289554369453,1,0.9579522593167777,0.41110877,"carboxy-lyase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9587081143100623,0.03072632,"magnesium chelatase activity"),
                     c("GO:0004088","carbamoyl-phosphate synthase (glutamine-hydrolyzing) activity",0.05676560958911491,1,0.9534244275930002,0.49740294,"magnesium chelatase activity"),
                     c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9467497923961009,0.43449047,"magnesium chelatase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,8,0.9576985621983991,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,7,0.9564055825885092,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,14,0.8623500997419943,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004427","inorganic diphosphate phosphatase activity",0.053847388478162325,1,0.9216085953602194,0.49693897,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,3,0.8986105214288743,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0004636","phosphoribosyl-ATP diphosphatase activity",0.02123693281124081,1,0.9263463719477826,0.45732449,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.903912906933785,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9041668247449802,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9019281482572177,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0017057","6-phosphogluconolactonase activity",0.03804774910465384,1,0.8969452217644882,0.49474271,"ATP hydrolysis activity"),
                     c("GO:0035529","NADH pyrophosphatase activity",0.024569825159084076,1,0.922203574775657,0.46310876,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,1,0.9261377070261394,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047714","galactolipase activity",0.0037120836928074678,1,0.9070072417778053,0.35781567,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9288473613191554,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9262021708941423,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9440469655202162,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,3,0.9774252496409356,0.04901613,"carbohydrate binding"),
                     c("GO:0030414","peptidase inhibitor activity",0.25608720744556174,2,0.9711923644749612,-0,"peptidase inhibitor activity"),
                     c("GO:0005094","Rho GDP-dissociation inhibitor activity",0.010056331867910313,1,0.9778961338048868,0.49485349,"peptidase inhibitor activity"),
                     c("GO:0043879","glycolate transmembrane transporter activity",0.00030379657462044384,1,0.9484034003314635,-0,"glycolate transmembrane transporter activity"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9456378334776795,0.2866085,"glycolate transmembrane transporter activity"),
                     c("GO:0009678","diphosphate hydrolysis-driven proton transmembrane transporter activity",0.01596151638042303,1,0.9316782497316297,0.42292952,"glycolate transmembrane transporter activity"),
                     c("GO:0015116","sulfate transmembrane transporter activity",0.08090967604508129,1,0.9275559067027012,0.28203837,"glycolate transmembrane transporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9161894325329503,0.35239237,"glycolate transmembrane transporter activity"),
                     c("GO:0015165","pyrimidine nucleotide-sugar transmembrane transporter activity",0.046425438585792796,1,0.9464063856271305,0.22592317,"glycolate transmembrane transporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9398362287280138,0.47124298,"glycolate transmembrane transporter activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9754100624782831,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9946516324784443,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8859491001696689,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9590108201390692,0.02641673,"carotenoid isomerase activity"),
                     c("GO:0004165","delta(3)-delta(2)-enoyl-CoA isomerase activity",0.021527424426388823,1,0.9516491001271186,0.49767449,"carotenoid isomerase activity"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9574994292018759,0.44355864,"carotenoid isomerase activity"),
                     c("GO:0009982","pseudouridine synthase activity",0.21083704027983347,1,0.9417551056910062,0.41191472,"carotenoid isomerase activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,1,0.9756286294269892,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9825458914650622,0.49249984,"ER retention sequence binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,5,0.925699410242818,0.0463004,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,3,0.9304151030336154,0.40485813,"unfolded protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,4,0.922196046603978,0.43412799,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9382433373927986,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,4,0.920971180283221,0.42444044,"unfolded protein binding"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.931962704354966,0.39461121,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,1,0.9345748502102412,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9338519918891887,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,2,0.9329831565563115,0.38789043,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9542223152176578,0.25391209,"unfolded protein binding"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,1,0.9310567427136419,0.4006017,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,2,0.9278226197158785,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,2,0.926708098819392,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,5,0.9328164464472001,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,2,0.9315624979994015,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,3,0.9379549488732737,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9211149610924724,0.4679269,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9505342250815753,0.26331899,"unfolded protein binding"),
                     c("GO:0070678","preprotein binding",3.326239868107049E-05,1,0.9592829924781565,0.22354599,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.9431309481722314,0.31312228,"unfolded protein binding"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,1,0.9378810432751246,0.35601456,"unfolded protein binding"),
                     c("GO:0140098","catalytic activity, acting on RNA",3.705433430564499,1,0.9134036898318094,0.06018061,"catalytic activity, acting on RNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9403686812373283,0.42852185,"catalytic activity, acting on RNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8788659907729516,0.07971789,"phosphatidic acid biosynthetic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9412418563513237,0.08924942,"biosynthetic process"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.8973081327946559,0.07593517,"L-ascorbic acid metabolic process"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8683715830175431,-0,"response to chemical"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,1,0.8684741202947966,0.3661968,"response to chemical"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8759140500361217,0.47274124,"response to chemical"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8227013881005394,0.46171648,"response to chemical"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8663195115509706,0.42391909,"response to chemical"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8440369243320586,0.36420456,"response to chemical"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8116479779945531,-0,"proteolysis involved in protein catabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8317136089843344,0.35126771,"proteolysis involved in protein catabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8332193434844108,0.38896703,"proteolysis involved in protein catabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,1,0.8316816435968096,0.47291952,"proteolysis involved in protein catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.917821196366319,0.12747018,"proteolysis involved in protein catabolic process"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9254120545412087,0.17097988,"proteolysis involved in protein catabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.8835428366968173,0.28805767,"proteolysis involved in protein catabolic process"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.983541002359366,-0,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.9651269603562943,0.12847604,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9612644040144976,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9124506377030752,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9573645219029069,0.35287398,"acquisition of seed longevity"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.8973712742484832,0,"proton transmembrane transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.89047721545884,0.4188028,"proton transmembrane transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8663934856093372,0.34990498,"proton transmembrane transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.8968533185711883,0.26647644,"proton transmembrane transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000502","proteasome complex",0.37652382533874423,1,0.7864587092017606,-0,"proteasome complex"),
                     c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
                     c("GO:0005737","cytoplasm",43.076660800468886,5,0.876265050477576,0.09354521,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,4,0.8431111418366491,0.17160779,"cytoplasm"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7902280744401448,0,"vacuole"),
                     c("GO:0005634","nucleus",16.5161752456724,3,0.7269108695242845,0.35145602,"vacuole"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7787985968197498,0.21968674,"vacuole"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7570212294648768,0.25147478,"vacuole"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7054070726892268,0.18309297,"vacuole"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7479319672799656,0.22703414,"vacuole"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.8041507518715392,0.18302188,"vacuole"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7761956513619133,0.21139748,"vacuole"),
                     c("GO:0005886","plasma membrane",17.177321395000487,2,0.9824062225654548,0.00015294,"plasma membrane"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.999929335396792,4.05E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,2,0.9998501802029992,7.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()


# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9203079142606666,0.03137215,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.8942300932197402,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9121037444404009,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9624906396952049,0.0386588,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,1,0.9507531512872165,0.05532665,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,3,0.9342685727484303,0.07494256,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9402160380017902,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,1,0.9349396667078802,0.12712323,"hydrolase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9694685659422415,0.02287106,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.950875340135772,0.03319032,"isomerase activity"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9074217519528892,-0,"glutathione binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9035224028271102,0.31015716,"glutathione binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,1,0.8641751877358295,0.26746377,"glutathione binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.904150293886289,0.12083997,"glutathione binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.8800172737207748,0.02352266,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0046961","proton-transporting ATPase activity, rotational mechanism",0.20284297713014954,1,0.9592398403830547,-0,"proton-transporting ATPase activity, rotational mechanism"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,1,0.8665610802126985,0,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.9143662058292141,0.42844795,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9001839323705105,0.30149286,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8878879162354781,0.34866016,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,1,0.8424252332778702,0.20765779,"N-acylphosphatidylethanolamine-specific phospholipase D activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,3,0.9998978355557501,0,"extracellular region"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,2,0.6437253130912475,5.394E-05,"Golgi apparatus"),
                     c("GO:0005634","nucleus",16.5161752456724,9,0.6698834627919237,0.35145602,"Golgi apparatus"),
                     c("GO:0005730","nucleolus",1.2114693153196516,1,0.7546055924229319,0.20773614,"Golgi apparatus"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8668703845762429,0.1047991,"Golgi apparatus"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.7104688145276024,0.2503098,"Golgi apparatus"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,1,0.5665642583874467,0.20731432,"Golgi apparatus"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.7086931775584476,0.19462727,"Golgi apparatus"),
                     c("GO:0009536","plastid",0.7624897559369845,1,0.7618898553001906,0.19658734,"Golgi apparatus"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.652862337588118,0.25147478,"Golgi apparatus"),
                     c("GO:0016020","membrane",49.2542153160787,8,0.9998355907512326,9.537E-05,"membrane"),
                     c("GO:0031012","extracellular matrix",0.6517835565339344,1,0.9626729933883422,4.559E-05,"extracellular matrix"),
                     c("GO:0005886","plasma membrane",17.177321395000487,5,0.90272072198916,0.39486578,"extracellular matrix"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,1,0.900045229888717,-0,"ribonucleoprotein complex"),
                     c("GO:0030117","membrane coat",0.29445950681101296,1,0.7080464488256948,0.29013538,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}





dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9423177302167721,-0,"DNA-binding transcription factor activity, RNA polymerase II-specific"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,3,0.7712757143447455,0,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8720970226076804,0.35580262,"protein kinase activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,6,0.939878503395764,0.06452231,"protein binding"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.939315351942652,0.07956199,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9445879769274047,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,3,0.9386239715185933,0.12712323,"transferase activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,1,0.9556232196159946,0.04919141,"carbohydrate binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9013056627527537,0.04828594,"protein dimerization activity"),
                     c("GO:0008083","growth factor activity",0.15856185451266308,1,0.8796806122877222,0.41153522,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9148935368882537,0.41374711,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9212121829196258,0.37611117,"protein dimerization activity"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8874789418037062,0.0474854,"dioxygenase activity"),
                     c("GO:0003954","NADH dehydrogenase activity",0.34065796483880617,1,0.8939506935615665,0.36657586,"dioxygenase activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8648493115254148,0.43460917,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.8821277040160563,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9255406726276421,0.45929075,"dioxygenase activity"),
                     c("GO:0051287","NAD binding",0.8468961503119814,1,0.8756421678359987,-0,"NAD binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,3,0.8164320127475632,0.42621458,"NAD binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,1,0.8745510615030264,0.31664803,"NAD binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.872334675051111,0.18774095,"NAD binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,5,0.8481322633103344,0.45223656,"NAD binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.8897315936146176,0.12297344,"NAD binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.9111581046147357,0.14439095,"NAD binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9075733758048868,0.21206178,"NAD binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,4,0.8758810110223696,0.45949964,"NAD binding"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9587271113193424,0.03057439,"pectin acetylesterase activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.937769329963123,0.17194065,"pectin acetylesterase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.9365366538460628,0.29096067,"pectin acetylesterase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}

# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,1,0.9802923737144542,0.05630604,"carbohydrate metabolic process"),
                     c("GO:0008152","metabolic process",57.597931274565454,1,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9924509524870899,0.00908442,"cell population proliferation"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,2,0.7343922169732996,0,"auxin-activated signaling pathway"),
                     c("GO:0001558","regulation of cell growth",0.15345225803226778,1,0.9035121795101533,0.16011872,"auxin-activated signaling pathway"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.1438718130046987,1,0.890996607463209,0.21941245,"auxin-activated signaling pathway"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.8192798869174048,0.24614308,"auxin-activated signaling pathway"),
                     c("GO:0007165","signal transduction",8.756627809544993,2,0.7069829307373783,0.47698259,"auxin-activated signaling pathway"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.7910388589017517,0.31987157,"auxin-activated signaling pathway"),
                     c("GO:0009643","photosynthetic acclimation",0.000325345701991027,1,0.8819377050944387,0.45106647,"auxin-activated signaling pathway"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8208711584340853,0.4259419,"auxin-activated signaling pathway"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9175885584765684,0.20908612,"auxin-activated signaling pathway"),
                     c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,1,0.8994986328820143,0.35474033,"auxin-activated signaling pathway"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,1,0.9124823657207489,0.28549545,"auxin-activated signaling pathway"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8411220649809261,0.47074713,"auxin-activated signaling pathway"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,1,0.8432302370486209,0.19569204,"auxin-activated signaling pathway"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,1,0.8757135091819515,0.42863919,"auxin-activated signaling pathway"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.9127106875120198,0.18626259,"auxin-activated signaling pathway"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9027616437444468,0.1517125,"auxin-activated signaling pathway"),
                     c("GO:0060918","auxin transport",0.015500750907739157,1,0.8740624142281647,0.11396033,"auxin-activated signaling pathway"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.8723847562786669,0.32241615,"auxin-activated signaling pathway"),
                     c("GO:1902882","regulation of response to oxidative stress",0.01609475328788944,1,0.884966155902162,0.12060378,"auxin-activated signaling pathway"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9723606319039315,0.03847046,"coumarin biosynthetic process"),
                     c("GO:0006564","L-serine biosynthetic process",0.06954757328091521,1,0.9625908650794829,0.11280602,"coumarin biosynthetic process"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",0.0004411884898211653,1,0.960418450021046,0.00539809,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,1,0.9416590292713564,0.36594593,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,3,0.9638208410633946,0.25284502,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0030154","cell differentiation",2.279586420543629,1,0.9743129543280651,0.01240155,"cell differentiation"),
                     c("GO:0010162","seed dormancy process",0.001212652161966555,1,0.9636117781828409,0.40373092,"cell differentiation"),
                     c("GO:0044843","cell cycle G1/S phase transition",0.05461124787132716,1,0.9928079783171261,0.00716697,"cell cycle G1/S phase transition"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9378672686661721,0.00884962,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,1,0.9510317837341355,0.4063955,"cell wall organization"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1,0.9103438370458188,0.37183399,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])


hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)


if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}







dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9717640153210997,0.04505158,"calmodulin binding"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.871105943940637,-0,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.9163165329214024,0.35458321,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9059564859780858,0.35290783,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,2,0.9583291816290221,0.0830445,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9660993454886375,0.05997119,"lyase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9665127633927455,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9443095294759732,0.04301538,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9559580965807031,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8937513058969788,0.05781294,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.9109092169090824,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,2,0.9568103026547466,0,"signaling receptor activity"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.9143895476121282,0.06306905,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,2,0.9025609564036702,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.8696048369605115,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,2,0.8401965582431901,0.26407481,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.890337407140178,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9092895267009665,0.01849091,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9228090868021512,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9870662401598954,-0,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

# 
# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_MF_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998642200963308,4.927E-05,"extracellular region"),
                     c("GO:0005730","nucleolus",1.2114693153196516,1,0.6971843266754033,3.325E-05,"nucleolus"),
                     c("GO:0000793","condensed chromosome",0.4223503377863461,1,0.7422472078416456,0.41377726,"nucleolus"),
                     c("GO:0005634","nucleus",16.5161752456724,3,0.7134999513385691,0.2729216,"nucleolus"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8347985501739508,0.09158706,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,1,0.8691419232895021,0.17160779,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,4,0.9998131780536398,0.00015294,"plasma membrane"),
                     c("GO:0016020","membrane",49.2542153160787,5,0.9997497606099126,9.537E-05,"membrane"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999181436783204,0,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 

# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_CC_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,4,0.8558294592345888,0,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,2,0.7502799574733543,0.44154323,"monoatomic ion transport"),
                     c("GO:0035672","oligopeptide transmembrane transport",0.23299928216907387,1,0.8709360234047956,0.30173276,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.7380855241749573,-0,"defense response"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,1,0.8967939814685623,0.38212221,"defense response"),
                     c("GO:0009416","response to light stimulus",0.2589924319660237,1,0.8123843145513597,0.27535725,"defense response"),
                     c("GO:0010039","response to iron ion",0.05924742412469929,1,0.8065962760979498,0.37669922,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.7571506143619126,0.46097507,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.6737446503385388,0.37522238,"defense response"),
                     c("GO:0071230","cellular response to amino acid stimulus",0.039006977876590854,1,0.7761054891889361,0.232997,"defense response"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,1,0.9670960296176608,0.07105297,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9626264511905962,0.12277613,"lipid catabolic process"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9913546418996747,0.00723774,"chromosome condensation"),
                     c("GO:0045739","positive regulation of DNA repair",0.0330817425160876,1,0.9351177414386274,-0,"positive regulation of DNA repair"),
                     c("GO:0048316","seed development",0.03280076213709535,1,0.968402770812103,-0,"seed development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.968402770812103,0.42998645,"seed development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9078452347744271,-0,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9431701279846126,0.15925058,"DNA demethylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# # 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_BP_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,2,1,-0,"DNA-binding transcription factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.9468857898693089,0.06915977,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,1,0.9464323281816921,0.12712323,"transferase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8754614116022313,0.04503758,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.875440697222909,0.32372808,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8706072740378815,0.18924409,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.8508700477574961,0.1977699,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0016301","kinase activity",5.872540794447113,1,0.8253480268700466,0.40437782,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,1,0.9067805247798328,-0,"phosphatase activity"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.8004191068519949,0,"sequence-specific DNA binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.818817523666658,0.2421934,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,3,0.7681033053829566,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.8336257217278011,0.33073906,"sequence-specific DNA binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,2,0.8662630963641703,0.09221581,"metal ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,1,0.7616494213246765,0.38577556,"metal ion binding"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.92645908734501,0.10689728,"metal ion binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AllZmDEG_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, All Zm DE, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005634","nucleus",16.5161752456724,4,1,0,"nucleus"),
                     c("GO:0005737","cytoplasm",43.076660800468886,3,1,0.18211679,"nucleus"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AllZmDEG_NoDEGC3_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, All Zm DE, No DEG in C3, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.8780550412304704,0.08992384,"lipid metabolic process"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,1,-0,"muscle organ development"),
                     c("GO:0009834","plant-type secondary cell wall biogenesis",0.02117458277124934,1,0.9935816576963363,0.00689751,"plant-type secondary cell wall biogenesis"),
                     c("GO:0010047","fruit dehiscence",0.002082705440775892,1,1,-0,"fruit dehiscence"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,1,0.9446492347469846,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.45210426038827783,-0,"phosphatidylinositol dephosphorylation"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,1,0.5429517765871298,0.49264795,"phosphatidylinositol dephosphorylation"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.44722560840293274,0.46733961,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.7031332987012088,0.39551818,"phosphatidylinositol dephosphorylation"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,1,0,"mRNA transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AllZmDEG_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, All Zm DE, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,9,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,6,0.9556687776402477,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,31,0.9300196903212681,0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,2,0.9476843341207317,0.42104634,"mRNA binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,5,0.9193094385579699,0.31664803,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,2,0.9517726307655306,0.19519762,"mRNA binding"),
                     c("GO:0003684","damaged DNA binding",0.33181460177613226,1,0.9403499025777665,0.40218543,"mRNA binding"),
                     c("GO:0003697","single-stranded DNA binding",0.40260363864918647,1,0.9393076541150771,0.41085789,"mRNA binding"),
                     c("GO:0003713","transcription coactivator activity",0.24395086691346182,3,0.9693393828269257,0.29461398,"mRNA binding"),
                     c("GO:0003723","RNA binding",6.099813894661886,48,0.927900773543416,0.40176997,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9401925612517594,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.939848294598848,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,8,0.9172452462340493,0.24940265,"mRNA binding"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8881524585665113,0.3505363,"mRNA binding"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8893166604050521,0.33108637,"mRNA binding"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8885916556311917,0.3394079,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9468722449079794,0.428381,"mRNA binding"),
                     c("GO:0008327","methyl-CpG binding",0.009309036644208929,1,0.9390418609120744,0.46750139,"mRNA binding"),
                     c("GO:0008878","glucose-1-phosphate adenylyltransferase activity",0.021398809818155354,1,0.8890101964432305,0.29083106,"mRNA binding"),
                     c("GO:0016780","phosphotransferase activity, for other substituted phosphate groups",0.32920461222629094,1,0.8716406785187323,0.35102513,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9578412310689256,0.13060477,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9488267691429012,0.41077037,"mRNA binding"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8701342346445498,0.27195941,"mRNA binding"),
                     c("GO:0044183","protein folding chaperone",0.39473153762799984,1,0.9779221791394449,0.32443613,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,6,0.9957656723800591,-0,"structural constituent of ribosome"),
                     c("GO:0003774","cytoskeletal motor activity",0.3926027441124113,1,1,-0,"cytoskeletal motor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,13,1,-0,"catalytic activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,2,0.9847701438688098,-0,"GTPase activator activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,64,0.9676666389088371,0.07766851,"protein binding"),
                     c("GO:0005543","phospholipid binding",0.6999096105403309,2,0.9593664947520258,0.04714527,"phospholipid binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,9,0.9407825034733968,0.05738816,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,48,0.9174794025435898,0.48160605,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,36,0.9411278764904909,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,4,0.9515337555035566,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,9,0.9294139947809361,0.19635374,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9433231958232332,0.29697106,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9588095165721036,0.26177282,"zinc ion binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9751546291666777,0.2693746,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9706463507141544,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.9348882157731002,0.33990845,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,5,0.9495984971076281,0.18168089,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,1,0.9466720577769359,0.29755636,"zinc ion binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,1,0.9411153726895474,0.30272555,"zinc ion binding"),
                     c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9568488657554564,0.25647876,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,4,0.9740360747584395,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0015288","porin activity",0.12362746592455742,2,0.93884781391447,-0,"porin activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9130671375083721,0.39044024,"porin activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9417239126056423,0.47124298,"porin activity"),
                     c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,1,0.9411586434365969,0.31858441,"porin activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,1,0.9623999851096069,0.18410598,"porin activity"),
                     c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9564962400079116,0.04318065,"CoA-ligase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,56,0.9472493489854454,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,19,0.9511880784179745,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,39,0.9467167859017094,0.12712323,"transferase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,8,0.874378058566644,0.05667901,"glycosyltransferase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9146616131447959,0.20205602,"glycosyltransferase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,8,0.8279886835043281,0.33973871,"glycosyltransferase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,3,0.8786654680684459,0.2917808,"glycosyltransferase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,2,0.8785936626683389,0.28934569,"glycosyltransferase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8690028745887369,0.36085078,"glycosyltransferase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8844175700784569,0.30154932,"glycosyltransferase activity"),
                     c("GO:0016767","geranylgeranyl-diphosphate geranylgeranyltransferase activity",0.019117009268633918,1,0.9020886228737552,0.20439581,"glycosyltransferase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.8864348573497074,0.49223587,"glycosyltransferase activity"),
                     c("GO:0047216","inositol 3-alpha-galactosyltransferase activity",0.0009202596968429504,1,0.9192503600880888,0.43743078,"glycosyltransferase activity"),
                     c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.912186342655241,0.17480439,"glycosyltransferase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.8656573719390799,0.45431535,"glycosyltransferase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8670665535523243,0.39355689,"glycosyltransferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,4,0.9571217737769077,0.05997119,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,6,0.958900996678351,0.05646078,"isomerase activity"),
                     c("GO:0016871","cycloartenol synthase activity",0.0006563780006397912,1,0.9622379275230969,0.02586138,"cycloartenol synthase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,3,0.9576146728115937,0.0589874,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,17,0.8653089897370425,-0,"ATP hydrolysis activity"),
                     c("GO:0000146","microfilament motor activity",0.092387421083296,1,0.9701164243091338,0.49829267,"ATP hydrolysis activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,2,0.8780590553391917,0.42844795,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.915502130688966,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9398070340024738,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,2,0.9046890764331185,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0004712","protein serine/threonine/tyrosine kinase activity",0.056286631048107494,1,0.8424023554143275,0.44384691,"ATP hydrolysis activity"),
                     c("GO:0004715","non-membrane spanning protein tyrosine kinase activity",0.060526478133321286,1,0.8395108893642167,0.44607451,"ATP hydrolysis activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,3,0.8410457538068864,0.42141607,"ATP hydrolysis activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,5,0.818388215743473,0.48775459,"ATP hydrolysis activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,14,0.8383411349272583,0.3656957,"ATP hydrolysis activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,4,0.8511527842038522,0.45595338,"ATP hydrolysis activity"),
                     c("GO:0008568","microtubule severing ATPase activity",0.013493446398287598,1,0.8880141756101257,0.41923201,"ATP hydrolysis activity"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9607133938713633,0.45552111,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8995729857130226,0.39048386,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.9103623179383683,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0033919","glucan 1,3-alpha-glucosidase activity",0.000915824710352141,1,0.940637920264108,0.4867618,"ATP hydrolysis activity"),
                     c("GO:0052742","phosphatidylinositol kinase activity",0.07485148449863564,1,0.8728322439305342,0.40778693,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9359184516152639,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,5,0.8377396328998755,0.37283533,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9351685479926555,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9371411775817016,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9746923776739819,0.04901613,"carbohydrate binding"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.972489785549527,0.06283846,"protein-containing complex binding"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9721699882542946,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.979940750299162,0.49249984,"ER retention sequence binding"),
                     c("GO:0047769","arogenate dehydratase activity",0.02556769711951619,1,0.9482814083541394,0.03410747,"arogenate dehydratase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.954541698643944,0.36924486,"arogenate dehydratase activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,4,0.9154767838810479,0.0463004,"unfolded protein binding"),
                     c("GO:0000149","SNARE binding",0.26943873427614345,1,0.9206649941369229,0.40559978,"unfolded protein binding"),
                     c("GO:0003779","actin binding",0.7421262469463454,3,0.9093692418236538,0.44654025,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9207908819958428,0.40485813,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9296308010906706,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,2,0.9080501720806812,0.42444044,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,2,0.9254854729266847,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,2,0.9280593508043852,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9236883761991721,0.38789043,"unfolded protein binding"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9284681538408958,0.36031388,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,2,0.9477390078220179,0.25391209,"unfolded protein binding"),
                     c("GO:0032051","clathrin light chain binding",0.02013705616152008,1,0.9308049285888775,0.3284969,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,3,0.9178683476412992,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9166127381504338,0.42967902,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9220851670854019,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,3,0.9293047565349698,0.3555383,"unfolded protein binding"),
                     c("GO:0046982","protein heterodimerization activity",0.29517274338906496,1,0.917948559934934,0.40897571,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9103192038989392,0.47948026,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9459807677703748,0.26331899,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.936850692254367,0.31312228,"unfolded protein binding"),
                     c("GO:1990935","splicing factor binding",0.0007095978385295039,1,0.9458890993090746,0.26381104,"unfolded protein binding"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9171112021374623,0.048252,"dioxygenase activity"),
                     c("GO:0004322","ferroxidase activity",0.06266857660838222,1,0.931427241325847,0.3160356,"dioxygenase activity"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,2,0.9129073152886723,0.37878622,"dioxygenase activity"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9332768155056413,0.2837958,"dioxygenase activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9043323631708091,0.44913707,"dioxygenase activity"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9054338817098204,0.38257989,"dioxygenase activity"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9150426384166358,0.41045836,"dioxygenase activity"),
                     c("GO:0016651","oxidoreductase activity, acting on NAD(P)H",0.7914588191768638,1,0.9166500814288117,0.39829104,"dioxygenase activity"),
                     c("GO:0016671","oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor",0.10547063123118373,1,0.9237378380452589,0.33002637,"dioxygenase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}


dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,5,1,-0,"cellular_component"),
c("GO:0005576","extracellular region",3.8687535442060237,6,0.9999257567330232,6.017E-05,"extracellular region"),
c("GO:0005739","mitochondrion",4.856981674016684,30,0.7570150031465227,0,"mitochondrion"),
c("GO:0000138","Golgi trans cisterna",0.007541260384887046,1,0.8004846910996123,0.37607379,"mitochondrion"),
c("GO:0000325","plant-type vacuole",0.04066068696484073,5,0.800433056704888,0.15656526,"mitochondrion"),
c("GO:0000785","chromatin",1.087181392930179,3,0.6297811278249654,0.45756389,"mitochondrion"),
c("GO:0005634","nucleus",16.5161752456724,122,0.7288138569252248,0.35145602,"mitochondrion"),
c("GO:0005730","nucleolus",1.2114693153196516,16,0.6654899531811092,0.22801288,"mitochondrion"),
c("GO:0005737","cytoplasm",43.076660800468886,92,0.8926316389378861,0.12447513,"mitochondrion"),
c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,9,0.6351655235600465,0.22750481,"mitochondrion"),
c("GO:0005773","vacuole",1.35235298008496,6,0.7859504209325622,0.23143591,"mitochondrion"),
c("GO:0005777","peroxisome",0.7476308639759881,4,0.7897848457014061,0.21411812,"mitochondrion"),
c("GO:0005794","Golgi apparatus",2.349961096041566,21,0.7093335411531074,0.2503098,"mitochondrion"),
c("GO:0005811","lipid droplet",0.20804685033482948,3,0.7886697141064198,0.38610723,"mitochondrion"),
c("GO:0005829","cytosol",14.048080989011824,39,0.8614320568974416,0.19224382,"mitochondrion"),
c("GO:0005905","clathrin-coated pit",0.09808111076528711,2,0.8969327790575943,0.46642691,"mitochondrion"),
c("GO:0009507","chloroplast",0.6990388085931706,29,0.7199715797459495,0.21231674,"mitochondrion"),
c("GO:0009536","plastid",0.7624897559369845,25,0.7967142509721892,0.21465144,"mitochondrion"),
c("GO:0009579","thylakoid",0.26089530737804617,4,0.7848965332746932,0.39454155,"mitochondrion"),
c("GO:0017177","glucosidase II complex",0.015944734065838607,2,0.6873665946024858,0.39861008,"mitochondrion"),
c("GO:0031410","cytoplasmic vesicle",2.424891655569294,6,0.6983167248005325,0.25147478,"mitochondrion"),
c("GO:0031982","vesicle",2.6690048632308567,3,0.7950210667277179,0.2461316,"mitochondrion"),
c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,1,0.7857697422896927,0.4052793,"mitochondrion"),
c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7422829307765494,0.486126,"mitochondrion"),
c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9068606864093187,0.10673752,"mitochondrion"),
c("GO:0070971","endoplasmic reticulum exit site",0.06951625759076932,1,0.7773003449898026,0.45185985,"mitochondrion"),
c("GO:0009506","plasmodesma",0.10230521048664064,8,0.9772103571276574,3.812E-05,"plasmodesma"),
c("GO:0009524","phragmoplast",0.002879842104146322,1,0.9315447531842158,0.0711282,"phragmoplast"),
c("GO:0009986","cell surface",0.6578761991908514,1,0.9999382845144928,4.693E-05,"cell surface"),
c("GO:0012505","endomembrane system",6.893425119210306,1,0.9999203715136444,0.00011166,"endomembrane system"),
c("GO:0016020","membrane",49.2542153160787,115,0.9998928023305879,0.00010118,"membrane"),
c("GO:0032991","protein-containing complex",19.75554538569037,3,1,-0,"protein-containing complex"),
c("GO:0044297","cell body",0.22864057885869896,1,0.9999438437429754,4.148E-05,"cell body"),
c("GO:0048046","apoplast",0.0960361495472436,2,0.9999476834174484,3.787E-05,"apoplast"),
c("GO:0070469","respirasome",0.5984977809313304,2,0.9690463907323354,4.638E-05,"respirasome"),
c("GO:0090406","pollen tube",0.001702063711251277,3,0.9968157547408188,2.697E-05,"pollen tube"),
c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.9968236570771475,0.27222822,"pollen tube"),
c("GO:0098552","side of membrane",0.7241349304670945,1,0.9685386317153775,0.0715415,"side of membrane"),
c("GO:0005886","plasma membrane",17.177321395000487,32,0.9489593005345854,0.10744228,"side of membrane"),
c("GO:1990904","ribonucleoprotein complex",4.315842197772249,15,0.785047308750822,-0,"ribonucleoprotein complex"),
c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7454113519299168,0.47495898,"ribonucleoprotein complex"),
c("GO:0000178","exosome (RNase complex)",0.10685978389207654,1,0.756601954765338,0.34132399,"ribonucleoprotein complex"),
c("GO:0000502","proteasome complex",0.37652382533874423,4,0.7489665034469939,0.49612005,"ribonucleoprotein complex"),
c("GO:0000974","Prp19 complex",0.0635900941581645,1,0.8446460881447212,0.24779784,"ribonucleoprotein complex"),
c("GO:0005681","spliceosomal complex",0.7626239332222511,10,0.5997265732831089,0.3245659,"ribonucleoprotein complex"),
c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8070784917476633,0.22864326,"ribonucleoprotein complex"),
c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7915274515176007,0.255018,"ribonucleoprotein complex"),
c("GO:0005942","phosphatidylinositol 3-kinase complex",0.07251039886384637,1,0.7405014019971716,0.46964777,"ribonucleoprotein complex"),
c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.7954019921746794,0.43971777,"ribonucleoprotein complex"),
c("GO:0016272","prefoldin complex",0.03795726314317447,1,0.849760484954309,0.23619435,"ribonucleoprotein complex"),
c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,2,0.7757571761165165,0.42147812,"ribonucleoprotein complex"),
c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7368856835896797,0.49563216,"ribonucleoprotein complex"),
c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7937655158975829,0.41997104,"ribonucleoprotein complex"),
c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7736831980262837,0.41725271,"ribonucleoprotein complex"),
c("GO:0045252","oxoglutarate dehydrogenase complex",0.07266196950090673,1,0.7234534952614657,0.46971574,"ribonucleoprotein complex"),
c("GO:0046930","pore complex",0.06221601936645362,2,0.7953428861431854,0.24728353,"ribonucleoprotein complex"),
c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7017571486001692,0.46572404,"ribonucleoprotein complex"),
c("GO:0071439","clathrin complex",0.011926869801468804,1,0.7618966594313018,0.36298129,"ribonucleoprotein complex"),
c("GO:0071819","DUBm complex",0.008179844872174023,1,0.71553362488728,0.43933875,"ribonucleoprotein complex"),
c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,6,0.7708417521135674,0.25161359,"ribonucleoprotein complex"),
c("GO:0089701","U2AF complex",0.020742814396387827,1,0.6995805439409414,0.46994592,"ribonucleoprotein complex"),
c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.78791741550791,0.40438888,"ribonucleoprotein complex"),
c("GO:0160064","multi-pass translocon complex",0.0006336149582030303,1,0.8428987684788185,0.17222646,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,5,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,6,0.9999257567330232,6.017E-05,"extracellular region"),
                     c("GO:0005739","mitochondrion",4.856981674016684,30,0.7570150031465227,0,"mitochondrion"),
                     c("GO:0000138","Golgi trans cisterna",0.007541260384887046,1,0.8004846910996123,0.37607379,"mitochondrion"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,5,0.800433056704888,0.15656526,"mitochondrion"),
                     c("GO:0000785","chromatin",1.087181392930179,3,0.6297811278249654,0.45756389,"mitochondrion"),
                     c("GO:0005634","nucleus",16.5161752456724,122,0.7288138569252248,0.35145602,"mitochondrion"),
                     c("GO:0005730","nucleolus",1.2114693153196516,16,0.6654899531811092,0.22801288,"mitochondrion"),
                     c("GO:0005737","cytoplasm",43.076660800468886,92,0.8926316389378861,0.12447513,"mitochondrion"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,9,0.6351655235600465,0.22750481,"mitochondrion"),
                     c("GO:0005773","vacuole",1.35235298008496,6,0.7859504209325622,0.23143591,"mitochondrion"),
                     c("GO:0005777","peroxisome",0.7476308639759881,4,0.7897848457014061,0.21411812,"mitochondrion"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,21,0.7093335411531074,0.2503098,"mitochondrion"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,3,0.7886697141064198,0.38610723,"mitochondrion"),
                     c("GO:0005829","cytosol",14.048080989011824,39,0.8614320568974416,0.19224382,"mitochondrion"),
                     c("GO:0005905","clathrin-coated pit",0.09808111076528711,2,0.8969327790575943,0.46642691,"mitochondrion"),
                     c("GO:0009507","chloroplast",0.6990388085931706,29,0.7199715797459495,0.21231674,"mitochondrion"),
                     c("GO:0009536","plastid",0.7624897559369845,25,0.7967142509721892,0.21465144,"mitochondrion"),
                     c("GO:0009579","thylakoid",0.26089530737804617,4,0.7848965332746932,0.39454155,"mitochondrion"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,2,0.6873665946024858,0.39861008,"mitochondrion"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,6,0.6983167248005325,0.25147478,"mitochondrion"),
                     c("GO:0031982","vesicle",2.6690048632308567,3,0.7950210667277179,0.2461316,"mitochondrion"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,1,0.7857697422896927,0.4052793,"mitochondrion"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7422829307765494,0.486126,"mitochondrion"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9068606864093187,0.10673752,"mitochondrion"),
                     c("GO:0070971","endoplasmic reticulum exit site",0.06951625759076932,1,0.7773003449898026,0.45185985,"mitochondrion"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,8,0.9772103571276574,3.812E-05,"plasmodesma"),
                     c("GO:0009524","phragmoplast",0.002879842104146322,1,0.9315447531842158,0.0711282,"phragmoplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999382845144928,4.693E-05,"cell surface"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.9999203715136444,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,115,0.9998928023305879,0.00010118,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,3,1,-0,"protein-containing complex"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999438437429754,4.148E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9999476834174484,3.787E-05,"apoplast"),
                     c("GO:0070469","respirasome",0.5984977809313304,2,0.9690463907323354,4.638E-05,"respirasome"),
                     c("GO:0090406","pollen tube",0.001702063711251277,3,0.9968157547408188,2.697E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.9968236570771475,0.27222822,"pollen tube"),
                     c("GO:0098552","side of membrane",0.7241349304670945,1,0.9685386317153775,0.0715415,"side of membrane"),
                     c("GO:0005886","plasma membrane",17.177321395000487,32,0.9489593005345854,0.10744228,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,15,0.785047308750822,-0,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7454113519299168,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000178","exosome (RNase complex)",0.10685978389207654,1,0.756601954765338,0.34132399,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,4,0.7489665034469939,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0000974","Prp19 complex",0.0635900941581645,1,0.8446460881447212,0.24779784,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,10,0.5997265732831089,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8070784917476633,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7915274515176007,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005942","phosphatidylinositol 3-kinase complex",0.07251039886384637,1,0.7405014019971716,0.46964777,"ribonucleoprotein complex"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.7954019921746794,0.43971777,"ribonucleoprotein complex"),
                     c("GO:0016272","prefoldin complex",0.03795726314317447,1,0.849760484954309,0.23619435,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,2,0.7757571761165165,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7368856835896797,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7937655158975829,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7736831980262837,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0045252","oxoglutarate dehydrogenase complex",0.07266196950090673,1,0.7234534952614657,0.46971574,"ribonucleoprotein complex"),
                     c("GO:0046930","pore complex",0.06221601936645362,2,0.7953428861431854,0.24728353,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7017571486001692,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071439","clathrin complex",0.011926869801468804,1,0.7618966594313018,0.36298129,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.71553362488728,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,6,0.7708417521135674,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.6995805439409414,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.78791741550791,0.40438888,"ribonucleoprotein complex"),
                     c("GO:0160064","multi-pass translocon complex",0.0006336149582030303,1,0.8428987684788185,0.17222646,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",0.9427113541805001,1,1,-0,"immune system process"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,8,0.8866644093239825,0.01352937,"chromatin organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9170930274600421,0.39503669,"chromatin organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9090306703822497,0.38546222,"chromatin organization"),
                     c("GO:0009658","chloroplast organization",0.0906284959258338,2,0.9171762363863234,0.32420414,"chromatin organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9487168285546999,0.21045842,"chromatin organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9114009532857971,0.45272646,"chromatin organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,3,0.8947949768582132,0.36480707,"chromatin organization"),
                     c("GO:0043622","cortical microtubule organization",0.008850881938255893,2,0.9250852165518301,0.26893696,"chromatin organization"),
                     c("GO:0051289","protein homotetramerization",0.01763275115184702,1,0.9066141946465268,0.48332609,"chromatin organization"),
                     c("GO:0061025","membrane fusion",0.31015797308444587,1,0.9134725287105906,0.36373751,"chromatin organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9062544766699169,0.48565697,"chromatin organization"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,3,0.9055404150789912,0.4063955,"chromatin organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9432061285641596,0.20937319,"chromatin organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9384183206495306,0.22611843,"chromatin organization"),
                     c("GO:0006397","mRNA processing",1.242261085587905,17,0.8085872367994686,0.01340109,"mRNA processing"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,2,0.8673451110596677,0.4354711,"mRNA processing"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8511791646442186,0.45078305,"mRNA processing"),
                     c("GO:0006012","galactose metabolic process",0.12319264300693568,1,0.9151486398359728,0.38165312,"mRNA processing"),
                     c("GO:0006412","translation",4.38869169324396,11,0.795526755579601,0.48035482,"mRNA processing"),
                     c("GO:0006457","protein folding",1.174377211919444,7,0.8399956489397382,0.40866092,"mRNA processing"),
                     c("GO:0006491","N-glycan processing",0.05083033645576476,2,0.9014579330039085,0.27495136,"mRNA processing"),
                     c("GO:0008652","amino acid biosynthetic process",2.679426429329935,4,0.8274795977949865,0.49906638,"mRNA processing"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8754694014898333,0.367811,"mRNA processing"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9347135507485762,0.32178408,"mRNA processing"),
                     c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8565169445990376,0.21599583,"mRNA processing"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,14,0.8678840100509494,0.19691418,"mRNA processing"),
                     c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,1,0.8562525213532106,0.35937708,"mRNA processing"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9285053680018841,0.33268781,"mRNA processing"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9095752142783583,0.27916289,"mRNA processing"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8847586056951727,0.25012422,"mRNA processing"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,1,0.8519189150891596,0.43903701,"mRNA processing"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8977734452569671,0.21843208,"mRNA processing"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.903313014554542,0.40604291,"mRNA processing"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8460257228812044,0.49920408,"mRNA processing"),
                     c("GO:0140040","mitochondrial polycistronic RNA processing",3.4506362332381646E-05,1,0.9003799201379421,0.41070687,"mRNA processing"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.804790687542519,0.46167258,"mRNA processing"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,2,0.863956424377826,0.33620138,"mRNA processing"),
                     c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,1,0.8640716719978603,0.41679514,"mRNA processing"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9519203852425954,0.10630662,"mRNA processing"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9269586764061225,0.17943885,"mRNA processing"),
                     c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,2,0.9906419939902001,-0,"intracellular iron ion homeostasis"),
                     c("GO:0006915","apoptotic process",0.3972840732335429,1,0.9865014295873916,0.01170435,"apoptotic process"),
                     c("GO:0007049","cell cycle",2.8073119376140743,2,0.988382107485483,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.982669832660182,0.01255057,"chromosome segregation"),
                     c("GO:0000919","cell plate assembly",0.001451731958126628,1,0.9222325983880028,0.47630759,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,2,0.9964178461450757,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,4,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,6,1,-0,"metabolic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,4,0.9505650520180743,0.09593543,"biosynthetic process"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,2,0.9467146471733705,0.23095606,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9186099655298551,0.16846091,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9360747008128004,0.17227516,"biosynthetic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9538084716959496,0.07301458,"toxin metabolic process"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,7,0.9113692898218173,-0,"response to salt stress"),
                     c("GO:0006283","transcription-coupled nucleotide-excision repair",0.06059070751549558,1,0.8230631407135823,0.33995763,"response to salt stress"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,4,0.9155380914298296,0.41360326,"response to salt stress"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.9253478567570326,0.44363439,"response to salt stress"),
                     c("GO:0009737","response to abscisic acid",0.05830589338105859,4,0.9014287750437808,0.1955272,"response to salt stress"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.9011652813787087,0.17275947,"response to salt stress"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9385510663190229,0.27664386,"response to salt stress"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8661801521465045,0.32283054,"response to salt stress"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8796247927802366,0.30241174,"response to salt stress"),
                     c("GO:0034059","response to anoxia",0.0003401341429906191,1,0.9333167070151557,0.45207873,"response to salt stress"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8431627364187204,0.46171648,"response to salt stress"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9234929233589982,0.27995591,"response to salt stress"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9839300014980518,0.48813486,"response to salt stress"),
                     c("GO:0051607","defense response to virus",0.17287441054506547,3,0.9056042757659867,0.44531989,"response to salt stress"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9437897153287094,0.12932521,"response to salt stress"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,3,0.8906728355846558,0.46817928,"response to salt stress"),
                     c("GO:1901562","response to paraquat",0.00022182661499388205,1,0.9360756072227834,0.38531986,"response to salt stress"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.9361097246591891,0.31787134,"response to salt stress"),
                     c("GO:0009853","photorespiration",0.014256057123606818,2,0.9594272139884765,0.06301944,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,10,0.9032255220494619,-0,"flower development"),
                     c("GO:0001824","blastocyst development",0.009126932836914946,1,0.9344416855485165,0.41203625,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,5,0.9288625525967364,0.4203896,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,2,0.9425335743316722,0.48835672,"flower development"),
                     c("GO:0009944","polarity specification of adaxial/abaxial axis",0.0024376280247661035,2,0.939214395325934,0.38300054,"flower development"),
                     c("GO:0010091","trichome branching",0.0006408324433156592,1,0.9407692954396792,0.31954118,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9232166220374267,0.29044182,"flower development"),
                     c("GO:0021987","cerebral cortex development",0.017519373104183483,1,0.9321073431898245,0.44335853,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.9383818769417046,0.39539149,"flower development"),
                     c("GO:0061137","bud dilation",9.612486649734888E-05,1,0.9441900191532959,0.29194579,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9492885830919227,0.28605225,"flower development"),
                     c("GO:1905177","tracheary element differentiation",0.0002661919379926584,1,0.9484684178267991,0.2910313,"flower development"),
                     c("GO:0015031","protein transport",3.093438694074183,21,0.8786053956854146,0,"protein transport"),
                     c("GO:0000055","ribosomal large subunit export from nucleus",0.042115015226671805,1,0.8427177124621217,0.47380685,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,9,0.9267581658958487,0.42145531,"protein transport"),
                     c("GO:0006833","water transport",0.07339010320064257,1,0.9491204356525017,0.25687781,"protein transport"),
                     c("GO:0006897","endocytosis",0.6399968963991822,3,0.9206566561305087,0.32211657,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9500934210874007,0.21356797,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9420683543832767,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,10,0.9218652825945151,0.38566904,"protein transport"),
                     c("GO:0045037","protein import into chloroplast stroma",0.006780500198312994,1,0.9169274294077039,0.49429905,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9266755636689286,0.49534489,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,3,0.9090932213665075,0.35175462,"protein transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1,0.9239553918946453,0.36406786,"protein transport"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,4,0.8711225736488413,0.0993369,"lipid catabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,4,0.9350954547555893,0.15163831,"lipid catabolic process"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.9127011287290062,0.49989178,"lipid catabolic process"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9350535605849191,0.15178462,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.9335782313079345,0.12277613,"lipid catabolic process"),
                     c("GO:0006914","autophagy",0.44134623319182764,2,0.9009225509892886,0.43293744,"lipid catabolic process"),
                     c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.863787864849355,0.44469555,"lipid catabolic process"),
                     c("GO:0033611","oxalate catabolic process",0.006792823899145988,1,0.9006615297708536,0.3986442,"lipid catabolic process"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,1,0.9063676543193916,0.45718665,"lipid catabolic process"),
                     c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1,0.9241632086625796,0.27800701,"lipid catabolic process"),
                     c("GO:0071040","nuclear polyadenylation-dependent antisense transcript catabolic process",0.008039982423444923,1,0.846140933066354,0.40606321,"lipid catabolic process"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8712411897571667,0.33283067,"lipid catabolic process"),
                     c("GO:0022900","electron transport chain",0.8207017864535318,3,0.9275162410801076,0.09080927,"electron transport chain"),
                     c("GO:0030048","actin filament-based movement",0.0797343443894676,1,0.9809486708225112,0.00993277,"actin filament-based movement"),
                     c("GO:0032259","methylation",2.6278542060840238,8,0.964081561638524,0.05828286,"methylation"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,8,0.8874822614448646,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.00031548674132463224,1,0.9431985304868815,0.44967237,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.937240340085944,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,5,0.9042985549496313,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.8694391801668475,0.48392513,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9269597144511158,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0009852","auxin catabolic process",3.4506362332381646E-05,1,0.9134896317867919,0.37814571,"positive regulation of DNA-templated transcription"),
                     c("GO:0009926","auxin polar transport",0.004261535748049133,1,0.9090351724516167,0.4567727,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9449242648127614,0.40013642,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,1,0.9527225817305781,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9449787012194721,0.2722228,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9146725916713815,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0032784","regulation of DNA-templated transcription elongation",0.17948484367188314,1,0.9268571465544592,0.37447493,"positive regulation of DNA-templated transcription"),
                     c("GO:0032876","negative regulation of DNA endoreduplication",0.0004264000488215732,1,0.9357402602755096,0.48325207,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9407954170213159,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042659","regulation of cell fate specification",0.0026569898995933866,1,0.9395438038696537,0.4268527,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,1,0.9498742269467771,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9388218162506468,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9369179017178604,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,6,0.8935182768487635,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,3,0.9307484475450518,0.14815537,"positive regulation of DNA-templated transcription"),
                     c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,2,0.939597542154595,0.42385341,"positive regulation of DNA-templated transcription"),
                     c("GO:0048586","regulation of long-day photoperiodism, flowering",0.0036552096670658557,3,0.9331763297600029,0.13376852,"positive regulation of DNA-templated transcription"),
                     c("GO:0050821","protein stabilization",0.10143638155636905,1,0.941029090531027,0.36692594,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9368560642337764,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0060195","negative regulation of antisense RNA transcription",0.00012077226816333578,1,0.9407265990034013,0.2205473,"positive regulation of DNA-templated transcription"),
                     c("GO:0061635","regulation of protein complex stability",0.0017622892191180627,1,0.9537144360669078,0.1205724,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9361167435317312,0.33352541,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.9535274617843031,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:1900057","positive regulation of leaf senescence",0.0003204162216578296,1,0.9325822322013414,0.38480392,"positive regulation of DNA-templated transcription"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,1,0.9280254393514106,0.47223923,"positive regulation of DNA-templated transcription"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9339593553083878,0.31791635,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9413740422442916,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:1901000","regulation of response to salt stress",0.0019915100546117406,2,0.9394272056531792,0.47172255,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9287448698467682,0.454995,"positive regulation of DNA-templated transcription"),
                     c("GO:1903329","regulation of iron-sulfur cluster assembly",0.0002464740166598689,1,0.9566642416934493,0.12070797,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,2,0.9417079737045797,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,1,0.9396185448834918,0.35770111,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,2,0.940687055674097,0.33609152,"positive regulation of DNA-templated transcription"),
                     c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,1,0.9409886171697722,0.2795365,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,3,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,1,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,2,0.9890436581646581,0.01381155,"cell division"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])


hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}



dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",0.9427113541805001,1,1,-0,"immune system process"),
c("GO:0006325","chromatin organization",1.3384303174082528,8,0.8866644093239825,0.01352937,"chromatin organization"),
c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9170930274600421,0.39503669,"chromatin organization"),
c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9090306703822497,0.38546222,"chromatin organization"),
c("GO:0009658","chloroplast organization",0.0906284959258338,2,0.9171762363863234,0.32420414,"chromatin organization"),
c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9487168285546999,0.21045842,"chromatin organization"),
c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9114009532857971,0.45272646,"chromatin organization"),
c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,3,0.8947949768582132,0.36480707,"chromatin organization"),
c("GO:0043622","cortical microtubule organization",0.008850881938255893,2,0.9250852165518301,0.26893696,"chromatin organization"),
c("GO:0051289","protein homotetramerization",0.01763275115184702,1,0.9066141946465268,0.48332609,"chromatin organization"),
c("GO:0061025","membrane fusion",0.31015797308444587,1,0.9134725287105906,0.36373751,"chromatin organization"),
c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9062544766699169,0.48565697,"chromatin organization"),
c("GO:0071555","cell wall organization",0.8943925879544993,3,0.9055404150789912,0.4063955,"chromatin organization"),
c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9432061285641596,0.20937319,"chromatin organization"),
c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9384183206495306,0.22611843,"chromatin organization"),
c("GO:0006397","mRNA processing",1.242261085587905,17,0.8085872367994686,0.01340109,"mRNA processing"),
c("GO:0000373","Group II intron splicing",0.009405448475740597,2,0.8673451110596677,0.4354711,"mRNA processing"),
c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8511791646442186,0.45078305,"mRNA processing"),
c("GO:0006012","galactose metabolic process",0.12319264300693568,1,0.9151486398359728,0.38165312,"mRNA processing"),
c("GO:0006412","translation",4.38869169324396,11,0.795526755579601,0.48035482,"mRNA processing"),
c("GO:0006457","protein folding",1.174377211919444,7,0.8399956489397382,0.40866092,"mRNA processing"),
c("GO:0006491","N-glycan processing",0.05083033645576476,2,0.9014579330039085,0.27495136,"mRNA processing"),
c("GO:0008652","amino acid biosynthetic process",2.679426429329935,4,0.8274795977949865,0.49906638,"mRNA processing"),
c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8754694014898333,0.367811,"mRNA processing"),
c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9347135507485762,0.32178408,"mRNA processing"),
c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8565169445990376,0.21599583,"mRNA processing"),
c("GO:0016567","protein ubiquitination",1.2678648064385323,14,0.8678840100509494,0.19691418,"mRNA processing"),
c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,1,0.8562525213532106,0.35937708,"mRNA processing"),
c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9285053680018841,0.33268781,"mRNA processing"),
c("GO:0030091","protein repair",0.06087415263465443,1,0.9095752142783583,0.27916289,"mRNA processing"),
c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8847586056951727,0.25012422,"mRNA processing"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,1,0.8519189150891596,0.43903701,"mRNA processing"),
c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8977734452569671,0.21843208,"mRNA processing"),
c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.903313014554542,0.40604291,"mRNA processing"),
c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8460257228812044,0.49920408,"mRNA processing"),
c("GO:0140040","mitochondrial polycistronic RNA processing",3.4506362332381646E-05,1,0.9003799201379421,0.41070687,"mRNA processing"),
c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.804790687542519,0.46167258,"mRNA processing"),
c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,2,0.863956424377826,0.33620138,"mRNA processing"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,1,0.8640716719978603,0.41679514,"mRNA processing"),
c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9519203852425954,0.10630662,"mRNA processing"),
c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9269586764061225,0.17943885,"mRNA processing"),
c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,2,0.9906419939902001,-0,"intracellular iron ion homeostasis"),
c("GO:0006915","apoptotic process",0.3972840732335429,1,0.9865014295873916,0.01170435,"apoptotic process"),
c("GO:0007049","cell cycle",2.8073119376140743,2,0.988382107485483,0.01495111,"cell cycle"),
c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.982669832660182,0.01255057,"chromosome segregation"),
c("GO:0000919","cell plate assembly",0.001451731958126628,1,0.9222325983880028,0.47630759,"chromosome segregation"),
c("GO:0007623","circadian rhythm",0.10049731555289494,2,0.9964178461450757,-0,"circadian rhythm"),
c("GO:0008150","biological_process",100,4,1,-0,"biological_process"),
c("GO:0008152","metabolic process",57.597931274565454,6,1,-0,"metabolic process"),
c("GO:0009058","biosynthetic process",29.004480601854056,4,0.9505650520180743,0.09593543,"biosynthetic process"),
c("GO:0044238","primary metabolic process",45.47569830278978,2,0.9467146471733705,0.23095606,"biosynthetic process"),
c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9186099655298551,0.16846091,"biosynthetic process"),
c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9360747008128004,0.17227516,"biosynthetic process"),
c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9538084716959496,0.07301458,"toxin metabolic process"),
c("GO:0009651","response to salt stress",0.07343446852364134,7,0.9113692898218173,-0,"response to salt stress"),
c("GO:0006283","transcription-coupled nucleotide-excision repair",0.06059070751549558,1,0.8230631407135823,0.33995763,"response to salt stress"),
c("GO:0006979","response to oxidative stress",0.8191810417707402,4,0.9155380914298296,0.41360326,"response to salt stress"),
c("GO:0009611","response to wounding",0.16569462243976346,1,0.9253478567570326,0.44363439,"response to salt stress"),
c("GO:0009737","response to abscisic acid",0.05830589338105859,4,0.9014287750437808,0.1955272,"response to salt stress"),
c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.9011652813787087,0.17275947,"response to salt stress"),
c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9385510663190229,0.27664386,"response to salt stress"),
c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8661801521465045,0.32283054,"response to salt stress"),
c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8796247927802366,0.30241174,"response to salt stress"),
c("GO:0034059","response to anoxia",0.0003401341429906191,1,0.9333167070151557,0.45207873,"response to salt stress"),
c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8431627364187204,0.46171648,"response to salt stress"),
c("GO:0042221","response to chemical",4.8560360057130705,1,0.9234929233589982,0.27995591,"response to salt stress"),
c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9839300014980518,0.48813486,"response to salt stress"),
c("GO:0051607","defense response to virus",0.17287441054506547,3,0.9056042757659867,0.44531989,"response to salt stress"),
c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9437897153287094,0.12932521,"response to salt stress"),
c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,3,0.8906728355846558,0.46817928,"response to salt stress"),
c("GO:1901562","response to paraquat",0.00022182661499388205,1,0.9360756072227834,0.38531986,"response to salt stress"),
c("GO:1902074","response to salt",0.002343967898435353,3,0.9361097246591891,0.31787134,"response to salt stress"),
c("GO:0009853","photorespiration",0.014256057123606818,2,0.9594272139884765,0.06301944,"photorespiration"),
c("GO:0009908","flower development",0.02997124042584006,10,0.9032255220494619,-0,"flower development"),
c("GO:0001824","blastocyst development",0.009126932836914946,1,0.9344416855485165,0.41203625,"flower development"),
c("GO:0009555","pollen development",0.012900450031977538,5,0.9288625525967364,0.4203896,"flower development"),
c("GO:0009846","pollen germination",0.003164726373912717,2,0.9425335743316722,0.48835672,"flower development"),
c("GO:0009944","polarity specification of adaxial/abaxial axis",0.0024376280247661035,2,0.939214395325934,0.38300054,"flower development"),
c("GO:0010091","trichome branching",0.0006408324433156592,1,0.9407692954396792,0.31954118,"flower development"),
c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9232166220374267,0.29044182,"flower development"),
c("GO:0021987","cerebral cortex development",0.017519373104183483,1,0.9321073431898245,0.44335853,"flower development"),
c("GO:0042335","cuticle development",0.004384772756379067,1,0.9383818769417046,0.39539149,"flower development"),
c("GO:0061137","bud dilation",9.612486649734888E-05,1,0.9441900191532959,0.29194579,"flower development"),
c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9492885830919227,0.28605225,"flower development"),
c("GO:1905177","tracheary element differentiation",0.0002661919379926584,1,0.9484684178267991,0.2910313,"flower development"),
c("GO:0015031","protein transport",3.093438694074183,21,0.8786053956854146,0,"protein transport"),
c("GO:0000055","ribosomal large subunit export from nucleus",0.042115015226671805,1,0.8427177124621217,0.47380685,"protein transport"),
c("GO:0006811","monoatomic ion transport",4.7761710300947735,9,0.9267581658958487,0.42145531,"protein transport"),
c("GO:0006833","water transport",0.07339010320064257,1,0.9491204356525017,0.25687781,"protein transport"),
c("GO:0006897","endocytosis",0.6399968963991822,3,0.9206566561305087,0.32211657,"protein transport"),
c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9500934210874007,0.21356797,"protein transport"),
c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9420683543832767,0.25083478,"protein transport"),
c("GO:0016192","vesicle-mediated transport",2.6087919056355493,10,0.9218652825945151,0.38566904,"protein transport"),
c("GO:0045037","protein import into chloroplast stroma",0.006780500198312994,1,0.9169274294077039,0.49429905,"protein transport"),
c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9266755636689286,0.49534489,"protein transport"),
c("GO:1902600","proton transmembrane transport",1.312853708699458,3,0.9090932213665075,0.35175462,"protein transport"),
c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1,0.9239553918946453,0.36406786,"protein transport"),
c("GO:0016042","lipid catabolic process",1.4093137798594644,4,0.8711225736488413,0.0993369,"lipid catabolic process"),
c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,4,0.9350954547555893,0.15163831,"lipid catabolic process"),
c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.9127011287290062,0.49989178,"lipid catabolic process"),
c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9350535605849191,0.15178462,"lipid catabolic process"),
c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.9335782313079345,0.12277613,"lipid catabolic process"),
c("GO:0006914","autophagy",0.44134623319182764,2,0.9009225509892886,0.43293744,"lipid catabolic process"),
c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.863787864849355,0.44469555,"lipid catabolic process"),
c("GO:0033611","oxalate catabolic process",0.006792823899145988,1,0.9006615297708536,0.3986442,"lipid catabolic process"),
c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,1,0.9063676543193916,0.45718665,"lipid catabolic process"),
c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1,0.9241632086625796,0.27800701,"lipid catabolic process"),
c("GO:0071040","nuclear polyadenylation-dependent antisense transcript catabolic process",0.008039982423444923,1,0.846140933066354,0.40606321,"lipid catabolic process"),
c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8712411897571667,0.33283067,"lipid catabolic process"),
c("GO:0022900","electron transport chain",0.8207017864535318,3,0.9275162410801076,0.09080927,"electron transport chain"),
c("GO:0030048","actin filament-based movement",0.0797343443894676,1,0.9809486708225112,0.00993277,"actin filament-based movement"),
c("GO:0032259","methylation",2.6278542060840238,8,0.964081561638524,0.05828286,"methylation"),
c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,8,0.8874822614448646,-0,"positive regulation of DNA-templated transcription"),
c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.00031548674132463224,1,0.9431985304868815,0.44967237,"positive regulation of DNA-templated transcription"),
c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.937240340085944,0.27139997,"positive regulation of DNA-templated transcription"),
c("GO:0006417","regulation of translation",1.3923070727099334,5,0.9042985549496313,0.40629436,"positive regulation of DNA-templated transcription"),
c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.8694391801668475,0.48392513,"positive regulation of DNA-templated transcription"),
c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9269597144511158,0.38387397,"positive regulation of DNA-templated transcription"),
c("GO:0009852","auxin catabolic process",3.4506362332381646E-05,1,0.9134896317867919,0.37814571,"positive regulation of DNA-templated transcription"),
c("GO:0009926","auxin polar transport",0.004261535748049133,1,0.9090351724516167,0.4567727,"positive regulation of DNA-templated transcription"),
c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9449242648127614,0.40013642,"positive regulation of DNA-templated transcription"),
c("GO:0010119","regulation of stomatal movement",0.016952482865865783,1,0.9527225817305781,0.16052978,"positive regulation of DNA-templated transcription"),
c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9449787012194721,0.2722228,"positive regulation of DNA-templated transcription"),
c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9146725916713815,0.41251019,"positive regulation of DNA-templated transcription"),
c("GO:0032784","regulation of DNA-templated transcription elongation",0.17948484367188314,1,0.9268571465544592,0.37447493,"positive regulation of DNA-templated transcription"),
c("GO:0032876","negative regulation of DNA endoreduplication",0.0004264000488215732,1,0.9357402602755096,0.48325207,"positive regulation of DNA-templated transcription"),
c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9407954170213159,0.21405542,"positive regulation of DNA-templated transcription"),
c("GO:0042659","regulation of cell fate specification",0.0026569898995933866,1,0.9395438038696537,0.4268527,"positive regulation of DNA-templated transcription"),
c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,1,0.9498742269467771,0.17055333,"positive regulation of DNA-templated transcription"),
c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9388218162506468,0.22344122,"positive regulation of DNA-templated transcription"),
c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9369179017178604,0.47414277,"positive regulation of DNA-templated transcription"),
c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,6,0.8935182768487635,0.47074713,"positive regulation of DNA-templated transcription"),
c("GO:0045995","regulation of embryonic development",0.01473914619626016,3,0.9307484475450518,0.14815537,"positive regulation of DNA-templated transcription"),
c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,2,0.939597542154595,0.42385341,"positive regulation of DNA-templated transcription"),
c("GO:0048586","regulation of long-day photoperiodism, flowering",0.0036552096670658557,3,0.9331763297600029,0.13376852,"positive regulation of DNA-templated transcription"),
c("GO:0050821","protein stabilization",0.10143638155636905,1,0.941029090531027,0.36692594,"positive regulation of DNA-templated transcription"),
c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9368560642337764,0.23293868,"positive regulation of DNA-templated transcription"),
c("GO:0060195","negative regulation of antisense RNA transcription",0.00012077226816333578,1,0.9407265990034013,0.2205473,"positive regulation of DNA-templated transcription"),
c("GO:0061635","regulation of protein complex stability",0.0017622892191180627,1,0.9537144360669078,0.1205724,"positive regulation of DNA-templated transcription"),
c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9361167435317312,0.33352541,"positive regulation of DNA-templated transcription"),
c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.9535274617843031,0.15711534,"positive regulation of DNA-templated transcription"),
c("GO:1900057","positive regulation of leaf senescence",0.0003204162216578296,1,0.9325822322013414,0.38480392,"positive regulation of DNA-templated transcription"),
c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,1,0.9280254393514106,0.47223923,"positive regulation of DNA-templated transcription"),
c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9339593553083878,0.31791635,"positive regulation of DNA-templated transcription"),
c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9413740422442916,0.42373749,"positive regulation of DNA-templated transcription"),
c("GO:1901000","regulation of response to salt stress",0.0019915100546117406,2,0.9394272056531792,0.47172255,"positive regulation of DNA-templated transcription"),
c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9287448698467682,0.454995,"positive regulation of DNA-templated transcription"),
c("GO:1903329","regulation of iron-sulfur cluster assembly",0.0002464740166598689,1,0.9566642416934493,0.12070797,"positive regulation of DNA-templated transcription"),
c("GO:2000024","regulation of leaf development",0.0023045320557697744,2,0.9417079737045797,0.42373749,"positive regulation of DNA-templated transcription"),
c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,1,0.9396185448834918,0.35770111,"positive regulation of DNA-templated transcription"),
c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,2,0.940687055674097,0.33609152,"positive regulation of DNA-templated transcription"),
c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,1,0.9409886171697722,0.2795365,"positive regulation of DNA-templated transcription"),
c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
c("GO:0050896","response to stimulus",17.567785530535815,3,1,-0,"response to stimulus"),
c("GO:0051179","localization",19.75810399172557,1,1,-0,"localization"),
c("GO:0051301","cell division",1.5693197819947182,2,0.9890436581646581,0.01381155,"cell division"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,9,1,-0,"molecular_function"),
c("GO:0003682","chromatin binding",0.6287702097345026,6,0.9556687776402477,0.04661592,"chromatin binding"),
c("GO:0003729","mRNA binding",1.136762432364793,31,0.9300196903212681,0,"mRNA binding"),
c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,2,0.9476843341207317,0.42104634,"mRNA binding"),
c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,5,0.9193094385579699,0.31664803,"mRNA binding"),
c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,2,0.9517726307655306,0.19519762,"mRNA binding"),
c("GO:0003684","damaged DNA binding",0.33181460177613226,1,0.9403499025777665,0.40218543,"mRNA binding"),
c("GO:0003697","single-stranded DNA binding",0.40260363864918647,1,0.9393076541150771,0.41085789,"mRNA binding"),
c("GO:0003713","transcription coactivator activity",0.24395086691346182,3,0.9693393828269257,0.29461398,"mRNA binding"),
c("GO:0003723","RNA binding",6.099813894661886,48,0.927900773543416,0.40176997,"mRNA binding"),
c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9401925612517594,0.48968021,"mRNA binding"),
c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.939848294598848,0.47960958,"mRNA binding"),
c("GO:0003743","translation initiation factor activity",0.37315089336372126,8,0.9172452462340493,0.24940265,"mRNA binding"),
c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8881524585665113,0.3505363,"mRNA binding"),
c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8893166604050521,0.33108637,"mRNA binding"),
c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8885916556311917,0.3394079,"mRNA binding"),
c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9468722449079794,0.428381,"mRNA binding"),
c("GO:0008327","methyl-CpG binding",0.009309036644208929,1,0.9390418609120744,0.46750139,"mRNA binding"),
c("GO:0008878","glucose-1-phosphate adenylyltransferase activity",0.021398809818155354,1,0.8890101964432305,0.29083106,"mRNA binding"),
c("GO:0016780","phosphotransferase activity, for other substituted phosphate groups",0.32920461222629094,1,0.8716406785187323,0.35102513,"mRNA binding"),
c("GO:0020037","heme binding",1.453841791525211,2,0.9578412310689256,0.13060477,"mRNA binding"),
c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9488267691429012,0.41077037,"mRNA binding"),
c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8701342346445498,0.27195941,"mRNA binding"),
c("GO:0044183","protein folding chaperone",0.39473153762799984,1,0.9779221791394449,0.32443613,"mRNA binding"),
c("GO:0003735","structural constituent of ribosome",2.128498588986873,6,0.9957656723800591,-0,"structural constituent of ribosome"),
c("GO:0003774","cytoskeletal motor activity",0.3926027441124113,1,1,-0,"cytoskeletal motor activity"),
c("GO:0003824","catalytic activity",60.727864217389204,13,1,-0,"catalytic activity"),
c("GO:0005096","GTPase activator activity",0.43493469016718705,2,0.9847701438688098,-0,"GTPase activator activity"),
c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
c("GO:0005515","protein binding",8.610051728351934,64,0.9676666389088371,0.07766851,"protein binding"),
c("GO:0005543","phospholipid binding",0.6999096105403309,2,0.9593664947520258,0.04714527,"phospholipid binding"),
c("GO:0008270","zinc ion binding",3.7731179768939866,9,0.9407825034733968,0.05738816,"zinc ion binding"),
c("GO:0000166","nucleotide binding",18.476222463002568,48,0.9174794025435898,0.48160605,"zinc ion binding"),
c("GO:0003676","nucleic acid binding",20.580503808256374,36,0.9411278764904909,0.20275684,"zinc ion binding"),
c("GO:0005509","calcium ion binding",1.3422486614434648,4,0.9515337555035566,0.3693843,"zinc ion binding"),
c("GO:0005525","GTP binding",1.780363236924562,9,0.9294139947809361,0.19635374,"zinc ion binding"),
c("GO:0016208","AMP binding",0.08038413014592037,1,0.9433231958232332,0.29697106,"zinc ion binding"),
c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9588095165721036,0.26177282,"zinc ion binding"),
c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9751546291666777,0.2693746,"zinc ion binding"),
c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9706463507141544,0.11877911,"zinc ion binding"),
c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.9348882157731002,0.33990845,"zinc ion binding"),
c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,5,0.9495984971076281,0.18168089,"zinc ion binding"),
c("GO:0070403","NAD+ binding",0.16643173804060435,1,0.9466720577769359,0.29755636,"zinc ion binding"),
c("GO:0071949","FAD binding",0.630202710371034,1,0.9411153726895474,0.30272555,"zinc ion binding"),
c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9568488657554564,0.25647876,"zinc ion binding"),
c("GO:0008289","lipid binding",1.291468066123697,4,0.9740360747584395,0.05834827,"lipid binding"),
c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
c("GO:0015288","porin activity",0.12362746592455742,2,0.93884781391447,-0,"porin activity"),
c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9130671375083721,0.39044024,"porin activity"),
c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9417239126056423,0.47124298,"porin activity"),
c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,1,0.9411586434365969,0.31858441,"porin activity"),
c("GO:0032977","membrane insertase activity",0.038238453523758646,1,0.9623999851096069,0.18410598,"porin activity"),
c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9564962400079116,0.04318065,"CoA-ligase activity"),
c("GO:0016740","transferase activity",20.627439270288612,56,0.9472493489854454,0.0817378,"transferase activity"),
c("GO:0016491","oxidoreductase activity",11.23609814678324,19,0.9511880784179745,0.10406278,"transferase activity"),
c("GO:0016787","hydrolase activity",22.243200201100027,39,0.9467167859017094,0.12712323,"transferase activity"),
c("GO:0016757","glycosyltransferase activity",2.478461155483397,8,0.874378058566644,0.05667901,"glycosyltransferase activity"),
c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9146616131447959,0.20205602,"glycosyltransferase activity"),
c("GO:0008168","methyltransferase activity",2.5264987116585993,8,0.8279886835043281,0.33973871,"glycosyltransferase activity"),
c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,3,0.8786654680684459,0.2917808,"glycosyltransferase activity"),
c("GO:0008483","transaminase activity",0.6992399275802186,2,0.8785936626683389,0.28934569,"glycosyltransferase activity"),
c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8690028745887369,0.36085078,"glycosyltransferase activity"),
c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8844175700784569,0.30154932,"glycosyltransferase activity"),
c("GO:0016767","geranylgeranyl-diphosphate geranylgeranyltransferase activity",0.019117009268633918,1,0.9020886228737552,0.20439581,"glycosyltransferase activity"),
c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.8864348573497074,0.49223587,"glycosyltransferase activity"),
c("GO:0047216","inositol 3-alpha-galactosyltransferase activity",0.0009202596968429504,1,0.9192503600880888,0.43743078,"glycosyltransferase activity"),
c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.912186342655241,0.17480439,"glycosyltransferase activity"),
c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.8656573719390799,0.45431535,"glycosyltransferase activity"),
c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8670665535523243,0.39355689,"glycosyltransferase activity"),
c("GO:0016829","lyase activity",3.6221510367468346,4,0.9571217737769077,0.05997119,"lyase activity"),
c("GO:0016853","isomerase activity",2.413124934500793,6,0.958900996678351,0.05646078,"isomerase activity"),
c("GO:0016871","cycloartenol synthase activity",0.0006563780006397912,1,0.9622379275230969,0.02586138,"cycloartenol synthase activity"),
c("GO:0016874","ligase activity",3.248266153118885,3,0.9576146728115937,0.0589874,"ligase activity"),
c("GO:0016887","ATP hydrolysis activity",4.018519084389943,17,0.8653089897370425,-0,"ATP hydrolysis activity"),
c("GO:0000146","microfilament motor activity",0.092387421083296,1,0.9701164243091338,0.49829267,"ATP hydrolysis activity"),
c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,2,0.8780590553391917,0.42844795,"ATP hydrolysis activity"),
c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.915502130688966,0.3692388,"ATP hydrolysis activity"),
c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9398070340024738,0.19229254,"ATP hydrolysis activity"),
c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,2,0.9046890764331185,0.31839638,"ATP hydrolysis activity"),
c("GO:0004712","protein serine/threonine/tyrosine kinase activity",0.056286631048107494,1,0.8424023554143275,0.44384691,"ATP hydrolysis activity"),
c("GO:0004715","non-membrane spanning protein tyrosine kinase activity",0.060526478133321286,1,0.8395108893642167,0.44607451,"ATP hydrolysis activity"),
c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,3,0.8410457538068864,0.42141607,"ATP hydrolysis activity"),
c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,5,0.818388215743473,0.48775459,"ATP hydrolysis activity"),
c("GO:0008233","peptidase activity",4.167383840940452,14,0.8383411349272583,0.3656957,"ATP hydrolysis activity"),
c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,4,0.8511527842038522,0.45595338,"ATP hydrolysis activity"),
c("GO:0008568","microtubule severing ATPase activity",0.013493446398287598,1,0.8880141756101257,0.41923201,"ATP hydrolysis activity"),
c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9607133938713633,0.45552111,"ATP hydrolysis activity"),
c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8995729857130226,0.39048386,"ATP hydrolysis activity"),
c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.9103623179383683,0.32827501,"ATP hydrolysis activity"),
c("GO:0033919","glucan 1,3-alpha-glucosidase activity",0.000915824710352141,1,0.940637920264108,0.4867618,"ATP hydrolysis activity"),
c("GO:0052742","phosphatidylinositol kinase activity",0.07485148449863564,1,0.8728322439305342,0.40778693,"ATP hydrolysis activity"),
c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9359184516152639,0.17583997,"ATP hydrolysis activity"),
c("GO:0106310","protein serine kinase activity",0.08584138102286135,5,0.8377396328998755,0.37283533,"ATP hydrolysis activity"),
c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9351685479926555,0.12676211,"ATP hydrolysis activity"),
c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9371411775817016,0.39401957,"ATP hydrolysis activity"),
c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9746923776739819,0.04901613,"carbohydrate binding"),
c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.972489785549527,0.06283846,"protein-containing complex binding"),
c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9721699882542946,0.03325955,"ER retention sequence binding"),
c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.979940750299162,0.49249984,"ER retention sequence binding"),
c("GO:0047769","arogenate dehydratase activity",0.02556769711951619,1,0.9482814083541394,0.03410747,"arogenate dehydratase activity"),
c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.954541698643944,0.36924486,"arogenate dehydratase activity"),
c("GO:0051082","unfolded protein binding",0.5891679978648201,4,0.9154767838810479,0.0463004,"unfolded protein binding"),
c("GO:0000149","SNARE binding",0.26943873427614345,1,0.9206649941369229,0.40559978,"unfolded protein binding"),
c("GO:0003779","actin binding",0.7421262469463454,3,0.9093692418236538,0.44654025,"unfolded protein binding"),
c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9207908819958428,0.40485813,"unfolded protein binding"),
c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9296308010906706,0.35368117,"unfolded protein binding"),
c("GO:0019900","kinase binding",0.4400460120978449,2,0.9080501720806812,0.42444044,"unfolded protein binding"),
c("GO:0030276","clathrin binding",0.11822565237875157,2,0.9254854729266847,0.37746303,"unfolded protein binding"),
c("GO:0030544","Hsp70 protein binding",0.07279586826014549,2,0.9280593508043852,0.36265291,"unfolded protein binding"),
c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9236883761991721,0.38789043,"unfolded protein binding"),
c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9284681538408958,0.36031388,"unfolded protein binding"),
c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,2,0.9477390078220179,0.25391209,"unfolded protein binding"),
c("GO:0032051","clathrin light chain binding",0.02013705616152008,1,0.9308049285888775,0.3284969,"unfolded protein binding"),
c("GO:0042393","histone binding",0.41579772345934446,3,0.9178683476412992,0.42217449,"unfolded protein binding"),
c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9166127381504338,0.42967902,"unfolded protein binding"),
c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9220851670854019,0.39725472,"unfolded protein binding"),
c("GO:0043621","protein self-association",0.05684543934594948,3,0.9293047565349698,0.3555383,"unfolded protein binding"),
c("GO:0046982","protein heterodimerization activity",0.29517274338906496,1,0.917948559934934,0.40897571,"unfolded protein binding"),
c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9103192038989392,0.47948026,"unfolded protein binding"),
c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9459807677703748,0.26331899,"unfolded protein binding"),
c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.936850692254367,0.31312228,"unfolded protein binding"),
c("GO:1990935","splicing factor binding",0.0007095978385295039,1,0.9458890993090746,0.26381104,"unfolded protein binding"),
c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9171112021374623,0.048252,"dioxygenase activity"),
c("GO:0004322","ferroxidase activity",0.06266857660838222,1,0.931427241325847,0.3160356,"dioxygenase activity"),
c("GO:0004601","peroxidase activity",0.4792135952914281,2,0.9129073152886723,0.37878622,"dioxygenase activity"),
c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9332768155056413,0.2837958,"dioxygenase activity"),
c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9043323631708091,0.44913707,"dioxygenase activity"),
c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9054338817098204,0.38257989,"dioxygenase activity"),
c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9150426384166358,0.41045836,"dioxygenase activity"),
c("GO:0016651","oxidoreductase activity, acting on NAD(P)H",0.7914588191768638,1,0.9166500814288117,0.39829104,"dioxygenase activity"),
c("GO:0016671","oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor",0.10547063123118373,1,0.9237378380452589,0.33002637,"dioxygenase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9295405986985653,0.04087,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.906244091234332,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9221916821559171,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9277573645136502,0.04086192,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.952334997282249,0.06542584,"protein binding"),
                     c("GO:0008233","peptidase activity",4.167383840940452,2,0.8083806670198322,0,"peptidase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.8782833034831994,0.42844795,"peptidase activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8929127249359601,0.30149286,"peptidase activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8343939322621955,0.46981862,"peptidase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8094514507876345,0.45595338,"peptidase activity"),
                     c("GO:0070139","SUMO-specific endopeptidase activity",0.006703482080858407,1,0.8342692672402504,0.29849456,"peptidase activity"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,1,0.8633346840746455,0.20765779,"peptidase activity"),
                     c("GO:0008289","lipid binding",1.291468066123697,1,0.9633534171978166,0.04613007,"lipid binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,4,0.935934694463014,0.08368982,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9415366565207313,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.9365708020414661,0.12712323,"hydrolase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9694178287693785,0.03203827,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9515729212267898,0.05675833,"isomerase activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9671609285137994,0.03203534,"identical protein binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9200546853179998,-0,"glutathione binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9125915317292687,0.31015716,"glutathione binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,1,0.8807569969649474,0.26746377,"glutathione binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.9148165787056828,0.12083997,"glutathione binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.8957533932258503,0.03333167,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0046961","proton-transporting ATPase activity, rotational mechanism",0.20284297713014954,1,0.9665876989887132,-0,"proton-transporting ATPase activity, rotational mechanism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
# 
# 
# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_MF_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8780001626499354,0.07971789,"phosphatidic acid biosynthetic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9398500057751906,0.08924942,"biosynthetic process"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.88685504952729,0,"response to salt stress"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,1,0.880047965361647,0.41360326,"response to salt stress"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8374780888310137,0.46171648,"response to salt stress"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8799541303160944,0.27995591,"response to salt stress"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8777411113302784,0.47274124,"response to salt stress"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8571310177279547,0.36420456,"response to salt stress"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.8662253906003201,0.07593517,"L-ascorbic acid metabolic process"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8030098185215787,-0,"proteolysis involved in protein catabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8158383632249968,0.35126771,"proteolysis involved in protein catabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.823447421984105,0.38896703,"proteolysis involved in protein catabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,2,0.8185128493832967,0.47291952,"proteolysis involved in protein catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9162140529896023,0.12747018,"proteolysis involved in protein catabolic process"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.911256018331753,0.17097988,"proteolysis involved in protein catabolic process"),
                     c("GO:0016926","protein desumoylation",0.02830261133305275,1,0.8730304986610063,0.2699367,"proteolysis involved in protein catabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.8738742220129585,0.28805767,"proteolysis involved in protein catabolic process"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9850316312530654,-0,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.968271662148011,0.12847604,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9299120410949966,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9041585510539594,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9353176309391981,0.35287398,"acquisition of seed longevity"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.9058910008046408,-0,"proton transmembrane transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.9001714182300902,0.4188028,"proton transmembrane transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8775240803746142,0.34990498,"proton transmembrane transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.9059983457890353,0.26647644,"proton transmembrane transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_BP_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
                     c("GO:0005737","cytoplasm",43.076660800468886,5,0.8667092671986336,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,4,0.8395125649390568,0.17160779,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,2,0.9858907243946746,0.00015294,"plasma membrane"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.8886866252011267,-0,"riboflavin synthase complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,1,0.8020284431559179,0.34749127,"riboflavin synthase complex"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7416531767183827,0,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,5,0.6957121104143723,0.35145602,"chloroplast"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7295874156668745,0.21968674,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7361187877932611,0.25147478,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7735171760825509,0.18302188,"chloroplast"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7101535089856432,0.18309297,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7357334093526413,0.22703414,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,1,0.7870061907625301,0.18475413,"chloroplast"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7578537489540528,0.21139748,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999303881506661,3.782E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,4,0.999850748184091,6.66E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
# 
# 
# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_CC_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0004672","protein kinase activity",3.5248673905776853,1,0.6979658095374469,0.0474854,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,1,0.9528504320620206,0.04547388,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,1,0.9096080929897793,0.08079952,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,1,0.9185086455647052,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,1,0.9106066208180881,0.12712323,"hydrolase activity"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.8988589484428764,-0,"L-ascorbic acid binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,1,0.8091385920543784,0.38577556,"L-ascorbic acid binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,1,0.8674928216631242,0.15860651,"L-ascorbic acid binding"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8797825248633305,0.02725726,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.8517859791708552,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.8941469198425223,0.45929075,"dioxygenase activity"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9619044822091802,0,"pectin acetylesterase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
n_hogs = length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, DEG in any C4, No DEG in C3, ReviGO 0.5, MF, # of HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006952","defense response",1.1604144588756624,1,0.9420584978843916,-0,"defense response"),
                     c("GO:0002238","response to molecule of fungal origin",0.007682595099288114,1,0.9420584978843916,0.20582114,"defense response"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9671470320818433,0.04068139,"coumarin biosynthetic process"),
                     c("GO:0031349","positive regulation of defense response",0.11892617777855335,1,0.8148885201816622,-0,"positive regulation of defense response"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9103045072845678,0.00664115,"protein autophosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.90100445252002,0.2739859,"protein autophosphorylation"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9920523869125935,0,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
n_hogs = length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, DEG in any C4, No DEG in C3, ReviGO 0.5, BP, # of HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998587407058247,5.614E-05,"extracellular region"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,1,0.645963007813005,0,"endoplasmic reticulum"),
                     c("GO:0005634","nucleus",16.5161752456724,1,0.770715974951891,0.32142968,"endoplasmic reticulum"),
                     c("GO:0005886","plasma membrane",17.177321395000487,1,0.9511769999912638,0.00015294,"plasma membrane"),
                     c("GO:0016020","membrane",49.2542153160787,1,0.9997727611556281,9.537E-05,"membrane"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
n_hogs = length(hogs$X)
# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, DEG in any C4, No DEG in C3, ReviGO 0.5, CC, # of HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,5,0.9999045347550348,4.647E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,8,0.8822387482013975,0.08417172,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,2,0.8355245771390575,0.17160779,"cytoplasm"),
                     c("GO:0009536","plastid",0.7624897559369845,6,0.7936096862849928,0,"plastid"),
                     c("GO:0000326","protein storage vacuole",0.0009591191132014498,1,0.875582815181251,0.10321831,"plastid"),
                     c("GO:0005634","nucleus",16.5161752456724,14,0.7116947111055937,0.35145602,"plastid"),
                     c("GO:0005730","nucleolus",1.2114693153196516,1,0.7907191209749264,0.20773614,"plastid"),
                     c("GO:0005739","mitochondrion",4.856981674016684,3,0.7465274424450951,0.2503098,"plastid"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,1,0.6335297289670339,0.20731432,"plastid"),
                     c("GO:0005773","vacuole",1.35235298008496,1,0.7810637377536104,0.21057364,"plastid"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,4,0.6453017114013119,0.19658734,"plastid"),
                     c("GO:0009507","chloroplast",0.6990388085931706,3,0.7655552769582101,0.17236349,"plastid"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.9998958906035237,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,17,0.9998477538758384,9.537E-05,"membrane"),
                     c("GO:0031012","extracellular matrix",0.6517835565339344,1,0.9257280167627275,3.812E-05,"extracellular matrix"),
                     c("GO:0005886","plasma membrane",17.177321395000487,8,0.8889188103189856,0.39486578,"extracellular matrix"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.999909750066312,5.465E-05,"cell projection"),
                     c("GO:0045202","synapse",1.0801172073373506,1,0.8640693522854003,4.855E-05,"synapse"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,1,0.9235279018696466,-0,"ribonucleoprotein complex"),
                     c("GO:0000015","phosphopyruvate hydratase complex",0.041644653723461905,1,0.8727375127462986,0.23819845,"ribonucleoprotein complex"),
                     c("GO:0030117","membrane coat",0.29445950681101296,1,0.7865848129568144,0.29013538,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_NoDEGC3_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)



# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE, No C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_NoDEGC3_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}







dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,1,1,-0,"intracellular iron ion homeostasis"),
                     c("GO:0007155","cell adhesion",1.1091577223710762,1,0.9916155097207778,0.01321071,"cell adhesion"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,2,0.8606460541627977,-0,"steroid metabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.949640318311361,0.15163831,"steroid metabolic process"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,1,0.8409969347756139,0.48590569,"steroid metabolic process"),
                     c("GO:0006564","L-serine biosynthetic process",0.06954757328091521,1,0.9040368338642071,0.39922301,"steroid metabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9485627677629662,0.10901633,"steroid metabolic process"),
                     c("GO:0008203","cholesterol metabolic process",0.2200421431132646,2,0.8431957573738564,0.45443823,"steroid metabolic process"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.8251796931729694,0.49264795,"steroid metabolic process"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9338566208781196,0.15725843,"steroid metabolic process"),
                     c("GO:0042759","long-chain fatty acid biosynthetic process",0.04469313344093403,1,0.8462410930269255,0.3981042,"steroid metabolic process"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.8365661311305623,0.42785472,"steroid metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9931813036421222,0.01017253,"cell population proliferation"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,3,0.8013592893626906,0.00987901,"auxin-activated signaling pathway"),
                     c("GO:0002238","response to molecule of fungal origin",0.007682595099288114,2,0.8594852318194514,0.47317835,"auxin-activated signaling pathway"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8787446648236106,0.4896406,"auxin-activated signaling pathway"),
                     c("GO:0006952","defense response",1.1604144588756624,3,0.8736176779792161,0.24614308,"auxin-activated signaling pathway"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8469156943498791,0.44951817,"auxin-activated signaling pathway"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.8910905201233614,0.45786068,"auxin-activated signaling pathway"),
                     c("GO:0009624","response to nematode",0.0011042035946362127,1,0.907037149048351,0.34321484,"auxin-activated signaling pathway"),
                     c("GO:0009643","photosynthetic acclimation",0.000325345701991027,1,0.9157020751352279,0.45106647,"auxin-activated signaling pathway"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.8741992433566023,0.4259419,"auxin-activated signaling pathway"),
                     c("GO:0010030","positive regulation of seed germination",0.0013630013121290752,1,0.9369115085247481,0.39302598,"auxin-activated signaling pathway"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.8970714416880673,0.17080919,"auxin-activated signaling pathway"),
                     c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,1,0.9283775213730041,0.33236762,"auxin-activated signaling pathway"),
                     c("GO:0032012","regulation of ARF protein signal transduction",0.05576721100946194,1,0.9172845681365718,0.46282204,"auxin-activated signaling pathway"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.7962498520573225,0.42828273,"auxin-activated signaling pathway"),
                     c("GO:0040008","regulation of growth",0.2209540969749061,1,0.9337256774892784,0.15374745,"auxin-activated signaling pathway"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.896899404294374,0.46472215,"auxin-activated signaling pathway"),
                     c("GO:0045927","positive regulation of growth",0.06802189911779062,1,0.9150735748532668,0.13260447,"auxin-activated signaling pathway"),
                     c("GO:0048209","regulation of vesicle targeting, to, from or within Golgi",2.4647401665986893E-06,1,0.9582528435805274,0.33066454,"auxin-activated signaling pathway"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,1,0.9181090724303207,0.42863919,"auxin-activated signaling pathway"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9321430218387391,0.1517125,"auxin-activated signaling pathway"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9213380167007119,0.3459764,"auxin-activated signaling pathway"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.9098844483473312,0.32241615,"auxin-activated signaling pathway"),
                     c("GO:1902882","regulation of response to oxidative stress",0.01609475328788944,1,0.9232242019189649,0.12060378,"auxin-activated signaling pathway"),
                     c("GO:1905421","regulation of plant organ morphogenesis",0.005429822587016912,1,0.9432231769503464,0.11290105,"auxin-activated signaling pathway"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9291044302566093,0.10716362,"auxin-activated signaling pathway"),
                     c("GO:0009793","embryo development ending in seed dormancy",0.0206816347379296,2,0.8674384876125509,-0,"embryo development ending in seed dormancy"),
                     c("GO:0010015","root morphogenesis",0.013467340270295237,1,0.8824575224266706,0.42719975,"embryo development ending in seed dormancy"),
                     c("GO:0010047","fruit dehiscence",0.002082705440775892,1,0.9427483459295852,0.46834396,"embryo development ending in seed dormancy"),
                     c("GO:0010087","phloem or xylem histogenesis",0.012814184126146584,1,0.906333814698279,0.36761229,"embryo development ending in seed dormancy"),
                     c("GO:0030154","cell differentiation",2.279586420543629,1,0.8786203199885205,0.4970999,"embryo development ending in seed dormancy"),
                     c("GO:0015031","protein transport",3.093438694074183,4,0.9071033113370945,0,"protein transport"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.1438718130046987,2,0.9156298282232382,0.2741348,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,2,0.9351114608168289,0.38566904,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9325408168304119,0.41259033,"protein transport"),
                     c("GO:0016135","saponin biosynthetic process",9.858960666394757E-06,1,0.9608118920285658,0.03120971,"saponin biosynthetic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,6,0.9038753398353797,0.05779371,"phosphorylation"),
                     c("GO:0044843","cell cycle G1/S phase transition",0.05461124787132716,1,0.9860181708808496,0.00959067,"cell cycle G1/S phase transition"),
                     c("GO:0000911","cytokinesis by cell plate formation",0.008712856488926366,1,0.9867939267757008,0.45372294,"cell cycle G1/S phase transition"),
                     c("GO:0050896","response to stimulus",17.567785530535815,2,1,-0,"response to stimulus"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,2,0.9688729196316782,0.01286368,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,1,0.9744268424189512,0.4063955,"cell wall organization"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9546893929296576,0.05291306,"macromolecule depalmitoylation"),
                     c("GO:0006013","mannose metabolic process",0.050640551462936674,1,0.9461162590969598,0.33438746,"macromolecule depalmitoylation"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1,0.8751974817615407,0.4946872,"macromolecule depalmitoylation"),
                     c("GO:0006412","translation",4.38869169324396,1,0.864790037921645,0.25211911,"macromolecule depalmitoylation"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,1,0.8788215882286222,0.43192082,"macromolecule depalmitoylation"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",0.0004411884898211653,1,0.9233227092369831,0.26905503,"macromolecule depalmitoylation"),
                     c("GO:0045489","pectin biosynthetic process",0.012338489273993038,1,0.9253633586050592,0.10510175,"macromolecule depalmitoylation"),
                     c("GO:1902066","regulation of cell wall pectin metabolic process",0.00012816648866313183,1,0.94974196728296,0.09252635,"regulation of cell wall pectin metabolic process"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1,0.9240142372941997,0.17609066,"regulation of cell wall pectin metabolic process"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.943261757992971,0.12668641,"regulation of cell wall pectin metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_NoDEGC3_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)



# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE, No C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)



# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_NoDEGC3_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}






dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100,3,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9674763781434593,0.03485749,"phosphopyruvate hydratase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,5,0.7495939076464454,0,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8800192775174226,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8806390324819089,0.48555173,"protein kinase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9064507700838336,0.2197454,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.7911728437529135,0.41376722,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8246484213155818,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8330343102594516,0.35123485,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8843305413378043,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8940567007377429,0.31489566,"protein kinase activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9220076937755323,0.21182638,"protein kinase activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,10,0.9513977847096342,0.06960742,"protein binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.8906488280094385,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0003954","NADH dehydrogenase activity",0.34065796483880617,1,0.9039063723566827,0.39393019,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.8921535204635909,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8819424727684569,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.918923140274019,0.32214763,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9352520491503108,0.19910628,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8974800889555598,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,9,0.9359688855353118,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,4,0.9419777446730677,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,8,0.9366806644692963,0.12712323,"hydrolase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8877803497394471,0.05402108,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8790333722593702,0.29096067,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8793658005888485,0.32827501,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9011474878550547,0.17478095,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9498762419375432,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9522034682976344,0.05541399,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9551689980615783,0.03465455,"intramolecular lyase activity"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9367313112785355,0.05445618,"heme binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,2,0.9014601427576004,0.31664803,"heme binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9624150224777499,0.35731485,"heme binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,2,0.9061840449122336,0.19568843,"heme binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,8,0.8858787347320689,0.45223656,"heme binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.9170075268241012,0.13060477,"heme binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8945491782583571,0.24450519,"heme binding"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9636134301546603,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9505756254153493,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9168542342812899,0.05306855,"protein dimerization activity"),
                     c("GO:0008083","growth factor activity",0.15856185451266308,1,0.8944806339717133,0.41153522,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.922727368133231,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.929258472289716,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9279111023338894,0.42311414,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9348016294787308,0.37611117,"protein dimerization activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9131244255135076,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,7,0.8582856764593443,0.42621458,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9100641587982604,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.912898418532941,0.3693843,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8933372358306384,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9114834216265212,0.21206178,"iron-sulfur cluster binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,1,0.9046146732125713,0.15052371,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)



# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE, No C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% topGO_ReviGO_Intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3),
                vp = "data")
      
    }
  )
}




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_NoDEGC3_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}







dev.off()


# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,6,0.9999072427077728,4.598E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,16,0.874050616778705,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,9,0.8327770134344908,0.17160779,"cytoplasm"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,6,0.9999378240664689,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,14,0.6925545081595136,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,2,0.7785261420993793,0.13280054,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,22,0.6830629422535041,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.7226568607680995,0.26090848,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,2,0.7726781380112001,0.17354881,"chloroplast"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,9,0.5870049357272622,0.2030659,"chloroplast"),
                     c("GO:0005788","endoplasmic reticulum lumen",0.11812570632538061,4,0.6719479964932582,0.48748182,"chloroplast"),
                     c("GO:0005797","Golgi medial cisterna",0.015012947362598856,1,0.6885502534137447,0.40558994,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,1,0.7649686227732585,0.43399423,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,13,0.7722484751190409,0.17236349,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,2,0.8129298426243511,0.12200202,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7273357215929946,0.30599076,"chloroplast"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9785429697301748,0.45176992,"chloroplast"),
                     c("GO:0070062","extracellular exosome",0.10153244871408715,1,0.7636087017257679,0.1380203,"chloroplast"),
                     c("GO:0012505","endomembrane system",6.893425119210306,2,0.9998990366728665,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,47,0.9998535964571377,9.537E-05,"membrane"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,1,0.8793774991137475,3.69E-05,"monoatomic ion channel complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,1,0.8390650576120395,0.21034994,"monoatomic ion channel complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7838973991835493,0.38761376,"monoatomic ion channel complex"),
                     c("GO:0005886","plasma membrane",17.177321395000487,18,0.9482739130660679,0.30272327,"monoatomic ion channel complex"),
                     c("GO:0019005","SCF ubiquitin ligase complex",0.14169866753507532,1,0.9369826166490048,0.21998964,"monoatomic ion channel complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.748904942431479,0.21281912,"monoatomic ion channel complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.07909597668631853,2,0.9233570780811821,0.05769476,"malate metabolic process"),
                     c("GO:0019752","carboxylic acid metabolic process",8.649401753337283,1,0.8886076936146385,0.46114007,"malate metabolic process"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,3,0.928749783285233,-0,"monoatomic ion transport"),
                     c("GO:0008643","carbohydrate transport",0.7503112720356397,1,0.9294830324797846,0.49791143,"monoatomic ion transport"),
                     c("GO:0009904","chloroplast accumulation movement",0.008823769796423306,1,0.9222736982608623,0.36807496,"monoatomic ion transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9452620568785269,0.42144071,"monoatomic ion transport"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,1,0.9336799109899854,0.46287171,"monoatomic ion transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9542245166376163,0.43140125,"monoatomic ion transport"),
                     c("GO:0070588","calcium ion transmembrane transport",0.4224811119566813,2,0.9073379524549021,0.3228373,"monoatomic ion transport"),
                     c("GO:0090150","establishment of protein localization to membrane",0.5833670263314108,1,0.9198806965903903,0.32758039,"monoatomic ion transport"),
                     c("GO:0120029","proton export across plasma membrane",0.0197376392541223,1,0.925506493023376,0.47511322,"monoatomic ion transport"),
                     c("GO:1905039","carboxylic acid transmembrane transport",1.2848271482650644,1,0.9022470319690351,0.37138613,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,10,0.863263325923334,0,"defense response"),
                     c("GO:0000165","MAPK cascade",0.14246691110973742,1,0.8460728647544693,0.41973941,"defense response"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.873976849600759,0.42761529,"defense response"),
                     c("GO:0006955","immune response",0.6433982378290884,1,0.8943936616876339,0.30172547,"defense response"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7699816731780171,2,0.8145752817286029,0.33768695,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8497324615346508,0.41395528,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,2,0.8798559923119899,0.22702886,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9069191088322223,0.20406574,"defense response"),
                     c("GO:0010332","response to gamma radiation",0.006090372951665361,1,0.8883604686898008,0.49899247,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8631862260538792,0.47637699,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,3,0.867231571110003,0.45104908,"defense response"),
                     c("GO:0070887","cellular response to chemical stimulus",2.703999888790924,1,0.8376580120760221,0.48850507,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9060059266188695,0.26100449,"defense response"),
                     c("GO:0071456","cellular response to hypoxia",0.026404761404771757,1,0.8333950216785795,0.40112537,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.886965027591401,0.26816179,"defense response"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,1,0.9916900306794391,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9533799742893749,0.08346192,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,1,0.9510640654740616,0.21768792,"biosynthetic process"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9563859857628859,0.08334549,"regulation of meristem structural organization"),
                     c("GO:0015979","photosynthesis",0.228607115192195,1,0.9550362872050443,0.07668617,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,6,0.8999483207890893,-0,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9385033876884507,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9166532622516966,0.43492456,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,2,0.8897313976828648,0.21165279,"protein ubiquitination"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8806067671586962,0.37832162,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.8827698762370719,0.37292242,"protein ubiquitination"),
                     c("GO:0006508","proteolysis",5.2622572267907,3,0.8941160077956305,0.47291952,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9384646078262686,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9371013497368873,0.12094833,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.8801258810506879,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,8,0.8965881077900265,0.44883397,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8784345080678305,0.14129671,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9349726911387148,0.45767193,"protein ubiquitination"),
                     c("GO:0042350","GDP-L-fucose biosynthetic process",0.02421607213683212,1,0.8664500406119435,0.2978867,"protein ubiquitination"),
                     c("GO:0042853","L-alanine catabolic process",0.0402368832197236,1,0.8981498040660442,0.47848651,"protein ubiquitination"),
                     c("GO:0046373","L-arabinose metabolic process",0.04285443727665141,1,0.906154900117013,0.41352792,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,1,0.8755268463641728,0.41108531,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9212692120426993,0.20437656,"protein ubiquitination"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8733651035091083,0.39252566,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9098462069731325,0.45115494,"protein ubiquitination"),
                     c("GO:0022900","electron transport chain",0.8207017864535318,1,0.9293086623013649,0.06987243,"electron transport chain"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9925821627158591,0.00919456,"establishment of cell polarity"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,2,0.8581895699032914,0.07681302,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006665","sphingolipid metabolic process",0.4116436494441469,1,0.8629569790534699,0.45795093,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016132","brassinosteroid biosynthetic process",0.019789398797620875,1,0.8247705280825638,0.39632639,"phosphatidylinositol dephosphorylation"),
                     c("GO:0019363","pyridine nucleotide biosynthetic process",0.26556096451000916,1,0.8322391304268336,0.39508505,"phosphatidylinositol dephosphorylation"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8568413888929544,0.39387188,"phosphatidylinositol dephosphorylation"),
                     c("GO:0048364","root development",0.0376070054619628,2,0.9018649371825134,-0,"root development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9285704594032439,0.38663663,"root development"),
                     c("GO:0009880","embryonic pattern specification",0.03627851051216611,1,0.9228303562653469,0.45392292,"root development"),
                     c("GO:0010228","vegetative to reproductive phase transition of meristem",0.011640967806845608,1,0.9100733787957433,0.43860092,"root development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.9280550655776705,0.43353892,"root development"),
                     c("GO:0055046","microgametogenesis",0.0014049018949612527,1,0.9325707904457887,0.37648736,"root development"),
                     c("GO:0050896","response to stimulus",17.567785530535815,3,1,-0,"response to stimulus"),
                     c("GO:0070814","hydrogen sulfide biosynthetic process",0.03377679924306844,1,0.9377314863417591,0.03784896,"hydrogen sulfide biosynthetic process"),
                     c("GO:0000103","sulfate assimilation",0.10156701278519879,1,0.9540059327190739,0.48363934,"hydrogen sulfide biosynthetic process"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,4,0.9448533068263878,-0,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9695139034288462,0.25762726,"cell wall organization"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9704702672104824,-0,"intracellular auxin homeostasis"),
                     c("GO:0098771","inorganic ion homeostasis",1.0570555799893464,1,0.9537452556217786,0.45565197,"intracellular auxin homeostasis"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,2,0.9075848927416297,-0,"negative regulation of defense response to bacterium"),
                     c("GO:0033314","mitotic DNA replication checkpoint signaling",0.042962885843981745,1,0.8283305829772818,0.48461656,"negative regulation of defense response to bacterium"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8909571607761799,0.47074713,"negative regulation of defense response to bacterium"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,2,0.9043708970367523,0.12505063,"negative regulation of defense response to bacterium"),
                     c("GO:0051513","regulation of monopolar cell growth",0.0030119124835835983,1,0.9353597894689105,0.41466435,"negative regulation of defense response to bacterium"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,1,0.9348753215790656,0.23293868,"negative regulation of defense response to bacterium"),
                     c("GO:1901001","negative regulation of response to salt stress",0.00037710524548959944,1,0.9213494043710891,0.43908873,"negative regulation of defense response to bacterium"),
                     c("GO:1901141","regulation of lignin biosynthetic process",0.00025140349699306625,1,0.9545320325178722,0.19785358,"negative regulation of defense response to bacterium"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9432801072575365,0.37546575,"negative regulation of defense response to bacterium"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,1,0.9244258774748108,0.37952976,"negative regulation of defense response to bacterium"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.9291745959451281,0.34953816,"negative regulation of defense response to bacterium"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9372987064370121,0.15519254,"negative regulation of defense response to bacterium"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9413769943469327,0.26299316,"negative regulation of defense response to bacterium"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_BP_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,3,0.9234219310168503,0.066121,"transcription cis-regulatory region binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9356732319507387,0.22593798,"transcription cis-regulatory region binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,11,0.9208157986383229,0.45223656,"transcription cis-regulatory region binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9558081254739234,0.31120183,"transcription cis-regulatory region binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,7,0.9892278958137051,0.40628023,"transcription cis-regulatory region binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.9419095650686188,0.31664803,"transcription cis-regulatory region binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9344428751834606,0.27764537,"transcription cis-regulatory region binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9559633361948003,0.14801235,"transcription cis-regulatory region binding"),
                     c("GO:0003674","molecular_function",100,5,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,1,0.9750857260813484,0.0533165,"chromatin binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,4,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,8,0.7712364958384572,0,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,2,0.8956709135511683,0.34856872,"protein kinase activity"),
                     c("GO:0004180","carboxypeptidase activity",0.4650748583587277,2,0.8278945384231705,0.43795071,"protein kinase activity"),
                     c("GO:0004781","sulfate adenylyltransferase (ATP) activity",0.04230755362907627,1,0.8917147765889218,0.41713542,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8730524232571509,0.30161168,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,1,0.8518403620480547,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8739778927909299,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,5,0.879331752206891,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,2,0.8896342995963697,0.21201958,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.9061901827145281,0.20092889,"protein kinase activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9606524111978103,0.38075906,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,1,0.880676507859515,0.48733868,"protein kinase activity"),
                     c("GO:0061630","ubiquitin protein ligase activity",0.6665452071699717,3,0.8153959195366739,0.45681379,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,20,0.9655182250160117,0.07766851,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,4,0.9259323685080195,-0,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,1,0.9597188024449858,0.21260669,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,2,0.94135628999035,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,2,0.9210495842930069,0.19635374,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,1,0.9537792124994369,0.24645234,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9511244949970609,0.26177282,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,1,0.9278312378715666,0.31800282,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9313622633952597,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.940630917527751,0.28523818,"zinc ion binding"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,1,0.9932715597699995,-0,"photoreceptor activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,1,0.9851425715980695,0.02993406,"auxin binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9091395406627808,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9103773871717047,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9339076030275826,0.31075808,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8921448713017424,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9225127918880266,0.35139584,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,20,0.9383201183959232,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,9,0.9431975421019412,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,14,0.9376609006114417,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,5,0.9504884412870329,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,4,0.9526540868602282,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,3,0.9510894178819743,0.05784577,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9491572006995631,0.04406353,"acid-amino acid ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9612669306650197,0.44358297,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9681016194434066,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,3,0.8546466589755487,0.05997119,"ATP hydrolysis activity"),
                     c("GO:0004566","beta-glucuronidase activity",0.008790143224784231,1,0.9039648109878481,0.4776459,"ATP hydrolysis activity"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,2,0.8256388430999533,0.32020834,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8803578753850313,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.8806899373247634,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.8794282273547864,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0046556","alpha-L-arabinofuranosidase activity",0.046915504593027235,1,0.8929193777485857,0.2148739,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9737960570931955,0.05647972,"carbohydrate binding"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.9657882125052901,-0,"membrane insertase activity"),
                     c("GO:0015085","calcium ion transmembrane transporter activity",0.312553455446547,1,0.9213806366238194,0.40515784,"membrane insertase activity"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,1,0.9248862352263971,0.43240123,"membrane insertase activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9433329921465046,0.40847859,"membrane insertase activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9264479109891621,0.21624994,"membrane insertase activity"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9943526289119874,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8916867778567433,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9279089282571554,0.05290415,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,2,0.9325759945518426,0.40485813,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9351237887764527,0.38789043,"unfolded protein binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9456650303771515,0.31955157,"unfolded protein binding"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,1,0.9372015203400825,0.37419503,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9300078306032469,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,3,0.928905480268335,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.9349583140706972,0.38898648,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9400609143876874,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9233978415106203,0.4679269,"unfolded protein binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.9317970257605913,0.41008684,"unfolded protein binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9332831307323307,0.0567159,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9528444094497709,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9653813869915613,0.01868194,"pterocarpan synthase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,1,0.9347660347962944,0.48814337,"pterocarpan synthase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9458411415056922,0.32867735,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)




# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_MF_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}

dev.off()


script for checking base mean
library("DESeq2")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\dea")

#load dea data



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


dea_Ta <- readRDS("dea_Triticum_aestivum")
dea_Os <- readRDS("dea_Oryza_sativa")
dea_Ms <- readRDS("dea_Miscanthus_sinensis")
dea_Sb <- readRDS("dea_Sorghum_bicolor")
dea_Zm <- readRDS("dea_Zea_mays")

dea_Ta_counts <- counts(dea_Ta)
dea_Os_counts <- counts(dea_Os)
dea_Ms_counts <- counts(dea_Ms)
dea_Sb_counts <- counts(dea_Sb)
dea_Zm_counts <- counts(dea_Zm)

reformat_col_names <- function(dea_table)
{
  names <- colnames(dea_table)
  names[grepl("control", names)] <- "control"
  names[grepl("drought", names)] <- "drought"
  colnames(dea_table) <- names
  
  return(dea_table)
}


dea_Ta_counts <- reformat_col_names(dea_Ta_counts)
dea_Os_counts <- reformat_col_names(dea_Os_counts)
dea_Ms_counts <- reformat_col_names(dea_Ms_counts)
dea_Sb_counts <- reformat_col_names(dea_Sb_counts)
dea_Zm_counts <- reformat_col_names(dea_Zm_counts)

Ms_genes <- rownames(dea_Ms_counts)
print("Misin03G292100" %in% Ms_genes)

get_baseMean <- function(HOG, add_genes = c())
{
  
  print(which_hypotheses(HOG))
  #get base means
  baseMeans = HOG_DE.a2tea[HOG_DE.a2tea$HOG == HOG,]
  gene_list = baseMeans$gene
  baseMean_list = baseMeans$baseMean
  control_m = c()
  control_sd = c()
  drought_m = c()
  drought_sd = c()
  for(g in add_genes)
  {
    gene_list <- append(gene_list, g)
    baseMean_list <- append(baseMean_list, HOG_DE.a2tea[HOG_DE.a2tea$gene == g,]$baseMean)
  }
  
  #get count means
  for(g in gene_list)
  {
    gene_calc = c()    
    if(g %in% rownames(dea_Ta_counts))gene_calc = dea_Ta_counts[g,]
    if(g %in% rownames(dea_Os_counts))gene_calc = dea_Os_counts[g,]
    if(g %in% rownames(dea_Ms_counts))gene_calc = dea_Ms_counts[g,]
    if(g %in% rownames(dea_Sb_counts))gene_calc = dea_Sb_counts[g,]
    if(g %in% rownames(dea_Zm_counts))gene_calc = dea_Zm_counts[g,]
    control_m <- append(control_m, mean(gene_calc[grepl("control", names(gene_calc))]))
    drought_m <- append(drought_m, mean(gene_calc[grepl("drought", names(gene_calc))]))
  }
  
  df <- data.frame(gene = gene_list, baseMeans = baseMean_list, control_m = control_m, drought_m = drought_m)
  
  return(df)  
}

#Msi Sb Zm vs Osj Ta

print(get_baseMean("N0.HOG0001331"))

print(get_baseMean("N0.HOG0002025"))

print(get_baseMean("N0.HOG0004017"))

print(get_baseMean("N0.HOG0004343", c("TraesCS2B02G215000", "TraesCS2B02G215100", "TraesCS2D02G195900")))

print(get_baseMean("N0.HOG0004368"))

#Zm vs Osj Ta

print(get_baseMean("N0.HOG0006824"))

print(get_baseMean("N0.HOG0017121", c("TraesCS1A02G209100")))

print(get_baseMean("N0.HOG0004652"))

print(get_baseMean("N0.HOG0000969"))

print(get_baseMean("N0.HOG0006817"))

print(get_baseMean("N0.HOG0008319"))

print(get_baseMean("N0.HOG0011584", c("TraesCS4A02G348100")))

print(get_baseMean("N0.HOG0000951"))

print(get_baseMean("N0.HOG0013739", c("TraesCS2B02G407100")))

print(get_baseMean("N0.HOG0005274"))

print(get_baseMean("N0.HOG0017487", c("TraesCS2D02G430700", "TraesCS2B02G454000", "Misin17G072400")))

print(get_baseMean("N0.HOG0006053"))

print(get_baseMean("N0.HOG0006529"))

print(get_baseMean("N0.HOG0009376"))

print(get_baseMean("N0.HOG0009828"))

print(get_baseMean("N0.HOG0010839", c("TraesCS3B02G469000")))

print(get_baseMean("N0.HOG0015605", c("TraesCS3A02G414400")))

print(get_baseMean("N0.HOG0007568"))

print(get_baseMean("N0.HOG0008285"))

print(get_baseMean("N0.HOG0009114"))

print(get_baseMean("N0.HOG0009637"))

print(get_baseMean("N0.HOG0014378"))

print(get_baseMean("N0.HOG0003297"))

print(get_baseMean("N0.HOG0008731"))

print(get_baseMean("N0.HOG0009150"))

print(get_baseMean("N0.HOG0011281"))

print(get_baseMean("N0.HOG0014437", c("TraesCS2A02G189600")))

print(get_baseMean("N0.HOG0006225"))

print(get_baseMean("N0.HOG0003381"))

print(get_baseMean("N0.HOG0004652"))

print(get_baseMean("N0.HOG0008146"))

print(get_baseMean("N0.HOG0008294"))

print(get_baseMean("N0.HOG0008863"))

print(get_baseMean("N0.HOG0009387"))

print(get_baseMean("N0.HOG0009416"))

print(get_baseMean("N0.HOG0008725"))

print(get_baseMean("N0.HOG0009290"))

print(get_baseMean("N0.HOG0004356"))

print(get_baseMean("N0.HOG0008535", c("TraesCS4B02G246900")))

print(get_baseMean("N0.HOG0005967"))

print(get_baseMean("N0.HOG0006609"))

print(get_baseMean("N0.HOG0016008", c("TraesCS1A02G083100", "TraesCS1B02G100600")))

print(get_baseMean("N0.HOG0006352"))

print(get_baseMean("N0.HOG0006529"))

print(get_baseMean("N0.HOG0007767"))

print(get_baseMean("N0.HOG0008263"))

print(get_baseMean("N0.HOG0008325"))

print(get_baseMean("N0.HOG0008383"))

print(get_baseMean("N0.HOG0008411"))

print(get_baseMean("N0.HOG0008493"))

print(get_baseMean("N0.HOG0009688"))

print(get_baseMean("N0.HOG0009883"))

print(get_baseMean("N0.HOG0009913"))

print(get_baseMean("N0.HOG0013207"))

print(get_baseMean("N0.HOG0014865"))

print(get_baseMean("N0.HOG0003181", c("TraesCS1A02G363700")))

print(get_baseMean("N0.HOG0006579"))

print(get_baseMean("N0.HOG0006996"))

print(get_baseMean("N0.HOG0001890"))

print(get_baseMean("N0.HOG0007568"))

print(get_baseMean("N0.HOG0008285"))

print(get_baseMean("N0.HOG0008440"))

print(get_baseMean("N0.HOG0009375"))

print(get_baseMean("N0.HOG0009637"))

print(get_baseMean("N0.HOG0009753"))

print(get_baseMean("N0.HOG0009914"))

print(get_baseMean("N0.HOG0011233"))

print(get_baseMean("N0.HOG0014378"))

print(get_baseMean("N0.HOG0014868"))

print(get_baseMean("N0.HOG0015144"))

print(get_baseMean("N0.HOG0016753", c("TraesCS4D02G191200", "TraesCS4B02G189800")))

print(get_baseMean("N0.HOG0004918"))

print(get_baseMean("N0.HOG0007725"))

print(get_baseMean("N0.HOG0012401"))

print(get_baseMean("N0.HOG0008768"))

print(get_baseMean("N0.HOG0008855"))

print(get_baseMean("N0.HOG0002218"))

print(get_baseMean("N0.HOG0008853"))

print(get_baseMean("N0.HOG0000491"))

print(get_baseMean("N0.HOG0004383"))

print(get_baseMean("N0.HOG0008974"))

print(get_baseMean("N0.HOG0013697"))

print(get_baseMean("N0.HOG0015040"))

print(get_baseMean("N0.HOG0002752"))

print(get_baseMean("N0.HOG0006546"))

print(get_baseMean("N0.HOG0007531"))

print(get_baseMean("N0.HOG0005085"))

print(get_baseMean("N0.HOG0005684"))

print(get_baseMean("N0.HOG0002005"))

print(get_baseMean("N0.HOG0009749"))  

print(get_baseMean("N0.HOG0002889", c("TraesCS5D02G052300", "TraesCS5B02G047200")))

print(get_baseMean("N0.HOG0007890"))

print(get_baseMean("N0.HOG0007925"))

print(get_baseMean("N0.HOG0008330"))

print(get_baseMean("N0.HOG0008381"))

print(get_baseMean("N0.HOG0008550"))

print(get_baseMean("N0.HOG0009033"))

print(get_baseMean("N0.HOG0008814"))



#Only all C3 upregulated all C4 expanded
print(get_baseMean("N0.HOG0004685"))

print(get_baseMean("N0.HOG0003448"))

print(get_baseMean("N0.HOG0003994"))

print(get_baseMean("N0.HOG0004908"))

print(get_baseMean("N0.HOG0003607"))

print(get_baseMean("N0.HOG0002654"))

print(get_baseMean("N0.HOG0003137"))

print(get_baseMean("N0.HOG0003710"))

print(get_baseMean("N0.HOG0003362"))

print(get_baseMean("N0.HOG0003870"))

print(get_baseMean("N0.HOG0004459"))

print(get_baseMean("N0.HOG0005632"))

print(get_baseMean("N0.HOG0001278"))

print(get_baseMean("N0.HOG0001572"))

print(get_baseMean("N0.HOG0003297"))

print(get_baseMean("N0.HOG0000135"))

print(get_baseMean("N0.HOG0000208"))

print(get_baseMean("N0.HOG0000340", c("TraesCS7D02G427900")))

print(get_baseMean("N0.HOG0000362"))

print(get_baseMean("N0.HOG0000585"))

print(get_baseMean("N0.HOG0000995"))

print(get_baseMean("N0.HOG0001045"))

print(get_baseMean("N0.HOG0001331"))

print(get_baseMean("N0.HOG0011860", c("Misin07G493400")))

print(get_baseMean("N0.HOG0008923"))

print(get_baseMean("N0.HOG0007510"))

print(get_baseMean("N0.HOG0007675"))

print(get_baseMean("N0.HOG0005262"))

print(get_baseMean("N0.HOG0014219", c("TraesCS7D02G384700")))

print(get_baseMean("N0.HOG0002869"))

print(get_baseMean("N0.HOG0004931", c("TraesCS1B02G141900")))

print(get_baseMean("N0.HOG0004908"))

print(get_baseMean("N0.HOG0000900"))

print(get_baseMean("N0.HOG0004512"))

print(get_baseMean("N0.HOG0003607"))

print(get_baseMean("N0.HOG0005093"))

print(get_baseMean("N0.HOG0016896", c("TraesCS4B02G052000", "TraesCS4A02G262900")))

print(get_baseMean("N0.HOG0002654"))

print(get_baseMean("N0.HOG0003137"))

print(get_baseMean("N0.HOG0005130", c("TraesCS1A02G291800")))

print(get_baseMean("N0.HOG0003710"))

print(get_baseMean("N0.HOG0003362"))

print(get_baseMean("N0.HOG0004560"))

print(get_baseMean("N0.HOG0003895"))

print(get_baseMean("N0.HOG0004459"))

print(get_baseMean("N0.HOG0004110"))

print(get_baseMean("N0.HOG0004099", c("TraesCS7D02G515300")))

print(get_baseMean("N0.HOG0005521"))

print(get_baseMean("N0.HOG0005249", c("TraesCS1A02G045700")))

print(get_baseMean("N0.HOG0004139"))

print(get_baseMean("N0.HOG0005093"))

print(get_baseMean("N0.HOG0005521"))


#Plot basemeans

library(ggplot2)
library(forcats)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\basemean")


HOG2869 <- get_baseMean("N0.HOG0002869")

HOG2869_order <- rev(c("TraesCS4A02G259700", "TraesCS4B02G055000","TraesCS4D02G055200",
                   "Os03g0724100","SORBI_3001G108700", "Misin01G098700", "Misin02G086400",
                   "SORBI_3001G108600", "Misin01G098600", "Misin02G086300", "Zm00001eb055720",
                   "Zm00001eb055730", "SORBI_3001G108500",
                   "Misin02G086200", "Misin01G098500"))

HOG2869 <- HOG2869[match(HOG2869_order, HOG2869$gene),]

pdf(".\\N0.HOG0002869_basemean.pdf", 4)

ggplot(HOG2869, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG4931 <- get_baseMean("N0.HOG0004931", c("TraesCS1B02G141900"))

HOG4931_order <- rev(c("TraesCS1D02G123600", "TraesCS1B02G141900","TraesCS1A02G122700",
                       "Os05g0187000","SORBI_3003G035601", "Zm00001eb304090", "Misin04G092000",
                       "SORBI_3002G087000", "Misin03G070000", "Zm00001eb268750", "Misin18G002100",
                       "SORBI_3010G007400", "Misin18G001300"))

HOG4931 <- HOG4931[match(HOG4931_order, HOG4931$gene),]

pdf(".\\N0.HOG0004931_basemean.pdf", 4)

ggplot(HOG4931, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG3137 <- get_baseMean("N0.HOG0003137")

HOG3137_order <- rev(c("TraesCS1A02G224800", "TraesCS1D02G226200","TraesCS1B02G238000",
                       "Os05g0349500","Misin06G248400", "Misin05G243200", "SORBI_3003G267800",
                       "Zm00001eb155830", "Misin06G274900", "Misin17G046300", "SORBI_3009G046700",
                       "Misin16G036400", "Zm00001eb360880", "Misin05G283600"))

HOG3137 <- HOG3137[match(HOG3137_order, HOG3137$gene),]

pdf(".\\N0.HOG0003137_basemean.pdf", 4)

ggplot(HOG3137, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG4139 <- get_baseMean("N0.HOG0004139")

HOG4139_order <- rev(c("SORBI_3010G185400", "Zm00001eb386290","Misin18G183800",
                       "Misin19G177000","Os09g0249700", "TraesCS5D02G172700", "TraesCS5B02G165200",
                       "TraesCS5A02G168400", "Zm00001eb269610", "Zm00001eb378360", "SORBI_3010G020400",
                       "Misin18G014700", "Misin19G020900"))

HOG4139 <- HOG4139[match(HOG4139_order, HOG4139$gene),]

pdf(".\\N0.HOG0004139_basemean.pdf", 4)

ggplot(HOG4139, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()





plotting the upset plots for HOGs containing DEG and HOGs species composition as well as HOGs where species are not DEG
# plot upset to show conserved HOGs, where species is DEG for each species


#load data ####################
library("DESeq2")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")


# conserved set HOGs that contain all 5 species
all_hogs <- unique(HOG_DE.a2tea$HOG)
conserved_hogs <- c()
x <- 0
for(hog in unique(HOG_DE.a2tea$HOG))
{
  x <- x + 1
  if(x %% 1000 == 0)
    cat(x, " of ", length(unique(HOG_DE.a2tea$HOG)), "\n")
  # you need to look at the gene ids and not the species column, because sometime
  # there is no entry for species for genes
  species <- unique(HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
  species <- species[!is.na(species)]
  species <- unique(substr(species, 1, 2))
  if(length(species) == 5)
  {
    conserved_hogs <- append(conserved_hogs, hog)
  }
  
}

conserved_set <- HOG_DE.a2tea[HOG_DE.a2tea$HOG %in% conserved_hogs, ]


# filter DE in HOGs per species


#filter HOG per species and DEG

Ta_conserved_DEG <- conserved_set[conserved_set$species == "Triticum_aestivum" & conserved_set$significant == "yes",]
Os_conserved_DEG <- conserved_set[conserved_set$species == "Oryza_sativa" & conserved_set$significant == "yes",]
Mi_conserved_DEG <- conserved_set[conserved_set$species == "Miscanthus_sinensis" & conserved_set$significant == "yes",]
Sb_conserved_DEG <- conserved_set[conserved_set$species == "Sorghum_bicolor" & conserved_set$significant == "yes",]
Zm_conserved_DEG <- conserved_set[conserved_set$species == "Zea_mays" & conserved_set$significant == "yes",]

#filter unique HOGs
Ta_conserved_HOG_DEG_unique <- unique(Ta_conserved_DEG$HOG)
Os_conserved_HOG_DEG_unique <- unique(Os_conserved_DEG$HOG)
Mi_conserved_HOG_DEG_unique <- unique(Mi_conserved_DEG$HOG)
Sb_conserved_HOG_DEG_unique <- unique(Sb_conserved_DEG$HOG)
Zm_conserved_HOG_DEG_unique <- unique(Zm_conserved_DEG$HOG)

All_conserved_HOG_DEG_unique <- unique(c(Ta_conserved_HOG_DEG_unique, Os_conserved_HOG_DEG_unique,
                                         Mi_conserved_HOG_DEG_unique, Sb_conserved_HOG_DEG_unique,
                                         Zm_conserved_HOG_DEG_unique))



library(UpSetR)

listInput <- list(Zea_mays = Zm_conserved_HOG_DEG_unique,
                  Sorghum_bicolor = Sb_conserved_HOG_DEG_unique, 
                  Miscanthus_sinensis = Mi_conserved_HOG_DEG_unique,
                  Oryza_sativa = Os_conserved_HOG_DEG_unique,
                  Triticum_aestivum = Ta_conserved_HOG_DEG_unique)

print(length(Ta_conserved_HOG_DEG_unique))
print(length(Os_conserved_HOG_DEG_unique))
print(length(Mi_conserved_HOG_DEG_unique))
print(length(Sb_conserved_HOG_DEG_unique))
print(length(Zm_conserved_HOG_DEG_unique))
print(length(All_conserved_HOG_DEG_unique))

sets <- c("Zea_mays", "Sorghum_bicolor", "Miscanthus_sinensis", "Oryza_sativa", "Triticum_aestivum")


pdf(".\\SetEnrichmentPlots\\Upset_HOGs_AtLeast1DEGfromSpecies.pdf", width = 9)

upset(fromList(listInput), nsets=5, nintersects=NA, sets.x.label = "Number of HOGs with at least \n1 DEG from SPECIES", mainbar.y.label = "Intersection of HOGs",
      keep.order = T, order.by = "freq",
      sets = sets,
      number.angles = 0,
      #show.numbers = T,
      point.size = 2, line.size = 1,
      text.scale = c(2,1.5,1,1.5,1.5,1.5)
      )

dev.off()
#load data ####################
library("DESeq2")
library(dplyr)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\OGSizes")

all_HOGs <- unique(HOG_DE.a2tea$HOG)
HOG_Ta <- unique(HOG_DE.a2tea$HOG[grepl("Traes", HOG_DE.a2tea$gene)])
HOG_Os <- unique(HOG_DE.a2tea$HOG[grepl("Os", HOG_DE.a2tea$gene)])
HOG_Ms <- unique(HOG_DE.a2tea$HOG[grepl("Misin", HOG_DE.a2tea$gene)])
HOG_Sb <- unique(HOG_DE.a2tea$HOG[grepl("SORBI", HOG_DE.a2tea$gene)])
HOG_Zm <- unique(HOG_DE.a2tea$HOG[grepl("Zm", HOG_DE.a2tea$gene)])
singletons <- sum(HOG_DE.a2tea$HOG == "singleton")




df <- data.frame(c(HOG_Ta, HOG_Os, HOG_Ms, HOG_Sb, HOG_Zm))


library(UpSetR)

listInput <- list(Zea_mays = HOG_Zm,  Sorghum_bicolor = HOG_Sb,Miscanthus_sinensis = HOG_Ms,
                  Oryza_sativa = HOG_Os, Triticum_aestivum = HOG_Ta)

sets <- c("Zea_mays", "Sorghum_bicolor", "Miscanthus_sinensis", "Oryza_sativa", "Triticum_aestivum")


pdf(".\\HOGs_Per_Species.pdf", width = 15)

upset(fromList(listInput), nsets=5, nintersects=NA, sets.x.label = "Number of HOGs\nfor each Species", mainbar.y.label = "Intersection of HOGs",
      keep.order = T, order.by = "freq",
      sets = sets,
      number.angles = 0,
      #show.numbers = T,
      point.size = 2, line.size = 1,
      text.scale = c(2,1.5,1,1.5,1.5,1.5)
)

dev.off()
# plot upset to show conserved HOGs, where species is DEG for each species

library(dplyr)
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


#load data ####################
library("DESeq2")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

# conserved set HOGs that contain all 5 species
all_hogs <- unique(HOG_DE.a2tea$HOG)
conserved_hogs <- c()
x <- 0
for(hog in unique(HOG_DE.a2tea$HOG))
{
  x <- x + 1
  if(x %% 1000 == 0)
    cat(x, " of ", length(unique(HOG_DE.a2tea$HOG)), "\n")
  # you need to look at the gene ids and not the species column, because sometime
  # there is no entry for species for genes
  species <- unique(HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]$gene)
  species <- species[!is.na(species)]
  species <- unique(substr(species, 1, 2))
  if(length(species) == 5)
  {
    conserved_hogs <- append(conserved_hogs, hog)
  }
  
}
conserved_set <- HOG_DE.a2tea[HOG_DE.a2tea$HOG %in% conserved_hogs, ]

# filter DE in HOGs per species


conserved_set$expansion <- rep("no", length(conserved_set$HOG))
conserved_set$Mi_expansion <- rep("no", length(conserved_set$HOG))
conserved_set$Sb_expansion <- rep("no", length(conserved_set$HOG))
conserved_set$Zm_expansion <- rep("no", length(conserved_set$HOG))

#mark where Mi Sb or Zm are expanded/at least one is expanded

for(hog in conserved_hogs)
{
  expanded <- c(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion,
                HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion,
                HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion)
  
  if("yes" %in% expanded)
    conserved_set[conserved_set$HOG == hog,]$expansion <- "yes"
  
  if(hog %in% HOG_level_list$hypothesis_3$HOG)
    if(HOG_level_list$hypothesis_3[HOG_level_list$hypothesis_3$HOG == hog,]$expansion == "yes")
      conserved_set[conserved_set$HOG == hog,]$Mi_expansion <- "yes"
  
  if(hog %in% HOG_level_list$hypothesis_6$HOG)
    if(HOG_level_list$hypothesis_6[HOG_level_list$hypothesis_6$HOG == hog,]$expansion == "yes")
      conserved_set[conserved_set$HOG == hog,]$Sb_expansion <- "yes"
  
  if(hog %in% HOG_level_list$hypothesis_9$HOG)
    if(HOG_level_list$hypothesis_9[HOG_level_list$hypothesis_9$HOG == hog,]$expansion == "yes")
      conserved_set[conserved_set$HOG == hog,]$Zm_expansion <- "yes"
}



#filter for Ta and Os which have at least one C4 expanded and at least one DE
Ta_conserved_expanded_DEG <- conserved_set[conserved_set$expansion == "yes" & conserved_set$species == "Triticum_aestivum" & conserved_set$significant == "yes",]
Os_conserved_expanded_DEG <- conserved_set[conserved_set$expansion == "yes" & conserved_set$species == "Oryza_sativa" & conserved_set$significant == "yes",]



#sum up and combine both HOG lists
Ta_conserved_expanded_HOG_DEG_unique <- unique(Ta_conserved_expanded_DEG$HOG)
Os_conserved_expanded_HOG_DEG_unique <- unique(Os_conserved_expanded_DEG$HOG)
Ta_Os_conserved_expanded_HOG_DEG_unique <- unique(intersect(Ta_conserved_expanded_HOG_DEG_unique, Os_conserved_expanded_HOG_DEG_unique))
#get all the info of every HOG which is relevant
Ta_Os_conserved_expanded_DEG <- conserved_set[conserved_set$HOG %in% Ta_Os_conserved_expanded_HOG_DEG_unique,]
Ta_Os_conserved_expanded_DEG <- Ta_Os_conserved_expanded_DEG[!is.na(Ta_Os_conserved_expanded_DEG$species),]



#go through and check every HOG in the C4 expanded C3 DE where C4 is expanded but not DE per species
Mi_conserved_expanded_HOG_DEG_unique <- c()
Sb_conserved_expanded_HOG_DEG_unique <- c()
Zm_conserved_expanded_HOG_DEG_unique <- c()

#go through every hog which as at least one C4 expanded and both C3 are DEG
for(hog in Ta_Os_conserved_expanded_HOG_DEG_unique)
{
  #is Mi expanded
  if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog, ]$Mi_expansion == "yes") > 0)
  {
    # is no Mi DEG
    if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog &
                                    Ta_Os_conserved_expanded_DEG$species == "Miscanthus_sinensis",]$significant == "yes") == 0)
      Mi_conserved_expanded_HOG_DEG_unique <- append(Mi_conserved_expanded_HOG_DEG_unique, hog)
  }
  
  #is Sb expanded
  if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog, ]$Sb_expansion == "yes") > 0)
  {
    # is no Sb DEG
    if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog &
                                        Ta_Os_conserved_expanded_DEG$species == "Sorghum_bicolor",]$significant == "yes") == 0)
      Sb_conserved_expanded_HOG_DEG_unique <- append(Sb_conserved_expanded_HOG_DEG_unique, hog)
  }
  
  #is Zm expanded
  if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog, ]$Zm_expansion == "yes") > 0 )
  {
    # is no Zm DEG
    if(sum(Ta_Os_conserved_expanded_DEG[Ta_Os_conserved_expanded_DEG$HOG == hog &
                                        Ta_Os_conserved_expanded_DEG$species == "Zea_mays",]$significant == "yes") == 0)
      Zm_conserved_expanded_HOG_DEG_unique <- append(Zm_conserved_expanded_HOG_DEG_unique, hog)
  }
  
  
  
}

library(UpSetR)

listInput <- list(Miscanthus_sinensis = Mi_conserved_expanded_HOG_DEG_unique, Sorghum_bicolor = Sb_conserved_expanded_HOG_DEG_unique,
                  Zea_mays = Zm_conserved_expanded_HOG_DEG_unique)

sets <- c("Miscanthus_sinensis", "Sorghum_bicolor", "Zea_mays")


pdf(".\\SetEnrichmentPlots\\Upset_HOGs_NoDEGInC4_AtLeastOneDEGPerC3.pdf", width = 9)

upset(fromList(listInput), nsets=5, nintersects=NA, sets.x.label = "Number of HOGs with\nExpansion and no DEG from C4\nAt Least One DEG per C3", mainbar.y.label = "Intersection of HOGs",
      keep.order = T, order.by = "freq",
      sets = sets,
      number.angles = 0,
      #show.numbers = T,
      point.size = 2, line.size = 1,
      text.scale = c(2,1.5,1,1.5,1.5,1.5)
      )

dev.off()

Script to caclulate parameters of hog sizes
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
script to calculate set enrichments
# pyhper("success in sample - 1",
#        "success in background",
#        "failure in background",
#        "sample size",
#        lower.tail = FALSE)


set <- c()
p <- c()

#All

{
#14771  Conserved
#3768   Conserved, No DEG in C4, At least One of Each C3 DEG
#183    Conserved, All C4 Expanded
#36    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG")
  
p <- append(p, phyper(36 -1,
       3768,
       14771 - 3768,
       183,
       lower.tail = FALSE))

#14771  Conserved
#1208   Conserved, No DEG in C4, At least One of Each C3 DEG, Only Upregulated
#183    Conserved, All C4 Expanded
#14    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Upregulated

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Upregulated")

p <- append(p, phyper(14 -1,
       1208,
       14771 - 1208,
       183,
       lower.tail = FALSE))

#14771  Conserved
#915   Conserved, No DEG in C4, At least One of Each C3 DEG, Only Downregulated
#183    Conserved, All C4 Expanded
#10    Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Downregulated

set <- append(set, "Conserved, All C4 Expanded, No DEG in C4, At least One of Each C3 DEG, Only Downregulated")

p <- append(p, phyper(10 -1,
       915,
       14771 - 915,
       183,
       lower.tail = FALSE))

}

#14771  Conserved
#180 Conserved At least one DEG in Zm, No C3 DE
#3786 Conserved Zm Expanded
#52 Conerved Zm Expanded, At least one DEG in Zm, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm, No C3 DE")

p <- append(p, phyper(52 - 1,
       180,
       14771 - 180,
       3786,
       lower.tail = FALSE))

#14771  Conserved
#90 Conserved At least one DEG in Zm Only Down, No C3 DE
#3786 Conserved Zm Expanded
#22 Conerved Zm Expanded, At least one DEG in Zm Only Down, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm Only Down, No C3 DE")

p <- append(p, phyper(22 - 1,
       90,
       14771 - 90,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#90 Conserved At least one DEG in Zm Only Up, No C3 DE
#3786 Conserved Zm Expanded
#30 Conerved Zm Expanded, At least one DEG in Zm Only Up, No C3 DE

set <- append(set, "Conerved Zm Expanded, At least one DEG in Zm Only Up, No C3 DE")

p <- append(p, phyper(30 - 1,
       90,
       14771 - 90,
       3786,
       lower.tail = FALSE))


#Constitutive Hypothesis

#Zm
{

#14771  Conserved
#3937 Conserved No DEG in Zm, At least one of Each C3 DE
#3786 Conserved Zm Expanded
#974 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 DE

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 DE")
  
  p <- append(p, phyper(974 - 1,
       3937,
       14771 - 3937,
       3786,
       lower.tail = FALSE))

#14771  Conserved
#964 Conserved No DEG in Zm, At least one of Each C3 Down
#3786 Conserved Zm Expanded
#236 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Down

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Down")
  
  p <- append(p, phyper(236 - 1,
       964,
       14771 - 964,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#1258 Conserved No DEG in Zm, At least one of Each C3 Up
#3786 Conserved Zm Expanded
#345 Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Up

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, At least one of Each C3 Up")
  
  p <- append(p, phyper(345 - 1,
       1258,
       14771 - 1258,
       3786,
       lower.tail = FALSE))


#*
#14771  Conserved
#1680 Conserved No DEG in Zm, All C3 DE
#3786 Conserved Zm Expanded
#513 Conserved Zm Expanded, No DEG in Zm, All C3 DE

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE")
  
  p <- append(p, phyper(513 - 1,
       1680,
       14771 - 1680,
       3786,
       lower.tail = FALSE))


#14771  Conserved
#382 Conserved No DEG in Zm, All C3 DE Only Down
#3786 Conserved Zm Expanded
#105 Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Down

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Down")
  
  p <- append(p, phyper(105 - 1,
       382,
       14771 - 382,
       3786,
       lower.tail = FALSE))

#*
#14771  Conserved
#610 Conserved No DEG in Zm, All C3 DE Only Up
#3786 Conserved Zm Expanded
#204 Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Up

  set <- append(set, "Conserved Zm Expanded, No DEG in Zm, All C3 DE Only Up")
  
  p <- append(p, phyper(204 - 1,
       610,
       14771 - 610,
       3786,
       lower.tail = FALSE))


}

#Ms
{
  
  #14771  Conserved
  #4812 Conserved No DEG in Ms, At least one of Each C3 DE
  #452 Conserved Ms Expanded
  #122 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 DE
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 DE")
  
  p <- append(p, phyper(122 - 1,
         4812,
         14771 - 4812,
         452,
         lower.tail = FALSE))
  
  #14771  Conserved
  #1192 Conserved No DEG in Ms, At least one of Each C3 Down
  #452 Conserved Ms Expanded
  #34 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Down
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Down")
  
  p <- append(p, phyper(34 - 1,
         1192,
         14771 - 1192,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #1560 Conserved No DEG in Ms, At least one of Each C3 Up
  #452 Conserved Ms Expanded
  #40 Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Up
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, At least one of Each C3 Up")
  
  p <- append(p, phyper(40 - 1,
         1560,
         14771 - 1560,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #2086 Conserved No DEG in Ms, All C3 DE
  #452 Conserved Ms Expanded
  #62 Conserved Ms Expanded, No DEG in Ms, All C3 DE
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE")
  
  p <- append(p, phyper(62 - 1,
         2086,
         14771 - 2086,
         452,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #473 Conserved No DEG in Ms, All C3 DE Only Down
  #452 Conserved Ms Expanded
  #16 Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Down
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Down")
  
  p <- append(p, phyper(16 - 1,
         473,
         14771 - 473,
         452,
         lower.tail = FALSE))
  
  #14771  Conserved
  #787 Conserved No DEG in Ms, All C3 DE Only Up
  #452 Conserved Ms Expanded
  #24 Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Up
  
  set <- append(set, "Conserved Ms Expanded, No DEG in Ms, All C3 DE Only Up")
  
  p <- append(p, phyper(24 - 1,
         787,
         14771 - 787,
         452,
         lower.tail = FALSE))
  
  
}

#Sb
{
  
  #14771  Conserved
  #4913 Conserved No DEG in Sb, At least one of Each C3 DE
  #981 Conserved Sb Expanded
  #292 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 DE
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 DE")
  
  p <- append(p, phyper(292 - 1,
         4913,
         14771 - 4913,
         981,
         lower.tail = FALSE))
  
  #14771  Conserved
  #1233 Conserved No DEG in Sb, At least one of Each C3 Down
  #981 Conserved Sb Expanded
  #80 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Down
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Down")
  
  p <- append(p, phyper(80 - 1,
         1233,
         14771 - 1233,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #1577 Conserved No DEG in Sb, At least one of Each C3 Up
  #981 Conserved Sb Expanded
  #104 Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Up
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, At least one of Each C3 Up")
  
  p <- append(p, phyper(104 - 1,
         1577,
         14771 - 1577,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #2146 Conserved No DEG in Sb, All C3 DE
  #981 Conserved Sb Expanded
  #139 Conserved Sb Expanded, No DEG in Sb, All C3 DE
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE")
  
  p <- append(p, phyper(139 - 1,
         2146,
         14771 - 2146,
         981,
         lower.tail = FALSE))
  
  
  #14771  Conserved
  #489 Conserved No DEG in Sb, All C3 DE Only Down
  #981 Conserved Sb Expanded
  #31 Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Down
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Down")
  
  p <- append(p, phyper(31 - 1,
         489,
         14771 - 489,
         981,
         lower.tail = FALSE))
  
  #14771  Conserved
  #814 Conserved No DEG in Sb, All C3 DE Only Up
  #981 Conserved Sb Expanded
  #55 Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Up
  
  set <- append(set, "Conserved Sb Expanded, No DEG in Sb, All C3 DE Only Up")
  
  p <- append(p, phyper(55 - 1,
         814,
         14771 - 814,
         981,
         lower.tail = FALSE))
  
  
}


results <- data.frame(set, p)
results$padj <- p.adjust(results$p, method="bonferroni")
script to compare the all vs the conserved set
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
script to dicern which hypotheses contain a certain HOG
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

Script to calculate GO-Term enrichment which are in the At Least One of each C3 DE sets
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")







# 14 HOGs
# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp ################################################
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs$x
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs)
plot_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")



# 10 HOGs
# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown ################################################
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs$x
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs)
plot_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
