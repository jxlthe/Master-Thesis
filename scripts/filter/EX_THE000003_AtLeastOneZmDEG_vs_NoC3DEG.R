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
write.csv(AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, ".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")