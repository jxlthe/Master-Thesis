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

