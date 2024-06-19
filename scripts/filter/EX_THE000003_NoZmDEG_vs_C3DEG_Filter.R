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

