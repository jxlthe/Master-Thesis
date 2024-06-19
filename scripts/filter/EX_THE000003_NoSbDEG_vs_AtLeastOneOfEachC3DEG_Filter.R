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


