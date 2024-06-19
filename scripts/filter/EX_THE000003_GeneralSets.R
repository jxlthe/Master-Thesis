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
