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


