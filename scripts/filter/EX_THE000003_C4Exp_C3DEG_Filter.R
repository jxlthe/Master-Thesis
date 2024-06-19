setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts\\filter")
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




