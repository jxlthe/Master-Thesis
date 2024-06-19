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



