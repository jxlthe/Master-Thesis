library("DESeq2")
library(dplyr)

#calculate GO Term enrichment by hand. Use topGO instead, discard this script

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

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




calc_HOG_enrich <- function(set)
{
  
  p <- rep(0, length(set$GOID))
  total <- rep(0, length(set$GOID))
  interest <- rep(0, length(set$GOID))
  set_total <- rep(0, length(set$GOID))
  set_interest <- rep(0, length(set$GOID))
  n <- 1
  
  #conserved set
  SFA_bound <- bind_rows(SFA$Oryza_sativa, SFA$Triticum_aestivum, SFA$Sorghum_bicolor, SFA$Miscanthus_sinensis, SFA$Zea_mays)
  SFA_bound <- SFA_bound[SFA_bound$HOG %in% conserved_HOGs,]
  
  for(go in set$GOID)
  {
    #bind all SFA sets
    hog_count_total <- unique(SFA_bound$HOG) # get unique HOGs
    hog_count_total <- length(hog_count_total[!is.na(hog_count_total)]) #filter NA
    hog_count_interest <- unique(SFA_bound[grepl(go, SFA_bound$'Gene-Ontology-Term'),]$HOG)
    hog_count_interest <- length(hog_count_interest[!is.na(hog_count_interest)])
    set_count_total <- length(unique(unlist(str_split(set$HOGs_with_GO, " "))))
    set_count_interest <- set[set$GOID == go,]$HOGs_with_GO_COUNT
    
    #cat(hog_count_total, " ", hog_count_interest, " ", set_count_total, " ", set_count_interest, "\n")
    
    p[n] <- phyper(set_count_interest - 1,
                     hog_count_interest,
                     hog_count_total - hog_count_interest,
                     set_count_total,
                     lower.tail= FALSE)
    total[n] <- hog_count_total
    interest[n] <- hog_count_interest
    set_total[n] <- set_count_total
    set_interest[n] <- set_count_interest
    n <- n + 1
  }
  
  p_adjust <- p.adjust(p, method = "BH")
  
  result <- data.frame(set$GOID, set$GOTERM, p_adjust, total, interest, set_total, set_interest)
  
  result <- result[order(result$p_adjust, decreasing = F),]
  
  return(result)
  
}

library(stringr)

MSZ_expanded_deg_any_exp_species_GO_data <- read.csv(".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_unique_GO.csv")
MSZ_expanded_deg_any_exp_species_enrich_result <- calc_HOG_enrich(MSZ_expanded_deg_any_exp_species_GO_data)
write.csv(MSZ_expanded_deg_any_exp_species_enrich_result, ".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_GO_HOG_count_pvalue.csv")

AllC4Exp_NoDEGC4_AllC3Up_GO_data <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_unique_GO.csv")
AllC4Exp_NoDEGC4_AllC3Up_enrich_result <- calc_HOG_enrich(AllC4Exp_NoDEGC4_AllC3Up_GO_data)
write.csv(AllC4Exp_NoDEGC4_AllC3Up_enrich_result, ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_GO_HOG_count_pvalue.csv")

AllC4Exp_DEGAnyC4_NoDEGC3_GO_data <- read.csv(".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_unique_GO.csv")
AllC4Exp_DEGAnyC4_NoDEGC3_enrich_result <- calc_HOG_enrich(AllC4Exp_DEGAnyC4_NoDEGC3_GO_data)
write.csv(AllC4Exp_DEGAnyC4_NoDEGC3_enrich_result, ".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_GO_HOG_count_pvalue.csv")

ZmExp_OnlyAnyZmDEG_NoDEGC3_GO_data <- read.csv(".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_unique_GO.csv")
ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich_result <- calc_HOG_enrich(ZmExp_OnlyAnyZmDEG_NoDEGC3_GO_data)
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich_result, ".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_GO_HOG_count_pvalue.csv")

ZmExp_AllZmDEG_NoDEGC3_GO_data <- read.csv(".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_unique_GO.csv")
ZmExp_AllZmDEG_NoDEGC3_enrich_result <- calc_HOG_enrich(ZmExp_AllZmDEG_NoDEGC3_GO_data)
write.csv(ZmExp_AllZmDEG_NoDEGC3_enrich_result, ".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_GO_HOG_count_pvalue.csv")
