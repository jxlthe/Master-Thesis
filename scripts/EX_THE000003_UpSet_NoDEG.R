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

