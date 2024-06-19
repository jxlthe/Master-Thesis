# plot upset to show conserved HOGs, where species is DEG for each species


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


#filter HOG per species and DEG

Ta_conserved_DEG <- conserved_set[conserved_set$species == "Triticum_aestivum" & conserved_set$significant == "yes",]
Os_conserved_DEG <- conserved_set[conserved_set$species == "Oryza_sativa" & conserved_set$significant == "yes",]
Mi_conserved_DEG <- conserved_set[conserved_set$species == "Miscanthus_sinensis" & conserved_set$significant == "yes",]
Sb_conserved_DEG <- conserved_set[conserved_set$species == "Sorghum_bicolor" & conserved_set$significant == "yes",]
Zm_conserved_DEG <- conserved_set[conserved_set$species == "Zea_mays" & conserved_set$significant == "yes",]

#filter unique HOGs
Ta_conserved_HOG_DEG_unique <- unique(Ta_conserved_DEG$HOG)
Os_conserved_HOG_DEG_unique <- unique(Os_conserved_DEG$HOG)
Mi_conserved_HOG_DEG_unique <- unique(Mi_conserved_DEG$HOG)
Sb_conserved_HOG_DEG_unique <- unique(Sb_conserved_DEG$HOG)
Zm_conserved_HOG_DEG_unique <- unique(Zm_conserved_DEG$HOG)

All_conserved_HOG_DEG_unique <- unique(c(Ta_conserved_HOG_DEG_unique, Os_conserved_HOG_DEG_unique,
                                         Mi_conserved_HOG_DEG_unique, Sb_conserved_HOG_DEG_unique,
                                         Zm_conserved_HOG_DEG_unique))



library(UpSetR)

listInput <- list(Zea_mays = Zm_conserved_HOG_DEG_unique,
                  Sorghum_bicolor = Sb_conserved_HOG_DEG_unique, 
                  Miscanthus_sinensis = Mi_conserved_HOG_DEG_unique,
                  Oryza_sativa = Os_conserved_HOG_DEG_unique,
                  Triticum_aestivum = Ta_conserved_HOG_DEG_unique)

print(length(Ta_conserved_HOG_DEG_unique))
print(length(Os_conserved_HOG_DEG_unique))
print(length(Mi_conserved_HOG_DEG_unique))
print(length(Sb_conserved_HOG_DEG_unique))
print(length(Zm_conserved_HOG_DEG_unique))
print(length(All_conserved_HOG_DEG_unique))

sets <- c("Zea_mays", "Sorghum_bicolor", "Miscanthus_sinensis", "Oryza_sativa", "Triticum_aestivum")


pdf(".\\SetEnrichmentPlots\\Upset_HOGs_AtLeast1DEGfromSpecies.pdf", width = 9)

upset(fromList(listInput), nsets=5, nintersects=NA, sets.x.label = "Number of HOGs with at least \n1 DEG from SPECIES", mainbar.y.label = "Intersection of HOGs",
      keep.order = T, order.by = "freq",
      sets = sets,
      number.angles = 0,
      #show.numbers = T,
      point.size = 2, line.size = 1,
      text.scale = c(2,1.5,1,1.5,1.5,1.5)
      )

dev.off()
