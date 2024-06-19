#load data ####################
library("DESeq2")
library(dplyr)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\OGSizes")

all_HOGs <- unique(HOG_DE.a2tea$HOG)
HOG_Ta <- unique(HOG_DE.a2tea$HOG[grepl("Traes", HOG_DE.a2tea$gene)])
HOG_Os <- unique(HOG_DE.a2tea$HOG[grepl("Os", HOG_DE.a2tea$gene)])
HOG_Ms <- unique(HOG_DE.a2tea$HOG[grepl("Misin", HOG_DE.a2tea$gene)])
HOG_Sb <- unique(HOG_DE.a2tea$HOG[grepl("SORBI", HOG_DE.a2tea$gene)])
HOG_Zm <- unique(HOG_DE.a2tea$HOG[grepl("Zm", HOG_DE.a2tea$gene)])
singletons <- sum(HOG_DE.a2tea$HOG == "singleton")




df <- data.frame(c(HOG_Ta, HOG_Os, HOG_Ms, HOG_Sb, HOG_Zm))


library(UpSetR)

listInput <- list(Zea_mays = HOG_Zm,  Sorghum_bicolor = HOG_Sb,Miscanthus_sinensis = HOG_Ms,
                  Oryza_sativa = HOG_Os, Triticum_aestivum = HOG_Ta)

sets <- c("Zea_mays", "Sorghum_bicolor", "Miscanthus_sinensis", "Oryza_sativa", "Triticum_aestivum")


pdf(".\\HOGs_Per_Species.pdf", width = 15)

upset(fromList(listInput), nsets=5, nintersects=NA, sets.x.label = "Number of HOGs\nfor each Species", mainbar.y.label = "Intersection of HOGs",
      keep.order = T, order.by = "freq",
      sets = sets,
      number.angles = 0,
      #show.numbers = T,
      point.size = 2, line.size = 1,
      text.scale = c(2,1.5,1,1.5,1.5,1.5)
)

dev.off()
