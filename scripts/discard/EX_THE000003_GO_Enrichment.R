setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\filtered_set")

#GO Term enrichment was done using the tool https://www.geneontology.org/
# this was deemed not usefull because Msi was not included and overall it was checked against a different BG
# this is discarded

#Ms Sb Zm vs Ta Os, Expanded, DEG in any expanded species

library(ggplot2)
library(tidyverse)
library(cowplot)

plot_GO_graph <- function(Terms)
{
  if(dim(Terms)[1] == 0)
  {
    ggplot() + ggtitle("\tNo Data")
  }
  else
  {
    maxn <- max(Terms[3])
    Terms_order <- Terms[order(Terms$upload_1..FDR., decreasing=T),]
    ggplot(Terms_order, ) + 
      geom_point(aes(x = upload_1..FDR., y = fct_inorder(Terms_order[,1]), size = Terms_order[,3])) +
      geom_vline(xintercept=0.05) +
      geom_vline(xintercept=0.01) +
      geom_vline(xintercept=0.001) +
      scale_size_binned_area(name = "Annotated", limits = c(0, maxn),
                             breaks = seq(0, maxn, 
                                          maxn/5)) +
      scale_x_reverse() +
      xlab("FDR") +
      ylab("GO Terms")
  }
  
}

plot_GO <- function(BP, MF, CC, title, pdf, width=10, height=10)
{
  bp <- plot_GO_graph(BP)
  mf <- plot_GO_graph(MF)
  cc <- plot_GO_graph(CC)
  
  grid <- plot_grid(bp, mf, cc, labels = c("BP", "MF", "CC"), ncol=1)
  t <- ggdraw() + draw_label(title)
  plot <- plot_grid(t, grid, ncol=1, rel_heights=c(0.1, 1))
  ggsave(filename = pdf, plot = plot, width=width, height=height)
}

#hypotheses Ms, Sb, Zm vs Ta, Os, Expanded, DEG in any expanded species (51 HOGs)

Zm_BP <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Zeamays_BP.txt", header = T, sep = "\t", skip = 11)
Zm_MF <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Zeamays_MF.txt", header = T, sep = "\t", skip = 11)
Zm_CC <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Zeamays_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Zm_BP, Zm_MF, Zm_CC, "Enriched GO-Terms for Zea mays\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_Zm_GO.pdf", width = 10, height = 10)


Sb_BP <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Sorghumbicolor_BP.txt", header = T, sep = "\t", skip = 11)
Sb_MF <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Sorghumbicolor_MF.txt", header = T, sep = "\t", skip = 11)
Sb_CC <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Sorghumbicolor_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Sb_BP, Sb_MF, Sb_CC, "Enriched GO-Terms for Sorghum bicolor\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_Sb_GO.pdf", width = 10, height = 20)


Os_BP <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Oryzasativa_BP.txt", header = T, sep = "\t", skip = 11)
Os_MF <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Oryzasativa_MF.txt", header = T, sep = "\t", skip = 11)
Os_CC <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Oryzasativa_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Os_BP, Os_MF, Os_CC, "Enriched GO-Terms for Oryza sativa\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_Os_GO.pdf", width = 10, height = 10)


Ta_BP <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Triticumaestivum_BP.txt", header = T, sep = "\t", skip = 11)
Ta_MF <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Triticumaestivum_MF.txt", header = T, sep = "\t", skip = 11)
Ta_CC <- read.table("MsSbZmvsTaOs_Expanded_DEGinAnyExpSpecies_GOEnrichment_Triticumaestivum_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Ta_BP, Ta_MF, Ta_CC, "Enriched GO-Terms for Triticum aestivum\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_Ta_GO.pdf", width = 10, height = 20)


# only C3 plants are upregulated and all C4 expanded (14 HOGs


Zm_BP <- read.table("OnlyC3UpAllC4Expanded_Zeamays_BP.txt", header = T, sep = "\t", skip = 11)
Zm_MF <- read.table("OnlyC3UpAllC4Expanded_Zeamays_MF.txt", header = T, sep = "\t", skip = 11)
Zm_CC <- read.table("OnlyC3UpAllC4Expanded_Zeamays_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Zm_BP, Zm_MF, Zm_CC, "Enriched GO-Terms for Zea mays\nOnly All C3 upregulated, All C4 expanded", "OnlyC3UpAllC4Expanded_Zm_GO.pdf", width = 10, height = 10)


Sb_BP <- read.table("OnlyC3UpAllC4Expanded_Sorghumbicolor_BP.txt", header = T, sep = "\t", skip = 11)
Sb_MF <- read.table("OnlyC3UpAllC4Expanded_Sorghumbicolor_MF.txt", header = T, sep = "\t", skip = 11)
Sb_CC <- read.table("OnlyC3UpAllC4Expanded_Sorghumbicolor_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Sb_BP, Sb_MF, Sb_CC, "Enriched GO-Terms for Sorghum bicolor\nOnly All C3 upregulated, All C4 expanded", "OnlyC3UpAllC4Expanded_Sb_GO.pdf", width = 10, height = 20)


Os_BP <- read.table("OnlyC3UpAllC4Expanded_Oryzasativa_BP.txt", header = T, sep = "\t", skip = 11)
Os_MF <- read.table("OnlyC3UpAllC4Expanded_Oryzasativa_MF.txt", header = T, sep = "\t", skip = 11)
Os_CC <- read.table("OnlyC3UpAllC4Expanded_Oryzasativa_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Os_BP, Os_MF, Os_CC, "Enriched GO-Terms for Oryza sativa\nOnly All C3 upregulated, All C4 expanded", "OnlyC3UpAllC4Expanded_Os_GO.pdf", width = 10, height = 10)


Ta_BP <- read.table("OnlyC3UpAllC4Expanded_Triticumaestivum_BP.txt", header = T, sep = "\t", skip = 11)
Ta_MF <- read.table("OnlyC3UpAllC4Expanded_Triticumaestivum_MF.txt", header = T, sep = "\t", skip = 11)
Ta_CC <- read.table("OnlyC3UpAllC4Expanded_Triticumaestivum_CC.txt", header = T, sep = "\t", skip = 11)
plot_GO(Ta_BP, Ta_MF, Ta_CC, "Enriched GO-Terms for Triticum aestivum\nOnly All C3 upregulated, All C4 expanded", "OnlyC3UpAllC4Expanded_Ta_GO.pdf", width = 10, height = 20)



