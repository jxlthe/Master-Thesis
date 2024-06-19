library(topGO)

#The analysis was done for sets with more than 20 HOGs

#load data ####################
library("DESeq2")
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(forcats)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

#make a list of conserved HOGs
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

#get all HOGs/GO list from conserved set
SFA_bound <- bind_rows(SFA$Oryza_sativa, SFA$Triticum_aestivum, SFA$Sorghum_bicolor, SFA$Miscanthus_sinensis, SFA$Zea_mays)
SFA_conserved <- SFA_bound[SFA_bound$HOG %in% conserved_HOGs,]


#make a HOG-GO annotation list

#split DF by HOG
HOG2GO <- rep("NA", length(unique(SFA_conserved$HOG[!is.na(SFA_conserved$HOG)])))
names(HOG2GO) <- unique(SFA_conserved$HOG[!is.na(SFA_conserved$HOG)])


#get HOG and GO pair
annotation <- SFA_conserved[,c("HOG", "Gene-Ontology-Term")]
# split the strings of the GO Terms
annotation$`Gene-Ontology-Term` <- str_split(SFA_conserved$`Gene-Ontology-Term`, ", ")
#split the dataframe by HOG number
annotation <- split(annotation, annotation$HOG)
#unlist each GOlist of each HOG independently
for(i in 1:length(annotation))
{
  HOG2GO[unique(annotation[[i]]$HOG)] <- list(unique(unlist(annotation[[i]]$`Gene-Ontology-Term`)))
}


#default background is the conserved HOGs
test_GO_enrich <- function(interest, background = conserved_HOGs, nTop = 200)
{
  #my HOGs are my genes
  
  #list of HOG Of Interest
  HOGs_interest <- interest
  
  #filter for background so you can have different backgrounds too
  SFA_background <- SFA_conserved[SFA_conserved$HOG %in% background,]
  
  #make the gene universe (HOGs_all) and mark everything that is of interest
  #init HOGs all with not of interest
  HOGs_all <- rep(0, length(unique(SFA_background$HOG[!is.na(SFA_background$HOG)])))
  #set HOGs as names
  names(HOGs_all) <- unique(SFA_background$HOG[!is.na(SFA_background$HOG)])
  #set values (0 = not interesting, 1 = interesting)
  HOGs_all[names(HOGs_all) %in% HOGs_interest] <- 1
  #convert to factor
  HOGs_all <- as.factor(HOGs_all)
  
  #object which includes enrichment test result of all three ontologies
  allRes <- list(BP = c(), MF = c(), CC = c())
  
  #go through all three ontologies
  for(ontology in c("BP", "MF", "CC"))
  {
    #make a topGodata object to analyse
    GOdata <- new("topGOdata", ontology = ontology, allGenes = HOGs_all,
                  annot = annFUN.gene2GO, gene2GO = HOG2GO)
    print(GOdata)
    #use the fishers weight01 function
    result_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    #make a results table
    res <- GenTable(GOdata, weight01 = result_weight01, orderBy = "weight01", ranksOf = "classic", topNodes = nTop)
    res <- res[res$Significant > 1,]
    res$Ontology <- rep(ontology, length(res$GO.ID))
    allRes[ontology] <- list(res)
  }
  
  allRes <- bind_rows(allRes$BP, allRes$MF, allRes$CC)
  
  allRes$weight01 <- as.numeric(allRes$weight01)
  
  return(allRes)
  
}


plot_GO_graph <- function(Terms, top)
{
  #combine ID with Term
  Terms$Term <- paste(Terms$GO.ID, Terms$Term)
  
  #only plot top n terms and filter out NA
  Terms <- Terms[1:top,]
  Terms <- Terms[!is.na(Terms$GO.ID),]
  
  if(dim(Terms)[1] == 0)
  {
    ggplot() + ggtitle("\tNo Data")
  }
  else
  {
    maxn <- max(Terms$Significant)
    Terms_order <- Terms[order(Terms$weight01, decreasing=T),]
    ggplot(Terms_order) + 
      geom_point(aes(x = Terms_order$weight01, y = fct_inorder(Terms_order$Term), size = Terms_order$Significant)) +
      geom_vline(xintercept=0.05) +
      geom_vline(xintercept=0.01) +
      geom_vline(xintercept=0.001) +
      scale_size_binned_area(name = "Significant\nAnnotated", limits = c(0, maxn),
                             breaks = seq(0, maxn, maxn/5)) +
      scale_x_reverse() +
      xlab("p-Value") +
      ylab("GO Term")
  }
  
}

plot_GO <- function(data, title, pdf, top = 10, width=10, height=10)
{
  bp <- plot_GO_graph(data[data$Ontology == "BP",], top = top)
  mf <- plot_GO_graph(data[data$Ontology == "MF",], top = top)
  cc <- plot_GO_graph(data[data$Ontology == "CC",], top = top)
  
  grid <- plot_grid(bp, mf, cc, labels = c("BP", "MF", "CC"), ncol=1)
  t <- ggdraw() + draw_label(title)
  plot <- plot_grid(t, grid, ncol=1, rel_heights=c(0.1, 1))
  ggsave(filename = pdf, plot = plot, width=width, height=height)
}

