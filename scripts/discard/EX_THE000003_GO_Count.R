setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

#plot the GO_Count by reduced GO Terms(Revigo)

#Ms Sb Zm vs Ta Os, Expanded, DEG in any expanded species

library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)

plot_GO_count_graph <- function(Terms)
{
  if(dim(Terms)[1] == 0)
  {
    ggplot() + ggtitle("\tNo Data")
  }
  else
  {
    #order the terms by Total count and clustered
    Terms_order <- Terms[order(Terms$TOTAL_COUNT.Freq, decreasing=F),]
  
    clustered_color <- ifelse(Terms_order$Clustered == TRUE, "black", "red")
  
    #combine ID with term
    Terms_order$GOTERM <- paste(Terms_order$GOTERM, Terms_order[]$GOID)
    #only get interesting bits
    Terms_order <- Terms_order[,c("GOTERM", "Ta_COUNT.Freq", "Os_COUNT.Freq", "Ms_COUNT.Freq", "Sb_COUNT.Freq", "Zm_COUNT.Freq")]
    
    # melt to plot a stacked bar plot
    Terms_order <- melt(Terms_order, id.vars = "GOTERM")
    
    ggplot(Terms_order) + 
      geom_bar(aes(x = value , y = fct_inorder(GOTERM), fill = variable),
               stat = "identity") +
      xlab("Count") +
      ylab("GO-Terms") +
      scale_fill_discrete(name="GO-Term Counts", labels=c("Triticum aestivum", "Oryza sativa", "Miscanthus sinensis", "Sorghum bicolor", "Zea mays")) +
      theme(axis.text.y = element_text(colour = clustered_color))
  }
  
}

plot_GO_count <- function(data, title, pdf, width=10, height=40, rel_h = c(1, 1, 1))
{
  bp <- plot_GO_count_graph(data[data$ONTOLOGY == "BP",])
  mf <- plot_GO_count_graph(data[data$ONTOLOGY == "MF",])
  cc <- plot_GO_count_graph(data[data$ONTOLOGY == "CC",])
  
  grid <- plot_grid(bp, mf, cc, labels = c("BP", "MF", "CC"), ncol=1, rel_heights=rel_h)
  t <- ggdraw() + draw_label(title)
  plot <- plot_grid(t, grid, ncol=1, rel_heights=c(0.1, 1))
  ggsave(filename = pdf, plot = plot, width=width, height=height)
}

#function to move count data over to the reduced set
#count_data unreduced GOTerms
#count_data_reuced <- raw output of the ReviGO reduction
#count_table <- the GOTerms that the GOTerms were reduced to
get_count_data <- function(count_table, count_data_reduced, count_data)
{
  count_table$TOTAL_COUNT.Freq <- rep(0, length(count_table$TermID))
  count_table$Ta_COUNT.Freq <- rep(0, length(count_table$TermID))
  count_table$Os_COUNT.Freq <- rep(0, length(count_table$TermID))
  count_table$Ms_COUNT.Freq <- rep(0, length(count_table$TermID))
  count_table$Sb_COUNT.Freq <- rep(0, length(count_table$TermID))
  count_table$Zm_COUNT.Freq <- rep(0, length(count_table$TermID))
  
  
  #mark if clustered
  count_table$Clustered <- rep(FALSE, length(count_table$TermID))
  
  for(GOTERM_parent in count_table$TermID)
  {

    #set count of parents to count in table
    count_table[count_table$TermID == GOTERM_parent,]$TOTAL_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$TOTAL_COUNT.Freq
    count_table[count_table$TermID == GOTERM_parent,]$Ta_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$Ta_COUNT.Freq
    count_table[count_table$TermID == GOTERM_parent,]$Os_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$Os_COUNT.Freq
    count_table[count_table$TermID == GOTERM_parent,]$Ms_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$Ms_COUNT.Freq
    count_table[count_table$TermID == GOTERM_parent,]$Sb_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$Sb_COUNT.Freq
    count_table[count_table$TermID == GOTERM_parent,]$Zm_COUNT.Freq <- count_data[count_data$GOID == GOTERM_parent,]$Zm_COUNT.Freq
    
    
    #if children exists add them to the count
    GOTERM_parent_short <- str_remove(tail(str_split(GOTERM_parent, ":")[[1]], n = 1), "^0+")
    GOTERM_children <- count_data_reduced[count_data_reduced$Representative == GOTERM_parent_short,]$TermID
    
    if(length(GOTERM_children) > 0)
    {
      count_table[count_table$TermID == GOTERM_parent,]$TOTAL_COUNT.Freq <- 
        count_table[count_table$TermID == GOTERM_parent,]$TOTAL_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$TOTAL_COUNT.Freq)
      count_table[count_table$TermID == GOTERM_parent,]$Ta_COUNT.Freq <-
        count_table[count_table$TermID == GOTERM_parent,]$Ta_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$Ta_COUNT.Freq)
      count_table[count_table$TermID == GOTERM_parent,]$Os_COUNT.Freq <- 
        count_table[count_table$TermID == GOTERM_parent,]$Os_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$Os_COUNT.Freq)
      count_table[count_table$TermID == GOTERM_parent,]$Ms_COUNT.Freq <- 
        count_table[count_table$TermID == GOTERM_parent,]$Ms_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$Ms_COUNT.Freq)
      count_table[count_table$TermID == GOTERM_parent,]$Sb_COUNT.Freq <- 
        count_table[count_table$TermID == GOTERM_parent,]$Sb_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$Sb_COUNT.Freq)
      count_table[count_table$TermID == GOTERM_parent,]$Zm_COUNT.Freq <- 
        count_table[count_table$TermID == GOTERM_parent,]$Zm_COUNT.Freq + sum(count_data[count_data$GOID %in% GOTERM_children,]$Zm_COUNT.Freq)
      # mark as clustered
      count_table[count_table$TermID == GOTERM_parent,]$Clustered = TRUE
    }
}
  
  return(count_table)
  
}

#reduced GO Terms with ReviGO 0.4
MSZ_expanded_deg_any_exp_species_GO_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_unique_GO.csv")

#read data and add ontology column
MSZ_expanded_deg_any_exp_species_reduced04_GO_BP_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced04_GO_BP.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced04_GO_BP_data$ONTOLOGY <- rep("BP", length(MSZ_expanded_deg_any_exp_species_reduced04_GO_BP_data$TermID))
MSZ_expanded_deg_any_exp_species_reduced04_GO_CC_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced04_GO_CC.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced04_GO_CC_data$ONTOLOGY <- rep("CC", length(MSZ_expanded_deg_any_exp_species_reduced04_GO_CC_data$TermID))
MSZ_expanded_deg_any_exp_species_reduced04_GO_MF_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced04_GO_MF.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced04_GO_MF_data$ONTOLOGY <- rep("MF", length(MSZ_expanded_deg_any_exp_species_reduced04_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(MSZ_expanded_deg_any_exp_species_reduced04_GO_BP_data, MSZ_expanded_deg_any_exp_species_reduced04_GO_CC_data, MSZ_expanded_deg_any_exp_species_reduced04_GO_MF_data)
MSZ_expanded_deg_any_exp_species_reduced04_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
MSZ_expanded_deg_any_exp_species_reduced04_GO_data$GOID <- MSZ_expanded_deg_any_exp_species_reduced04_GO_data$TermID
MSZ_expanded_deg_any_exp_species_reduced04_GO_data$GOTERM <- MSZ_expanded_deg_any_exp_species_reduced04_GO_data$Name

#make a count table

MSZ_expanded_deg_any_exp_species_reduced04_GO_count_table <- MSZ_expanded_deg_any_exp_species_reduced04_GO_data[MSZ_expanded_deg_any_exp_species_reduced04_GO_data$Representative == "null",]


MSZ_expanded_deg_any_exp_species_reduced04_GO_count_table <- get_count_data(MSZ_expanded_deg_any_exp_species_reduced04_GO_count_table,
                                                                            MSZ_expanded_deg_any_exp_species_reduced04_GO_data,
                                                                            MSZ_expanded_deg_any_exp_species_GO_data)
plot_GO_count(MSZ_expanded_deg_any_exp_species_reduced04_GO_count_table, "Overall Reduced (0.4) GO-Term Count\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\GOCount_Plots\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_reduced04_GO_count.pdf", width = 20, height = 40, rel_h = c(1, 2, 1))



#reduced GO Terms with ReviGO 0.5

#read data and add ontology column
MSZ_expanded_deg_any_exp_species_reduced05_GO_BP_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced05_GO_BP.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced05_GO_BP_data$ONTOLOGY <- rep("BP", length(MSZ_expanded_deg_any_exp_species_reduced05_GO_BP_data$TermID))
MSZ_expanded_deg_any_exp_species_reduced05_GO_CC_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced05_GO_CC.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced05_GO_CC_data$ONTOLOGY <- rep("CC", length(MSZ_expanded_deg_any_exp_species_reduced05_GO_CC_data$TermID))
MSZ_expanded_deg_any_exp_species_reduced05_GO_MF_data <- read.csv(".\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\raw\\MSZ_expanded_deg_any_exp_species_reduced05_GO_MF.tsv", sep = "\t")
MSZ_expanded_deg_any_exp_species_reduced05_GO_MF_data$ONTOLOGY <- rep("MF", length(MSZ_expanded_deg_any_exp_species_reduced05_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(MSZ_expanded_deg_any_exp_species_reduced05_GO_BP_data, MSZ_expanded_deg_any_exp_species_reduced05_GO_CC_data, MSZ_expanded_deg_any_exp_species_reduced05_GO_MF_data)
MSZ_expanded_deg_any_exp_species_reduced05_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
MSZ_expanded_deg_any_exp_species_reduced05_GO_data$GOID <- MSZ_expanded_deg_any_exp_species_reduced05_GO_data$TermID
MSZ_expanded_deg_any_exp_species_reduced05_GO_data$GOTERM <- MSZ_expanded_deg_any_exp_species_reduced05_GO_data$Name

#make a count table

MSZ_expanded_deg_any_exp_species_reduced05_GO_count_table <- MSZ_expanded_deg_any_exp_species_reduced05_GO_data[MSZ_expanded_deg_any_exp_species_reduced05_GO_data$Representative == "null",]

MSZ_expanded_deg_any_exp_species_reduced05_GO_count_table <- get_count_data(MSZ_expanded_deg_any_exp_species_reduced05_GO_count_table,
                                                                            MSZ_expanded_deg_any_exp_species_reduced05_GO_data,
                                                                            MSZ_expanded_deg_any_exp_species_GO_data)
plot_GO_count(MSZ_expanded_deg_any_exp_species_reduced05_GO_count_table, "Overall Reduced (0.5) GO-Term Count\nC4 vs C3, Expanded, DEG in any expanded species", "MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\GOCount_Plots\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_reduced05_GO_count.pdf", width = 20, height = 40, rel_h = c(1, 2, 1))






#All C4 Exp No C4 DEG All C3 Up
AllC4Exp_NoDEGC4_AllC3Up_GO_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_unique_GO.csv")


#reduced GO Terms with ReviGO 0.5

#read data and add ontology column
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_BP_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_BP.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_BP_data$ONTOLOGY <- rep("BP", length(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_BP_data$TermID))
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_CC_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_CC.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_CC_data$ONTOLOGY <- rep("CC", length(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_CC_data$TermID))
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_MF_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_MF.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_MF_data$ONTOLOGY <- rep("MF", length(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_BP_data, AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_CC_data, AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_MF_data)
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data$GOID <- AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data$TermID
AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data$GOTERM <- AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data$Name

#make a count table

AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_count_table <- AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data[AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data$Representative == "null",]

AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_count_table <- get_count_data(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_count_table,
                                                                    AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_data,
                                                                    AllC4Exp_NoDEGC4_AllC3Up_GO_data)
plot_GO_count(AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_count_table, "Overall Reduced (0.5) GO-Term Count\nAll C4 expanded, No C4 DEG, All C3 upregulated", "AllC4Exp_NoDEGC4_AllC3Up\\GOCount_Plots\\AllC4Exp_NoDEGC4_AllC3Up_reduced05_GO_count.pdf", width = 20, height = 40)


#reduced GO Terms with ReviGO 0.7

#read data and add ontology column
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_BP_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_BP.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_BP_data$ONTOLOGY <- rep("BP", length(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_BP_data$TermID))
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_CC_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_CC.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_CC_data$ONTOLOGY <- rep("CC", length(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_CC_data$TermID))
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_MF_data <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3Up\\raw\\AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_MF.tsv", sep = "\t")
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_MF_data$ONTOLOGY <- rep("MF", length(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_BP_data, AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_CC_data, AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_MF_data)
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data$GOID <- AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data$TermID
AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data$GOTERM <- AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data$Name

#make a count table

AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_count_table <- AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data[AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data$Representative == "null",]

AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_count_table <- get_count_data(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_count_table,
                                                                    AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_data,
                                                                    AllC4Exp_NoDEGC4_AllC3Up_GO_data)
plot_GO_count(AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_count_table, "Overall Reduced (0.7) GO-Term Count\nAll C4 expanded, No C4 DEG, All C3 upregulated", "AllC4Exp_NoDEGC4_AllC3Up\\GOCount_Plots\\AllC4Exp_NoDEGC4_AllC3Up_reduced07_GO_count.pdf", width = 20, height = 40)





#NoDEGC3 All C4 Exp DEGAnyC4


AllC4Exp_DEGAnyC4_NoDEGC3_GO_data <- read.csv(".\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_unique_GO.csv")

#reduced GO Terms with ReviGO 0.5

#read data and add ontology column
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_BP_data <- read.csv(".\\AllC4Exp_DEGAnyC4_NoDEGC3\\raw\\AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_BP.tsv", sep = "\t")
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_BP_data$ONTOLOGY <- rep("BP", length(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_BP_data$TermID))
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_CC_data <- read.csv(".\\AllC4Exp_DEGAnyC4_NoDEGC3\\raw\\AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_CC.tsv", sep = "\t")
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_CC_data$ONTOLOGY <- rep("CC", length(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_CC_data$TermID))
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_MF_data <- read.csv(".\\AllC4Exp_DEGAnyC4_NoDEGC3\\raw\\AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_MF.tsv", sep = "\t")
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_MF_data$ONTOLOGY <- rep("MF", length(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_BP_data, AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_CC_data, AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_MF_data)
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data$GOID <- AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data$TermID
AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data$GOTERM <- AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data$Name

#make a count table

AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_count_table <- AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data[AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data$Representative == "null",]

AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_count_table <- get_count_data(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_count_table,
                                                                    AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_data,
                                                                    AllC4Exp_DEGAnyC4_NoDEGC3_GO_data)
plot_GO_count(AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_count_table, "Overall Reduced (0.5) GO-Term Count\nAll C4 Expanded, Any C4 DEG, No DEG C3", "AllC4Exp_DEGAnyC4_NoDEGC3\\GOCount_Plots\\AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_GO_count.pdf", width = 20, height = 40)




#Zm expanded, no DEG in C3, only any Zm DEG


ZmExp_OnlyAnyZmDEG_NoDEGC3_GO_data <- read.csv("ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_unique_GO.csv")
#reduced GO Terms with ReviGO 0.5

#read data and add ontology column
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_BP_data <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\raw\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_BP.tsv", sep = "\t")
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_BP_data$ONTOLOGY <- rep("BP", length(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_BP_data$TermID))
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_CC_data <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\raw\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_CC.tsv", sep = "\t")
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_CC_data$ONTOLOGY <- rep("CC", length(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_CC_data$TermID))
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_MF_data <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\raw\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_MF.tsv", sep = "\t")
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_MF_data$ONTOLOGY <- rep("MF", length(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_BP_data, ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_CC_data, ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_MF_data)
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data$GOID <- ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data$TermID
ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data$GOTERM <- ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data$Name

#make a count table

ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_count_table <- ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data[ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data$Representative == "null",]

ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_count_table <- get_count_data(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_count_table,
                                                                     ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_data,
                                                                     ZmExp_OnlyAnyZmDEG_NoDEGC3_GO_data)
plot_GO_count(ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_count_table, "Overall Reduced (0.5) GO-Term Count\nZm expanded, Only Any Zm DEG, No C3 DEG", "ZmExp_OnlyAnyZmDEG_NoDEGC3\\GOCount_Plots\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_GO_count.pdf", width = 20, height = 40)


#Zm expanded, no DEG in C3, All Zm DEG


ZmExp_AllZmDEG_NoDEGC3_GO_data <- read.csv("ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_unique_GO.csv")
#reduced GO Terms with ReviGO 0.5

#read data and add ontology column
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_BP_data <- read.csv(".\\ZmExp_AllZmDEG_NoDEGC3\\raw\\ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_BP.tsv", sep = "\t")
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_BP_data$ONTOLOGY <- rep("BP", length(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_BP_data$TermID))
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_CC_data <- read.csv(".\\ZmExp_AllZmDEG_NoDEGC3\\raw\\ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_CC.tsv", sep = "\t")
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_CC_data$ONTOLOGY <- rep("CC", length(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_CC_data$TermID))
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_MF_data <- read.csv(".\\ZmExp_AllZmDEG_NoDEGC3\\raw\\ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_MF.tsv", sep = "\t")
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_MF_data$ONTOLOGY <- rep("MF", length(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_MF_data$TermID))
#and combine into one single dataframe
df_list <- list(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_BP_data, ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_CC_data, ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_MF_data)
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
#make new renamed columns. You could also rename 
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data$GOID <- ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data$TermID
ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data$GOTERM <- ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data$Name

#make a count table

ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_count_table <- ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data[ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data$Representative == "null",]

ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_count_table <- get_count_data(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_count_table,
                                                                      ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_data,
                                                                      ZmExp_AllZmDEG_NoDEGC3_GO_data)
plot_GO_count(ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_count_table, "Overall Reduced (0.5) GO-Term Count\nZm expanded, All Zm DEG, No C3 DEG", "ZmExp_AllZmDEG_NoDEGC3\\GOCount_Plots\\ZmExp_AllZmDEG_NoDEGC3_reduced05_GO_count.pdf", width = 20, height = 40)



