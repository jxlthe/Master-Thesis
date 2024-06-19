#preparation
setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\Exp_DEG_GO_Comparison")



#get unique GOIDs
calc_GO_ID_diff <- function(exp, deg, exp_deg)
{
  SFA_bound <- bind_rows(SFA$Oryza_sativa, SFA$Triticum_aestivum, SFA$Sorghum_bicolor, SFA$Miscanthus_sinensis, SFA$Zea_mays)
  
  exp_id <- SFA_bound[SFA_bound$HOG %in% exp,]$'Gene-Ontology-Term'
  exp_id <- unique(exp_id[!is.na(exp_id)])
  deg_id <- SFA_bound[SFA_bound$HOG %in% deg,]$'Gene-Ontology-Term'
  deg_id <- unique(deg_id[!is.na(deg_id)])
  exp_deg_id <- SFA_bound[SFA_bound$HOG %in% exp_deg,]$'Gene-Ontology-Term'
  exp_deg_id <- unique(exp_deg_id[!is.na(exp_deg_id)])
  diff <- setdiff(exp_deg_id, exp_id)
  diff <- setdiff(diff, deg_id)
  return(diff)
}

#actual data

#compare GO-Terms which are included in Zm expanded, Zm DEG and ZM Expanded and DEG
#in the conserved set

# Zm expanded, All Zm DEG, no DEG in C3

ZmExp_AllZmDEG_hogs = c()
ZmExp_hogs = c()
AllZmDEG_hogs = c()

x = 0

for(hog in conserved_HOGs)
{
  
  x = x + 1
  if(x %% 1000 == 1)
  {
    cat(x, "of", length(conserved_HOGs), "\n")
  }
  genes <- HOG_DE.a2tea[HOG_DE.a2tea$HOG == hog,]
  genes$gene <- substr(genes$gene, 1, 2)
  count <- table(genes$gene)
  
  expanded <- FALSE
  DEG <- FALSE
  
  #check for HOGs, that have Zm expanded
  if(count["Zm"] > 1)
  {
    expanded <- TRUE
  }
  
  #check for NA in Zm Log2FoldChanges
  if(sum(is.na(genes[genes$gene == "Zm",]$log2FoldChange)) == 0)
  {  #check for Zm that are all not nonsignificant
    if(sum(genes[genes$gene == "Zm",]$significant == "no") == 0)
    {
      DEG <- TRUE
    }
  }
  if(expanded)ZmExp_hogs <- append(ZmExp_hogs, hog)
  if(DEG)AllZmDEG_hogs <- append(AllZmDEG_hogs, hog)
  if(expanded && DEG)ZmExp_AllZmDEG_hogs <- append(ZmExp_AllZmDEG_hogs, hog)
}
#the Background is set to be the set ZmExp or AllZmDEG
ZmExp_AllZmDEG_GO_diff <- calc_GO_ID_diff(ZmExp_hogs, AllZmDEG_hogs, ZmExp_AllZmDEG_hogs)
ZmExp_AllZmDEG_BG <- unique(c(ZmExp_hogs, AllZmDEG_hogs))

#test for enrichment
ZmExp_vs_Conserved <- test_GO_enrich(ZmExp_hogs)
AllZmDEG_vs_Conserved <- test_GO_enrich(AllZmDEG_hogs)
ZmExp_AllZmDEG_vs_Conserved <- test_GO_enrich(ZmExp_AllZmDEG_hogs)
ZmExp_AllZmDEG_vs_BG <- test_GO_enrich(ZmExp_AllZmDEG_hogs, background = ZmExp_AllZmDEG_BG)
#write and plot
write.csv(ZmExp_AllZmDEG_GO_diff, file = ".\\ZmExp_AllZmDEG\\ZmExp_AllZmDEG_Diff.csv")
write.csv(ZmExp_vs_Conserved, file = ".\\ZmExp_AllZmDEG\\ZmExp_vs_Conserved.csv")
write.csv(AllZmDEG_vs_Conserved, file = ".\\ZmExp_AllZmDEG\\AllZmDEG_vs_Conserved.csv")
write.csv(ZmExp_AllZmDEG_vs_Conserved, file = ".\\ZmExp_AllZmDEG\\ZmExp_AllZmDEG_vs_Conserved.csv")
write.csv(ZmExp_AllZmDEG_vs_BG, file = ".\\ZmExp_AllZmDEG\\ZmExp_AllZmDEG_vs_BG.csv")

plot_GO(ZmExp_vs_Conserved, "Zm Expanded vs Conserved Set", ".\\ZmExp_AllZmDEG\\ZmExp_vs_Conserved.pdf")
plot_GO(AllZmDEG_vs_Conserved, "DEG in Zm vs conserved Set", ".\\ZmExp_AllZmDEG\\AllZmDEG_vs_Conserved.pdf")
plot_GO(ZmExp_AllZmDEG_vs_Conserved, "Zm Expanded and DEG in Zm vs Conserved Set", ".\\ZmExp_AllZmDEG\\ZmExp_AllZmDEG_vs_Conserved.pdf")
plot_GO(ZmExp_AllZmDEG_vs_BG, "Zm Expanded AND DEG in Zm vs Zm Expanded OR DEG in Zm", ".\\ZmExp_AllZmDEG\\ZmExp_AllZmDEG_vs_BG.pdf")



