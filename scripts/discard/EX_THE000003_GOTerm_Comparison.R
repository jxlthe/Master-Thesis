library(VennDiagram)

#plot the # of GO-Terms that are enriched when you filter for DEG, Expansion and both
plot_venn <- function(title, deg_exp, exp, deg)
{
  #pass GO_ID list
  overlap_deg_vs_deg_exp <- table(deg_exp %in% deg)["TRUE"]
  overlap_exp_vs_deg_exp <- table(deg_exp %in% exp)["TRUE"]
  overlap_deg_vs_exp <- table(deg %in% exp)["TRUE"]
  overlap_three <- table(deg[deg %in% exp] %in% deg_exp)["TRUE"]
  
  if(is.na(overlap_deg_vs_deg_exp))overlap_deg_vs_deg_exp = 0
  
  if(is.na(overlap_exp_vs_deg_exp))overlap_exp_vs_deg_exp = 0
  
  if(is.na(overlap_deg_vs_exp))overlap_deg_vs_exp = 0
  
  if(is.na(overlap_three))overlap_three = 0
  
  plot.new();
  title(title);
  venn_sub1 <- draw.triple.venn(area1 = length(deg_exp),
                                area2 = length(deg),
                                area3 = length(exp),
                                n12 = overlap_deg_vs_deg_exp,
                                n23 = overlap_deg_vs_exp,
                                n13 = overlap_exp_vs_deg_exp,
                                n123 = overlap_three,
                                cex = 2,
                                #cat.pos = c(0, 180, 0),
                                #cat.dist = c(0.04, 0.11, 0.03),
                                cat.cex = 1.3,
                                cat.fontfamily = 1,
                                fill = c("red", "yellow", "blue"),
                                category = c("DEG\nand\nGene Family Expansion", "DEG", "Gene Family Expansion"));
  
}


#Zm vs Osj, Ta = a
{
#BP
a_BP_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_BP_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
a_BP_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_BP_Exp_GOEnrichment.tsv", sep="\t")
a_BP_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_BP_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
plot_venn("Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology BP",
          a_BP_exp_deg$GO.ID,
          a_BP_exp$GO.ID,
          a_BP_deg$GO.ID)
#MF
a_MF_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_MF_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
a_MF_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_MF_Exp_GOEnrichment.tsv", sep="\t")
a_MF_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_MF_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
plot_venn("Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology MF",
          a_MF_exp_deg$GO.ID,
          a_MF_exp$GO.ID,
          a_MF_deg$GO.ID)
#CC
a_CC_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_CC_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
a_CC_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_CC_Exp_GOEnrichment.tsv", sep="\t")
a_CC_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Zm_vs_Osj_Ta_CC_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
plot_venn("Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology CC",
          a_CC_exp_deg$GO.ID,
          a_CC_exp$GO.ID,
          a_CC_deg$GO.ID)
#BP + MF + CC
a_BP_MF_CC_exp_deg <- c(a_BP_exp_deg$GO.ID, a_MF_exp_deg$GO.ID, a_CC_exp_deg$GO.ID)
a_BP_MF_CC_exp <- c(a_BP_exp$GO.ID, a_MF_exp$GO.ID, a_CC_exp$GO.ID)
a_BP_MF_CC_deg <- c(a_BP_deg$GO.ID, a_MF_deg$GO.ID, a_CC_deg$GO.ID)
plot_venn("Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology BP + MF + CC",
          a_BP_MF_CC_exp_deg,
          a_BP_MF_CC_exp,
          a_BP_MF_CC_deg)
}


#Msi, Sb, Zm vs Osj, Ta = b
{
  #BP
  b_BP_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_BP_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
  b_BP_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_BP_Exp_GOEnrichment.tsv", sep="\t")
  b_BP_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_BP_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
  plot_venn("Msi Sb Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology BP",
            b_BP_exp_deg$GO.ID,
            b_BP_exp$GO.ID,
            b_BP_deg$GO.ID)
  #MF
  b_MF_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_MF_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
  b_MF_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_MF_Exp_GOEnrichment.tsv", sep="\t")
  b_MF_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_MF_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
  plot_venn("Msi Sb Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology MF",
            b_MF_exp_deg$GO.ID,
            b_MF_exp$GO.ID,
            b_MF_deg$GO.ID)
  #CC
  b_CC_exp_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_CC_Exp_DEG_AnySpecies_GOEnrichment.tsv", sep = "\t")
  b_CC_exp <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_CC_Exp_GOEnrichment.tsv", sep="\t")
  b_CC_deg <- read.csv("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\GoTermEnrichment\\Msi_Sb_Zm_vs_Osj_Ta_CC_DEG_AnySpecies_GOEnrichment.tsv", sep="\t")
  plot_venn("Msi Sb Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology CC",
            b_CC_exp_deg$GO.ID,
            b_CC_exp$GO.ID,
            b_CC_deg$GO.ID)
  #BP + MF + CC
  b_BP_MF_CC_exp_deg <- c(b_BP_exp_deg$GO.ID, b_MF_exp_deg$GO.ID, b_CC_exp_deg$GO.ID)
  b_BP_MF_CC_exp <- c(b_BP_exp$GO.ID, b_MF_exp$GO.ID, b_CC_exp$GO.ID)
  b_BP_MF_CC_deg <- c(b_BP_deg$GO.ID, b_MF_deg$GO.ID, b_CC_deg$GO.ID)
  plot_venn("Msi Sb Zm vs Osj, Ta\nGO-Term Enrichment in the sets Expanded + DEG, Expanded and DEG\nfor Ontology BP + MF + CC",
            b_BP_MF_CC_exp_deg,
            b_BP_MF_CC_exp,
            b_BP_MF_CC_deg)
  
}
