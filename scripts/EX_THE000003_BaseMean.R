library("DESeq2")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

load("EX_THE000003_A2TEA_finished.RData")

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\dea")

#load dea data



which_hypotheses <- function(hog)
{
  
  print(paste(hog, "is in:"))
  
  if(HOG_level_list$hypothesis_1[which(HOG_level_list$hypothesis_1$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 1")
  
  if(HOG_level_list$hypothesis_2[which(HOG_level_list$hypothesis_2$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 2")
  
  if(HOG_level_list$hypothesis_3[which(HOG_level_list$hypothesis_3$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 3")
  
  if(HOG_level_list$hypothesis_4[which(HOG_level_list$hypothesis_4$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 4")
  
  if(HOG_level_list$hypothesis_5[which(HOG_level_list$hypothesis_5$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 5")
  
  if(HOG_level_list$hypothesis_6[which(HOG_level_list$hypothesis_6$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 6")
  
  if(HOG_level_list$hypothesis_7[which(HOG_level_list$hypothesis_7$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 7")
  
  if(HOG_level_list$hypothesis_8[which(HOG_level_list$hypothesis_8$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 8")
  
  if(HOG_level_list$hypothesis_9[which(HOG_level_list$hypothesis_9$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 9")
  
  if(HOG_level_list$hypothesis_10[which(HOG_level_list$hypothesis_10$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 10")
  
  if(HOG_level_list$hypothesis_11[which(HOG_level_list$hypothesis_11$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 11")
  
  if(HOG_level_list$hypothesis_12[which(HOG_level_list$hypothesis_12$HOG %in% hog),]$expansion == "yes")
    print("hypotheses 12")
  
}


dea_Ta <- readRDS("dea_Triticum_aestivum")
dea_Os <- readRDS("dea_Oryza_sativa")
dea_Ms <- readRDS("dea_Miscanthus_sinensis")
dea_Sb <- readRDS("dea_Sorghum_bicolor")
dea_Zm <- readRDS("dea_Zea_mays")

dea_Ta_counts <- counts(dea_Ta)
dea_Os_counts <- counts(dea_Os)
dea_Ms_counts <- counts(dea_Ms)
dea_Sb_counts <- counts(dea_Sb)
dea_Zm_counts <- counts(dea_Zm)

reformat_col_names <- function(dea_table)
{
  names <- colnames(dea_table)
  names[grepl("control", names)] <- "control"
  names[grepl("drought", names)] <- "drought"
  colnames(dea_table) <- names
  
  return(dea_table)
}


dea_Ta_counts <- reformat_col_names(dea_Ta_counts)
dea_Os_counts <- reformat_col_names(dea_Os_counts)
dea_Ms_counts <- reformat_col_names(dea_Ms_counts)
dea_Sb_counts <- reformat_col_names(dea_Sb_counts)
dea_Zm_counts <- reformat_col_names(dea_Zm_counts)

Ms_genes <- rownames(dea_Ms_counts)
print("Misin03G292100" %in% Ms_genes)

get_baseMean <- function(HOG, add_genes = c())
{
  
  print(which_hypotheses(HOG))
  #get base means
  baseMeans = HOG_DE.a2tea[HOG_DE.a2tea$HOG == HOG,]
  gene_list = baseMeans$gene
  baseMean_list = baseMeans$baseMean
  control_m = c()
  control_sd = c()
  drought_m = c()
  drought_sd = c()
  for(g in add_genes)
  {
    gene_list <- append(gene_list, g)
    baseMean_list <- append(baseMean_list, HOG_DE.a2tea[HOG_DE.a2tea$gene == g,]$baseMean)
  }
  
  #get count means
  for(g in gene_list)
  {
    gene_calc = c()    
    if(g %in% rownames(dea_Ta_counts))gene_calc = dea_Ta_counts[g,]
    if(g %in% rownames(dea_Os_counts))gene_calc = dea_Os_counts[g,]
    if(g %in% rownames(dea_Ms_counts))gene_calc = dea_Ms_counts[g,]
    if(g %in% rownames(dea_Sb_counts))gene_calc = dea_Sb_counts[g,]
    if(g %in% rownames(dea_Zm_counts))gene_calc = dea_Zm_counts[g,]
    control_m <- append(control_m, mean(gene_calc[grepl("control", names(gene_calc))]))
    drought_m <- append(drought_m, mean(gene_calc[grepl("drought", names(gene_calc))]))
  }
  
  df <- data.frame(gene = gene_list, baseMeans = baseMean_list, control_m = control_m, drought_m = drought_m)
  
  return(df)  
}

#Msi Sb Zm vs Osj Ta

print(get_baseMean("N0.HOG0001331"))

print(get_baseMean("N0.HOG0002025"))

print(get_baseMean("N0.HOG0004017"))

print(get_baseMean("N0.HOG0004343", c("TraesCS2B02G215000", "TraesCS2B02G215100", "TraesCS2D02G195900")))

print(get_baseMean("N0.HOG0004368"))

#Zm vs Osj Ta

print(get_baseMean("N0.HOG0006824"))

print(get_baseMean("N0.HOG0017121", c("TraesCS1A02G209100")))

print(get_baseMean("N0.HOG0004652"))

print(get_baseMean("N0.HOG0000969"))

print(get_baseMean("N0.HOG0006817"))

print(get_baseMean("N0.HOG0008319"))

print(get_baseMean("N0.HOG0011584", c("TraesCS4A02G348100")))

print(get_baseMean("N0.HOG0000951"))

print(get_baseMean("N0.HOG0013739", c("TraesCS2B02G407100")))

print(get_baseMean("N0.HOG0005274"))

print(get_baseMean("N0.HOG0017487", c("TraesCS2D02G430700", "TraesCS2B02G454000", "Misin17G072400")))

print(get_baseMean("N0.HOG0006053"))

print(get_baseMean("N0.HOG0006529"))

print(get_baseMean("N0.HOG0009376"))

print(get_baseMean("N0.HOG0009828"))

print(get_baseMean("N0.HOG0010839", c("TraesCS3B02G469000")))

print(get_baseMean("N0.HOG0015605", c("TraesCS3A02G414400")))

print(get_baseMean("N0.HOG0007568"))

print(get_baseMean("N0.HOG0008285"))

print(get_baseMean("N0.HOG0009114"))

print(get_baseMean("N0.HOG0009637"))

print(get_baseMean("N0.HOG0014378"))

print(get_baseMean("N0.HOG0003297"))

print(get_baseMean("N0.HOG0008731"))

print(get_baseMean("N0.HOG0009150"))

print(get_baseMean("N0.HOG0011281"))

print(get_baseMean("N0.HOG0014437", c("TraesCS2A02G189600")))

print(get_baseMean("N0.HOG0006225"))

print(get_baseMean("N0.HOG0003381"))

print(get_baseMean("N0.HOG0004652"))

print(get_baseMean("N0.HOG0008146"))

print(get_baseMean("N0.HOG0008294"))

print(get_baseMean("N0.HOG0008863"))

print(get_baseMean("N0.HOG0009387"))

print(get_baseMean("N0.HOG0009416"))

print(get_baseMean("N0.HOG0008725"))

print(get_baseMean("N0.HOG0009290"))

print(get_baseMean("N0.HOG0004356"))

print(get_baseMean("N0.HOG0008535", c("TraesCS4B02G246900")))

print(get_baseMean("N0.HOG0005967"))

print(get_baseMean("N0.HOG0006609"))

print(get_baseMean("N0.HOG0016008", c("TraesCS1A02G083100", "TraesCS1B02G100600")))

print(get_baseMean("N0.HOG0006352"))

print(get_baseMean("N0.HOG0006529"))

print(get_baseMean("N0.HOG0007767"))

print(get_baseMean("N0.HOG0008263"))

print(get_baseMean("N0.HOG0008325"))

print(get_baseMean("N0.HOG0008383"))

print(get_baseMean("N0.HOG0008411"))

print(get_baseMean("N0.HOG0008493"))

print(get_baseMean("N0.HOG0009688"))

print(get_baseMean("N0.HOG0009883"))

print(get_baseMean("N0.HOG0009913"))

print(get_baseMean("N0.HOG0013207"))

print(get_baseMean("N0.HOG0014865"))

print(get_baseMean("N0.HOG0003181", c("TraesCS1A02G363700")))

print(get_baseMean("N0.HOG0006579"))

print(get_baseMean("N0.HOG0006996"))

print(get_baseMean("N0.HOG0001890"))

print(get_baseMean("N0.HOG0007568"))

print(get_baseMean("N0.HOG0008285"))

print(get_baseMean("N0.HOG0008440"))

print(get_baseMean("N0.HOG0009375"))

print(get_baseMean("N0.HOG0009637"))

print(get_baseMean("N0.HOG0009753"))

print(get_baseMean("N0.HOG0009914"))

print(get_baseMean("N0.HOG0011233"))

print(get_baseMean("N0.HOG0014378"))

print(get_baseMean("N0.HOG0014868"))

print(get_baseMean("N0.HOG0015144"))

print(get_baseMean("N0.HOG0016753", c("TraesCS4D02G191200", "TraesCS4B02G189800")))

print(get_baseMean("N0.HOG0004918"))

print(get_baseMean("N0.HOG0007725"))

print(get_baseMean("N0.HOG0012401"))

print(get_baseMean("N0.HOG0008768"))

print(get_baseMean("N0.HOG0008855"))

print(get_baseMean("N0.HOG0002218"))

print(get_baseMean("N0.HOG0008853"))

print(get_baseMean("N0.HOG0000491"))

print(get_baseMean("N0.HOG0004383"))

print(get_baseMean("N0.HOG0008974"))

print(get_baseMean("N0.HOG0013697"))

print(get_baseMean("N0.HOG0015040"))

print(get_baseMean("N0.HOG0002752"))

print(get_baseMean("N0.HOG0006546"))

print(get_baseMean("N0.HOG0007531"))

print(get_baseMean("N0.HOG0005085"))

print(get_baseMean("N0.HOG0005684"))

print(get_baseMean("N0.HOG0002005"))

print(get_baseMean("N0.HOG0009749"))  

print(get_baseMean("N0.HOG0002889", c("TraesCS5D02G052300", "TraesCS5B02G047200")))

print(get_baseMean("N0.HOG0007890"))

print(get_baseMean("N0.HOG0007925"))

print(get_baseMean("N0.HOG0008330"))

print(get_baseMean("N0.HOG0008381"))

print(get_baseMean("N0.HOG0008550"))

print(get_baseMean("N0.HOG0009033"))

print(get_baseMean("N0.HOG0008814"))



#Only all C3 upregulated all C4 expanded
print(get_baseMean("N0.HOG0004685"))

print(get_baseMean("N0.HOG0003448"))

print(get_baseMean("N0.HOG0003994"))

print(get_baseMean("N0.HOG0004908"))

print(get_baseMean("N0.HOG0003607"))

print(get_baseMean("N0.HOG0002654"))

print(get_baseMean("N0.HOG0003137"))

print(get_baseMean("N0.HOG0003710"))

print(get_baseMean("N0.HOG0003362"))

print(get_baseMean("N0.HOG0003870"))

print(get_baseMean("N0.HOG0004459"))

print(get_baseMean("N0.HOG0005632"))

print(get_baseMean("N0.HOG0001278"))

print(get_baseMean("N0.HOG0001572"))

print(get_baseMean("N0.HOG0003297"))

print(get_baseMean("N0.HOG0000135"))

print(get_baseMean("N0.HOG0000208"))

print(get_baseMean("N0.HOG0000340", c("TraesCS7D02G427900")))

print(get_baseMean("N0.HOG0000362"))

print(get_baseMean("N0.HOG0000585"))

print(get_baseMean("N0.HOG0000995"))

print(get_baseMean("N0.HOG0001045"))

print(get_baseMean("N0.HOG0001331"))

print(get_baseMean("N0.HOG0011860", c("Misin07G493400")))

print(get_baseMean("N0.HOG0008923"))

print(get_baseMean("N0.HOG0007510"))

print(get_baseMean("N0.HOG0007675"))

print(get_baseMean("N0.HOG0005262"))

print(get_baseMean("N0.HOG0014219", c("TraesCS7D02G384700")))

print(get_baseMean("N0.HOG0002869"))

print(get_baseMean("N0.HOG0004931", c("TraesCS1B02G141900")))

print(get_baseMean("N0.HOG0004908"))

print(get_baseMean("N0.HOG0000900"))

print(get_baseMean("N0.HOG0004512"))

print(get_baseMean("N0.HOG0003607"))

print(get_baseMean("N0.HOG0005093"))

print(get_baseMean("N0.HOG0016896", c("TraesCS4B02G052000", "TraesCS4A02G262900")))

print(get_baseMean("N0.HOG0002654"))

print(get_baseMean("N0.HOG0003137"))

print(get_baseMean("N0.HOG0005130", c("TraesCS1A02G291800")))

print(get_baseMean("N0.HOG0003710"))

print(get_baseMean("N0.HOG0003362"))

print(get_baseMean("N0.HOG0004560"))

print(get_baseMean("N0.HOG0003895"))

print(get_baseMean("N0.HOG0004459"))

print(get_baseMean("N0.HOG0004110"))

print(get_baseMean("N0.HOG0004099", c("TraesCS7D02G515300")))

print(get_baseMean("N0.HOG0005521"))

print(get_baseMean("N0.HOG0005249", c("TraesCS1A02G045700")))

print(get_baseMean("N0.HOG0004139"))

print(get_baseMean("N0.HOG0005093"))

print(get_baseMean("N0.HOG0005521"))


#Plot basemeans

library(ggplot2)
library(forcats)

setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\basemean")


HOG2869 <- get_baseMean("N0.HOG0002869")

HOG2869_order <- rev(c("TraesCS4A02G259700", "TraesCS4B02G055000","TraesCS4D02G055200",
                   "Os03g0724100","SORBI_3001G108700", "Misin01G098700", "Misin02G086400",
                   "SORBI_3001G108600", "Misin01G098600", "Misin02G086300", "Zm00001eb055720",
                   "Zm00001eb055730", "SORBI_3001G108500",
                   "Misin02G086200", "Misin01G098500"))

HOG2869 <- HOG2869[match(HOG2869_order, HOG2869$gene),]

pdf(".\\N0.HOG0002869_basemean.pdf", 4)

ggplot(HOG2869, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG4931 <- get_baseMean("N0.HOG0004931", c("TraesCS1B02G141900"))

HOG4931_order <- rev(c("TraesCS1D02G123600", "TraesCS1B02G141900","TraesCS1A02G122700",
                       "Os05g0187000","SORBI_3003G035601", "Zm00001eb304090", "Misin04G092000",
                       "SORBI_3002G087000", "Misin03G070000", "Zm00001eb268750", "Misin18G002100",
                       "SORBI_3010G007400", "Misin18G001300"))

HOG4931 <- HOG4931[match(HOG4931_order, HOG4931$gene),]

pdf(".\\N0.HOG0004931_basemean.pdf", 4)

ggplot(HOG4931, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG3137 <- get_baseMean("N0.HOG0003137")

HOG3137_order <- rev(c("TraesCS1A02G224800", "TraesCS1D02G226200","TraesCS1B02G238000",
                       "Os05g0349500","Misin06G248400", "Misin05G243200", "SORBI_3003G267800",
                       "Zm00001eb155830", "Misin06G274900", "Misin17G046300", "SORBI_3009G046700",
                       "Misin16G036400", "Zm00001eb360880", "Misin05G283600"))

HOG3137 <- HOG3137[match(HOG3137_order, HOG3137$gene),]

pdf(".\\N0.HOG0003137_basemean.pdf", 4)

ggplot(HOG3137, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG4139 <- get_baseMean("N0.HOG0004139")

HOG4139_order <- rev(c("SORBI_3010G185400", "Zm00001eb386290","Misin18G183800",
                       "Misin19G177000","Os09g0249700", "TraesCS5D02G172700", "TraesCS5B02G165200",
                       "TraesCS5A02G168400", "Zm00001eb269610", "Zm00001eb378360", "SORBI_3010G020400",
                       "Misin18G014700", "Misin19G020900"))

HOG4139 <- HOG4139[match(HOG4139_order, HOG4139$gene),]

pdf(".\\N0.HOG0004139_basemean.pdf", 4)

ggplot(HOG4139, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()





