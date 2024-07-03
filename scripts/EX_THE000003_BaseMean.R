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


# the other hogs


HOG0900 <- get_baseMean("N0.HOG0000900")

HOG0900_order <- c("Zm00001eb007650", "Misin01G446600","Misin02G422800",
                       "SORBI_3001G464900","Os03g0203700", "TraesCS4A02G046800", "TraesCS4B02G258000",
                       "TraesCS4D02G257900", "SORBI_3007G215200", "Zm00001eb199600", "Zm00001eb034200",
                       "Misin13G189000", "Misin07G224500")

HOG0900 <- HOG0900[match(HOG0900_order, HOG0900$gene),]

pdf(".\\N0.HOG0000900_basemean.pdf", 4)

ggplot(HOG0900, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG2654 <- get_baseMean("N0.HOG0002654")

HOG2654_order <- rev(c("Misin17G100300", "SORBI_3009G109925","Misin16G126600",
                       "SORBI_3009G110000","Misin17G100400", "Zm00001eb285800", "SORBI_3003G147000",
                       "Misin06G127400", "Misin05G138300", "Misin06G127600", "Zm00001eb339870",
                       "Os05g0358700", "TraesCS1B02G242500", "TraesCS1D02G230100", "TraesCS1A02G228300"))

HOG2654 <- HOG2654[match(HOG2654_order, HOG2654$gene),]

pdf(".\\N0.HOG0002654_basemean.pdf", 4)

ggplot(HOG2654, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG3362 <- get_baseMean("N0.HOG0003362")

HOG3362_order <- rev(c("TraesCS6A02G068900", "TraesCS6D02G066700","TraesCS6B02G093100",
                       "Os02g0121300","Misin07G021100", "SORBI_3004G018400", "Misin07G019900",
                       "SORBI_3007G164300", "MisinT275300", "Zm00001eb181260", "Misin08G149800",
                       "SORBI_3004G162000", "Zm00001eb242390", "Zm00001eb242400"))

HOG3362 <- HOG3362[match(HOG3362_order, HOG3362$gene),]

pdf(".\\N0.HOG0003362_basemean.pdf", 4)

ggplot(HOG3362, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()


HOG3607 <- get_baseMean("N0.HOG0003607")

HOG3607_order <- rev(c("Misin04G348100", "Misin03G331200","SORBI_3002G373100",
                       "Zm00001eb092910","Zm00001eb109110", "Zm00001eb113400", "TraesCS2A02G170200",
                       "TraesCS2B02G196500", "TraesCS2D02G177700", "Os07g0614000", "Misin04G348000",
                       "Misin03G331500", "Misin04G347900", "SORBI_3002G373200"))

HOG3607 <- HOG3607[match(HOG3607_order, HOG3607$gene),]

pdf(".\\N0.HOG0003607_basemean.pdf", 4)

ggplot(HOG3607, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG3895 <- get_baseMean("N0.HOG0003895")

HOG3895_order <- rev(c("TraesCS5D02G173300", "TraesCS5A02G169000","TraesCS5B02G165800",
                       "Os01g0218032","Zm00001eb241240", "Misin14G079500", "Misin15G097700",
                       "Misin15G097800", "SORBI_3008G085300", "Zm00001eb202980", "Misin07G374100",
                       "SORBI_3004G149800", "Misin08G138700"))

HOG3895 <- HOG3895[match(HOG3895_order, HOG3895$gene),]

pdf(".\\N0.HOG0003895_basemean.pdf", 4)

ggplot(HOG3895, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()


HOG4099 <- get_baseMean("N0.HOG0004099", c("TraesCS7D02G515300"))

HOG4099_order <- rev(c("TraesCS7A02G526800", "TraesCS7D02G515300","TraesCS7B02G444100",
                       "Zm00001eb067090","SORBI_3006G257100", "Misin11G258100", "Os04g0666800",
                       "Misin12G265100", "Zm00001eb433080", "Misin12G264800", "Misin11G258300",
                       "Misin12G264600", "SORBI_3005G159100", "SORBI_3006G256900"))

HOG4099 <- HOG4099[match(HOG4099_order, HOG4099$gene),]

pdf(".\\N0.HOG0004099_basemean.pdf", 4)

ggplot(HOG4099, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()


HOG4110 <- get_baseMean("N0.HOG0004110")

HOG4110_order <- rev(c("TraesCS2D02G316500", "TraesCS2A02G318800","TraesCS2B02G337100",
                       "Os04g0450600","Misin12G099800", "Zm00001eb423830", "Zm00001eb081290",
                       "Misin11G103600", "SORBI_3006G100332", "Misin12G100000", "Misin12G099900",
                       "SORBI_3006G100400", "Zm00001eb081280"))

HOG4110 <- HOG4110[match(HOG4110_order, HOG4110$gene),]

pdf(".\\N0.HOG0004110_basemean.pdf", 4)

ggplot(HOG4110, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG4459 <- get_baseMean("N0.HOG0004459")

HOG4459_order <- rev(c("Misin15G077700", "SORBI_3008G042800","Misin14G028300",
                       "Zm00001eb164080","TraesCS4D02G114100", "TraesCS4A02G199000", "TraesCS4B02G116400",
                       "Os11g0156800", "Zm00001eb197530", "Misin10G026800", "SORBI_3005G043000",
                       "Misin09G039700"))

HOG4459 <- HOG4459[match(HOG4459_order, HOG4459$gene),]

pdf(".\\N0.HOG0004459_basemean.pdf", 4)

ggplot(HOG4459, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG4512 <- get_baseMean("N0.HOG0004512")

HOG4512_order <- rev(c("TraesCS1D02G443100", "TraesCS1A02G433500","TraesCS1B02G469400",
                       "Os10g0497600","Misin16G254100", "Misin17G264000", "SORBI_3009G254000",
                       "Zm00001eb297460", "Zm00001eb095070", "Misin03G074000", "Misin04G086700",
                       "SORBI_3002G082100"))

HOG4512 <- HOG4512[match(HOG4512_order, HOG4512$gene),]

pdf(".\\N0.HOG0004512_basemean.pdf", 4)

ggplot(HOG4512, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG4560 <- get_baseMean("N0.HOG0004560")

HOG4560_order <- rev(c("Misin12G192300", "Misin11G184300","SORBI_3006G192200",
                       "Zm00001eb072840","TraesCS2A02G456300", "TraesCS2D02G456600", "TraesCS2B02G478400",
                       "Os02g0685600", "Zm00001eb250950", "SORBI_3004G279200", "Misin08G243600",
                       "Misin07G448100"))

HOG4560 <- HOG4560[match(HOG4560_order, HOG4560$gene),]

pdf(".\\N0.HOG0004560_basemean.pdf", 4)

ggplot(HOG4560, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG4908 <- get_baseMean("N0.HOG0004908")

HOG4908_order <- rev(c("Zm00001eb106630", "Zm00001eb323620","SORBI_3002G332750",
                       "Misin04G316100","Misin03G290500", "Misin03G290400", "TraesCS2B02G237200",
                       "TraesCS2A02G212100", "TraesCS2D02G218000", "Os07g0549700", "Misin07G544600",
                       "SORBI_3004G347600"))

HOG4908 <- HOG4908[match(HOG4908_order, HOG4908$gene),]

pdf(".\\N0.HOG0004908_basemean.pdf", 4)

ggplot(HOG4908, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG5130 <- get_baseMean("N0.HOG0005130", c("TraesCS1A02G291800"))

HOG5130 <- HOG5130[!is.na(HOG5130$gene),]

HOG5130_order <- rev(c("Misin16G168700", "Misin17G177900","SORBI_3009G167600",
                       "Zm00001eb289950","TraesCS1A02G291800", "TraesCS1B02G301200", "TraesCS1D02G290300",
                       "Os01g0833400", "Misin05G319300", "Zm00001eb368310", "SORBI_3003G345600",
                       "Misin06G312700"))

HOG5130 <- HOG5130[match(HOG5130_order, HOG5130$gene),]

pdf(".\\N0.HOG0005130_basemean.pdf", 4)

ggplot(HOG5130, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG5249 <- get_baseMean("N0.HOG0005249", c("TraesCS1A02G045700"))

HOG5249_order <- rev(c("Misin16G013900", "SORBI_3009G017732","Misin17G016900",
                       "Zm00001eb355550","TraesCS1B02G059100", "TraesCS1A02G045700", "Os05g0116100",
                       "Misin16G014000", "Misin17G016800", "SORBI_3009G017800", "Zm00001eb355540",
                       "Zm00001eb266260"))

HOG5249 <- HOG5249[match(HOG5249_order, HOG5249$gene),]

pdf(".\\N0.HOG0005249_basemean.pdf", 4)

ggplot(HOG5249, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()




HOG5521 <- get_baseMean("N0.HOG0005521")

HOG5521_order <- rev(c("MisinT054900", "Misin14G001900","MisinT103100",
                       "Misin14G086900","SORBI_3005G002000", "Os12g0562400", "TraesCS5D02G160300",
                       "Zm00001eb229100", "Zm00001eb229070", "Misin14G002300", "SORBI_3008G002000"))

HOG5521 <- HOG5521[match(HOG5521_order, HOG5521$gene),]

pdf(".\\N0.HOG0005521_basemean.pdf", 4)

ggplot(HOG5521, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()



HOG3710 <- get_baseMean("N0.HOG0003710")

HOG3710_order <- rev(c("TraesCS7D02G317400", "TraesCS7A02G320500","TraesCS7B02G221500",
                       "Os08g0110800","Zm00001eb260490", "Zm00001eb363030", "Zm00001eb174080",
                       "Misin06G230000", "SORBI_3003G276901", "Misin05G256600", "Misin07G211100",
                       "SORBI_3007G011600", "Misin13G005800"))

HOG3710 <- HOG3710[match(HOG3710_order, HOG3710$gene),]

pdf(".\\N0.HOG0003710_basemean.pdf", 4)

ggplot(HOG3710, aes(baseMeans, fct_inorder(gene))) +
  geom_bar(stat="identity", ) +
  xlab("Base Mean") +
  ylab("Gene")

dev.off()


