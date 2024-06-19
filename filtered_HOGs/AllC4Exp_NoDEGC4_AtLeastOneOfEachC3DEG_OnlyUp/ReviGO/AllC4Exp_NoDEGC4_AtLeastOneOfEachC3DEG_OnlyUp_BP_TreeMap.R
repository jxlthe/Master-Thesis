# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8780001626499354,0.07971789,"phosphatidic acid biosynthetic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9398500057751906,0.08924942,"biosynthetic process"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.88685504952729,0,"response to salt stress"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,1,0.880047965361647,0.41360326,"response to salt stress"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8374780888310137,0.46171648,"response to salt stress"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8799541303160944,0.27995591,"response to salt stress"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8777411113302784,0.47274124,"response to salt stress"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8571310177279547,0.36420456,"response to salt stress"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.8662253906003201,0.07593517,"L-ascorbic acid metabolic process"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8030098185215787,-0,"proteolysis involved in protein catabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8158383632249968,0.35126771,"proteolysis involved in protein catabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.823447421984105,0.38896703,"proteolysis involved in protein catabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,2,0.8185128493832967,0.47291952,"proteolysis involved in protein catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9162140529896023,0.12747018,"proteolysis involved in protein catabolic process"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.911256018331753,0.17097988,"proteolysis involved in protein catabolic process"),
                     c("GO:0016926","protein desumoylation",0.02830261133305275,1,0.8730304986610063,0.2699367,"proteolysis involved in protein catabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.8738742220129585,0.28805767,"proteolysis involved in protein catabolic process"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9850316312530654,-0,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.968271662148011,0.12847604,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9299120410949966,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9041585510539594,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9353176309391981,0.35287398,"acquisition of seed longevity"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.9058910008046408,-0,"proton transmembrane transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.9001714182300902,0.4188028,"proton transmembrane transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8775240803746142,0.34990498,"proton transmembrane transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.9059983457890353,0.26647644,"proton transmembrane transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_BP_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

