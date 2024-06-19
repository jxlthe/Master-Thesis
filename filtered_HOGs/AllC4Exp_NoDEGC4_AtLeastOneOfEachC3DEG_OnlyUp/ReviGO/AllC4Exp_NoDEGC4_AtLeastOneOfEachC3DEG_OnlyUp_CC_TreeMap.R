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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
                     c("GO:0005737","cytoplasm",43.076660800468886,5,0.8667092671986336,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,4,0.8395125649390568,0.17160779,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,2,0.9858907243946746,0.00015294,"plasma membrane"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.8886866252011267,-0,"riboflavin synthase complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,1,0.8020284431559179,0.34749127,"riboflavin synthase complex"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7416531767183827,0,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,5,0.6957121104143723,0.35145602,"chloroplast"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7295874156668745,0.21968674,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7361187877932611,0.25147478,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7735171760825509,0.18302188,"chloroplast"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7101535089856432,0.18309297,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7357334093526413,0.22703414,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,1,0.7870061907625301,0.18475413,"chloroplast"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7578537489540528,0.21139748,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999303881506661,3.782E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,4,0.999850748184091,6.66E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
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
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
# 
# 
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
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_CC_Table.tsv", sep = "\t")
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

