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
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998642200963308,4.927E-05,"extracellular region"),
                     c("GO:0005730","nucleolus",1.2114693153196516,1,0.6971843266754033,3.325E-05,"nucleolus"),
                     c("GO:0000793","condensed chromosome",0.4223503377863461,1,0.7422472078416456,0.41377726,"nucleolus"),
                     c("GO:0005634","nucleus",16.5161752456724,3,0.7134999513385691,0.2729216,"nucleolus"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8347985501739508,0.09158706,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,1,0.8691419232895021,0.17160779,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,4,0.9998131780536398,0.00015294,"plasma membrane"),
                     c("GO:0016020","membrane",49.2542153160787,5,0.9997497606099126,9.537E-05,"membrane"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999181436783204,0,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
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
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_CC_Table.tsv", sep = "\t")
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

