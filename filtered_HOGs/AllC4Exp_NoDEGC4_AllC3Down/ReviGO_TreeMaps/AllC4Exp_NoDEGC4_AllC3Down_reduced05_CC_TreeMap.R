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
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,1,0.9998480528425548,0,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8388123989769932,0.00024363,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,1,0.8751259862390451,0.17160779,"cytoplasm"),
                     c("GO:0043231","intracellular membrane-bounded organelle",29.853178741885195,1,0.7517918558885712,0.23467685,"cytoplasm"),
                     c("GO:0005886","plasma membrane",17.177321395000487,1,0.9997853522349387,0.00015294,"plasma membrane"),
                     c("GO:0016020","membrane",49.2542153160787,2,0.9997056601254689,9.537E-05,"membrane"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999114077929783,3.703E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()