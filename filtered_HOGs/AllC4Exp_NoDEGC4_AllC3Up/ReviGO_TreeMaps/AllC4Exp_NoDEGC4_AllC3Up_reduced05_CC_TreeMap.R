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
revigo.data <- rbind(c("GO:0000502","proteasome complex",0.37652382533874423,1,0.7864587092017606,-0,"proteasome complex"),
                     c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
                     c("GO:0005737","cytoplasm",43.076660800468886,5,0.876265050477576,0.09354521,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,4,0.8431111418366491,0.17160779,"cytoplasm"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7902280744401448,0,"vacuole"),
                     c("GO:0005634","nucleus",16.5161752456724,3,0.7269108695242845,0.35145602,"vacuole"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7787985968197498,0.21968674,"vacuole"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7570212294648768,0.25147478,"vacuole"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7054070726892268,0.18309297,"vacuole"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7479319672799656,0.22703414,"vacuole"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.8041507518715392,0.18302188,"vacuole"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7761956513619133,0.21139748,"vacuole"),
                     c("GO:0005886","plasma membrane",17.177321395000487,2,0.9824062225654548,0.00015294,"plasma membrane"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.999929335396792,4.05E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,2,0.9998501802029992,7.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()


