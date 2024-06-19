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
revigo.data <- rbind(c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.8780550412304704,0.08992384,"lipid metabolic process"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,1,-0,"muscle organ development"),
                     c("GO:0009834","plant-type secondary cell wall biogenesis",0.02117458277124934,1,0.9935816576963363,0.00689751,"plant-type secondary cell wall biogenesis"),
                     c("GO:0010047","fruit dehiscence",0.002082705440775892,1,1,-0,"fruit dehiscence"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,1,0.9446492347469846,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.45210426038827783,-0,"phosphatidylinositol dephosphorylation"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,1,0.5429517765871298,0.49264795,"phosphatidylinositol dephosphorylation"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.44722560840293274,0.46733961,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.7031332987012088,0.39551818,"phosphatidylinositol dephosphorylation"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,1,0,"mRNA transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AllZmDEG_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, All Zm DE, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

