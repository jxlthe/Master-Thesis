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
revigo.data <- rbind(c("GO:0006952","defense response",1.1604144588756624,1,0.9420584978843916,-0,"defense response"),
                     c("GO:0002238","response to molecule of fungal origin",0.007682595099288114,1,0.9420584978843916,0.20582114,"defense response"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9671470320818433,0.04068139,"coumarin biosynthetic process"),
                     c("GO:0031349","positive regulation of defense response",0.11892617777855335,1,0.8148885201816622,-0,"positive regulation of defense response"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9103045072845678,0.00664115,"protein autophosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.90100445252002,0.2739859,"protein autophosphorylation"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9920523869125935,0,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
n_hogs = length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, DEG in any C4, No DEG in C3, ReviGO 0.5, BP, # of HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

