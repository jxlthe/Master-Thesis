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
revigo.data <- rbind(c("GO:0004672","protein kinase activity",3.5248673905776853,1,0.6979658095374469,0.0474854,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,1,0.9528504320620206,0.04547388,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,1,0.9096080929897793,0.08079952,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,1,0.9185086455647052,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,1,0.9106066208180881,0.12712323,"hydrolase activity"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.8988589484428764,-0,"L-ascorbic acid binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,1,0.8091385920543784,0.38577556,"L-ascorbic acid binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,1,0.8674928216631242,0.15860651,"L-ascorbic acid binding"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8797825248633305,0.02725726,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.8517859791708552,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.8941469198425223,0.45929075,"dioxygenase activity"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9619044822091802,0,"pectin acetylesterase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_DEGAnyC4_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
n_hogs = length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, DEG in any C4, No DEG in C3, ReviGO 0.5, MF, # of HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

