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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003677","DNA binding",12.252016067246458,1,0.9175114432932019,0.07367259,"DNA binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9702347536865078,0.04505158,"calmodulin binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,2,0.9426878477320969,0.0830445,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9533385771708953,0.05997119,"lyase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.953905728294811,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.923513779183175,0.04301538,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9394411123873799,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8548448400924897,0.05781294,"ATP hydrolysis activity"),
                     c("GO:0005388","P-type calcium transporter activity",0.05516679695917812,1,0.9553370293207659,0.4743178,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.8780630243161062,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.8435989880930058,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,2,0.8257511538368354,0.38577556,"iron-sulfur cluster binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,2,0.872864697242199,0.23204379,"iron-sulfur cluster binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.8758674032257879,0.01849091,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.8942184621241172,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9821612378176074,0,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()