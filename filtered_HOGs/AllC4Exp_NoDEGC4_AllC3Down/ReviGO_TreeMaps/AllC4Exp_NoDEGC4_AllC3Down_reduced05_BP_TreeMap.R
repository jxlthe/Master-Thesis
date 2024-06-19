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
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.9154451308330057,-0,"monoatomic ion transport"),
                     c("GO:0070588","calcium ion transmembrane transport",0.4224811119566813,1,0.8666153222276416,0.3228373,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.7015775537782857,0,"defense response"),
                     c("GO:0009416","response to light stimulus",0.2589924319660237,1,0.7854965966676558,0.27535725,"defense response"),
                     c("GO:0009733","response to auxin",0.1256056236300358,1,0.7978122908701059,0.25746396,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.675481479067557,0.37522238,"defense response"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,1,0.9411691623050732,0.07105297,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9332400755998402,0.12277613,"lipid catabolic process"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,1,-0,"gametophyte development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8441020893355486,-0,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9062375099417653,0.15925058,"DNA demethylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Down_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

