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
revigo.data <- rbind(c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8788659907729516,0.07971789,"phosphatidic acid biosynthetic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9412418563513237,0.08924942,"biosynthetic process"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.8973081327946559,0.07593517,"L-ascorbic acid metabolic process"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8683715830175431,-0,"response to chemical"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,1,0.8684741202947966,0.3661968,"response to chemical"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8759140500361217,0.47274124,"response to chemical"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8227013881005394,0.46171648,"response to chemical"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8663195115509706,0.42391909,"response to chemical"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8440369243320586,0.36420456,"response to chemical"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8116479779945531,-0,"proteolysis involved in protein catabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8317136089843344,0.35126771,"proteolysis involved in protein catabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8332193434844108,0.38896703,"proteolysis involved in protein catabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,1,0.8316816435968096,0.47291952,"proteolysis involved in protein catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.917821196366319,0.12747018,"proteolysis involved in protein catabolic process"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9254120545412087,0.17097988,"proteolysis involved in protein catabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.8835428366968173,0.28805767,"proteolysis involved in protein catabolic process"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.983541002359366,-0,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.9651269603562943,0.12847604,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9612644040144976,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9124506377030752,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9573645219029069,0.35287398,"acquisition of seed longevity"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.8973712742484832,0,"proton transmembrane transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,1,0.89047721545884,0.4188028,"proton transmembrane transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8663934856093372,0.34990498,"proton transmembrane transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.8968533185711883,0.26647644,"proton transmembrane transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

