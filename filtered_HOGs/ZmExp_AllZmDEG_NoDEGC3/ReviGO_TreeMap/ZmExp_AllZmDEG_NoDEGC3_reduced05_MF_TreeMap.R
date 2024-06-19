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
revigo.data <- rbind(c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,2,1,-0,"DNA-binding transcription factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.9468857898693089,0.06915977,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,1,0.9464323281816921,0.12712323,"transferase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8754614116022313,0.04503758,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.875440697222909,0.32372808,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8706072740378815,0.18924409,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.8508700477574961,0.1977699,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0016301","kinase activity",5.872540794447113,1,0.8253480268700466,0.40437782,"transferase activity, transferring alkyl or aryl (other than methyl) groups"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,1,0.9067805247798328,-0,"phosphatase activity"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.8004191068519949,0,"sequence-specific DNA binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.818817523666658,0.2421934,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,3,0.7681033053829566,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.8336257217278011,0.33073906,"sequence-specific DNA binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,2,0.8662630963641703,0.09221581,"metal ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,1,0.7616494213246765,0.38577556,"metal ion binding"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.92645908734501,0.10689728,"metal ion binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AllZmDEG_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, All Zm DE, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

