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
revigo.data <- rbind(c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9203079142606666,0.03137215,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.8942300932197402,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9121037444404009,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9624906396952049,0.0386588,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,1,0.9507531512872165,0.05532665,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,3,0.9342685727484303,0.07494256,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9402160380017902,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,1,0.9349396667078802,0.12712323,"hydrolase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9694685659422415,0.02287106,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.950875340135772,0.03319032,"isomerase activity"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9074217519528892,-0,"glutathione binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9035224028271102,0.31015716,"glutathione binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,1,0.8641751877358295,0.26746377,"glutathione binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.904150293886289,0.12083997,"glutathione binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.8800172737207748,0.02352266,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0046961","proton-transporting ATPase activity, rotational mechanism",0.20284297713014954,1,0.9592398403830547,-0,"proton-transporting ATPase activity, rotational mechanism"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,1,0.8665610802126985,0,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.9143662058292141,0.42844795,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9001839323705105,0.30149286,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8878879162354781,0.34866016,"N-acylphosphatidylethanolamine-specific phospholipase D activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,1,0.8424252332778702,0.20765779,"N-acylphosphatidylethanolamine-specific phospholipase D activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3Up_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

