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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,1,1,-0,"cellular_component"),
c("GO:0005576","extracellular region",3.0008202208427486,4,0.9998992004886821,5.282E-05,"extracellular region"),
c("GO:0009505","plant-type cell wall",0.013780776369110043,1,0.9872638457957056,3.275E-05,"plant-type cell wall"),
c("GO:0005886","plasma membrane",14.184672303876209,6,0.9722381493674755,0.29250802,"plant-type cell wall"),
c("GO:0009506","plasmodesma",0.02618690220789438,1,0.9999359590919521,3.43E-05,"plasmodesma"),
c("GO:0009536","plastid",0.5104294469183077,7,0.6823513310670152,0,"plastid"),
c("GO:0005634","nucleus",12.339605699963672,13,0.578123451304178,0.44231542,"plastid"),
c("GO:0005737","cytoplasm",25.080076249313628,5,0.7543622687518001,0.22725924,"plastid"),
c("GO:0005739","mitochondrion",2.7614025547937224,3,0.6289338728049949,0.33771423,"plastid"),
c("GO:0005764","lysosome",0.2904358594155184,1,0.6782503799265944,0.22601794,"plastid"),
c("GO:0005773","vacuole",0.6250852016776067,1,0.6767722261747079,0.28151423,"plastid"),
c("GO:0005783","endoplasmic reticulum",2.123313387573007,4,0.5768419972738681,0.27526608,"plastid"),
c("GO:0005811","lipid droplet",0.07544965542347233,1,0.7019048550948056,0.4118361,"plastid"),
c("GO:0005829","cytosol",1.696303522837096,3,0.7330563829423168,0.22154616,"plastid"),
c("GO:0005840","ribosome",3.532158825827589,2,0.591537393537252,0.22161604,"plastid"),
c("GO:0009507","chloroplast",0.4142191413803319,5,0.6059529774655197,0.2334722,"plastid"),
c("GO:0010369","chromocenter",0.0029358879747399397,1,0.7574317221125776,0.31480507,"plastid"),
c("GO:0012511","monolayer-surrounded lipid storage body",0.00680851841612842,1,0.7451383986254116,0.33527212,"plastid"),
c("GO:0022625","cytosolic large ribosomal subunit",0.11811713241045162,2,0.6404092045880214,0.43014036,"plastid"),
c("GO:0030018","Z disc",0.03448430804052515,1,0.6884421737859988,0.3833401,"plastid"),
c("GO:0031966","mitochondrial membrane",1.347511654066338,1,0.55442047181033,0.30801698,"plastid"),
c("GO:0043661","peribacteroid membrane",3.807896205888378E-05,1,0.7718186759491326,0.12528944,"plastid"),
c("GO:0016020","membrane",61.57163680895251,19,0.9998372235275705,0.00011624,"membrane"),
c("GO:1990904","ribonucleoprotein complex",4.368586074828206,2,0.9163188984692716,-0,"ribonucleoprotein complex"),
c("GO:0005854","nascent polypeptide-associated complex",0.01439765555446396,1,0.7808724333896949,0.26990927,"ribonucleoprotein complex"),
c("GO:0009317","acetyl-CoA carboxylase complex",0.08877729214408166,1,0.7457109324693386,0.3836939,"ribonucleoprotein complex"),
c("GO:0019005","SCF ubiquitin ligase complex",0.043566140491568935,1,0.9321001607798933,0.29740373,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])
hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)



library(grid)
library(dplyr)

with(
  # t$tm is a data frame with one row per rectangle.  Filter to the group we
  # want to highlight.
  t$tm %>%
    filter(description %in% topGO_ReviGO_Intersect),
  {
    # Use grid.rect to add a rectangle on top of the treemap.
    grid.rect(x = x0 + (w / 2),
              y = y0 + (h / 2),
              width = w,
              height = h,
              gp = gpar(col = "yellow", fill = NA, lwd = 3),
              vp = "data")
    
  }
)
dev.off()

