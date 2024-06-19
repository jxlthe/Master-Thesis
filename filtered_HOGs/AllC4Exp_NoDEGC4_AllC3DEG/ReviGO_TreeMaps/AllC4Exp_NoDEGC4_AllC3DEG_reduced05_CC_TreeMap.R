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
revigo.data <- rbind(c("GO:0000502","proteasome complex",0.37652382533874423,1,0.8271175970868048,-0,"proteasome complex"),
                     c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,1,0.9999075555427188,4.999E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9878249124521837,2.721E-05,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,8,0.876325831792117,0.09354521,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,6,0.8517450403470601,0.17160779,"cytoplasm"),
                     c("GO:0005773","vacuole",1.35235298008496,2,0.7954947398274242,0,"vacuole"),
                     c("GO:0005634","nucleus",16.5161752456724,5,0.722821401188212,0.35145602,"vacuole"),
                     c("GO:0005654","nucleoplasm",1.830446525605921,1,0.7846459032454033,0.21968674,"vacuole"),
                     c("GO:0005739","mitochondrion",4.856981674016684,1,0.7604285311621104,0.25147478,"vacuole"),
                     c("GO:0005774","vacuolar membrane",0.7015583598387308,1,0.7266424122076163,0.18309297,"vacuole"),
                     c("GO:0005777","peroxisome",0.7476308639759881,1,0.8084648576917414,0.18435889,"vacuole"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,1,0.7583222273090043,0.22703414,"vacuole"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.8098307472795678,0.18302188,"vacuole"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.7808471979237221,0.21139748,"vacuole"),
                     c("GO:0005886","plasma membrane",17.177321395000487,5,0.9808601938612675,0.06338235,"plasma membrane"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999273299581207,4.05E-05,"cell surface"),
                     c("GO:0016020","membrane",49.2542153160787,7,0.9998444774268712,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,1,1,-0,"protein-containing complex"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9999407554378136,3.357E-05,"apoplast"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, " Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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





# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_Table.tsv", sep = "\t")
representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
sig_reduced <- c()
for(sig in representatives)
{
  sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
}

sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])

if(length(sig_reduced_ReviGO_intersect) > 0)
{
  with(
    # t$tm is a data frame with one row per rectangle.  Filter to the group we
    # want to highlight.
    t$tm %>%
      filter(description %in% sig_reduced_ReviGO_intersect),
    {
      # Use grid.rect to add a rectangle on top of the treemap.
      grid.rect(x = x0 + (w / 2),
                y = y0 + (h / 2),
                width = w,
                height = h,
                gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
                vp = "data")
      
    }
  )
}




dev.off()

