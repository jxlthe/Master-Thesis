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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,2,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,4,0.9999045028833746,4.647E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,8,0.8807663144139185,0.08417172,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,2,0.8308974626973006,0.17160779,"cytoplasm"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,1,0.9653046351644721,3.018E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,8,0.9028962103285528,0.28538862,"plant-type cell wall"),
                     c("GO:0009536","plastid",0.7624897559369845,5,0.7919367026122143,0,"plastid"),
                     c("GO:0000326","protein storage vacuole",0.0009591191132014498,1,0.8744842036391205,0.10321831,"plastid"),
                     c("GO:0005634","nucleus",16.5161752456724,13,0.7099293946394818,0.35145602,"plastid"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.744532265502735,0.2503098,"plastid"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,1,0.6245014656964339,0.20731432,"plastid"),
                     c("GO:0005773","vacuole",1.35235298008496,1,0.7793048162302162,0.21057364,"plastid"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,4,0.6374782178787609,0.19658734,"plastid"),
                     c("GO:0009507","chloroplast",0.6990388085931706,2,0.7623733963062687,0.17236349,"plastid"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.9998957949175861,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,16,0.9998470848141064,9.537E-05,"membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.9999097497667101,5.465E-05,"cell projection"),
                     c("GO:0045202","synapse",1.0801172073373506,1,0.8571940349807317,4.855E-05,"synapse"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,1,0.9195873050105273,-0,"ribonucleoprotein complex"),
                     c("GO:0000015","phosphopyruvate hydratase complex",0.041644653723461905,1,0.8682150166326739,0.23819845,"ribonucleoprotein complex"),
                     c("GO:0030117","membrane coat",0.29445950681101296,1,0.7782435040005375,0.29013538,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_CC_Table.tsv", sep = "\t")
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



