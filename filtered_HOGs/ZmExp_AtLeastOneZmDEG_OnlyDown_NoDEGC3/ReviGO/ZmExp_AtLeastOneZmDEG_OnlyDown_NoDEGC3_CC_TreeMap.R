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
revigo.data <- rbind(c("GO:0005576","extracellular region",3.8687535442060237,3,0.9998978355557501,0,"extracellular region"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,2,0.6437253130912475,5.394E-05,"Golgi apparatus"),
                     c("GO:0005634","nucleus",16.5161752456724,9,0.6698834627919237,0.35145602,"Golgi apparatus"),
                     c("GO:0005730","nucleolus",1.2114693153196516,1,0.7546055924229319,0.20773614,"Golgi apparatus"),
                     c("GO:0005737","cytoplasm",43.076660800468886,2,0.8668703845762429,0.1047991,"Golgi apparatus"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.7104688145276024,0.2503098,"Golgi apparatus"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,1,0.5665642583874467,0.20731432,"Golgi apparatus"),
                     c("GO:0009507","chloroplast",0.6990388085931706,1,0.7086931775584476,0.19462727,"Golgi apparatus"),
                     c("GO:0009536","plastid",0.7624897559369845,1,0.7618898553001906,0.19658734,"Golgi apparatus"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,1,0.652862337588118,0.25147478,"Golgi apparatus"),
                     c("GO:0016020","membrane",49.2542153160787,8,0.9998355907512326,9.537E-05,"membrane"),
                     c("GO:0031012","extracellular matrix",0.6517835565339344,1,0.9626729933883422,4.559E-05,"extracellular matrix"),
                     c("GO:0005886","plasma membrane",17.177321395000487,5,0.90272072198916,0.39486578,"extracellular matrix"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,1,0.900045229888717,-0,"ribonucleoprotein complex"),
                     c("GO:0030117","membrane coat",0.29445950681101296,1,0.7080464488256948,0.29013538,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

library(grid)
library(dplyr)

if(length(topGO_ReviGO_Intersect) > 0)
{
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
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_CC_Table.tsv", sep = "\t")
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

