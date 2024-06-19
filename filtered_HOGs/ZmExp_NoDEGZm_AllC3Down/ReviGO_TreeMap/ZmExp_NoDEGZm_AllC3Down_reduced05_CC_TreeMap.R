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
                     c("GO:0005576","extracellular region",3.8687535442060237,6,0.9999072427077728,4.598E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,16,0.874050616778705,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,9,0.8327770134344908,0.17160779,"cytoplasm"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,6,0.9999378240664689,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,14,0.6925545081595136,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,2,0.7785261420993793,0.13280054,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,22,0.6830629422535041,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,2,0.7226568607680995,0.26090848,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,2,0.7726781380112001,0.17354881,"chloroplast"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,9,0.5870049357272622,0.2030659,"chloroplast"),
                     c("GO:0005788","endoplasmic reticulum lumen",0.11812570632538061,4,0.6719479964932582,0.48748182,"chloroplast"),
                     c("GO:0005797","Golgi medial cisterna",0.015012947362598856,1,0.6885502534137447,0.40558994,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,1,0.7649686227732585,0.43399423,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,13,0.7722484751190409,0.17236349,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,2,0.8129298426243511,0.12200202,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7273357215929946,0.30599076,"chloroplast"),
                     c("GO:0048046","apoplast",0.0960361495472436,1,0.9785429697301748,0.45176992,"chloroplast"),
                     c("GO:0070062","extracellular exosome",0.10153244871408715,1,0.7636087017257679,0.1380203,"chloroplast"),
                     c("GO:0012505","endomembrane system",6.893425119210306,2,0.9998990366728665,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,47,0.9998535964571377,9.537E-05,"membrane"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,1,0.8793774991137475,3.69E-05,"monoatomic ion channel complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,1,0.8390650576120395,0.21034994,"monoatomic ion channel complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7838973991835493,0.38761376,"monoatomic ion channel complex"),
                     c("GO:0005886","plasma membrane",17.177321395000487,18,0.9482739130660679,0.30272327,"monoatomic ion channel complex"),
                     c("GO:0019005","SCF ubiquitin ligase complex",0.14169866753507532,1,0.9369826166490048,0.21998964,"monoatomic ion channel complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.748904942431479,0.21281912,"monoatomic ion channel complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_CC_Table.tsv", sep = "\t")
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

