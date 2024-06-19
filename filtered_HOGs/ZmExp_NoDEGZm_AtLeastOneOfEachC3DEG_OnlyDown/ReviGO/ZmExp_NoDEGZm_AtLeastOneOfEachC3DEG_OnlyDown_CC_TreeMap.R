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
                     c("GO:0005576","extracellular region",3.8687535442060237,9,0.9999164107002521,4.598E-05,"extracellular region"),
                     c("GO:0005737","cytoplasm",43.076660800468886,37,0.888698863151965,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,18,0.8576699668760611,0.17160779,"cytoplasm"),
                     c("GO:0009504","cell plate",0.0007255512462560189,1,0.9373247065460408,0.05708197,"cell plate"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,1,0.9709662236773764,2.997E-05,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,9,0.9501447151073252,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,22,0.7106448688607407,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,2,0.804497965675289,0.13280054,"chloroplast"),
                     c("GO:0001673","male germ cell nucleus",0.006684016617906476,1,0.8360702504197768,0.1163509,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,59,0.7074440060782706,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,3,0.7432982904889285,0.2503098,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,2,0.7873356189967988,0.17354881,"chloroplast"),
                     c("GO:0005788","endoplasmic reticulum lumen",0.11812570632538061,4,0.7056582812642246,0.47469365,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,17,0.649438687571789,0.19658734,"chloroplast"),
                     c("GO:0005797","Golgi medial cisterna",0.015012947362598856,2,0.738297847695811,0.39669825,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,5,0.700122434902463,0.15783714,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,3,0.6941870352481858,0.49095344,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,18,0.7869536915557038,0.17236349,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,3,0.7599599827038843,0.43399423,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.7759755773508121,0.22189912,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7519763330461235,0.30090248,"chloroplast"),
                     c("GO:0012505","endomembrane system",6.893425119210306,5,0.9999099064894741,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,116,0.9998758904921088,9.537E-05,"membrane"),
                     c("GO:0019005","SCF ubiquitin ligase complex",0.14169866753507532,3,0.9091563726356336,-0,"SCF ubiquitin ligase complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,2,0.8257573386651481,0.19052614,"SCF ubiquitin ligase complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7852224787846288,0.38761376,"SCF ubiquitin ligase complex"),
                     c("GO:0005886","plasma membrane",17.177321395000487,44,0.9129233654725646,0.34881718,"SCF ubiquitin ligase complex"),
                     c("GO:0009898","cytoplasmic side of plasma membrane",0.2657902935257323,1,0.9378338399477214,0.20588623,"SCF ubiquitin ligase complex"),
                     c("GO:0010330","cellulose synthase complex",0.0019629639881584074,1,0.9333903733494184,0.38954564,"SCF ubiquitin ligase complex"),
                     c("GO:0016602","CCAAT-binding factor complex",0.03082101937862897,1,0.7691215119887675,0.41728323,"SCF ubiquitin ligase complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.7508330451999053,0.19254961,"SCF ubiquitin ligase complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,2,0.8719318836064335,0.21998964,"SCF ubiquitin ligase complex"),
                     c("GO:0045273","respiratory chain complex II",0.05444864540824706,1,0.8811186657595379,0.46599171,"SCF ubiquitin ligase complex"),
                     c("GO:0046930","pore complex",0.06221601936645362,1,0.8913857435781342,0.47090253,"SCF ubiquitin ligase complex"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.99992041171445,5.465E-05,"cell projection"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9906499036338272,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9064298792259897,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9913127437403347,2.369E-05,"pollen tube"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_CC_Table.tsv", sep = "\t")
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

