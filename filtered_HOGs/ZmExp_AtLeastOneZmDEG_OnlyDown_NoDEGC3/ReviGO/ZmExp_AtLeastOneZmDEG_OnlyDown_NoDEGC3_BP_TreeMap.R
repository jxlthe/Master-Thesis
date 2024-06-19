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
revigo.data <- rbind(c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,1,0.9802923737144542,0.05630604,"carbohydrate metabolic process"),
                     c("GO:0008152","metabolic process",57.597931274565454,1,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9924509524870899,0.00908442,"cell population proliferation"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,2,0.7343922169732996,0,"auxin-activated signaling pathway"),
                     c("GO:0001558","regulation of cell growth",0.15345225803226778,1,0.9035121795101533,0.16011872,"auxin-activated signaling pathway"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.1438718130046987,1,0.890996607463209,0.21941245,"auxin-activated signaling pathway"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.8192798869174048,0.24614308,"auxin-activated signaling pathway"),
                     c("GO:0007165","signal transduction",8.756627809544993,2,0.7069829307373783,0.47698259,"auxin-activated signaling pathway"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.7910388589017517,0.31987157,"auxin-activated signaling pathway"),
                     c("GO:0009643","photosynthetic acclimation",0.000325345701991027,1,0.8819377050944387,0.45106647,"auxin-activated signaling pathway"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8208711584340853,0.4259419,"auxin-activated signaling pathway"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9175885584765684,0.20908612,"auxin-activated signaling pathway"),
                     c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,1,0.8994986328820143,0.35474033,"auxin-activated signaling pathway"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,1,0.9124823657207489,0.28549545,"auxin-activated signaling pathway"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8411220649809261,0.47074713,"auxin-activated signaling pathway"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,1,0.8432302370486209,0.19569204,"auxin-activated signaling pathway"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,1,0.8757135091819515,0.42863919,"auxin-activated signaling pathway"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.9127106875120198,0.18626259,"auxin-activated signaling pathway"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9027616437444468,0.1517125,"auxin-activated signaling pathway"),
                     c("GO:0060918","auxin transport",0.015500750907739157,1,0.8740624142281647,0.11396033,"auxin-activated signaling pathway"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.8723847562786669,0.32241615,"auxin-activated signaling pathway"),
                     c("GO:1902882","regulation of response to oxidative stress",0.01609475328788944,1,0.884966155902162,0.12060378,"auxin-activated signaling pathway"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9723606319039315,0.03847046,"coumarin biosynthetic process"),
                     c("GO:0006564","L-serine biosynthetic process",0.06954757328091521,1,0.9625908650794829,0.11280602,"coumarin biosynthetic process"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",0.0004411884898211653,1,0.960418450021046,0.00539809,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,1,0.9416590292713564,0.36594593,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,3,0.9638208410633946,0.25284502,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0030154","cell differentiation",2.279586420543629,1,0.9743129543280651,0.01240155,"cell differentiation"),
                     c("GO:0010162","seed dormancy process",0.001212652161966555,1,0.9636117781828409,0.40373092,"cell differentiation"),
                     c("GO:0044843","cell cycle G1/S phase transition",0.05461124787132716,1,0.9928079783171261,0.00716697,"cell cycle G1/S phase transition"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9378672686661721,0.00884962,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,1,0.9510317837341355,0.4063955,"cell wall organization"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1,0.9103438370458188,0.37183399,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_BP_Table.tsv", sep = "\t")
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

