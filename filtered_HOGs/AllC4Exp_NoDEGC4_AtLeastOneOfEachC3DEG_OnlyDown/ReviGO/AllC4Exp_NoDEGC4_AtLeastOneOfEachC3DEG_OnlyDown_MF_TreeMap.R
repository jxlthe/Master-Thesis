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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,1,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9717640153210997,0.04505158,"calmodulin binding"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.871105943940637,-0,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.9163165329214024,0.35458321,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9059564859780858,0.35290783,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,2,0.9583291816290221,0.0830445,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9660993454886375,0.05997119,"lyase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9665127633927455,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9443095294759732,0.04301538,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9559580965807031,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8937513058969788,0.05781294,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.9109092169090824,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,2,0.9568103026547466,0,"signaling receptor activity"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.9143895476121282,0.06306905,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,2,0.9025609564036702,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.8696048369605115,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,2,0.8401965582431901,0.26407481,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.890337407140178,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9092895267009665,0.01849091,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9228090868021512,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9870662401598954,-0,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

# 
# 
# library(grid)
# library(dplyr)
# 
# if(length(topGO_ReviGO_Intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% topGO_ReviGO_Intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3),
#                 vp = "data")
#       
#     }
#   )
# }
# 
# 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_MF_Table.tsv", sep = "\t")
# representatives <- unique(ReviGO_table[ReviGO_table$TermID %in% not_plotted,]$Representative)
# sig_reduced <- c()
# for(sig in representatives)
# {
#   sig_reduced <- append(sig_reduced, ReviGO_table[grepl(sig, ReviGO_table$TermID),]$Name)
# }
# 
# sig_reduced_ReviGO_intersect <- intersect(sig_reduced, revigo.data[,2])
# 
# if(length(sig_reduced_ReviGO_intersect) > 0)
# {
#   with(
#     # t$tm is a data frame with one row per rectangle.  Filter to the group we
#     # want to highlight.
#     t$tm %>%
#       filter(description %in% sig_reduced_ReviGO_intersect),
#     {
#       # Use grid.rect to add a rectangle on top of the treemap.
#       grid.rect(x = x0 + (w / 2),
#                 y = y0 + (h / 2),
#                 width = w,
#                 height = h,
#                 gp = gpar(col = "yellow", fill = NA, lwd = 3, lty = 2),
#                 vp = "data")
#       
#     }
#   )
# }


dev.off()

