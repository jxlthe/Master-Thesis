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
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,4,0.8558294592345888,0,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,2,0.7502799574733543,0.44154323,"monoatomic ion transport"),
                     c("GO:0035672","oligopeptide transmembrane transport",0.23299928216907387,1,0.8709360234047956,0.30173276,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.7380855241749573,-0,"defense response"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,1,0.8967939814685623,0.38212221,"defense response"),
                     c("GO:0009416","response to light stimulus",0.2589924319660237,1,0.8123843145513597,0.27535725,"defense response"),
                     c("GO:0010039","response to iron ion",0.05924742412469929,1,0.8065962760979498,0.37669922,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.7571506143619126,0.46097507,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.6737446503385388,0.37522238,"defense response"),
                     c("GO:0071230","cellular response to amino acid stimulus",0.039006977876590854,1,0.7761054891889361,0.232997,"defense response"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,1,0.9670960296176608,0.07105297,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,1,0.9626264511905962,0.12277613,"lipid catabolic process"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9913546418996747,0.00723774,"chromosome condensation"),
                     c("GO:0045739","positive regulation of DNA repair",0.0330817425160876,1,0.9351177414386274,-0,"positive regulation of DNA repair"),
                     c("GO:0048316","seed development",0.03280076213709535,1,0.968402770812103,-0,"seed development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.968402770812103,0.42998645,"seed development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9078452347744271,-0,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9431701279846126,0.15925058,"DNA demethylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)


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
# # 
# 
# # also find the reduced form of topGO sig Terms and mark them too
# not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_BP_Table.tsv", sep = "\t")
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

