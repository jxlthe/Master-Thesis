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
revigo.data <- rbind(c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9295405986985653,0.04087,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.906244091234332,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9221916821559171,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,1,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9277573645136502,0.04086192,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.952334997282249,0.06542584,"protein binding"),
                     c("GO:0008233","peptidase activity",4.167383840940452,2,0.8083806670198322,0,"peptidase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.8782833034831994,0.42844795,"peptidase activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8929127249359601,0.30149286,"peptidase activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8343939322621955,0.46981862,"peptidase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8094514507876345,0.45595338,"peptidase activity"),
                     c("GO:0070139","SUMO-specific endopeptidase activity",0.006703482080858407,1,0.8342692672402504,0.29849456,"peptidase activity"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,1,0.8633346840746455,0.20765779,"peptidase activity"),
                     c("GO:0008289","lipid binding",1.291468066123697,1,0.9633534171978166,0.04613007,"lipid binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,4,0.935934694463014,0.08368982,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9415366565207313,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.9365708020414661,0.12712323,"hydrolase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9694178287693785,0.03203827,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9515729212267898,0.05675833,"isomerase activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9671609285137994,0.03203534,"identical protein binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9200546853179998,-0,"glutathione binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9125915317292687,0.31015716,"glutathione binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,1,0.8807569969649474,0.26746377,"glutathione binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.9148165787056828,0.12083997,"glutathione binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.8957533932258503,0.03333167,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0046961","proton-transporting ATPase activity, rotational mechanism",0.20284297713014954,1,0.9665876989887132,-0,"proton-transporting ATPase activity, rotational mechanism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
# 
# topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
# topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, Only Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
# 
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
# ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_MF_Table.tsv", sep = "\t")
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

