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
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9337110928954263,0.04668641,"translation initiation factor activity"),
                     c("GO:0003723","RNA binding",6.099813894661886,1,0.9074954800738403,0.34520254,"translation initiation factor activity"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9260545876803968,0.24940265,"translation initiation factor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0005515","protein binding",8.610051728351934,2,0.9562618187433756,0.06960742,"protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9728461942621843,0.04505158,"calmodulin binding"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9791953685095711,-0,"phosphatidylinositol transfer activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,6,0.9303677761411089,0.0830445,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9369434002734686,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,2,0.931141723353588,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9456232781319497,0.05997119,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9481756669577021,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9463330863851006,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9414959979979491,0.04472286,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9603176088388498,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,2,0.8295433660701551,0,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8716740809455901,0.35450436,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.7962202662229924,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9205756783951948,0.03322883,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.905059530511471,-0,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9136241846650163,0.31015716,"iron-sulfur cluster binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,4,0.8634289110966092,0.26746377,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,1,0.9129423235538945,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.8989871378015214,0.37942352,"iron-sulfur cluster binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9323615679242189,0.10462209,"iron-sulfur cluster binding"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,1,0.8447059239259296,0.03812803,"protein serine kinase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.866922162628014,0.42844795,"protein serine kinase activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8530913035056045,0.25738438,"protein serine kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9308067313434633,0.18603615,"protein serine kinase activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,1,0.7838849590751893,0.37283533,"protein serine kinase activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9157122115697741,0.05781294,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.935748566292739,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9712716402102628,0.01879944,"pterocarpan synthase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9593242515245328,0.23977693,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
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
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, " Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_Table.tsv", sep = "\t")
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

