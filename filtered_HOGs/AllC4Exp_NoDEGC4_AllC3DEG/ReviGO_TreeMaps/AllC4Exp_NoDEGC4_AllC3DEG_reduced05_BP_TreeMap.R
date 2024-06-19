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
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,2,0.859657561830051,-0,"monoatomic ion transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8452765335777909,0.4188028,"monoatomic ion transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.846647839368461,0.30464429,"monoatomic ion transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,1,0.840916109634737,0.37247286,"monoatomic ion transport"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,2,0.8772005106481118,-0,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.9292980500454672,0.12277613,"lipid catabolic process"),
                     c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.8925541854984664,0.46381088,"lipid catabolic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.9254223629521837,0.07741388,"phosphorylation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9474388241629303,0.11172681,"phosphorylation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,1,-0,"carbohydrate homeostasis"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,2,0.7760089427909355,0,"intracellular signal transduction"),
                     c("GO:0006952","defense response",1.1604144588756624,2,0.8131845409871619,0.37522238,"intracellular signal transduction"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,1,0.8258431389896448,0.47274124,"intracellular signal transduction"),
                     c("GO:0009733","response to auxin",0.1256056236300358,1,0.8526513918430121,0.49934907,"intracellular signal transduction"),
                     c("GO:0010468","regulation of gene expression",13.923871767653479,1,0.9428240885966082,0.39926133,"intracellular signal transduction"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8276119997819956,0.46171648,"intracellular signal transduction"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9605487851922724,0.21912257,"intracellular signal transduction"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,1,0.8305467899052755,0.42391909,"intracellular signal transduction"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9803322186234948,0.11943324,"intracellular signal transduction"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8115566104700873,0.42151038,"intracellular signal transduction"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8876779238953187,0.07105297,"DNA demethylation"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8589276766539276,0.36848795,"DNA demethylation"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8610960350815848,0.27735821,"DNA demethylation"),
                     c("GO:0006508","proteolysis",5.2622572267907,1,0.8629236559306371,0.44475214,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9180290003370447,0.15925058,"DNA demethylation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9344204935933396,0.23039539,"DNA demethylation"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.9016549919916592,0.15970956,"DNA demethylation"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9067097943832549,0.1131676,"DNA demethylation"),
                     c("GO:0140547","acquisition of seed longevity",0.00018978499282809906,1,0.9622253596069233,-0,"acquisition of seed longevity"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.8988350893712908,0.23389832,"acquisition of seed longevity"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9569296911885713,0.43353892,"acquisition of seed longevity"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9555047392627222,0.35287398,"acquisition of seed longevity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("All C4 Expanded, No DEG in C4, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_Table.tsv", sep = "\t")
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

