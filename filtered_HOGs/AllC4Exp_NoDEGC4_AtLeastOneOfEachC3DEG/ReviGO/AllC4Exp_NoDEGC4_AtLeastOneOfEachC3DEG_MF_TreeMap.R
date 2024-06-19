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
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9277321826146957,0.04070747,"glutathione transferase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,3,0.9637918578215985,0.07993224,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,3,0.934805058457749,0,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9300274339252134,0.31015716,"zinc ion binding"),
                     c("GO:0005524","ATP binding",12.418006524131227,4,0.8963589951387093,0.26746377,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.9237444654140653,0.37942352,"zinc ion binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9504269344941726,0.12083997,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,1,0.9298237490314022,0.18168089,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,1,0.9726340124201915,0.05834827,"lipid binding"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.9292742341183012,-0,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9578984371331922,0.29459182,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.951676000791839,0.35458321,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9497163330768981,0.35290783,"ligand-gated monoatomic ion channel activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,9,0.9359613624514527,0.0830445,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,4,0.9417497750842001,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.936645291325742,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9493991577398884,0.05997119,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9516616283484537,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,1,0.9500276011425831,0.0589874,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9490108836223107,0.04472286,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9662945343000091,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,2,0.8210042296808621,-0,"ATP hydrolysis activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,1,0.8792213438440544,0.42844795,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.8730328155525962,0.35450436,"ATP hydrolysis activity"),
                     c("GO:0004298","threonine-type endopeptidase activity",0.042207766433033055,1,0.8368571495135634,0.46981862,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,2,0.788574523776987,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008081","phosphoric diester hydrolase activity",0.424918273177694,1,0.8094718727251017,0.26945252,"ATP hydrolysis activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,2,0.7802928649580687,0.3656957,"ATP hydrolysis activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8054879923335817,0.45595338,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.854669789283242,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,1,0.866099721971034,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0047782","coniferin beta-glucosidase activity",2.66099189448564E-05,1,0.8784272638229184,0.3973634,"ATP hydrolysis activity"),
                     c("GO:0070139","SUMO-specific endopeptidase activity",0.006703482080858407,1,0.8474368226275719,0.18227358,"ATP hydrolysis activity"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,1,0.8634077729955099,0.37283533,"ATP hydrolysis activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,2,0.983830955048966,-0,"signaling receptor activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9578211223943313,0.05189795,"identical protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9599078328645053,0.39898091,"identical protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9612831540639439,0.38355794,"identical protein binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.925570731080364,0.06961566,"sequence-specific DNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,2,0.9133923051843015,0.48153711,"sequence-specific DNA binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9377979962929239,0.33073906,"sequence-specific DNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9435955158068591,0.28841983,"sequence-specific DNA binding"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9220990254326026,0.03322883,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9127134136818348,0.30799911,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0047037","salutaridine reductase (NADPH) activity",1.77399459632376E-05,1,0.9517261712274524,0.38914098,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0047501","(+)-neomenthol dehydrogenase activity",2.66099189448564E-05,1,0.9509287428710922,0.3973742,"glutathione dehydrogenase (ascorbate) activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9292986410001567,0.05781294,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.947165212862748,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9746685872764549,0.01879944,"pterocarpan synthase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9635594141554411,0.23977693,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_MF_Table.tsv", sep = "\t")
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

