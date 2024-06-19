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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,7,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,10,0.9999257530914285,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9777373502571378,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,93,0.8953425503676568,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,41,0.8659805935478573,0.17160779,"cytoplasm"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,11,0.9749079730975561,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,46,0.7318287361495825,0,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,3,0.8115344680290907,0.13280054,"chloroplast"),
                     c("GO:0000786","nucleosome",0.17397327417075828,3,0.6654943895171002,0.41609548,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,126,0.7380939282232801,0.35145602,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,17,0.7647495322700621,0.26090848,"chloroplast"),
                     c("GO:0005742","mitochondrial outer membrane translocase complex",0.03996992242217233,1,0.6251152373967991,0.45499251,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,6,0.792921125079359,0.21802429,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,9,0.7956693369993979,0.17354881,"chloroplast"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,19,0.7078727052934302,0.2030659,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,1,0.7855235907068862,0.42381054,"chloroplast"),
                     c("GO:0005856","cytoskeleton",3.1074887771882316,7,0.7293296072016058,0.15783714,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,3,0.7418478568689448,0.49095344,"chloroplast"),
                     c("GO:0008278","cohesin complex",0.03346232408674592,1,0.7116417682771674,0.4995934,"chloroplast"),
                     c("GO:0009522","photosystem I",0.033362933505067006,2,0.6778045196054654,0.35622153,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,47,0.8033821939658199,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,1,0.8308007132050742,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,5,0.78174159519074,0.43399423,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6933246414373261,0.40758865,"chloroplast"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,5,0.7125904709330525,0.23571918,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.8024356179131102,0.23050241,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.7995452269663076,0.4052793,"chloroplast"),
                     c("GO:0032586","protein storage vacuole membrane",0.00013169252072455135,1,0.8005647844849978,0.39770966,"chloroplast"),
                     c("GO:0033097","amyloplast membrane",7.4542936259180024E-06,1,0.7956395472568772,0.45278409,"chloroplast"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,1,0.7625147011790306,0.44880082,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999383493815162,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8871293438012026,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,4,0.9999203225276694,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,135,0.9998923458366435,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999439252066005,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9999477722348837,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,3,0.9101336775332723,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0062023","collagen-containing extracellular matrix",0.3197519250837527,1,0.9736314571567926,3.527E-05,"collagen-containing extracellular matrix"),
                     c("GO:0005886","plasma membrane",17.177321395000487,50,0.9340226189443323,0.35740582,"collagen-containing extracellular matrix"),
                     c("GO:0019897","extrinsic component of plasma membrane",0.23434808301161014,1,0.9444410808444494,0.22734667,"collagen-containing extracellular matrix"),
                     c("GO:0090404","pollen tube tip",0.00040998614942549014,1,0.9999633644803905,2.175E-05,"pollen tube tip"),
                     c("GO:0098552","side of membrane",0.7241349304670945,1,0.9689450210341369,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,11,0.785911280813043,-0,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7508173070633546,0.42313164,"ribonucleoprotein complex"),
                     c("GO:0005667","transcription regulator complex",0.7880231963702957,1,0.8151175856854667,0.32589738,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7250098958054277,0.38761376,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8106290470898265,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7951516237973681,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005853","eukaryotic translation elongation factor 1 complex",0.01002354016231774,1,0.8216842407429149,0.21073066,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.7697752543865839,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.7715584015735832,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0031519","PcG protein complex",0.09511430190217174,1,0.6829025228810525,0.25767523,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7937233990166496,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7507656679905836,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7743188138755523,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,5,0.7722807812800847,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7888424149240479,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
