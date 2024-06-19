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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,5,1,-0,"cellular_component"),
c("GO:0005576","extracellular region",3.8687535442060237,6,0.9999257567330232,6.017E-05,"extracellular region"),
c("GO:0005739","mitochondrion",4.856981674016684,30,0.7570150031465227,0,"mitochondrion"),
c("GO:0000138","Golgi trans cisterna",0.007541260384887046,1,0.8004846910996123,0.37607379,"mitochondrion"),
c("GO:0000325","plant-type vacuole",0.04066068696484073,5,0.800433056704888,0.15656526,"mitochondrion"),
c("GO:0000785","chromatin",1.087181392930179,3,0.6297811278249654,0.45756389,"mitochondrion"),
c("GO:0005634","nucleus",16.5161752456724,122,0.7288138569252248,0.35145602,"mitochondrion"),
c("GO:0005730","nucleolus",1.2114693153196516,16,0.6654899531811092,0.22801288,"mitochondrion"),
c("GO:0005737","cytoplasm",43.076660800468886,92,0.8926316389378861,0.12447513,"mitochondrion"),
c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,9,0.6351655235600465,0.22750481,"mitochondrion"),
c("GO:0005773","vacuole",1.35235298008496,6,0.7859504209325622,0.23143591,"mitochondrion"),
c("GO:0005777","peroxisome",0.7476308639759881,4,0.7897848457014061,0.21411812,"mitochondrion"),
c("GO:0005794","Golgi apparatus",2.349961096041566,21,0.7093335411531074,0.2503098,"mitochondrion"),
c("GO:0005811","lipid droplet",0.20804685033482948,3,0.7886697141064198,0.38610723,"mitochondrion"),
c("GO:0005829","cytosol",14.048080989011824,39,0.8614320568974416,0.19224382,"mitochondrion"),
c("GO:0005905","clathrin-coated pit",0.09808111076528711,2,0.8969327790575943,0.46642691,"mitochondrion"),
c("GO:0009507","chloroplast",0.6990388085931706,29,0.7199715797459495,0.21231674,"mitochondrion"),
c("GO:0009536","plastid",0.7624897559369845,25,0.7967142509721892,0.21465144,"mitochondrion"),
c("GO:0009579","thylakoid",0.26089530737804617,4,0.7848965332746932,0.39454155,"mitochondrion"),
c("GO:0017177","glucosidase II complex",0.015944734065838607,2,0.6873665946024858,0.39861008,"mitochondrion"),
c("GO:0031410","cytoplasmic vesicle",2.424891655569294,6,0.6983167248005325,0.25147478,"mitochondrion"),
c("GO:0031982","vesicle",2.6690048632308567,3,0.7950210667277179,0.2461316,"mitochondrion"),
c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,1,0.7857697422896927,0.4052793,"mitochondrion"),
c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7422829307765494,0.486126,"mitochondrion"),
c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9068606864093187,0.10673752,"mitochondrion"),
c("GO:0070971","endoplasmic reticulum exit site",0.06951625759076932,1,0.7773003449898026,0.45185985,"mitochondrion"),
c("GO:0009506","plasmodesma",0.10230521048664064,8,0.9772103571276574,3.812E-05,"plasmodesma"),
c("GO:0009524","phragmoplast",0.002879842104146322,1,0.9315447531842158,0.0711282,"phragmoplast"),
c("GO:0009986","cell surface",0.6578761991908514,1,0.9999382845144928,4.693E-05,"cell surface"),
c("GO:0012505","endomembrane system",6.893425119210306,1,0.9999203715136444,0.00011166,"endomembrane system"),
c("GO:0016020","membrane",49.2542153160787,115,0.9998928023305879,0.00010118,"membrane"),
c("GO:0032991","protein-containing complex",19.75554538569037,3,1,-0,"protein-containing complex"),
c("GO:0044297","cell body",0.22864057885869896,1,0.9999438437429754,4.148E-05,"cell body"),
c("GO:0048046","apoplast",0.0960361495472436,2,0.9999476834174484,3.787E-05,"apoplast"),
c("GO:0070469","respirasome",0.5984977809313304,2,0.9690463907323354,4.638E-05,"respirasome"),
c("GO:0090406","pollen tube",0.001702063711251277,3,0.9968157547408188,2.697E-05,"pollen tube"),
c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.9968236570771475,0.27222822,"pollen tube"),
c("GO:0098552","side of membrane",0.7241349304670945,1,0.9685386317153775,0.0715415,"side of membrane"),
c("GO:0005886","plasma membrane",17.177321395000487,32,0.9489593005345854,0.10744228,"side of membrane"),
c("GO:1990904","ribonucleoprotein complex",4.315842197772249,15,0.785047308750822,-0,"ribonucleoprotein complex"),
c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7454113519299168,0.47495898,"ribonucleoprotein complex"),
c("GO:0000178","exosome (RNase complex)",0.10685978389207654,1,0.756601954765338,0.34132399,"ribonucleoprotein complex"),
c("GO:0000502","proteasome complex",0.37652382533874423,4,0.7489665034469939,0.49612005,"ribonucleoprotein complex"),
c("GO:0000974","Prp19 complex",0.0635900941581645,1,0.8446460881447212,0.24779784,"ribonucleoprotein complex"),
c("GO:0005681","spliceosomal complex",0.7626239332222511,10,0.5997265732831089,0.3245659,"ribonucleoprotein complex"),
c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8070784917476633,0.22864326,"ribonucleoprotein complex"),
c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7915274515176007,0.255018,"ribonucleoprotein complex"),
c("GO:0005942","phosphatidylinositol 3-kinase complex",0.07251039886384637,1,0.7405014019971716,0.46964777,"ribonucleoprotein complex"),
c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.7954019921746794,0.43971777,"ribonucleoprotein complex"),
c("GO:0016272","prefoldin complex",0.03795726314317447,1,0.849760484954309,0.23619435,"ribonucleoprotein complex"),
c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,2,0.7757571761165165,0.42147812,"ribonucleoprotein complex"),
c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7368856835896797,0.49563216,"ribonucleoprotein complex"),
c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7937655158975829,0.41997104,"ribonucleoprotein complex"),
c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7736831980262837,0.41725271,"ribonucleoprotein complex"),
c("GO:0045252","oxoglutarate dehydrogenase complex",0.07266196950090673,1,0.7234534952614657,0.46971574,"ribonucleoprotein complex"),
c("GO:0046930","pore complex",0.06221601936645362,2,0.7953428861431854,0.24728353,"ribonucleoprotein complex"),
c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7017571486001692,0.46572404,"ribonucleoprotein complex"),
c("GO:0071439","clathrin complex",0.011926869801468804,1,0.7618966594313018,0.36298129,"ribonucleoprotein complex"),
c("GO:0071819","DUBm complex",0.008179844872174023,1,0.71553362488728,0.43933875,"ribonucleoprotein complex"),
c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,6,0.7708417521135674,0.25161359,"ribonucleoprotein complex"),
c("GO:0089701","U2AF complex",0.020742814396387827,1,0.6995805439409414,0.46994592,"ribonucleoprotein complex"),
c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.78791741550791,0.40438888,"ribonucleoprotein complex"),
c("GO:0160064","multi-pass translocon complex",0.0006336149582030303,1,0.8428987684788185,0.17222646,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

