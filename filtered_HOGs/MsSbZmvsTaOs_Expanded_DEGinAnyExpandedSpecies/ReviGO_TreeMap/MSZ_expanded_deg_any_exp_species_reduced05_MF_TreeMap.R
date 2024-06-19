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
revigo.data <- rbind(c("GO:0001671","ATPase activator activity",0.0341266115541096,1,1,-0,"ATPase activator activity"),
c("GO:0003700","DNA-binding transcription factor activity",4.532603540907692,4,0.974517641048897,-0,"DNA-binding transcription factor activity"),
c("GO:0003735","structural constituent of ribosome",2.2675249562595297,2,0.9876859144497072,-0,"structural constituent of ribosome"),
c("GO:0003989","acetyl-CoA carboxylase activity",0.09778675834076506,1,0.9680756504817798,0.03241463,"acetyl-CoA carboxylase activity"),
c("GO:0005515","protein binding",5.073796515930947,5,0.9671970182257469,0.07857259,"protein binding"),
c("GO:0008270","zinc ion binding",3.776138036817838,2,0.9118799298452052,0.06655862,"zinc ion binding"),
c("GO:0000287","magnesium ion binding",1.7504660637586555,1,0.9330218898914651,0.37813207,"zinc ion binding"),
c("GO:0005509","calcium ion binding",1.4167951420399318,1,0.9346426100333597,0.36751657,"zinc ion binding"),
c("GO:0008198","ferrous iron binding",0.04409462238483419,1,0.9422771150598274,0.47821211,"zinc ion binding"),
c("GO:0043169","cation binding",18.51293079024014,1,0.9318775708616098,0.27606386,"zinc ion binding"),
c("GO:0016491","oxidoreductase activity",11.721830439531795,7,0.948754590514374,0.09355129,"oxidoreductase activity"),
c("GO:0016740","transferase activity",21.216634542374862,13,0.9444826589956185,0.07040944,"transferase activity"),
c("GO:0016787","hydrolase activity",21.340557702986732,6,0.9444375928461459,0.11166207,"transferase activity"),
c("GO:0016757","glycosyltransferase activity",2.3960575386602994,5,0.8658684166356257,0,"glycosyltransferase activity"),
c("GO:0004185","serine-type carboxypeptidase activity",0.2130250273078274,1,0.8738354742290061,0.37435024,"glycosyltransferase activity"),
c("GO:0004478","methionine adenosyltransferase activity",0.0400000828170345,1,0.9046567569670109,0.21140702,"glycosyltransferase activity"),
c("GO:0004674","protein serine/threonine kinase activity",1.268552650275868,3,0.8008849369801089,0.30164209,"glycosyltransferase activity"),
c("GO:0004721","phosphoprotein phosphatase activity",0.5987290200478247,1,0.8497554998303689,0.41513428,"glycosyltransferase activity"),
c("GO:0008168","methyltransferase activity",2.684244473270584,1,0.8475540584513593,0.34937253,"glycosyltransferase activity"),
c("GO:0016746","acyltransferase activity",3.4236888973740025,1,0.8609030877987566,0.34378837,"glycosyltransferase activity"),
c("GO:0016763","pentosyltransferase activity",0.5429555612510428,1,0.8178482066765651,0.4643706,"glycosyltransferase activity"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.40741643979158,1,0.8467183847227898,0.41755188,"glycosyltransferase activity"),
c("GO:0016776","phosphotransferase activity, phosphate group as acceptor",0.3442981997320215,1,0.8730271763372275,0.45182766,"glycosyltransferase activity"),
c("GO:0050505","hydroquinone glucosyltransferase activity",0.00010352129313132951,1,0.8705521572792748,0.38255741,"glycosyltransferase activity"),
c("GO:0050734","hydroxycinnamoyltransferase activity",0.004307030643174526,1,0.9138157539858599,0.3128147,"glycosyltransferase activity"),
c("GO:0051753","mannan synthase activity",0.002838662827443299,1,0.859925513385278,0.34619078,"glycosyltransferase activity"),
c("GO:0052636","arabinosyltransferase activity",0.0025662383718345365,1,0.8681219548068662,0.17080964,"glycosyltransferase activity"),
c("GO:0102406","omega-hydroxypalmitate O-sinapoyl transferase activity",9.534855946306665E-05,1,0.9291050517417418,0.13884408,"glycosyltransferase activity"),
c("GO:0016829","lyase activity",3.4809062057855114,1,0.9557116593592981,0.05103459,"lyase activity"),
c("GO:0016838","carbon-oxygen lyase activity, acting on phosphates",0.12872327951969606,1,0.9306632613945891,0.03327247,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0008948","oxaloacetate decarboxylase activity",0.024681655678153825,1,0.9436743477893873,0.43322627,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0050552","(4S)-limonene synthase activity",1.3621222780438093E-05,1,0.9436599512926328,0.40007135,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0102903","gamma-terpinene synthase activity",2.7242445560876187E-06,1,0.9484835211799921,0.35927317,"carbon-oxygen lyase activity, acting on phosphates"),
c("GO:0016887","ATP hydrolysis activity",3.291701972876113,3,0.8903472934655432,0.04836706,"ATP hydrolysis activity"),
c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.7350856328136782,1,0.9135470321751943,0.29181845,"ATP hydrolysis activity"),
c("GO:0052793","pectin acetylesterase activity",0.0008881037252845636,1,0.9327902018153947,0.16150403,"ATP hydrolysis activity"),
c("GO:0032453","histone H3K4 demethylase activity",0.010000701765397649,1,0.9133527628479083,0.02670362,"histone H3K4 demethylase activity"),
c("GO:0042910","xenobiotic transmembrane transporter activity",0.1732101931206069,1,0.967766881976088,0.01034163,"xenobiotic transmembrane transporter activity"),
c("GO:0015171","amino acid transmembrane transporter activity",0.2472497116659562,1,0.9670156354346987,0.35414066,"xenobiotic transmembrane transporter activity"),
c("GO:0043565","sequence-specific DNA binding",1.6398235953577127,3,0.9421291030608252,-0,"sequence-specific DNA binding"),
c("GO:0000166","nucleotide binding",19.572068037245,9,0.9207263260943546,0.47112053,"sequence-specific DNA binding"),
c("GO:0000976","transcription cis-regulatory region binding",0.37429485653910227,2,0.9502676073579898,0.38815494,"sequence-specific DNA binding"),
c("GO:0003677","DNA binding",11.827531128307996,6,0.9284325863898248,0.41192319,"sequence-specific DNA binding"),
c("GO:0003729","mRNA binding",0.29463794571910035,2,0.9535506502136918,0.25887051,"sequence-specific DNA binding"),
c("GO:0005525","GTP binding",1.9488618978002925,1,0.923165453996413,0.34132961,"sequence-specific DNA binding"),
c("GO:0020037","heme binding",1.466676059861896,3,0.9579581837897007,0.13335963,"sequence-specific DNA binding"),
c("GO:0031418","L-ascorbic acid binding",0.07127168607636428,1,0.947832325379137,0.20191161,"sequence-specific DNA binding"),
c("GO:0051287","NAD binding",0.9436456232940781,1,0.9445296330886084,0.12665645,"sequence-specific DNA binding"),
c("GO:0051087","protein-folding chaperone binding",0.14845770708399478,2,0.9298292504190665,0.04628281,"protein-folding chaperone binding"),
c("GO:0030544","Hsp70 protein binding",0.04413276180861943,1,0.9354711628598626,0.39749369,"protein-folding chaperone binding"),
c("GO:0031072","heat shock protein binding",0.10834320599560461,1,0.9313805060735177,0.42425405,"protein-folding chaperone binding"),
c("GO:0043621","protein self-association",0.01613842475026305,1,0.9395344681131841,0.37126244,"protein-folding chaperone binding"),
c("GO:0051082","unfolded protein binding",0.4145292243879603,2,0.9243039485814422,0.47170041,"protein-folding chaperone binding"),
c("GO:0140297","DNA-binding transcription factor binding",0.06685568565094625,1,0.9336380105890657,0.40943701,"protein-folding chaperone binding"),
c("GO:0051213","dioxygenase activity",0.6317060003992653,2,0.881211831865022,0.03928991,"dioxygenase activity"),
c("GO:0004497","monooxygenase activity",1.2382263598775005,1,0.8749229607046218,0.43585558,"dioxygenase activity"),
c("GO:0004601","peroxidase activity",0.4509741680593005,2,0.8613915519395966,0.36163411,"dioxygenase activity"),
c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.3673746215479463,2,0.8542871256794812,0.40450551,"dioxygenase activity"),
c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.4015829604387384,1,0.8736960818213074,0.44216553,"dioxygenase activity"),
c("GO:0016712","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",0.05945663743661228,1,0.8691652305490726,0.37156298,"dioxygenase activity"),
c("GO:0033771","licodione synthase activity",1.3621222780438093E-05,1,0.9075974940010506,0.38603396,"dioxygenase activity"),
c("GO:0047890","flavanone 4-reductase activity",1.906971189261333E-05,1,0.9195822764721437,0.41114057,"dioxygenase activity"),
c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.451820100478857E-05,1,0.9055300158092532,0.25783003,"dioxygenase activity"),
c("GO:0102469","naringenin 2-hydroxylase activity",2.7242445560876187E-06,1,0.912862992625522,0.36819953,"dioxygenase activity"),
c("GO:0102472","eriodictyol 2-hydroxylase activity",2.7242445560876187E-06,1,0.912862992625522,0.16834365,"dioxygenase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );



# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])
hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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

dev.off()

