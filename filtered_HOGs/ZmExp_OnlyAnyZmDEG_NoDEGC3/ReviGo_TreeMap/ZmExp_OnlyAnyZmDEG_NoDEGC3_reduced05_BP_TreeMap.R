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
revigo.data <- rbind(c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,1,1,-0,"intracellular iron ion homeostasis"),
                     c("GO:0007155","cell adhesion",1.1091577223710762,1,0.9918552211007083,0.01321071,"cell adhesion"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,2,0.8534247524356472,-0,"steroid metabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9480000801318849,0.15163831,"steroid metabolic process"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,1,0.8377826242161245,0.48590569,"steroid metabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9468844395095692,0.10901633,"steroid metabolic process"),
                     c("GO:0008203","cholesterol metabolic process",0.2200421431132646,2,0.837604392634995,0.45443823,"steroid metabolic process"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.8177359333941377,0.49264795,"steroid metabolic process"),
                     c("GO:0009805","coumarin biosynthetic process",0.00857483103959684,1,0.9321240528154487,0.15725843,"steroid metabolic process"),
                     c("GO:0042759","long-chain fatty acid biosynthetic process",0.04469313344093403,1,0.8440047120986782,0.3981042,"steroid metabolic process"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.8282662942985284,0.42785472,"steroid metabolic process"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,3,0.7958593720914472,0.00987901,"auxin-activated signaling pathway"),
                     c("GO:0002238","response to molecule of fungal origin",0.007682595099288114,2,0.8561750466668394,0.47317835,"auxin-activated signaling pathway"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8798364605318134,0.4896406,"auxin-activated signaling pathway"),
                     c("GO:0006952","defense response",1.1604144588756624,3,0.876429535144067,0.24614308,"auxin-activated signaling pathway"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8443530521025387,0.44951817,"auxin-activated signaling pathway"),
                     c("GO:0009624","response to nematode",0.0011042035946362127,1,0.9061332936456427,0.34321484,"auxin-activated signaling pathway"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.8802376082010779,0.4259419,"auxin-activated signaling pathway"),
                     c("GO:0010030","positive regulation of seed germination",0.0013630013121290752,1,0.9330367249858501,0.39302598,"auxin-activated signaling pathway"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.8960262652098409,0.17080919,"auxin-activated signaling pathway"),
                     c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,1,0.9239894723365195,0.33236762,"auxin-activated signaling pathway"),
                     c("GO:0032012","regulation of ARF protein signal transduction",0.05576721100946194,1,0.9122333943999571,0.46282204,"auxin-activated signaling pathway"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.7926893926416296,0.42828273,"auxin-activated signaling pathway"),
                     c("GO:0040008","regulation of growth",0.2209540969749061,1,0.9296589724366165,0.15374745,"auxin-activated signaling pathway"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8906418163906689,0.46472215,"auxin-activated signaling pathway"),
                     c("GO:0045927","positive regulation of growth",0.06802189911779062,1,0.9098907835228487,0.13260447,"auxin-activated signaling pathway"),
                     c("GO:0048209","regulation of vesicle targeting, to, from or within Golgi",2.4647401665986893E-06,1,0.9556732001053965,0.33066454,"auxin-activated signaling pathway"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,1,0.9131070277074601,0.42863919,"auxin-activated signaling pathway"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9279811128506251,0.1517125,"auxin-activated signaling pathway"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9165286189212768,0.3459764,"auxin-activated signaling pathway"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.9093436558702916,0.32241615,"auxin-activated signaling pathway"),
                     c("GO:1902882","regulation of response to oxidative stress",0.01609475328788944,1,0.9185275199482019,0.12060378,"auxin-activated signaling pathway"),
                     c("GO:1905421","regulation of plant organ morphogenesis",0.005429822587016912,1,0.9397297178437349,0.11290105,"auxin-activated signaling pathway"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9247599954019007,0.10716362,"auxin-activated signaling pathway"),
                     c("GO:0009793","embryo development ending in seed dormancy",0.0206816347379296,2,0.8812401528239566,-0,"embryo development ending in seed dormancy"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.9164516963298328,0.40853039,"embryo development ending in seed dormancy"),
                     c("GO:0010015","root morphogenesis",0.013467340270295237,1,0.8895735484931658,0.43918043,"embryo development ending in seed dormancy"),
                     c("GO:0010047","fruit dehiscence",0.002082705440775892,1,0.9484688509873693,0.46834396,"embryo development ending in seed dormancy"),
                     c("GO:0010090","trichome morphogenesis",0.006709022733481632,1,0.897637974036138,0.35446469,"embryo development ending in seed dormancy"),
                     c("GO:0015031","protein transport",3.093438694074183,4,0.9022621150389559,0,"protein transport"),
                     c("GO:0006891","intra-Golgi vesicle-mediated transport",0.1438718130046987,2,0.9110832468027599,0.2741348,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,2,0.9319397887570026,0.38566904,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9292971562616306,0.41259033,"protein transport"),
                     c("GO:0016135","saponin biosynthetic process",9.858960666394757E-06,1,0.9597459419647616,0.03120971,"saponin biosynthetic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,6,0.9001037157078688,0.05779371,"phosphorylation"),
                     c("GO:0044843","cell cycle G1/S phase transition",0.05461124787132716,1,0.9857270870712774,0.00959067,"cell cycle G1/S phase transition"),
                     c("GO:0000911","cytokinesis by cell plate formation",0.008712856488926366,1,0.9864782001597643,0.45372294,"cell cycle G1/S phase transition"),
                     c("GO:0050896","response to stimulus",17.567785530535815,2,1,-0,"response to stimulus"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,2,0.9746267827770626,0.01286368,"cell wall organization"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9526808770000171,0.05291306,"macromolecule depalmitoylation"),
                     c("GO:0006013","mannose metabolic process",0.050640551462936674,1,0.9464052548391848,0.33438746,"macromolecule depalmitoylation"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1,0.8772859143564911,0.4946872,"macromolecule depalmitoylation"),
                     c("GO:0006412","translation",4.38869169324396,1,0.8626879215497344,0.25211911,"macromolecule depalmitoylation"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,1,0.8766733169127588,0.43192082,"macromolecule depalmitoylation"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",0.0004411884898211653,1,0.9205940057980978,0.26905503,"macromolecule depalmitoylation"),
                     c("GO:0045489","pectin biosynthetic process",0.012338489273993038,1,0.9231630072637175,0.10510175,"macromolecule depalmitoylation"),
                     c("GO:1902066","regulation of cell wall pectin metabolic process",0.00012816648866313183,1,0.9466438616787388,0.09252635,"regulation of cell wall pectin metabolic process"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1,0.9193648063327486,0.17609066,"regulation of cell wall pectin metabolic process"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9397706342454366,0.12668641,"regulation of cell wall pectin metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_BP_Table.tsv", sep = "\t")
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

