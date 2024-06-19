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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,4,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,5,0.9999224840330628,6.017E-05,"extracellular region"),
                     c("GO:0005739","mitochondrion",4.856981674016684,20,0.7429688158951828,0,"mitochondrion"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,4,0.7849669552020982,0.15656526,"mitochondrion"),
                     c("GO:0000785","chromatin",1.087181392930179,2,0.6267578627158689,0.45756389,"mitochondrion"),
                     c("GO:0005634","nucleus",16.5161752456724,75,0.7108148502548343,0.35145602,"mitochondrion"),
                     c("GO:0005730","nucleolus",1.2114693153196516,8,0.6525206369991348,0.22801288,"mitochondrion"),
                     c("GO:0005737","cytoplasm",43.076660800468886,64,0.8844138649186254,0.12447513,"mitochondrion"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,5,0.625087928442393,0.22750481,"mitochondrion"),
                     c("GO:0005773","vacuole",1.35235298008496,4,0.7746528732709345,0.23143591,"mitochondrion"),
                     c("GO:0005777","peroxisome",0.7476308639759881,3,0.7769793003190435,0.21411812,"mitochondrion"),
                     c("GO:0005783","endoplasmic reticulum",3.094528245337302,13,0.692503713783265,0.26090848,"mitochondrion"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,2,0.7819153382655257,0.38610723,"mitochondrion"),
                     c("GO:0005829","cytosol",14.048080989011824,27,0.8532490738086973,0.19224382,"mitochondrion"),
                     c("GO:0009507","chloroplast",0.6990388085931706,14,0.7294336041196315,0.21231674,"mitochondrion"),
                     c("GO:0009536","plastid",0.7624897559369845,10,0.7863656568074683,0.21465144,"mitochondrion"),
                     c("GO:0009579","thylakoid",0.26089530737804617,1,0.7779106746861505,0.39454155,"mitochondrion"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6909463303728424,0.40758865,"mitochondrion"),
                     c("GO:0031410","cytoplasmic vesicle",2.424891655569294,4,0.6894561784958844,0.25147478,"mitochondrion"),
                     c("GO:0031982","vesicle",2.6690048632308567,1,0.7825608346576599,0.2461316,"mitochondrion"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,1,0.7832725491827927,0.4052793,"mitochondrion"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7358403967214099,0.486126,"mitochondrion"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,1,0.9030105125468458,0.10673752,"mitochondrion"),
                     c("GO:0005886","plasma membrane",17.177321395000487,23,0.959322598108123,0.00015294,"plasma membrane"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,5,0.9795750356957036,3.812E-05,"plasmodesma"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999360816469345,4.693E-05,"cell surface"),
                     c("GO:0012505","endomembrane system",6.893425119210306,1,0.999916561578383,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,69,0.9998855005769441,0.00010118,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999420397195233,4.148E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,2,0.9999461291247308,3.787E-05,"apoplast"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9957105861594396,2.697E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.9957188065510704,0.27222822,"pollen tube"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,10,0.7917625287215776,-0,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7433502552256501,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,3,0.7515619719814391,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,5,0.5914476605766589,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8088404021989811,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7931079870659599,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005942","phosphatidylinositol 3-kinase complex",0.07251039886384637,1,0.7616223441935326,0.46964777,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.772242956869616,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0030131","clathrin adaptor complex",0.06838817448871372,1,0.7570390531220923,0.40377927,"ribonucleoprotein complex"),
                     c("GO:0033179","proton-transporting V-type ATPase, V0 domain",0.05699304429922707,1,0.7650237145721609,0.32669992,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7278739987734285,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.7991344382389276,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7445202547292793,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.778691591901441,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7012162813786493,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.7152465682937388,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,5,0.7773723415349952,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.6989975611723888,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7945554103286121,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_CC_Table.tsv", sep = "\t")
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
