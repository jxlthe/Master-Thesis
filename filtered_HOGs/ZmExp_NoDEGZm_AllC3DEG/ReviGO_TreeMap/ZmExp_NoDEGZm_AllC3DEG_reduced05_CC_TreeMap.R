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
                     c("GO:0005576","extracellular region",3.8687535442060237,20,0.9999272053216985,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9770013299882434,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,120,0.897241834188223,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,54,0.8676759780409261,0.17160779,"cytoplasm"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,2,0.9871015694677739,2.997E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,65,0.9402851631532837,0.3208057,"plant-type cell wall"),
                     c("GO:0009925","basal plasma membrane",0.13578741268972233,1,0.9499567674091494,0.18839421,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,14,0.9670987319565381,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,66,0.7267748257471457,0,"chloroplast"),
                     c("GO:0000138","Golgi trans cisterna",0.007541260384887046,1,0.7847126229203457,0.37607379,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,8,0.8097345480148721,0.13280054,"chloroplast"),
                     c("GO:0000785","chromatin",1.087181392930179,2,0.6310218292716503,0.45756389,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,147,0.7351725457812942,0.35145602,"chloroplast"),
                     c("GO:0005730","nucleolus",1.2114693153196516,11,0.6593840512562242,0.22801288,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,32,0.7631637114037019,0.21465144,"chloroplast"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,5,0.6487140970015214,0.22750481,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,8,0.7910318848950383,0.23143591,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,7,0.7950367744269272,0.17354881,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,30,0.7058158314534835,0.2503098,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,2,0.7805719258777527,0.38610723,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,5,0.7346995442253187,0.44106058,"chloroplast"),
                     c("GO:0005905","clathrin-coated pit",0.09808111076528711,1,0.8843034674173224,0.46642691,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,62,0.801418427855384,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,1,0.8262819239280932,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,9,0.7767429599132744,0.39454155,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,1,0.6864100932783908,0.39861008,"chloroplast"),
                     c("GO:0031977","thylakoid lumen",0.009951481990600534,2,0.7872051218152879,0.29998046,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,3,0.7976194431299694,0.2461316,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.7927351528465045,0.28709921,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7823312993438182,0.31004802,"chloroplast"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.737053290975465,0.486126,"chloroplast"),
                     c("GO:0043661","peribacteroid membrane",1.987811633578134E-05,1,0.8037970927187932,0.30368923,"chloroplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999393415696792,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8842947842593032,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,5,0.9999220057920939,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,189,0.9998955392513816,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,2,1,-0,"protein-containing complex"),
                     c("GO:0034045","phagophore assembly site membrane",0.08896450966078939,1,0.9106158180456937,0.07998094,"phagophore assembly site membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,1,0.9999304352273388,5.465E-05,"cell projection"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999447452210687,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,3,0.9951486090250801,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,2,0.9107301932146584,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0090406","pollen tube",0.001702063711251277,1,0.9926089325091942,2.369E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.993655243626757,0.27222822,"pollen tube"),
                     c("GO:0098552","side of membrane",0.7241349304670945,1,0.9680390950268467,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,12,0.7977674519751797,-0,"ribonucleoprotein complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,1,0.7854921537208012,0.25406433,"ribonucleoprotein complex"),
                     c("GO:0000159","protein phosphatase type 2A complex",0.06136622989309897,1,0.812960587135565,0.32835456,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7621299398938972,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,3,0.7627026089730917,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7270950971474648,0.4447901,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,5,0.6140553343058295,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.8160882818678357,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.8011306332882022,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.7764207224685161,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,1,0.7822423676723407,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0030131","clathrin adaptor complex",0.06838817448871372,1,0.7599273560528329,0.47444895,"ribonucleoprotein complex"),
                     c("GO:0033179","proton-transporting V-type ATPase, V0 domain",0.05699304429922707,1,0.7570345167166008,0.46766199,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7521800807563347,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.8054714902444956,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,2,0.7711540233433036,0.30930488,"ribonucleoprotein complex"),
                     c("GO:0034719","SMN-Sm protein complex",0.04507362879138419,1,0.7700242866731067,0.23993593,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.787512812863161,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7166713683204003,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.7298382771511723,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,6,0.7869422711005678,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.7145901294300672,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.8005408034760669,0.40438888,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_CC_Table.tsv", sep = "\t")
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

