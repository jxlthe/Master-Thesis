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
revigo.data <- rbind(c("GO:0005575","cellular_component",100,9,1,-0,"cellular_component"),
                     c("GO:0005576","extracellular region",3.8687535442060237,30,0.9999298415024575,4.598E-05,"extracellular region"),
                     c("GO:0005628","prospore membrane",0.0069225540139358525,1,0.9766709014171154,0.04954734,"prospore membrane"),
                     c("GO:0005737","cytoplasm",43.076660800468886,204,0.9006467084488148,0.08291202,"cytoplasm"),
                     c("GO:0005829","cytosol",14.048080989011824,86,0.8711955549649424,0.17160779,"cytoplasm"),
                     c("GO:0009504","cell plate",0.0007255512462560189,2,0.9395571191954046,0.05708197,"cell plate"),
                     c("GO:0009505","plant-type cell wall",0.04809510247442295,4,0.9802540780053186,2.997E-05,"plant-type cell wall"),
                     c("GO:0005886","plasma membrane",17.177321395000487,128,0.9365244522321454,0.34881718,"plant-type cell wall"),
                     c("GO:0009898","cytoplasmic side of plasma membrane",0.2657902935257323,1,0.9550365138301414,0.23284101,"plant-type cell wall"),
                     c("GO:0009925","basal plasma membrane",0.13578741268972233,1,0.9463298922396219,0.18839421,"plant-type cell wall"),
                     c("GO:0009506","plasmodesma",0.10230521048664064,25,0.9764118381717231,3.188E-05,"plasmodesma"),
                     c("GO:0009507","chloroplast",0.6990388085931706,106,0.7353647898965974,0,"chloroplast"),
                     c("GO:0000138","Golgi trans cisterna",0.007541260384887046,2,0.7868433828079426,0.37607379,"chloroplast"),
                     c("GO:0000325","plant-type vacuole",0.04066068696484073,11,0.8161229772690657,0.13280054,"chloroplast"),
                     c("GO:0000785","chromatin",1.087181392930179,4,0.6154547100122836,0.45756389,"chloroplast"),
                     c("GO:0001673","male germ cell nucleus",0.006684016617906476,1,0.8400958319844106,0.21436526,"chloroplast"),
                     c("GO:0005634","nucleus",16.5161752456724,279,0.7434541920863925,0.35145602,"chloroplast"),
                     c("GO:0005730","nucleolus",1.2114693153196516,25,0.6574244508381973,0.20773614,"chloroplast"),
                     c("GO:0005739","mitochondrion",4.856981674016684,49,0.7694728451028692,0.2503098,"chloroplast"),
                     c("GO:0005743","mitochondrial inner membrane",1.1915116865185273,12,0.6508642725908665,0.20731432,"chloroplast"),
                     c("GO:0005773","vacuole",1.35235298008496,12,0.7958642379536632,0.21057364,"chloroplast"),
                     c("GO:0005777","peroxisome",0.7476308639759881,9,0.8012526236944911,0.17354881,"chloroplast"),
                     c("GO:0005794","Golgi apparatus",2.349961096041566,59,0.707932123912852,0.19658734,"chloroplast"),
                     c("GO:0005811","lipid droplet",0.20804685033482948,4,0.7748705589208366,0.38610723,"chloroplast"),
                     c("GO:0005874","microtubule",0.7782307393103814,10,0.7214407684248566,0.44106058,"chloroplast"),
                     c("GO:0005905","clathrin-coated pit",0.09808111076528711,3,0.8806413166374408,0.46642691,"chloroplast"),
                     c("GO:0009536","plastid",0.7624897559369845,100,0.8057497688689464,0.17236349,"chloroplast"),
                     c("GO:0009574","preprophase band",0.001202626038314771,2,0.8190122800115598,0.49340982,"chloroplast"),
                     c("GO:0009579","thylakoid",0.26089530737804617,18,0.7710315941262023,0.39454155,"chloroplast"),
                     c("GO:0012511","monolayer-surrounded lipid storage body",0.0046142077544432435,1,0.8242011923160282,0.28396619,"chloroplast"),
                     c("GO:0017177","glucosidase II complex",0.015944734065838607,2,0.6825635969117441,0.39861008,"chloroplast"),
                     c("GO:0031977","thylakoid lumen",0.009951481990600534,4,0.7865019059678459,0.29998046,"chloroplast"),
                     c("GO:0031982","vesicle",2.6690048632308567,6,0.8028487298659475,0.22189912,"chloroplast"),
                     c("GO:0032044","DSIF complex",0.014580598332295611,1,0.6977513011508626,0.49997569,"chloroplast"),
                     c("GO:0032578","aleurone grain membrane",6.95734071752347E-05,2,0.788446534024034,0.28709921,"chloroplast"),
                     c("GO:0042406","extrinsic component of endoplasmic reticulum membrane",0.00027580886415896606,1,0.7810150470710108,0.31004802,"chloroplast"),
                     c("GO:0042470","melanosome",0.049131249288425556,1,0.7757444321211127,0.3330034,"chloroplast"),
                     c("GO:0042645","mitochondrial nucleoid",0.035025240983646726,1,0.7274865822451015,0.486126,"chloroplast"),
                     c("GO:0043661","peribacteroid membrane",1.987811633578134E-05,1,0.7980277845638879,0.30368923,"chloroplast"),
                     c("GO:0070971","endoplasmic reticulum exit site",0.06951625759076932,1,0.77292085047214,0.45185985,"chloroplast"),
                     c("GO:0009524","phragmoplast",0.002879842104146322,3,0.9348019928088872,0.06218567,"phragmoplast"),
                     c("GO:0009986","cell surface",0.6578761991908514,1,0.9999411529088681,3.782E-05,"cell surface"),
                     c("GO:0010168","ER body",0.00030314127412066545,1,0.8857840483760228,0.09598239,"ER body"),
                     c("GO:0012505","endomembrane system",6.893425119210306,9,0.9999250488252842,0.00011166,"endomembrane system"),
                     c("GO:0016020","membrane",49.2542153160787,378,0.9999011439048419,9.537E-05,"membrane"),
                     c("GO:0032991","protein-containing complex",19.75554538569037,6,1,-0,"protein-containing complex"),
                     c("GO:0034045","phagophore assembly site membrane",0.08896450966078939,1,0.9108104635559916,0.07998094,"phagophore assembly site membrane"),
                     c("GO:0042995","cell projection",2.575965339721232,3,0.999932835127444,5.465E-05,"cell projection"),
                     c("GO:0044297","cell body",0.22864057885869896,1,0.9999462429024586,3.42E-05,"cell body"),
                     c("GO:0048046","apoplast",0.0960361495472436,5,0.9965175146401201,3.171E-05,"apoplast"),
                     c("GO:0048471","perinuclear region of cytoplasm",0.2584900553014166,3,0.9120657479855075,0.08779222,"perinuclear region of cytoplasm"),
                     c("GO:0070469","respirasome",0.5984977809313304,3,0.968203139496543,0.0715415,"respirasome"),
                     c("GO:0090406","pollen tube",0.001702063711251277,4,0.9947034317355141,2.369E-05,"pollen tube"),
                     c("GO:0090395","plant cell papilla",2.4847645419726674E-05,1,0.995454037968868,0.27222822,"pollen tube"),
                     c("GO:0098552","side of membrane",0.7241349304670945,3,0.9677007152814887,4.617E-05,"side of membrane"),
                     c("GO:1990904","ribonucleoprotein complex",4.315842197772249,17,0.7938877775225776,-0,"ribonucleoprotein complex"),
                     c("GO:0000145","exocyst",0.08239479221181366,2,0.7836928687770224,0.25406433,"ribonucleoprotein complex"),
                     c("GO:0000148","1,3-beta-D-glucan synthase complex",0.012428792238947283,2,0.7603198759962925,0.41849798,"ribonucleoprotein complex"),
                     c("GO:0000172","ribonuclease MRP complex",0.04151047643819539,1,0.7658118738451125,0.47495898,"ribonucleoprotein complex"),
                     c("GO:0000178","exosome (RNase complex)",0.10685978389207654,1,0.7658931000593635,0.34132399,"ribonucleoprotein complex"),
                     c("GO:0000502","proteasome complex",0.37652382533874423,4,0.7616170306174925,0.49612005,"ribonucleoprotein complex"),
                     c("GO:0000974","Prp19 complex",0.0635900941581645,1,0.8509571458457826,0.24779784,"ribonucleoprotein complex"),
                     c("GO:0005667","transcription regulator complex",0.7880231963702957,1,0.8214255995273086,0.32589738,"ribonucleoprotein complex"),
                     c("GO:0005677","chromatin silencing complex",0.009745246533616803,1,0.7261137000342258,0.4447901,"ribonucleoprotein complex"),
                     c("GO:0005681","spliceosomal complex",0.7626239332222511,12,0.6207102252827627,0.3245659,"ribonucleoprotein complex"),
                     c("GO:0005851","eukaryotic translation initiation factor 2B complex",0.02637826037758184,1,0.813311779148287,0.22864326,"ribonucleoprotein complex"),
                     c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.08561256229366826,2,0.7984871088787376,0.255018,"ribonucleoprotein complex"),
                     c("GO:0005950","anthranilate synthase complex",0.00014660110797638738,1,0.8267219358756437,0.15701804,"ribonucleoprotein complex"),
                     c("GO:0005951","carbamoyl-phosphate synthase complex",0.027533675889599128,1,0.771298490982265,0.31126482,"ribonucleoprotein complex"),
                     c("GO:0009349","riboflavin synthase complex",0.02715350691467731,1,0.801646627727481,0.43971777,"ribonucleoprotein complex"),
                     c("GO:0010330","cellulose synthase complex",0.0019629639881584074,1,0.8288151557394036,0.37568323,"ribonucleoprotein complex"),
                     c("GO:0016272","prefoldin complex",0.03795726314317447,1,0.8558414471368166,0.23619435,"ribonucleoprotein complex"),
                     c("GO:0016281","eukaryotic translation initiation factor 4F complex",0.039140011065153454,2,0.7844583658805595,0.42147812,"ribonucleoprotein complex"),
                     c("GO:0016602","CCAAT-binding factor complex",0.03082101937862897,2,0.7042485204988549,0.48430429,"ribonucleoprotein complex"),
                     c("GO:0030289","protein phosphatase 4 complex",0.015278817168589934,1,0.818783800666579,0.29979865,"ribonucleoprotein complex"),
                     c("GO:0033290","eukaryotic 48S preinitiation complex",0.07039586423862765,2,0.7595578125779482,0.49563216,"ribonucleoprotein complex"),
                     c("GO:0033588","elongator holoenzyme complex",0.03708262602440009,1,0.8016531148164839,0.41997104,"ribonucleoprotein complex"),
                     c("GO:0034702","monoatomic ion channel complex",0.5134268973078122,3,0.7618393545642829,0.47090253,"ribonucleoprotein complex"),
                     c("GO:0035657","eRF1 methyltransferase complex",0.011842387807041733,1,0.7817374230316501,0.41725271,"ribonucleoprotein complex"),
                     c("GO:0045252","oxoglutarate dehydrogenase complex",0.07266196950090673,1,0.7324805062638388,0.46971574,"ribonucleoprotein complex"),
                     c("GO:0045273","respiratory chain complex II",0.05444864540824706,1,0.7586647868595982,0.39997785,"ribonucleoprotein complex"),
                     c("GO:0046930","pore complex",0.06221601936645362,3,0.7967552314432623,0.24728353,"ribonucleoprotein complex"),
                     c("GO:0070390","transcription export complex 2",0.01837731855242985,1,0.7158452399269217,0.46572404,"ribonucleoprotein complex"),
                     c("GO:0071439","clathrin complex",0.011926869801468804,1,0.7703665716373392,0.36298129,"ribonucleoprotein complex"),
                     c("GO:0071819","DUBm complex",0.008179844872174023,1,0.728818219916503,0.43933875,"ribonucleoprotein complex"),
                     c("GO:0080008","Cul4-RING E3 ubiquitin ligase complex",0.07457026866914172,9,0.7813093676340951,0.25161359,"ribonucleoprotein complex"),
                     c("GO:0089701","U2AF complex",0.020742814396387827,1,0.7137965779612403,0.46994592,"ribonucleoprotein complex"),
                     c("GO:0140535","intracellular protein-containing complex",3.7132246772303286,1,0.7966516257371146,0.40438888,"ribonucleoprotein complex"),
                     c("GO:0160064","multi-pass translocon complex",0.0006336149582030303,1,0.8478591733073518,0.17222646,"ribonucleoprotein complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_CC_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "CC" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, CC, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_CC_Table.tsv", sep = "\t")
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

