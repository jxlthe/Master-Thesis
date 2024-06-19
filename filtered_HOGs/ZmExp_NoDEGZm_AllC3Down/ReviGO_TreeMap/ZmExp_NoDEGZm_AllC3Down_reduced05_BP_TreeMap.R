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
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.07909597668631853,2,0.9233570780811821,0.05769476,"malate metabolic process"),
                     c("GO:0019752","carboxylic acid metabolic process",8.649401753337283,1,0.8886076936146385,0.46114007,"malate metabolic process"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,3,0.928749783285233,-0,"monoatomic ion transport"),
                     c("GO:0008643","carbohydrate transport",0.7503112720356397,1,0.9294830324797846,0.49791143,"monoatomic ion transport"),
                     c("GO:0009904","chloroplast accumulation movement",0.008823769796423306,1,0.9222736982608623,0.36807496,"monoatomic ion transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9452620568785269,0.42144071,"monoatomic ion transport"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,1,0.9336799109899854,0.46287171,"monoatomic ion transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9542245166376163,0.43140125,"monoatomic ion transport"),
                     c("GO:0070588","calcium ion transmembrane transport",0.4224811119566813,2,0.9073379524549021,0.3228373,"monoatomic ion transport"),
                     c("GO:0090150","establishment of protein localization to membrane",0.5833670263314108,1,0.9198806965903903,0.32758039,"monoatomic ion transport"),
                     c("GO:0120029","proton export across plasma membrane",0.0197376392541223,1,0.925506493023376,0.47511322,"monoatomic ion transport"),
                     c("GO:1905039","carboxylic acid transmembrane transport",1.2848271482650644,1,0.9022470319690351,0.37138613,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,10,0.863263325923334,0,"defense response"),
                     c("GO:0000165","MAPK cascade",0.14246691110973742,1,0.8460728647544693,0.41973941,"defense response"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.873976849600759,0.42761529,"defense response"),
                     c("GO:0006955","immune response",0.6433982378290884,1,0.8943936616876339,0.30172547,"defense response"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7699816731780171,2,0.8145752817286029,0.33768695,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8497324615346508,0.41395528,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,2,0.8798559923119899,0.22702886,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9069191088322223,0.20406574,"defense response"),
                     c("GO:0010332","response to gamma radiation",0.006090372951665361,1,0.8883604686898008,0.49899247,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8631862260538792,0.47637699,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,3,0.867231571110003,0.45104908,"defense response"),
                     c("GO:0070887","cellular response to chemical stimulus",2.703999888790924,1,0.8376580120760221,0.48850507,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9060059266188695,0.26100449,"defense response"),
                     c("GO:0071456","cellular response to hypoxia",0.026404761404771757,1,0.8333950216785795,0.40112537,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.886965027591401,0.26816179,"defense response"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,1,0.9916900306794391,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,1,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9533799742893749,0.08346192,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,1,0.9510640654740616,0.21768792,"biosynthetic process"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9563859857628859,0.08334549,"regulation of meristem structural organization"),
                     c("GO:0015979","photosynthesis",0.228607115192195,1,0.9550362872050443,0.07668617,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,6,0.8999483207890893,-0,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9385033876884507,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9166532622516966,0.43492456,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,2,0.8897313976828648,0.21165279,"protein ubiquitination"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8806067671586962,0.37832162,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.8827698762370719,0.37292242,"protein ubiquitination"),
                     c("GO:0006508","proteolysis",5.2622572267907,3,0.8941160077956305,0.47291952,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9384646078262686,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9371013497368873,0.12094833,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.8801258810506879,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,8,0.8965881077900265,0.44883397,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8784345080678305,0.14129671,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9349726911387148,0.45767193,"protein ubiquitination"),
                     c("GO:0042350","GDP-L-fucose biosynthetic process",0.02421607213683212,1,0.8664500406119435,0.2978867,"protein ubiquitination"),
                     c("GO:0042853","L-alanine catabolic process",0.0402368832197236,1,0.8981498040660442,0.47848651,"protein ubiquitination"),
                     c("GO:0046373","L-arabinose metabolic process",0.04285443727665141,1,0.906154900117013,0.41352792,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,1,0.8755268463641728,0.41108531,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9212692120426993,0.20437656,"protein ubiquitination"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.8236390666048707,1,0.8733651035091083,0.39252566,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9098462069731325,0.45115494,"protein ubiquitination"),
                     c("GO:0022900","electron transport chain",0.8207017864535318,1,0.9293086623013649,0.06987243,"electron transport chain"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9925821627158591,0.00919456,"establishment of cell polarity"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,2,0.8581895699032914,0.07681302,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006665","sphingolipid metabolic process",0.4116436494441469,1,0.8629569790534699,0.45795093,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016132","brassinosteroid biosynthetic process",0.019789398797620875,1,0.8247705280825638,0.39632639,"phosphatidylinositol dephosphorylation"),
                     c("GO:0019363","pyridine nucleotide biosynthetic process",0.26556096451000916,1,0.8322391304268336,0.39508505,"phosphatidylinositol dephosphorylation"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8568413888929544,0.39387188,"phosphatidylinositol dephosphorylation"),
                     c("GO:0048364","root development",0.0376070054619628,2,0.9018649371825134,-0,"root development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9285704594032439,0.38663663,"root development"),
                     c("GO:0009880","embryonic pattern specification",0.03627851051216611,1,0.9228303562653469,0.45392292,"root development"),
                     c("GO:0010228","vegetative to reproductive phase transition of meristem",0.011640967806845608,1,0.9100733787957433,0.43860092,"root development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.9280550655776705,0.43353892,"root development"),
                     c("GO:0055046","microgametogenesis",0.0014049018949612527,1,0.9325707904457887,0.37648736,"root development"),
                     c("GO:0050896","response to stimulus",17.567785530535815,3,1,-0,"response to stimulus"),
                     c("GO:0070814","hydrogen sulfide biosynthetic process",0.03377679924306844,1,0.9377314863417591,0.03784896,"hydrogen sulfide biosynthetic process"),
                     c("GO:0000103","sulfate assimilation",0.10156701278519879,1,0.9540059327190739,0.48363934,"hydrogen sulfide biosynthetic process"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,4,0.9448533068263878,-0,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9695139034288462,0.25762726,"cell wall organization"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9704702672104824,-0,"intracellular auxin homeostasis"),
                     c("GO:0098771","inorganic ion homeostasis",1.0570555799893464,1,0.9537452556217786,0.45565197,"intracellular auxin homeostasis"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,2,0.9075848927416297,-0,"negative regulation of defense response to bacterium"),
                     c("GO:0033314","mitotic DNA replication checkpoint signaling",0.042962885843981745,1,0.8283305829772818,0.48461656,"negative regulation of defense response to bacterium"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,1,0.8909571607761799,0.47074713,"negative regulation of defense response to bacterium"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,2,0.9043708970367523,0.12505063,"negative regulation of defense response to bacterium"),
                     c("GO:0051513","regulation of monopolar cell growth",0.0030119124835835983,1,0.9353597894689105,0.41466435,"negative regulation of defense response to bacterium"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,1,0.9348753215790656,0.23293868,"negative regulation of defense response to bacterium"),
                     c("GO:1901001","negative regulation of response to salt stress",0.00037710524548959944,1,0.9213494043710891,0.43908873,"negative regulation of defense response to bacterium"),
                     c("GO:1901141","regulation of lignin biosynthetic process",0.00025140349699306625,1,0.9545320325178722,0.19785358,"negative regulation of defense response to bacterium"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9432801072575365,0.37546575,"negative regulation of defense response to bacterium"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,1,0.9244258774748108,0.37952976,"negative regulation of defense response to bacterium"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.9291745959451281,0.34953816,"negative regulation of defense response to bacterium"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9372987064370121,0.15519254,"negative regulation of defense response to bacterium"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9413769943469327,0.26299316,"negative regulation of defense response to bacterium"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_BP_Table.tsv", sep = "\t")
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

