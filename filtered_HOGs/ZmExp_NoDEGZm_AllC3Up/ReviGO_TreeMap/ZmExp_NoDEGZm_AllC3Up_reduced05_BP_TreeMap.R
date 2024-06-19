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
revigo.data <- rbind(c("GO:0006325","chromatin organization",1.3384303174082528,6,0.9037388951885091,0.01352937,"chromatin organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9333411730790718,0.39503669,"chromatin organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9279374065611924,0.35552134,"chromatin organization"),
                     c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.8820687665564247,0.33525093,"chromatin organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9579440959878853,0.21045842,"chromatin organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9256312676607006,0.45272646,"chromatin organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,2,0.9115965269804759,0.36480707,"chromatin organization"),
                     c("GO:0051017","actin filament bundle assembly",0.08024701034412013,1,0.9174737417503959,0.46988587,"chromatin organization"),
                     c("GO:0051289","protein homotetramerization",0.01763275115184702,1,0.9215957100689234,0.48332609,"chromatin organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.921292474100556,0.48565697,"chromatin organization"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,2,0.9264926397767697,0.4063955,"chromatin organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9543856736003223,0.20937319,"chromatin organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9499810451975432,0.22611843,"chromatin organization"),
                     c("GO:0006412","translation",4.38869169324396,8,0.7926785473312071,0.01596295,"translation"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,1,0.8677237896585783,0.33478559,"translation"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8520901314386792,0.39091196,"translation"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.934003880398652,0.15163831,"translation"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.9172294311001228,0.49989178,"translation"),
                     c("GO:0006368","transcription elongation by RNA polymerase II",0.1500706345236944,1,0.8393995899047425,0.41096565,"translation"),
                     c("GO:0006397","mRNA processing",1.242261085587905,8,0.8084440392980892,0.48035482,"translation"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.835936321514754,0.47678922,"translation"),
                     c("GO:0006572","tyrosine catabolic process",0.043362173750970734,1,0.8371953264870066,0.47697306,"translation"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.9324316490173363,0.14656365,"translation"),
                     c("GO:0006744","ubiquinone biosynthetic process",0.2844852395091539,1,0.8859151989131273,0.49014462,"translation"),
                     c("GO:0006747","FAD biosynthetic process",0.030434611577160615,1,0.8516574668367464,0.49383913,"translation"),
                     c("GO:0006914","autophagy",0.44134623319182764,1,0.8875365731855603,0.45718665,"translation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9486550285633618,0.23095606,"translation"),
                     c("GO:0009072","aromatic amino acid metabolic process",0.7377977862098182,1,0.8480798733104356,0.32628769,"translation"),
                     c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.8559943480922595,0.44469555,"translation"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8758672995845684,0.36563143,"translation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9357301655135775,0.23899188,"translation"),
                     c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8582844465109717,0.49160318,"translation"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,4,0.8540507583586597,0.11630623,"translation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,7,0.9024381441247201,0.41915608,"translation"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,7,0.8716226622696978,0.4385246,"translation"),
                     c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,1,0.8549385071091905,0.39483568,"translation"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9307485795768821,0.24495457,"translation"),
                     c("GO:0022900","electron transport chain",0.8207017864535318,2,0.9333472366556608,0.105264,"translation"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9078036745139673,0.31208096,"translation"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8912983901142879,0.27892474,"translation"),
                     c("GO:0032544","plastid translation",0.006610433126817685,1,0.8419312670230298,0.46103483,"translation"),
                     c("GO:0033611","oxalate catabolic process",0.006792823899145988,1,0.8951568888908242,0.3986442,"translation"),
                     c("GO:0042364","water-soluble vitamin biosynthetic process",1.0506004254930246,1,0.8683508000478061,0.2442778,"translation"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,1,0.8974335419432055,0.39728749,"translation"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8558057009916047,0.435373,"translation"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,2,0.9445158923116669,0.11952856,"translation"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,1,0.8494579766589408,0.43903701,"translation"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8945424620674443,0.25604839,"translation"),
                     c("GO:0046473","phosphatidic acid metabolic process",0.11459562930583944,1,0.8854244492520411,0.46395125,"translation"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9169646294953643,0.20101185,"translation"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9040388038942404,0.23016052,"translation"),
                     c("GO:0046940","nucleoside monophosphate phosphorylation",0.18958288413443797,1,0.8374780159780952,0.37868979,"translation"),
                     c("GO:0071040","nuclear polyadenylation-dependent antisense transcript catabolic process",0.008039982423444923,1,0.8356113046566447,0.40606321,"translation"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8447728017865375,0.49920408,"translation"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.814945114628871,0.31056631,"translation"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9340960408288801,0.20646638,"translation"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,1,0.8616523601909517,0.41679514,"translation"),
                     c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,1,0.859860881669243,0.37088281,"translation"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.951378202843326,0.11478523,"translation"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8638776036655083,0.33283067,"translation"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9295141548961214,0.11660613,"translation"),
                     c("GO:0006915","apoptotic process",0.3972840732335429,1,0.9854642682152334,0.01170435,"apoptotic process"),
                     c("GO:0007049","cell cycle",2.8073119376140743,2,0.9886804497897039,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9829888947796039,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,2,0.9951870710351908,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,4,1,-0,"metabolic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9533893830047817,0.08207669,"toxin metabolic process"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,4,0.8994341575880751,-0,"response to salt stress"),
                     c("GO:0006979","response to oxidative stress",0.8191810417707402,3,0.9046004999747203,0.44531989,"response to salt stress"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.9156091405660307,0.38955838,"response to salt stress"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,3,0.8908940499207098,0.49159695,"response to salt stress"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,3,0.8393855848938463,0.31063873,"response to salt stress"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8921612280569129,0.17275947,"response to salt stress"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9304287573707649,0.23963034,"response to salt stress"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8561489198440102,0.38492887,"response to salt stress"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8683695017144614,0.42625684,"response to salt stress"),
                     c("GO:0034059","response to anoxia",0.0003401341429906191,1,0.9239237437324571,0.45207873,"response to salt stress"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,1,0.9003240931273356,0.46099612,"response to salt stress"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8286460161569248,0.46171648,"response to salt stress"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9116752148000451,0.27995591,"response to salt stress"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9784320132590988,0.48813486,"response to salt stress"),
                     c("GO:0051607","defense response to virus",0.17287441054506547,3,0.8905814624060662,0.36620957,"response to salt stress"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9346344614831905,0.43018112,"response to salt stress"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.9262375341932602,0.16037728,"response to salt stress"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9591183995889442,0.06965752,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,5,0.9065471644188652,-0,"flower development"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.9430738021597689,0.4229884,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,4,0.9356400235231461,0.4203896,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9464499525346355,0.48835672,"flower development"),
                     c("GO:0009944","polarity specification of adaxial/abaxial axis",0.0024376280247661035,2,0.9434678225980172,0.38300054,"flower development"),
                     c("GO:0010073","meristem maintenance",0.03536902139069119,1,0.9377404342222467,0.39934449,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9293437414221889,0.29044182,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.941796478243767,0.39539149,"flower development"),
                     c("GO:0048367","shoot system development",0.08373708242002387,1,0.9303744015693652,0.48899493,"flower development"),
                     c("GO:0061137","bud dilation",9.612486649734888E-05,1,0.9538573069389013,0.29194579,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9539044179192605,0.28605225,"flower development"),
                     c("GO:0015031","protein transport",3.093438694074183,12,0.8911444677877459,0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,3,0.9371998749958568,0.42145531,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9565182844571332,0.21356797,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9416038378199687,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,6,0.9317059101071705,0.38566904,"protein transport"),
                     c("GO:0048250","iron import into the mitochondrion",0.008732574410259155,1,0.9473693396780738,0.4861705,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9386774339596206,0.49534489,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,3,0.9277070548179531,0.35175462,"protein transport"),
                     c("GO:0032259","methylation",2.6278542060840238,4,0.9630381097319232,0.06915605,"methylation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,0.9890168330900395,-0,"carbohydrate homeostasis"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,4,0.8817156359903194,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.00031548674132463224,1,0.9386595886188672,0.44967237,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9327569894298242,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,4,0.89929908295257,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9283701550229393,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0009926","auxin polar transport",0.004261535748049133,1,0.9151122521821295,0.4567727,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9409409043106762,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,1,0.9499354952025029,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9441846315281429,0.19548819,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9090817532027623,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0032876","negative regulation of DNA endoreduplication",0.0004264000488215732,1,0.9325263834384296,0.48325207,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.937124668130686,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,1,0.9467761303273207,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9349975528089349,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9340602359715098,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,3,0.8886162915959869,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,2,0.924959055124064,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,1,0.932876854775429,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.935688272889682,0.33352541,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.9507973047224163,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,1,0.9278657196743837,0.47223923,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9370307576793555,0.38691897,"positive regulation of DNA-templated transcription"),
                     c("GO:1901000","regulation of response to salt stress",0.0019915100546117406,1,0.9391213556476719,0.47172255,"positive regulation of DNA-templated transcription"),
                     c("GO:1902183","regulation of shoot apical meristem development",0.0016217990296219374,1,0.9357692310781289,0.38065657,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9278599682057307,0.454995,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,2,0.9374785388067236,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,1,0.9395637552639283,0.3403934,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9403832387750986,0.12485439,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0051301","cell division",1.5693197819947182,1,0.9893374613953471,0.01381155,"cell division"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_BP_Table.tsv", sep = "\t")
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

