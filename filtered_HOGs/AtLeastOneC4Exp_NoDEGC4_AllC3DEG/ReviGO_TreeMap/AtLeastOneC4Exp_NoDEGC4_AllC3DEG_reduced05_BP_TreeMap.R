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
revigo.data <- rbind(c("GO:0002376","immune system process",0.9427113541805001,2,1,-0,"immune system process"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,7,0.9186356768617372,0.01352937,"chromatin organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,2,0.9390813746700919,0.36382358,"chromatin organization"),
                     c("GO:0007032","endosome organization",0.07057290519022028,1,0.9401718046072549,0.34948919,"chromatin organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,1,0.9350750787398051,0.41402298,"chromatin organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9372819581737211,0.38649952,"chromatin organization"),
                     c("GO:0010027","thylakoid membrane organization",0.08263041408522105,3,0.929376668550639,0.32157946,"chromatin organization"),
                     c("GO:0010405","arabinogalactan protein metabolic process",0.00041161160782198105,1,0.9116607404941047,0.48519003,"chromatin organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9384416533666503,0.37383592,"chromatin organization"),
                     c("GO:0042254","ribosome biogenesis",2.136224808753416,5,0.9175878132704369,0.40498561,"chromatin organization"),
                     c("GO:0051017","actin filament bundle assembly",0.08024701034412013,1,0.9292945375911447,0.45180791,"chromatin organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9401720574956033,0.4524177,"chromatin organization"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,5,0.9262728972101378,0.4063955,"chromatin organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9584379981320287,0.20937319,"chromatin organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9543707006702298,0.22611843,"chromatin organization"),
                     c("GO:0006952","defense response",1.1604144588756624,20,0.8957636459317844,-0,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,3,0.8903484054146903,0.47274124,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,4,0.8943590811635584,0.49159695,"defense response"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,6,0.835229802251021,0.23262547,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8867081186019515,0.25765651,"defense response"),
                     c("GO:0010188","response to microbial phytotoxin",0.00018485551249490168,1,0.9336715752723834,0.27319058,"defense response"),
                     c("GO:0010196","nonphotochemical quenching",0.0005693549784842973,1,0.9231741544045722,0.37622087,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9226374812104545,0.20406574,"defense response"),
                     c("GO:0010343","singlet oxygen-mediated programmed cell death",0.0019717921332789512,1,0.8852394962381888,0.44410187,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8480469496818611,0.38492887,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,2,0.8595221172259188,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8927239928546409,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,5,0.8212463874064555,0.40205681,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9073299609359704,0.46171648,"defense response"),
                     c("GO:0042631","cellular response to water deprivation",0.00080843477464437,3,0.8633771504728337,0.35855228,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,7,0.8864965580958846,0.45104908,"defense response"),
                     c("GO:0051103","DNA ligation involved in DNA repair",0.04254388001565997,1,0.8144245704915166,0.43114998,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9236295310126758,0.4243543,"defense response"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.7957322799946543,0.41244237,"defense response"),
                     c("GO:0071461","cellular response to redox state",0.0002957688199918427,1,0.9516903714411864,0.16681679,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,1,0.9260116244231166,0.31063873,"defense response"),
                     c("GO:1990414","replication-born double-strand break repair via sister chromatid exchange",0.027676567330736677,1,0.8186751991451701,0.39308191,"defense response"),
                     c("GO:0007049","cell cycle",2.8073119376140743,3,0.9887995062983475,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9821067489482855,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,3,0.9968868658947184,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,7,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,8,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9916705180602594,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9901453047902254,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,3,0.9485961331730296,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,2,0.9455384672030147,0.21768792,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.916406665515398,0.16846091,"biosynthetic process"),
                     c("GO:0071704","organic substance metabolic process",50.871943734517124,1,0.9430877798011215,0.29474266,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9338959248200909,0.17227516,"biosynthetic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9530544271030384,0.06359719,"toxin metabolic process"),
                     c("GO:0009813","flavonoid biosynthetic process",0.06199314467029023,3,0.9280455846162399,0.05649883,"flavonoid biosynthetic process"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,2,0.9003447070443624,0.20397686,"flavonoid biosynthetic process"),
                     c("GO:0009693","ethylene biosynthetic process",0.004658358914871523,1,0.9362773595837025,0.10463112,"flavonoid biosynthetic process"),
                     c("GO:0009697","salicylic acid biosynthetic process",0.021366832504244038,2,0.8527037207100914,0.31042729,"flavonoid biosynthetic process"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,1,0.8928899941250885,0.34289816,"flavonoid biosynthetic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,1,0.9006858251553228,0.3067044,"flavonoid biosynthetic process"),
                     c("GO:0042372","phylloquinone biosynthetic process",0.033601802691239926,2,0.8853206771987169,0.11806969,"flavonoid biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9568121087752541,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,9,0.8908392724438149,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,3,0.9230354927540166,0.4203896,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9421605704579668,0.48835672,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,1,0.9381000734135579,0.48720779,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.852418032883764,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9294022137406913,0.29044182,"flower development"),
                     c("GO:0022038","corpus callosum development",0.002092564401442287,1,0.9302073632677357,0.39347824,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.928925869834166,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9236609156308953,0.42767441,"flower development"),
                     c("GO:0048262","determination of dorsal/ventral asymmetry",0.0029848003417510126,2,0.9271829185931671,0.38718574,"flower development"),
                     c("GO:0090626","plant epidermis morphogenesis",0.016804598455869863,1,0.9183486446093306,0.47435352,"flower development"),
                     c("GO:1905392","plant organ morphogenesis",0.021231271795081108,1,0.9034616319762211,0.44850026,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,1,0.9321525162187911,0.34979162,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.993597558808241,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9549422791824785,0.04613723,"chlororespiration"),
                     c("GO:0015031","protein transport",3.093438694074183,14,0.9103422405483595,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,7,0.9460843053691552,0.42145531,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9507515579714438,0.35562027,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,7,0.9398786887052012,0.38566904,"protein transport"),
                     c("GO:0097339","glycolate transmembrane transport",0.0003376694028240204,1,0.951599441018805,0.41757087,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9485023667597943,0.49534489,"protein transport"),
                     c("GO:1901975","glycerate transmembrane transport",0.00010351908699714494,1,0.952832221721286,0.26680356,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,4,0.9283663106447313,0.35175462,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,3,0.9492701579522043,0.04477571,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,20,0.8787817727983146,0,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8608063849588752,0.40302693,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,5,0.9350254592557123,0.15163831,"protein ubiquitination"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",18.82728500884891,2,0.8436208381558299,0.27206623,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,8,0.8195192732168812,0.48035482,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,12,0.8074755633777055,0.4385246,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,7,0.8493320361736643,0.47678922,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,7,0.8528425168304801,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9349826835902781,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,9,0.9334748490371114,0.12094833,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,2,0.907561429820102,0.45718665,"protein ubiquitination"),
                     c("GO:0008652","amino acid biosynthetic process",2.679426429329935,4,0.8187448752720609,0.49906638,"protein ubiquitination"),
                     c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.8641141306992849,0.46840296,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.8636457055124153,0.16277389,"protein ubiquitination"),
                     c("GO:0018258","protein O-linked glycosylation via hydroxyproline",0.0002957688199918427,1,0.8965164020303701,0.36854741,"protein ubiquitination"),
                     c("GO:0018364","peptidyl-glutamine methylation",1.9717921332789515E-05,1,0.9407597825211361,0.31815889,"protein ubiquitination"),
                     c("GO:0030212","hyaluronan metabolic process",0.023962203899672456,1,0.916362091979296,0.14846587,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,3,0.8955529307937089,0.14129671,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,1,0.8509575240184966,0.28833637,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,1,0.9126631931540495,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,2,0.9141606037501719,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.8740122109151391,0.26614879,"protein ubiquitination"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8937726619685422,0.16428511,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,4,0.9073240796496012,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,1,0.9074103442869037,0.41108531,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8553631478800025,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.88211423643748,0.43527038,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,1,0.8748603118957694,0.32141303,"protein ubiquitination"),
                     c("GO:0032259","methylation",2.6278542060840238,4,0.9630452431500669,0.05843136,"methylation"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9743453680080014,0.03897815,"collagen metabolic process"),
                     c("GO:0045454","cell redox homeostasis",0.337647220162521,2,0.9782047884873202,-0,"cell redox homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9847712188746586,0.42883749,"cell redox homeostasis"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,9,0.8795963347957488,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9328147254894259,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,5,0.8997569399837453,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0006450","regulation of translational fidelity",0.29561600610151356,1,0.9347178755935446,0.18254504,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9451145884433402,0.38949219,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,2,0.9481091027917375,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9453653482900994,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9080353016145205,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0030512","negative regulation of transforming growth factor beta receptor signaling pathway",0.013122276646971421,1,0.9070792414433485,0.42782958,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.935007086825738,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9208022567925225,0.31416992,"positive regulation of DNA-templated transcription"),
                     c("GO:0042548","regulation of photosynthesis, light reaction",0.0070861279789712316,1,0.9351907319277292,0.23746411,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,2,0.9449907287239626,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0042981","regulation of apoptotic process",0.6083841390223874,1,0.9291369098141573,0.22270387,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9328435266119133,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0043484","regulation of RNA splicing",0.23003419974865566,1,0.9199613319251218,0.37266406,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9322503081007453,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,4,0.8854071590941985,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,1,0.9418206371631411,0.41192877,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9393485456220058,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9306901176114808,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080038","positive regulation of cytokinin-activated signaling pathway",0.0001996439534944938,1,0.9225478867658441,0.35146486,"positive regulation of DNA-templated transcription"),
                     c("GO:0090070","positive regulation of ribosome biogenesis",0.00022675609532707942,1,0.9433122010325001,0.35407617,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.946765199196527,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,3,0.9118699682618242,0.37396531,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9412071089323403,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,2,0.912096614000361,0.48792077,"positive regulation of DNA-templated transcription"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9439336984968786,0.49983935,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,3,0.9143430387136154,0.13990144,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,2,0.9270459558020184,0.36922831,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.9224841775726509,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000032","regulation of secondary shoot formation",0.007307954593965113,1,0.93652840146289,0.46431004,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9403446562038746,0.41762823,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9225784256688021,0.49708545,"positive regulation of DNA-templated transcription"),
                     c("GO:2000280","regulation of root development",0.008382581306602141,1,0.93764660082577,0.1419752,"positive regulation of DNA-templated transcription"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9397099562114293,0.15519254,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9340937151495133,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.926396799635331,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,2,0.8901036103505938,0.07681302,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,2,0.9113751604809034,0.3693131,"phosphatidylinositol dephosphorylation"),
                     c("GO:0006221","pyrimidine nucleotide biosynthetic process",0.5101395959817636,2,0.8020872555179069,0.41640622,"phosphatidylinositol dephosphorylation"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9355979539909652,0.31764343,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016310","phosphorylation",5.235381700014107,35,0.896548239281793,0.47540491,"phosphatidylinositol dephosphorylation"),
                     c("GO:0016311","dephosphorylation",0.5642283189377719,1,0.9168936756356868,0.32207069,"phosphatidylinositol dephosphorylation"),
                     c("GO:1901137","carbohydrate derivative biosynthetic process",5.186311188037295,1,0.8715682360898135,0.47912132,"phosphatidylinositol dephosphorylation"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,5,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,1,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,2,0.9894456293343421,0.01381155,"cell division"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9518651679962883,0.09043985,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9263950558989176,0.08958079,"intrachromosomal DNA recombination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8815458283936144,0.15020607,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

