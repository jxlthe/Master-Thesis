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
revigo.data <- rbind(c("GO:0002376","immune system process",0.9427113541805001,1,1,-0,"immune system process"),
c("GO:0006325","chromatin organization",1.3384303174082528,8,0.8866644093239825,0.01352937,"chromatin organization"),
c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9170930274600421,0.39503669,"chromatin organization"),
c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9090306703822497,0.38546222,"chromatin organization"),
c("GO:0009658","chloroplast organization",0.0906284959258338,2,0.9171762363863234,0.32420414,"chromatin organization"),
c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9487168285546999,0.21045842,"chromatin organization"),
c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9114009532857971,0.45272646,"chromatin organization"),
c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,3,0.8947949768582132,0.36480707,"chromatin organization"),
c("GO:0043622","cortical microtubule organization",0.008850881938255893,2,0.9250852165518301,0.26893696,"chromatin organization"),
c("GO:0051289","protein homotetramerization",0.01763275115184702,1,0.9066141946465268,0.48332609,"chromatin organization"),
c("GO:0061025","membrane fusion",0.31015797308444587,1,0.9134725287105906,0.36373751,"chromatin organization"),
c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9062544766699169,0.48565697,"chromatin organization"),
c("GO:0071555","cell wall organization",0.8943925879544993,3,0.9055404150789912,0.4063955,"chromatin organization"),
c("GO:0080119","ER body organization",0.00018239077232830298,1,0.9432061285641596,0.20937319,"chromatin organization"),
c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9384183206495306,0.22611843,"chromatin organization"),
c("GO:0006397","mRNA processing",1.242261085587905,17,0.8085872367994686,0.01340109,"mRNA processing"),
c("GO:0000373","Group II intron splicing",0.009405448475740597,2,0.8673451110596677,0.4354711,"mRNA processing"),
c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8511791646442186,0.45078305,"mRNA processing"),
c("GO:0006012","galactose metabolic process",0.12319264300693568,1,0.9151486398359728,0.38165312,"mRNA processing"),
c("GO:0006412","translation",4.38869169324396,11,0.795526755579601,0.48035482,"mRNA processing"),
c("GO:0006457","protein folding",1.174377211919444,7,0.8399956489397382,0.40866092,"mRNA processing"),
c("GO:0006491","N-glycan processing",0.05083033645576476,2,0.9014579330039085,0.27495136,"mRNA processing"),
c("GO:0008652","amino acid biosynthetic process",2.679426429329935,4,0.8274795977949865,0.49906638,"mRNA processing"),
c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8754694014898333,0.367811,"mRNA processing"),
c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9347135507485762,0.32178408,"mRNA processing"),
c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8565169445990376,0.21599583,"mRNA processing"),
c("GO:0016567","protein ubiquitination",1.2678648064385323,14,0.8678840100509494,0.19691418,"mRNA processing"),
c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,1,0.8562525213532106,0.35937708,"mRNA processing"),
c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9285053680018841,0.33268781,"mRNA processing"),
c("GO:0030091","protein repair",0.06087415263465443,1,0.9095752142783583,0.27916289,"mRNA processing"),
c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,2,0.8847586056951727,0.25012422,"mRNA processing"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,1,0.8519189150891596,0.43903701,"mRNA processing"),
c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.8977734452569671,0.21843208,"mRNA processing"),
c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.903313014554542,0.40604291,"mRNA processing"),
c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8460257228812044,0.49920408,"mRNA processing"),
c("GO:0140040","mitochondrial polycistronic RNA processing",3.4506362332381646E-05,1,0.9003799201379421,0.41070687,"mRNA processing"),
c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.804790687542519,0.46167258,"mRNA processing"),
c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,2,0.863956424377826,0.33620138,"mRNA processing"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,1,0.8640716719978603,0.41679514,"mRNA processing"),
c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9519203852425954,0.10630662,"mRNA processing"),
c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9269586764061225,0.17943885,"mRNA processing"),
c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,2,0.9906419939902001,-0,"intracellular iron ion homeostasis"),
c("GO:0006915","apoptotic process",0.3972840732335429,1,0.9865014295873916,0.01170435,"apoptotic process"),
c("GO:0007049","cell cycle",2.8073119376140743,2,0.988382107485483,0.01495111,"cell cycle"),
c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.982669832660182,0.01255057,"chromosome segregation"),
c("GO:0000919","cell plate assembly",0.001451731958126628,1,0.9222325983880028,0.47630759,"chromosome segregation"),
c("GO:0007623","circadian rhythm",0.10049731555289494,2,0.9964178461450757,-0,"circadian rhythm"),
c("GO:0008150","biological_process",100,4,1,-0,"biological_process"),
c("GO:0008152","metabolic process",57.597931274565454,6,1,-0,"metabolic process"),
c("GO:0009058","biosynthetic process",29.004480601854056,4,0.9505650520180743,0.09593543,"biosynthetic process"),
c("GO:0044238","primary metabolic process",45.47569830278978,2,0.9467146471733705,0.23095606,"biosynthetic process"),
c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9186099655298551,0.16846091,"biosynthetic process"),
c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9360747008128004,0.17227516,"biosynthetic process"),
c("GO:0009404","toxin metabolic process",0.08737257416575693,1,0.9538084716959496,0.07301458,"toxin metabolic process"),
c("GO:0009651","response to salt stress",0.07343446852364134,7,0.9113692898218173,-0,"response to salt stress"),
c("GO:0006283","transcription-coupled nucleotide-excision repair",0.06059070751549558,1,0.8230631407135823,0.33995763,"response to salt stress"),
c("GO:0006979","response to oxidative stress",0.8191810417707402,4,0.9155380914298296,0.41360326,"response to salt stress"),
c("GO:0009611","response to wounding",0.16569462243976346,1,0.9253478567570326,0.44363439,"response to salt stress"),
c("GO:0009737","response to abscisic acid",0.05830589338105859,4,0.9014287750437808,0.1955272,"response to salt stress"),
c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.9011652813787087,0.17275947,"response to salt stress"),
c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9385510663190229,0.27664386,"response to salt stress"),
c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8661801521465045,0.32283054,"response to salt stress"),
c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8796247927802366,0.30241174,"response to salt stress"),
c("GO:0034059","response to anoxia",0.0003401341429906191,1,0.9333167070151557,0.45207873,"response to salt stress"),
c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8431627364187204,0.46171648,"response to salt stress"),
c("GO:0042221","response to chemical",4.8560360057130705,1,0.9234929233589982,0.27995591,"response to salt stress"),
c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9839300014980518,0.48813486,"response to salt stress"),
c("GO:0051607","defense response to virus",0.17287441054506547,3,0.9056042757659867,0.44531989,"response to salt stress"),
c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9437897153287094,0.12932521,"response to salt stress"),
c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,3,0.8906728355846558,0.46817928,"response to salt stress"),
c("GO:1901562","response to paraquat",0.00022182661499388205,1,0.9360756072227834,0.38531986,"response to salt stress"),
c("GO:1902074","response to salt",0.002343967898435353,3,0.9361097246591891,0.31787134,"response to salt stress"),
c("GO:0009853","photorespiration",0.014256057123606818,2,0.9594272139884765,0.06301944,"photorespiration"),
c("GO:0009908","flower development",0.02997124042584006,10,0.9032255220494619,-0,"flower development"),
c("GO:0001824","blastocyst development",0.009126932836914946,1,0.9344416855485165,0.41203625,"flower development"),
c("GO:0009555","pollen development",0.012900450031977538,5,0.9288625525967364,0.4203896,"flower development"),
c("GO:0009846","pollen germination",0.003164726373912717,2,0.9425335743316722,0.48835672,"flower development"),
c("GO:0009944","polarity specification of adaxial/abaxial axis",0.0024376280247661035,2,0.939214395325934,0.38300054,"flower development"),
c("GO:0010091","trichome branching",0.0006408324433156592,1,0.9407692954396792,0.31954118,"flower development"),
c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9232166220374267,0.29044182,"flower development"),
c("GO:0021987","cerebral cortex development",0.017519373104183483,1,0.9321073431898245,0.44335853,"flower development"),
c("GO:0042335","cuticle development",0.004384772756379067,1,0.9383818769417046,0.39539149,"flower development"),
c("GO:0061137","bud dilation",9.612486649734888E-05,1,0.9441900191532959,0.29194579,"flower development"),
c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9492885830919227,0.28605225,"flower development"),
c("GO:1905177","tracheary element differentiation",0.0002661919379926584,1,0.9484684178267991,0.2910313,"flower development"),
c("GO:0015031","protein transport",3.093438694074183,21,0.8786053956854146,0,"protein transport"),
c("GO:0000055","ribosomal large subunit export from nucleus",0.042115015226671805,1,0.8427177124621217,0.47380685,"protein transport"),
c("GO:0006811","monoatomic ion transport",4.7761710300947735,9,0.9267581658958487,0.42145531,"protein transport"),
c("GO:0006833","water transport",0.07339010320064257,1,0.9491204356525017,0.25687781,"protein transport"),
c("GO:0006897","endocytosis",0.6399968963991822,3,0.9206566561305087,0.32211657,"protein transport"),
c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9500934210874007,0.21356797,"protein transport"),
c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9420683543832767,0.25083478,"protein transport"),
c("GO:0016192","vesicle-mediated transport",2.6087919056355493,10,0.9218652825945151,0.38566904,"protein transport"),
c("GO:0045037","protein import into chloroplast stroma",0.006780500198312994,1,0.9169274294077039,0.49429905,"protein transport"),
c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9266755636689286,0.49534489,"protein transport"),
c("GO:1902600","proton transmembrane transport",1.312853708699458,3,0.9090932213665075,0.35175462,"protein transport"),
c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1,0.9239553918946453,0.36406786,"protein transport"),
c("GO:0016042","lipid catabolic process",1.4093137798594644,4,0.8711225736488413,0.0993369,"lipid catabolic process"),
c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,4,0.9350954547555893,0.15163831,"lipid catabolic process"),
c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.9127011287290062,0.49989178,"lipid catabolic process"),
c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9350535605849191,0.15178462,"lipid catabolic process"),
c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.9335782313079345,0.12277613,"lipid catabolic process"),
c("GO:0006914","autophagy",0.44134623319182764,2,0.9009225509892886,0.43293744,"lipid catabolic process"),
c("GO:0009450","gamma-aminobutyric acid catabolic process",0.04595015092589936,1,0.863787864849355,0.44469555,"lipid catabolic process"),
c("GO:0033611","oxalate catabolic process",0.006792823899145988,1,0.9006615297708536,0.3986442,"lipid catabolic process"),
c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,1,0.9063676543193916,0.45718665,"lipid catabolic process"),
c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1,0.9241632086625796,0.27800701,"lipid catabolic process"),
c("GO:0071040","nuclear polyadenylation-dependent antisense transcript catabolic process",0.008039982423444923,1,0.846140933066354,0.40606321,"lipid catabolic process"),
c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8712411897571667,0.33283067,"lipid catabolic process"),
c("GO:0022900","electron transport chain",0.8207017864535318,3,0.9275162410801076,0.09080927,"electron transport chain"),
c("GO:0030048","actin filament-based movement",0.0797343443894676,1,0.9809486708225112,0.00993277,"actin filament-based movement"),
c("GO:0032259","methylation",2.6278542060840238,8,0.964081561638524,0.05828286,"methylation"),
c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,8,0.8874822614448646,-0,"positive regulation of DNA-templated transcription"),
c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.00031548674132463224,1,0.9431985304868815,0.44967237,"positive regulation of DNA-templated transcription"),
c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.937240340085944,0.27139997,"positive regulation of DNA-templated transcription"),
c("GO:0006417","regulation of translation",1.3923070727099334,5,0.9042985549496313,0.40629436,"positive regulation of DNA-templated transcription"),
c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.8694391801668475,0.48392513,"positive regulation of DNA-templated transcription"),
c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9269597144511158,0.38387397,"positive regulation of DNA-templated transcription"),
c("GO:0009852","auxin catabolic process",3.4506362332381646E-05,1,0.9134896317867919,0.37814571,"positive regulation of DNA-templated transcription"),
c("GO:0009926","auxin polar transport",0.004261535748049133,1,0.9090351724516167,0.4567727,"positive regulation of DNA-templated transcription"),
c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9449242648127614,0.40013642,"positive regulation of DNA-templated transcription"),
c("GO:0010119","regulation of stomatal movement",0.016952482865865783,1,0.9527225817305781,0.16052978,"positive regulation of DNA-templated transcription"),
c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9449787012194721,0.2722228,"positive regulation of DNA-templated transcription"),
c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9146725916713815,0.41251019,"positive regulation of DNA-templated transcription"),
c("GO:0032784","regulation of DNA-templated transcription elongation",0.17948484367188314,1,0.9268571465544592,0.37447493,"positive regulation of DNA-templated transcription"),
c("GO:0032876","negative regulation of DNA endoreduplication",0.0004264000488215732,1,0.9357402602755096,0.48325207,"positive regulation of DNA-templated transcription"),
c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9407954170213159,0.21405542,"positive regulation of DNA-templated transcription"),
c("GO:0042659","regulation of cell fate specification",0.0026569898995933866,1,0.9395438038696537,0.4268527,"positive regulation of DNA-templated transcription"),
c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,1,0.9498742269467771,0.17055333,"positive regulation of DNA-templated transcription"),
c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9388218162506468,0.22344122,"positive regulation of DNA-templated transcription"),
c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9369179017178604,0.47414277,"positive regulation of DNA-templated transcription"),
c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,6,0.8935182768487635,0.47074713,"positive regulation of DNA-templated transcription"),
c("GO:0045995","regulation of embryonic development",0.01473914619626016,3,0.9307484475450518,0.14815537,"positive regulation of DNA-templated transcription"),
c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,2,0.939597542154595,0.42385341,"positive regulation of DNA-templated transcription"),
c("GO:0048586","regulation of long-day photoperiodism, flowering",0.0036552096670658557,3,0.9331763297600029,0.13376852,"positive regulation of DNA-templated transcription"),
c("GO:0050821","protein stabilization",0.10143638155636905,1,0.941029090531027,0.36692594,"positive regulation of DNA-templated transcription"),
c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9368560642337764,0.23293868,"positive regulation of DNA-templated transcription"),
c("GO:0060195","negative regulation of antisense RNA transcription",0.00012077226816333578,1,0.9407265990034013,0.2205473,"positive regulation of DNA-templated transcription"),
c("GO:0061635","regulation of protein complex stability",0.0017622892191180627,1,0.9537144360669078,0.1205724,"positive regulation of DNA-templated transcription"),
c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9361167435317312,0.33352541,"positive regulation of DNA-templated transcription"),
c("GO:0090333","regulation of stomatal closure",0.012828972567146176,2,0.9535274617843031,0.15711534,"positive regulation of DNA-templated transcription"),
c("GO:1900057","positive regulation of leaf senescence",0.0003204162216578296,1,0.9325822322013414,0.38480392,"positive regulation of DNA-templated transcription"),
c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,1,0.9280254393514106,0.47223923,"positive regulation of DNA-templated transcription"),
c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9339593553083878,0.31791635,"positive regulation of DNA-templated transcription"),
c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9413740422442916,0.42373749,"positive regulation of DNA-templated transcription"),
c("GO:1901000","regulation of response to salt stress",0.0019915100546117406,2,0.9394272056531792,0.47172255,"positive regulation of DNA-templated transcription"),
c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9287448698467682,0.454995,"positive regulation of DNA-templated transcription"),
c("GO:1903329","regulation of iron-sulfur cluster assembly",0.0002464740166598689,1,0.9566642416934493,0.12070797,"positive regulation of DNA-templated transcription"),
c("GO:2000024","regulation of leaf development",0.0023045320557697744,2,0.9417079737045797,0.42373749,"positive regulation of DNA-templated transcription"),
c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,1,0.9396185448834918,0.35770111,"positive regulation of DNA-templated transcription"),
c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,2,0.940687055674097,0.33609152,"positive regulation of DNA-templated transcription"),
c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,1,0.9409886171697722,0.2795365,"positive regulation of DNA-templated transcription"),
c("GO:0048511","rhythmic process",0.1547413171393989,2,1,-0,"rhythmic process"),
c("GO:0050896","response to stimulus",17.567785530535815,3,1,-0,"response to stimulus"),
c("GO:0051179","localization",19.75810399172557,1,1,-0,"localization"),
c("GO:0051301","cell division",1.5693197819947182,2,0.9890436581646581,0.01381155,"cell division"));

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

