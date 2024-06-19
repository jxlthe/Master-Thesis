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
revigo.data <- rbind(c("GO:0000003","reproduction",1.1154428097959028,1,1,-0,"reproduction"),
                     c("GO:0002376","immune system process",0.9427113541805001,2,1,-0,"immune system process"),
                     c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,3,0.9686358216302681,-0,"intracellular iron ion homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.981020579886502,0.46430576,"intracellular iron ion homeostasis"),
                     c("GO:0006952","defense response",1.1604144588756624,44,0.9134724519311507,-0,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.879827777764831,0.44951817,"defense response"),
                     c("GO:0009414","response to water deprivation",0.019128848432972426,16,0.8785314214004789,0.38192823,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,8,0.9254516172618769,0.45786068,"defense response"),
                     c("GO:0009625","response to insect",0.000539778096485113,2,0.9378941685822639,0.49043501,"defense response"),
                     c("GO:0009646","response to absence of light",0.0015108857221249963,3,0.9253577229366186,0.4554337,"defense response"),
                     c("GO:0009873","ethylene-activated signaling pathway",0.038247837905278456,11,0.8531535145513461,0.23262547,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8982731987736398,0.25765651,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,2,0.933346065833084,0.4902216,"defense response"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9350745269183776,0.4570295,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9329884633724373,0.29372182,"defense response"),
                     c("GO:0010447","response to acidic pH",0.005513623752681268,1,0.9263956476869073,0.48753241,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8634333996745125,0.38492887,"defense response"),
                     c("GO:0019236","response to pheromone",0.03625386311050012,1,0.9148117746363615,0.4993889,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,3,0.8757630073551622,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,4,0.9108509718162283,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,6,0.8413285625884777,0.40205681,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9203174924393598,0.46171648,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,13,0.9025456246148672,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,2,0.9820782649233186,0.48259105,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9298786822575796,0.49002119,"defense response"),
                     c("GO:0071000","response to magnetism",7.887168533115806E-05,1,0.9402222825420283,0.39597089,"defense response"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9390096720600496,0.14771546,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9284492141524632,0.44641368,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9248163774630318,0.39541577,"defense response"),
                     c("GO:1901562","response to paraquat",0.00022182661499388205,1,0.9351793434848085,0.37754219,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.934434494507892,0.31063873,"defense response"),
                     c("GO:0007017","microtubule-based process",1.4522840599239462,1,0.990143487369369,0.01367267,"microtubule-based process"),
                     c("GO:0007049","cell cycle",2.8073119376140743,9,0.989492087424824,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,3,0.9757722735243382,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,9,0.9959005435852354,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,12,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,19,1,-0,"metabolic process"),
                     c("GO:0008219","cell death",0.4651211168388386,1,0.9910898639992425,0.01191292,"cell death"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,2,0.9920904124444596,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9901838028127955,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,8,0.9563835628896364,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,7,0.9540238308470911,0.21768792,"biosynthetic process"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,3,0.9530524036319646,0.27529491,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9280437098504655,0.16846091,"biosynthetic process"),
                     c("GO:0071704","organic substance metabolic process",50.871943734517124,1,0.9521326283355345,0.31960258,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9431382619255397,0.17227516,"biosynthetic process"),
                     c("GO:0009695","jasmonic acid biosynthetic process",0.005493905831348478,4,0.8820553239475333,0.06061134,"jasmonic acid biosynthetic process"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,4,0.9107687165294583,0.26380719,"jasmonic acid biosynthetic process"),
                     c("GO:0006694","steroid biosynthetic process",0.2949431320360321,3,0.8643322450747801,0.40616444,"jasmonic acid biosynthetic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,1,0.8889571051236331,0.46657264,"jasmonic acid biosynthetic process"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9322091152525966,0.38978709,"jasmonic acid biosynthetic process"),
                     c("GO:0019432","triglyceride biosynthetic process",0.09997725537774264,2,0.8930989767583173,0.49870367,"jasmonic acid biosynthetic process"),
                     c("GO:0019760","glucosinolate metabolic process",0.005523482713347663,2,0.8984957731210338,0.22547604,"jasmonic acid biosynthetic process"),
                     c("GO:0031407","oxylipin metabolic process",0.012508556345488349,1,0.9241968951629272,0.27726473,"jasmonic acid biosynthetic process"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,2,0.9032366017649239,0.31813021,"jasmonic acid biosynthetic process"),
                     c("GO:0033387","putrescine biosynthetic process from ornithine",0.01868765994315126,1,0.8776589989161953,0.2843421,"jasmonic acid biosynthetic process"),
                     c("GO:0033491","coniferin metabolic process",1.7253181166190823E-05,1,0.9325722453908004,0.3924192,"jasmonic acid biosynthetic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,2,0.9100973444117998,0.3067044,"jasmonic acid biosynthetic process"),
                     c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1,0.9270117887483821,0.26504175,"jasmonic acid biosynthetic process"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9080846791615016,0.48869939,"jasmonic acid biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,2,0.962069788115057,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,19,0.8822204902848627,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,7,0.9155596328357238,0.4203896,"flower development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,4,0.9119539320899355,0.38141807,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,2,0.9356805668650227,0.48720779,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.8583572508885404,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9190179371314678,0.29044182,"flower development"),
                     c("GO:0019827","stem cell population maintenance",0.017275363827690213,2,0.9245142188070515,0.36342025,"flower development"),
                     c("GO:0021987","cerebral cortex development",0.017519373104183483,1,0.9190617509063173,0.44335853,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.9261666268625263,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.920758101491622,0.42767441,"flower development"),
                     c("GO:0055047","generative cell mitosis",7.147746483136198E-05,1,0.8937433537239626,0.40521407,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.9379660468836518,0.28605225,"flower development"),
                     c("GO:0090708","specification of plant organ axis polarity",0.002785156388256519,2,0.9048766078118994,0.42388435,"flower development"),
                     c("GO:1905177","tracheary element differentiation",0.0002661919379926584,1,0.9383884191347898,0.2910313,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,2,0.9295369247237566,0.4298524,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.9938703007851393,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9593105824517872,0.04613723,"chlororespiration"),
                     c("GO:0015031","protein transport",3.093438694074183,39,0.8856881222322535,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,28,0.9299296200807329,0.42145531,"protein transport"),
                     c("GO:0006833","water transport",0.07339010320064257,2,0.9512306966604464,0.25687781,"protein transport"),
                     c("GO:0006897","endocytosis",0.6399968963991822,5,0.9251490553342954,0.32211657,"protein transport"),
                     c("GO:0009852","auxin catabolic process",3.4506362332381646E-05,1,0.918609329134268,0.44733804,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,2,0.9523856179247027,0.21356797,"protein transport"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9433319652038283,0.48373949,"protein transport"),
                     c("GO:0010966","regulation of phosphate transport",0.0010204024289718573,1,0.9540785580000554,0.48139262,"protein transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.951568450969917,0.18579698,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9498793467403328,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,14,0.9257846983818245,0.38566904,"protein transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,8,0.9031255984668667,0.44154323,"protein transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9561906223935729,0.30124551,"protein transport"),
                     c("GO:0045037","protein import into chloroplast stroma",0.006780500198312994,1,0.9230660938123257,0.49429905,"protein transport"),
                     c("GO:0050821","protein stabilization",0.10143638155636905,1,0.9390658629772833,0.43402974,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9224341182431762,0.44310395,"protein transport"),
                     c("GO:0051646","mitochondrion localization",0.04394385243028803,1,0.9470244054068786,0.2322362,"protein transport"),
                     c("GO:0060918","auxin transport",0.015500750907739157,3,0.9005987910456148,0.22426661,"protein transport"),
                     c("GO:0061635","regulation of protein complex stability",0.0017622892191180627,1,0.9519447151888122,0.33196934,"protein transport"),
                     c("GO:0072699","protein localization to cortical microtubule cytoskeleton",0.00015774337066231612,1,0.9467319327152923,0.36351639,"protein transport"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,2,0.9278062822324343,0.3859558,"protein transport"),
                     c("GO:0120010","intermembrane phospholipid transfer",0.007554428610624982,1,0.8902307770751422,0.42470759,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.9314889060544345,0.49534489,"protein transport"),
                     c("GO:1904823","purine nucleobase transmembrane transport",0.0731485586643159,1,0.9290749747825838,0.49367033,"protein transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1,0.9297548608646681,0.30569128,"protein transport"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9337466975912142,0.46789196,"protein transport"),
                     c("GO:2000694","regulation of phragmoplast microtubule organization",0.0019791863537787476,1,0.9443902205137444,0.42590434,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,9,0.9550108711717666,0.04477571,"photosynthesis"),
                     c("GO:0016477","cell migration",0.49093927008395993,1,0.991049257897926,0.01198612,"cell migration"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,47,0.8875503224800294,0,"protein ubiquitination"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,3,0.889910754013842,0.4354711,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8756768221539741,0.45078305,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,20,0.9424227992072893,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9378176150799266,0.43492456,"protein ubiquitination"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,1,0.9207936575917155,0.34119038,"protein ubiquitination"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,2,0.8961847229904666,0.46608863,"protein ubiquitination"),
                     c("GO:0006076","(1->3)-beta-D-glucan catabolic process",2.4647401665986893E-06,1,0.9369468412709423,0.49567686,"protein ubiquitination"),
                     c("GO:0006182","cGMP biosynthetic process",0.05174475505757287,1,0.8495420820947562,0.45290106,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,3,0.8740544592960273,0.4022316,"protein ubiquitination"),
                     c("GO:0006283","transcription-coupled nucleotide-excision repair",0.06059070751549558,1,0.8366645313379071,0.46757848,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,24,0.8417882101772437,0.19691418,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,18,0.8296522788341744,0.48035482,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,20,0.8674236250172993,0.40866092,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,10,0.8693454815537701,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,4,0.9423874205568638,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,32,0.9411424466613101,0.12094833,"protein ubiquitination"),
                     c("GO:0006747","FAD biosynthetic process",0.030434611577160615,1,0.8585945486528581,0.43629805,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,3,0.9201546126255483,0.45718665,"protein ubiquitination"),
                     c("GO:0008652","amino acid biosynthetic process",2.679426429329935,12,0.8381579891861963,0.49906638,"protein ubiquitination"),
                     c("GO:0009820","alkaloid metabolic process",0.01601588160255828,2,0.9485491751247405,0.14390223,"protein ubiquitination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8965124392057823,0.367811,"protein ubiquitination"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9437965846692097,0.32178408,"protein ubiquitination"),
                     c("GO:0015966","diadenosine tetraphosphate biosynthetic process",0.028509649507047038,1,0.8647900867621974,0.36610043,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,4,0.8811583262165019,0.20835326,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,77,0.9098459003605229,0.44883397,"protein ubiquitination"),
                     c("GO:0016554","cytidine to uridine editing",0.012560315888986921,1,0.8967884109527346,0.3977518,"protein ubiquitination"),
                     c("GO:0018130","heterocycle biosynthetic process",8.277253100322714,2,0.8745940089841773,0.46779521,"protein ubiquitination"),
                     c("GO:0018345","protein palmitoylation",0.0003918936864891916,1,0.9136811946983828,0.37471442,"protein ubiquitination"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9399102664551204,0.33268781,"protein ubiquitination"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9241297699418284,0.27916289,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,8,0.8875206936688738,0.25012422,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,2,0.8710247194655414,0.28833637,"protein ubiquitination"),
                     c("GO:0033358","UDP-L-arabinose biosynthetic process",2.21826614993882E-05,1,0.9080094382781345,0.46513452,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9496007741312238,0.45767193,"protein ubiquitination"),
                     c("GO:0042245","RNA repair",0.022010129687726296,1,0.9033176951860613,0.35967549,"protein ubiquitination"),
                     c("GO:0042726","flavin-containing compound metabolic process",0.22289631222618586,1,0.9078809166818013,0.23046044,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,2,0.9061961605967513,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,3,0.9234150366642541,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.888352038284189,0.31626431,"protein ubiquitination"),
                     c("GO:0046166","glyceraldehyde-3-phosphate biosynthetic process",0.03550458209985412,1,0.9017590795849828,0.33432601,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,9,0.9142523755930573,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,2,0.9025197296743325,0.41108531,"protein ubiquitination"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.8229889350754062,0.42789517,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.8736344757012352,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.895328418162861,0.45115494,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,3,0.9364364843856474,0.45105203,"protein ubiquitination"),
                     c("GO:0140040","mitochondrial polycistronic RNA processing",3.4506362332381646E-05,1,0.9158728863742371,0.41070687,"protein ubiquitination"),
                     c("GO:1900871","chloroplast mRNA modification",0.00043625900948796796,1,0.9104502567020882,0.47298773,"protein ubiquitination"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.8391362254428162,0.46167258,"protein ubiquitination"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",9.084396351119786,2,0.8812813452757842,0.47967538,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,4,0.8947764310896186,0.31086297,"protein ubiquitination"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8768921732476757,0.49989178,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9920083214141256,0.01031083,"establishment of cell polarity"),
                     c("GO:0032259","methylation",2.6278542060840238,15,0.9678323230608921,0.05843136,"methylation"),
                     c("GO:0032502","developmental process",4.156005433076044,1,1,-0,"developmental process"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9772237518509119,0.03897815,"collagen metabolic process"),
                     c("GO:0040011","locomotion",0.5148398554794674,1,1,-0,"locomotion"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,16,0.8915756535066078,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,5,0.9097910614586415,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009786","regulation of asymmetric cell division",0.004426673339211246,1,0.955108677528184,0.14531513,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,5,0.924008625745951,0.38244242,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,2,0.9428171668433051,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010075","regulation of meristem growth",0.0026200187970944065,2,0.8605703832674918,0.39444345,"positive regulation of DNA-templated transcription"),
                     c("GO:0010089","xylem development",0.003164726373912717,1,0.9289878202466676,0.47349401,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,3,0.9515244205675838,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9498995339125141,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010366","negative regulation of ethylene biosynthetic process",0.0002045734338276912,1,0.9350142978377944,0.45396999,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9178782014905182,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0010928","regulation of auxin mediated signaling pathway",0.004086539196220627,1,0.9331502998788527,0.40226117,"positive regulation of DNA-templated transcription"),
                     c("GO:0032784","regulation of DNA-templated transcription elongation",0.17948484367188314,2,0.9295724228994686,0.37447493,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9399437687794638,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042659","regulation of cell fate specification",0.0026569898995933866,2,0.9372355320335743,0.39470747,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,4,0.948889587906444,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9380556035294847,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9355785026747696,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,10,0.8940374526132623,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,3,0.9293109081356303,0.42989004,"positive regulation of DNA-templated transcription"),
                     c("GO:0048509","regulation of meristem development",0.004372449055546074,1,0.9374968356560581,0.4381267,"positive regulation of DNA-templated transcription"),
                     c("GO:0048510","regulation of timing of transition from vegetative to reproductive phase",0.002316855756602768,2,0.9385062593712483,0.39214155,"positive regulation of DNA-templated transcription"),
                     c("GO:0050777","negative regulation of immune response",0.030624396569988714,2,0.9105462894204667,0.45733164,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,3,0.9437530455211556,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,5,0.9361829156562242,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0060147","regulation of post-transcriptional gene silencing",0.021610841780737307,1,0.9400609881569504,0.26975857,"positive regulation of DNA-templated transcription"),
                     c("GO:0060195","negative regulation of antisense RNA transcription",0.00012077226816333578,1,0.9398665110039087,0.26135971,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9324036814431228,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,6,0.9497056054756448,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:0141005","retrotransposon silencing by heterochromatin formation",0.001355607091629279,1,0.9428929429040672,0.29342052,"positive regulation of DNA-templated transcription"),
                     c("GO:1900057","positive regulation of leaf senescence",0.0003204162216578296,1,0.9301175786325859,0.35419608,"positive regulation of DNA-templated transcription"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,5,0.9255845936231454,0.138781,"positive regulation of DNA-templated transcription"),
                     c("GO:1900457","regulation of brassinosteroid mediated signaling pathway",0.001054908791304239,1,0.9375717953283268,0.34087448,"positive regulation of DNA-templated transcription"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,3,0.9277197793630935,0.32535433,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9387344175913066,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1902326","positive regulation of chlorophyll biosynthetic process",0.00019717921332789513,1,0.9378805987248041,0.43671847,"positive regulation of DNA-templated transcription"),
                     c("GO:1902448","positive regulation of shade avoidance",4.190058283217772E-05,1,0.9348739791521216,0.30806662,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.9232814504088713,0.48792077,"positive regulation of DNA-templated transcription"),
                     c("GO:1903329","regulation of iron-sulfur cluster assembly",0.0002464740166598689,1,0.949749293047066,0.40066627,"positive regulation of DNA-templated transcription"),
                     c("GO:1903527","positive regulation of membrane tubulation",0.0019298915504467734,1,0.9308269458228253,0.46246533,"positive regulation of DNA-templated transcription"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9341913688746851,0.40819627,"positive regulation of DNA-templated transcription"),
                     c("GO:1904966","positive regulation of vitamin E biosynthetic process",1.9717921332789515E-05,1,0.9401339564910491,0.38536455,"positive regulation of DNA-templated transcription"),
                     c("GO:1905038","regulation of membrane lipid metabolic process",0.009506502822571143,1,0.9404618961423398,0.25241702,"positive regulation of DNA-templated transcription"),
                     c("GO:1905157","positive regulation of photosynthesis",0.00010598382716374363,1,0.9382733536899409,0.39735289,"positive regulation of DNA-templated transcription"),
                     c("GO:1905183","negative regulation of protein serine/threonine phosphatase activity",1.7253181166190823E-05,1,0.9393315556699844,0.24023597,"positive regulation of DNA-templated transcription"),
                     c("GO:2000022","regulation of jasmonic acid mediated signaling pathway",0.011707515791343773,1,0.9292449834339215,0.42518681,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,5,0.9395119566225884,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,4,0.9278259249908258,0.37396531,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,3,0.9351126910533618,0.36714483,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,2,0.9306424609003232,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000049","positive regulation of cell-cell adhesion mediated by cadherin",0.0027210731439249528,1,0.9389922460317134,0.41411697,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,3,0.9369698607582899,0.39204232,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,2,0.9337478031918454,0.49708545,"positive regulation of DNA-templated transcription"),
                     c("GO:2000134","negative regulation of G1/S transition of mitotic cell cycle",0.02027248787027422,1,0.9295991590062869,0.16280015,"positive regulation of DNA-templated transcription"),
                     c("GO:2000185","regulation of phosphate transmembrane transport",0.00017253181166190822,1,0.9556565290595118,0.11823538,"positive regulation of DNA-templated transcription"),
                     c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,1,0.9425165946111713,0.2795365,"positive regulation of DNA-templated transcription"),
                     c("GO:2000280","regulation of root development",0.008382581306602141,1,0.93530346561912,0.45378622,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,2,0.9400621691337849,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.9320504846101171,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,2,0.9564176021664041,0.12578285,"positive regulation of DNA-templated transcription"),
                     c("GO:2000762","regulation of phenylpropanoid metabolic process",0.005555524335513445,1,0.9435717501109856,0.23331889,"positive regulation of DNA-templated transcription"),
                     c("GO:2001006","regulation of cellulose biosynthetic process",0.002356291599268347,2,0.9450457632488612,0.22862385,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,3,1,-0,"rhythmic process"),
                     c("GO:0048544","recognition of pollen",0.04932191547380636,1,0.9924948346594075,0.00950258,"recognition of pollen"),
                     c("GO:0050896","response to stimulus",17.567785530535815,12,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,3,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,6,0.9900713261441032,0.01381155,"cell division"),
                     c("GO:0060361","flight",9.612486649734888E-05,1,0.9992041163488156,-0,"flight"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,24,0.9115489367113405,0.01286368,"cell wall organization"),
                     c("GO:0000741","karyogamy",0.004604134631206351,1,0.9409217300776628,0.24960876,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,18,0.909556730491422,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,2,0.9283815069843394,0.39503669,"cell wall organization"),
                     c("GO:0007005","mitochondrion organization",0.8276597479438399,1,0.9175331228084536,0.4723023,"cell wall organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,3,0.9205428040419101,0.38546222,"cell wall organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,2,0.9238270729386129,0.42383478,"cell wall organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9263348218637661,0.40792275,"cell wall organization"),
                     c("GO:0009658","chloroplast organization",0.0906284959258338,7,0.92573827825937,0.31305597,"cell wall organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9565340798418022,0.2057032,"cell wall organization"),
                     c("GO:0010275","NAD(P)H dehydrogenase complex assembly",0.0022404488114382086,1,0.9342787929768167,0.41693818,"cell wall organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.929799954548471,0.45272646,"cell wall organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,5,0.9160784415903175,0.35075221,"cell wall organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9274077633240769,0.39384196,"cell wall organization"),
                     c("GO:0034462","small-subunit processome assembly",0.005688620304509775,1,0.9171810923498835,0.44451454,"cell wall organization"),
                     c("GO:0043622","cortical microtubule organization",0.008850881938255893,3,0.9340381291657605,0.26122046,"cell wall organization"),
                     c("GO:0048564","photosystem I assembly",0.005678761343843379,1,0.8947685122755664,0.44445982,"cell wall organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9258000800482549,0.48565697,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9293688456269509,0.4524177,"cell wall organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,2,0.9505782803397702,0.20466634,"cell wall organization"),
                     c("GO:1905691","lipid droplet disassembly",0.0006679445851482447,1,0.9461820487839615,0.22063843,"cell wall organization"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.957105431386808,0.08070223,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9352773509874961,0.08958079,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches
#pdf( file="test_BP.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_BP_Table.tsv", sep = "\t")
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

# 
# 
# topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
# topGO_data <- topGO_data[topGO_data$Ontology == "BP",]
# 
# terms <- intersect(t$tm$description, topGO_data$Term)
# 
# 
# 
# p <- topGO_data[topGO_data$Term %in% terms,]$weight01
# l <- log10(p)
# 
# 
# colfunc<-colorRampPalette(c("blue", "yellow"))
# 
# colors <- colfunc(100)
# col_index <- round(l / log10(0.01) * 100)
# col_plot <- colors[col_index]
# 
# 
# with(
#   # t$tm is a data frame with one row per rectangle.  Filter to the group we
#   # want to highlight.
#   t$tm %>%
#     filter(description %in% terms),
#   {
#     # Use grid.rect to add a rectangle on top of the treemap.
#     grid.rect(x = x0 + (w / 2),
#               y = y0 + (h / 2),
#               width = w,
#               height = h,
#               gp = gpar(col = col_plot, fill = NA, lwd = 3),
#               vp = "data")
#     
#   }
# )




dev.off()

