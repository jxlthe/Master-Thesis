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
revigo.data <- rbind(c("GO:0006952","defense response",1.1604144588756624,21,0.9014415233402141,-0,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8752042541957586,0.44951817,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,2,0.9150452051149945,0.45786068,"defense response"),
                     c("GO:0009625","response to insect",0.000539778096485113,1,0.9291554169291139,0.49043501,"defense response"),
                     c("GO:0009646","response to absence of light",0.0015108857221249963,2,0.9120105516585164,0.48888107,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,7,0.8891454207552074,0.4259419,"defense response"),
                     c("GO:0009734","auxin-activated signaling pathway",0.07526084098709097,7,0.8436645479161755,0.24614308,"defense response"),
                     c("GO:0009751","response to salicylic acid",0.010802956150202055,2,0.9066274414130375,0.48305816,"defense response"),
                     c("GO:0010019","chloroplast-nucleus signaling pathway",0.00844173507060051,1,0.8945114961461071,0.26818052,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,1,0.9216148131432202,0.4902216,"defense response"),
                     c("GO:0010272","response to silver ion",0.00015774337066231612,1,0.9286548257423337,0.41832444,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9225360330506323,0.29372182,"defense response"),
                     c("GO:0010343","singlet oxygen-mediated programmed cell death",0.0019717921332789512,3,0.8891696093474535,0.38975058,"defense response"),
                     c("GO:0018872","arsonoacetate metabolic process",0.003095713649247954,1,0.8541037605927543,0.39975389,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,2,0.8704782250522954,0.46097507,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,4,0.8995467814880478,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,4,0.8347303631233572,0.42828273,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.9093418202069201,0.46171648,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,5,0.8888638884058782,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9800865143688429,0.48259105,"defense response"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,1,0.9232817742053951,0.44100881,"defense response"),
                     c("GO:0070987","error-free translesion synthesis",0.020632339934597628,1,0.8100738023559206,0.384161,"defense response"),
                     c("GO:0071000","response to magnetism",7.887168533115806E-05,1,0.9303898364080667,0.31824892,"defense response"),
                     c("GO:0071217","cellular response to external biotic stimulus",3.204162216578296E-05,1,0.9304188450813394,0.14771546,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9219237996016767,0.40537115,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9133752749105614,0.26816179,"defense response"),
                     c("GO:1902074","response to salt",0.002343967898435353,3,0.9277672192478896,0.32241615,"defense response"),
                     c("GO:0007017","microtubule-based process",1.4522840599239462,1,0.9899491638152037,0.01367267,"microtubule-based process"),
                     c("GO:0007018","microtubule-based movement",0.6202173565622278,2,0.9888196900259673,0.01231347,"microtubule-based movement"),
                     c("GO:0007049","cell cycle",2.8073119376140743,4,0.9892705821896848,0.01495111,"cell cycle"),
                     c("GO:0007059","chromosome segregation",0.7290602823192258,2,0.9790446456918097,0.01255057,"chromosome segregation"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,6,0.9959533166598166,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,7,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,9,1,-0,"metabolic process"),
                     c("GO:0008283","cell population proliferation",0.10238777126067615,1,0.9919661296944029,0.01017253,"cell population proliferation"),
                     c("GO:0008356","asymmetric cell division",0.009826919044228973,1,0.9910094743951549,0.00829585,"asymmetric cell division"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,6,0.9528457169537842,0.09593543,"biosynthetic process"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,4,0.9502312409889643,0.21768792,"biosynthetic process"),
                     c("GO:0044238","primary metabolic process",45.47569830278978,3,0.9491604339430597,0.27529491,"biosynthetic process"),
                     c("GO:0046483","heterocycle metabolic process",21.449154825808435,1,0.9229285056035131,0.16846091,"biosynthetic process"),
                     c("GO:1901360","organic cyclic compound metabolic process",22.80942274109275,1,0.9387905680641739,0.17227516,"biosynthetic process"),
                     c("GO:0009853","photorespiration",0.014256057123606818,1,0.9595194773116903,0.05587782,"photorespiration"),
                     c("GO:0009908","flower development",0.02997124042584006,10,0.8825629182428915,-0,"flower development"),
                     c("GO:0009555","pollen development",0.012900450031977538,4,0.9126127560427858,0.4203896,"flower development"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9137796479654485,0.3897097,"flower development"),
                     c("GO:0009846","pollen germination",0.003164726373912717,1,0.9372901018299922,0.48835672,"flower development"),
                     c("GO:0009877","nodulation",0.0005175954349857247,2,0.9342080187589882,0.3161767,"flower development"),
                     c("GO:0010077","maintenance of inflorescence meristem identity",0.0001355607091629279,2,0.9279314546127306,0.29658706,"flower development"),
                     c("GO:0010143","cutin biosynthetic process",0.008845952457922695,1,0.8517435466529647,0.36762071,"flower development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9240734370776816,0.29044182,"flower development"),
                     c("GO:0042335","cuticle development",0.004384772756379067,1,0.9235019774120778,0.39539149,"flower development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.917819379933899,0.42767441,"flower development"),
                     c("GO:0048262","determination of dorsal/ventral asymmetry",0.0029848003417510126,2,0.9201281255410524,0.38718574,"flower development"),
                     c("GO:0090392","sepal giant cell differentiation",0.00018485551249490168,1,0.937038912692902,0.28605225,"flower development"),
                     c("GO:0090626","plant epidermis morphogenesis",0.016804598455869863,1,0.9153740668740785,0.46768707,"flower development"),
                     c("GO:1905393","plant organ formation",0.0036354917457330667,1,0.9278642103639095,0.48720779,"flower development"),
                     c("GO:0010118","stomatal movement",0.0019570036922793594,1,0.9937958592974828,0.00736082,"stomatal movement"),
                     c("GO:0010478","chlororespiration",0.0006087908211498762,1,0.9585695496235732,0.04613723,"chlororespiration"),
                     c("GO:0010960","magnesium ion homeostasis",0.024344238625495253,2,0.9731856934415024,-0,"magnesium ion homeostasis"),
                     c("GO:0010268","brassinosteroid homeostasis",0.011392029050019141,1,0.9758982367300477,0.46517239,"magnesium ion homeostasis"),
                     c("GO:0060586","multicellular organismal-level iron ion homeostasis",0.00353443739890252,1,0.9233250804375102,0.48479779,"magnesium ion homeostasis"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9798449387920019,0.37586433,"magnesium ion homeostasis"),
                     c("GO:0015031","protein transport",3.093438694074183,20,0.9020310464464982,-0,"protein transport"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,10,0.9408560440643576,0.42145531,"protein transport"),
                     c("GO:0008361","regulation of cell size",0.13159987171520382,1,0.889803447393566,0.42491424,"protein transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9587691104465781,0.21356797,"protein transport"),
                     c("GO:0010600","regulation of auxin biosynthetic process",0.0002045734338276912,1,0.9479669069238335,0.48373949,"protein transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9606534356467149,0.39622712,"protein transport"),
                     c("GO:0015795","sorbitol transmembrane transport",3.4506362332381646E-05,1,0.9524081710689201,0.25083478,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,8,0.9355613033604673,0.38566904,"protein transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9617817562793647,0.30124551,"protein transport"),
                     c("GO:0048250","iron import into the mitochondrion",0.008732574410259155,1,0.9460649327887085,0.4861705,"protein transport"),
                     c("GO:0060918","auxin transport",0.015500750907739157,3,0.9119859198185138,0.22426661,"protein transport"),
                     c("GO:0120010","intermembrane phospholipid transfer",0.007554428610624982,1,0.9010444626280496,0.42470759,"protein transport"),
                     c("GO:0140570","extraction of mislocalized protein from mitochondrial outer membrane",0.006970285191141092,1,0.942342927030281,0.49534489,"protein transport"),
                     c("GO:1902600","proton transmembrane transport",1.312853708699458,4,0.9248036519119001,0.35175462,"protein transport"),
                     c("GO:0015979","photosynthesis",0.228607115192195,7,0.9523107492733556,0.04477571,"photosynthesis"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,21,0.8867418503322062,0,"protein ubiquitination"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,2,0.8857069838859444,0.17872097,"protein ubiquitination"),
                     c("GO:0002926","tRNA wobble base 5-methoxycarbonylmethyl-2-thiouridinylation",0.014958508071087444,1,0.8719000511074585,0.40302693,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,8,0.938475339077507,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.9378623204198039,0.43492456,"protein ubiquitination"),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.10908200555315818,1,0.9162879551515956,0.41634329,"protein ubiquitination"),
                     c("GO:0006108","malate metabolic process",0.07909597668631853,3,0.9076858446580361,0.40337462,"protein ubiquitination"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",18.82728500884891,2,0.8598648721274008,0.31102513,"protein ubiquitination"),
                     c("GO:0006221","pyrimidine nucleotide biosynthetic process",0.5101395959817636,2,0.8177561276725704,0.39188311,"protein ubiquitination"),
                     c("GO:0006397","mRNA processing",1.242261085587905,9,0.8362351537008291,0.48035482,"protein ubiquitination"),
                     c("GO:0006412","translation",4.38869169324396,15,0.8179568353632322,0.4385246,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,14,0.8564446603318003,0.47678922,"protein ubiquitination"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,5,0.8552032043842711,0.37517093,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,3,0.9384366369838889,0.15178462,"protein ubiquitination"),
                     c("GO:0006541","glutamine metabolic process",0.5478057552077248,3,0.8657742027539862,0.48867994,"protein ubiquitination"),
                     c("GO:0006572","tyrosine catabolic process",0.043362173750970734,2,0.8529629947027482,0.46616537,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,19,0.9370739768470074,0.12862028,"protein ubiquitination"),
                     c("GO:0006631","fatty acid metabolic process",1.936945636803579,6,0.8375113762255482,0.1033813,"protein ubiquitination"),
                     c("GO:0006744","ubiquinone biosynthetic process",0.2844852395091539,2,0.8806614315335874,0.31493359,"protein ubiquitination"),
                     c("GO:0006914","autophagy",0.44134623319182764,2,0.9105803245525947,0.45718665,"protein ubiquitination"),
                     c("GO:0009695","jasmonic acid biosynthetic process",0.005493905831348478,1,0.8734302267163684,0.46296398,"protein ubiquitination"),
                     c("GO:0010028","xanthophyll cycle",0.0017894013609506482,1,0.9271848026810958,0.37903497,"protein ubiquitination"),
                     c("GO:0010586","miRNA metabolic process",0.02921949467502746,1,0.8926605170285802,0.36563143,"protein ubiquitination"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9433089891320287,0.32178408,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,4,0.8758501588520232,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,37,0.9022626348237903,0.46958224,"protein ubiquitination"),
                     c("GO:0019919","peptidyl-arginine methylation, to asymmetrical-dimethyl arginine",4.6830063165375095E-05,1,0.9400872737254753,0.33268781,"protein ubiquitination"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.9190334711922511,0.27916289,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,4,0.887881288036158,0.14129671,"protein ubiquitination"),
                     c("GO:0031407","oxylipin metabolic process",0.012508556345488349,1,0.92142645830669,0.34584249,"protein ubiquitination"),
                     c("GO:0031408","oxylipin biosynthetic process",0.012412431478991,2,0.8988273425569252,0.34563628,"protein ubiquitination"),
                     c("GO:0033320","UDP-D-xylose biosynthetic process",0.014470489518100906,2,0.8658033640727484,0.28833637,"protein ubiquitination"),
                     c("GO:0033358","UDP-L-arabinose biosynthetic process",2.21826614993882E-05,1,0.9044380961315857,0.46513452,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9461418414882614,0.45767193,"protein ubiquitination"),
                     c("GO:0042245","RNA repair",0.022010129687726296,1,0.9008423249166922,0.26040271,"protein ubiquitination"),
                     c("GO:0042726","flavin-containing compound metabolic process",0.22289631222618586,1,0.9038603215380578,0.18149234,"protein ubiquitination"),
                     c("GO:0042732","D-xylose metabolic process",0.08287442336171433,2,0.906545601370787,0.43491138,"protein ubiquitination"),
                     c("GO:0042744","hydrogen peroxide catabolic process",0.18507240962956237,2,0.915799458724189,0.39247003,"protein ubiquitination"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,1,0.8857711662541896,0.26614879,"protein ubiquitination"),
                     c("GO:0044272","sulfur compound biosynthetic process",1.4955624325092522,1,0.9013344167394665,0.16428511,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,6,0.9116976020514135,0.40604291,"protein ubiquitination"),
                     c("GO:0046835","carbohydrate phosphorylation",0.3487410156523817,2,0.8999673429304007,0.41108531,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9041969061294882,0.44155167,"protein ubiquitination"),
                     c("GO:0071051","polyadenylation-dependent snoRNA 3'-end processing",0.053797883616349594,1,0.86697030996335,0.49920408,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8960007461728071,0.43527038,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9380390776392794,0.45105203,"protein ubiquitination"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,1,0.8377417472723272,0.33478559,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,3,0.8855880125377597,0.32141303,"protein ubiquitination"),
                     c("GO:1902000","homogentisate catabolic process",0.014150073296443074,1,0.8702808049653483,0.49716341,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9918814181577328,0.01031083,"establishment of cell polarity"),
                     c("GO:0032259","methylation",2.6278542060840238,7,0.9655103582987943,0.05843136,"methylation"),
                     c("GO:0032963","collagen metabolic process",0.0483212309661673,1,0.9757454380382784,0.03897815,"collagen metabolic process"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,9,0.8991566082954726,-0,"positive regulation of DNA-templated transcription"),
                     c("GO:0006111","regulation of gluconeogenesis",0.02427522590083049,1,0.9428039972514977,0.27139997,"positive regulation of DNA-templated transcription"),
                     c("GO:0006417","regulation of translation",1.3923070727099334,4,0.9135546939424173,0.40629436,"positive regulation of DNA-templated transcription"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9311564383206755,0.38387397,"positive regulation of DNA-templated transcription"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,2,0.9474527886383333,0.36714542,"positive regulation of DNA-templated transcription"),
                     c("GO:0010119","regulation of stomatal movement",0.016952482865865783,3,0.9547254021333066,0.16052978,"positive regulation of DNA-templated transcription"),
                     c("GO:0010310","regulation of hydrogen peroxide metabolic process",0.0011042035946362127,1,0.9529935368317014,0.20908612,"positive regulation of DNA-templated transcription"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,1,0.9219415586561889,0.41251019,"positive regulation of DNA-templated transcription"),
                     c("GO:0010928","regulation of auxin mediated signaling pathway",0.004086539196220627,1,0.9378189029717181,0.40226117,"positive regulation of DNA-templated transcription"),
                     c("GO:0042127","regulation of cell population proliferation",0.41876428378545044,1,0.9436478537127534,0.21405542,"positive regulation of DNA-templated transcription"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.9326499349548192,0.31416992,"positive regulation of DNA-templated transcription"),
                     c("GO:0042548","regulation of photosynthesis, light reaction",0.0070861279789712316,1,0.9443667897618462,0.23746411,"positive regulation of DNA-templated transcription"),
                     c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,3,0.9521274585792716,0.17055333,"positive regulation of DNA-templated transcription"),
                     c("GO:0043067","regulation of programmed cell death",0.627229542336201,1,0.9418312780824202,0.22344122,"positive regulation of DNA-templated transcription"),
                     c("GO:0045019","negative regulation of nitric oxide biosynthetic process",0.000704915687647225,1,0.9412598393586022,0.47414277,"positive regulation of DNA-templated transcription"),
                     c("GO:0045892","negative regulation of DNA-templated transcription",1.5258343712354174,6,0.9022931592291114,0.47074713,"positive regulation of DNA-templated transcription"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,2,0.9345363630100864,0.42373749,"positive regulation of DNA-templated transcription"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.947303756171555,0.19569204,"positive regulation of DNA-templated transcription"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9400267042149186,0.23293868,"positive regulation of DNA-templated transcription"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9380809659054875,0.33352541,"positive regulation of DNA-templated transcription"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,4,0.953751044933043,0.15711534,"positive regulation of DNA-templated transcription"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9373416751386906,0.46789196,"positive regulation of DNA-templated transcription"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,3,0.9237448337854424,0.34078171,"positive regulation of DNA-templated transcription"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9343407088175688,0.39696279,"positive regulation of DNA-templated transcription"),
                     c("GO:1900618","regulation of shoot system morphogenesis",0.0023045320557697744,1,0.9435583442134363,0.43888643,"positive regulation of DNA-templated transcription"),
                     c("GO:1901529","positive regulation of anion channel activity",1.2323700832993446E-05,1,0.9516654506605777,0.45462725,"positive regulation of DNA-templated transcription"),
                     c("GO:1902183","regulation of shoot apical meristem development",0.0016217990296219374,1,0.9439574409215425,0.38065657,"positive regulation of DNA-templated transcription"),
                     c("GO:1902448","positive regulation of shade avoidance",4.190058283217772E-05,1,0.9403999884252786,0.30806662,"positive regulation of DNA-templated transcription"),
                     c("GO:1902584","positive regulation of response to water deprivation",0.0008848417198089295,1,0.929191830567516,0.45547566,"positive regulation of DNA-templated transcription"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9434116934817975,0.40819627,"positive regulation of DNA-templated transcription"),
                     c("GO:2000024","regulation of leaf development",0.0023045320557697744,3,0.9443940753676973,0.12960503,"positive regulation of DNA-templated transcription"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,3,0.9324739861219885,0.36922831,"positive regulation of DNA-templated transcription"),
                     c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,3,0.9394989877255865,0.13634361,"positive regulation of DNA-templated transcription"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,1,0.93548146543556,0.37952976,"positive regulation of DNA-templated transcription"),
                     c("GO:2000032","regulation of secondary shoot formation",0.007307954593965113,1,0.9380226816663005,0.49488849,"positive regulation of DNA-templated transcription"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9414798963672286,0.42989004,"positive regulation of DNA-templated transcription"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1,0.9379112523451213,0.46345207,"positive regulation of DNA-templated transcription"),
                     c("GO:2000185","regulation of phosphate transmembrane transport",0.00017253181166190822,1,0.959784708305648,0.44318481,"positive regulation of DNA-templated transcription"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9431411244988532,0.15519254,"positive regulation of DNA-templated transcription"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,2,0.9435362881050269,0.26299316,"positive regulation of DNA-templated transcription"),
                     c("GO:2000436","positive regulation of protein neddylation",0.0005570312776513037,1,0.9371495805560641,0.46528664,"positive regulation of DNA-templated transcription"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,1,0.9618554096420213,0.12578285,"positive regulation of DNA-templated transcription"),
                     c("GO:2000762","regulation of phenylpropanoid metabolic process",0.005555524335513445,1,0.9473593492243972,0.23331889,"positive regulation of DNA-templated transcription"),
                     c("GO:0048511","rhythmic process",0.1547413171393989,3,1,-0,"rhythmic process"),
                     c("GO:0050896","response to stimulus",17.567785530535815,7,1,-0,"response to stimulus"),
                     c("GO:0051179","localization",19.75810399172557,2,1,-0,"localization"),
                     c("GO:0051301","cell division",1.5693197819947182,2,0.9898740850883301,0.01381155,"cell division"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,10,0.9223546881051858,0.01286368,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,8,0.9160040513190786,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.9355201785624562,0.4007549,"cell wall organization"),
                     c("GO:0007005","mitochondrion organization",0.8276597479438399,1,0.9257385419234894,0.4723023,"cell wall organization"),
                     c("GO:0007030","Golgi organization",0.24565079344422494,2,0.9287212961507355,0.42383478,"cell wall organization"),
                     c("GO:0007033","vacuole organization",0.29445264874287896,2,0.9314090845663998,0.39090469,"cell wall organization"),
                     c("GO:0009657","plastid organization",0.1855628929227155,1,0.9336721204720532,0.41402298,"cell wall organization"),
                     c("GO:0009658","chloroplast organization",0.0906284959258338,3,0.9309591048054646,0.31305597,"cell wall organization"),
                     c("GO:0009663","plasmodesma organization",0.0001996439534944938,1,0.9604528369759061,0.2057032,"cell wall organization"),
                     c("GO:0009920","cell plate formation involved in plant-type cell wall biogenesis",0.00015527863049571742,1,0.9296977845087128,0.45904668,"cell wall organization"),
                     c("GO:0010207","photosystem II assembly",0.0070861279789712316,1,0.8968887601852811,0.45155549,"cell wall organization"),
                     c("GO:0010275","NAD(P)H dehydrogenase complex assembly",0.0022404488114382086,1,0.9392468742119799,0.41693818,"cell wall organization"),
                     c("GO:0010387","COP9 signalosome assembly",0.007344925696464094,1,0.9350898198040947,0.45272646,"cell wall organization"),
                     c("GO:0016226","iron-sulfur cluster assembly",0.3194697614338557,2,0.9225334499289222,0.35075221,"cell wall organization"),
                     c("GO:0034462","small-subunit processome assembly",0.005688620304509775,1,0.9233327415361822,0.44451454,"cell wall organization"),
                     c("GO:0070072","vacuolar proton-transporting V-type ATPase complex assembly",0.01876406688831582,1,0.9313780013526536,0.48565697,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9344573376323081,0.4524177,"cell wall organization"),
                     c("GO:0080119","ER body organization",0.00018239077232830298,2,0.9555833813546034,0.20466634,"cell wall organization"),
                     c("GO:1901696","cannabinoid biosynthetic process",0.0003056277806582375,1,0.9541501166999655,0.03891167,"cannabinoid biosynthetic process"),
                     c("GO:1990067","intrachromosomal DNA recombination",3.4506362332381646E-05,1,0.9346695357018706,0.08958079,"intrachromosomal DNA recombination"),
                     c("GO:0006270","DNA replication initiation",0.18677061560434888,1,0.8915880117660419,0.25665205,"intrachromosomal DNA recombination"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_BP_Table.tsv", sep = "\t")
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

