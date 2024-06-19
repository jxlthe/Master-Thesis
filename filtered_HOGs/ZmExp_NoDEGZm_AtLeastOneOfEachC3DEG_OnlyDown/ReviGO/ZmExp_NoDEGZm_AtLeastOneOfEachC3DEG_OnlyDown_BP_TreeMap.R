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
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.07909597668631853,3,0.9217149565546587,0.05769476,"malate metabolic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,1,0.9229464063040455,0.3067044,"malate metabolic process"),
                     c("GO:0045487","gibberellin catabolic process",0.00595974172283563,2,0.8762033265144406,0.26507989,"malate metabolic process"),
                     c("GO:0006811","monoatomic ion transport",4.7761710300947735,10,0.9122748819072783,-0,"monoatomic ion transport"),
                     c("GO:0010497","plasmodesmata-mediated intercellular transport",0.008392440267268536,1,0.9432813994604033,0.22103265,"monoatomic ion transport"),
                     c("GO:0015031","protein transport",3.093438694074183,4,0.8694349327399881,0.42145531,"monoatomic ion transport"),
                     c("GO:0015692","lead ion transport",0.001227440602966147,1,0.9349354269439446,0.19142099,"monoatomic ion transport"),
                     c("GO:0015783","GDP-fucose transmembrane transport",0.013410651246463469,2,0.916147426399219,0.42609188,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,4,0.8762118526965018,0.44154323,"monoatomic ion transport"),
                     c("GO:0043090","amino acid import",7.887168533115806E-05,1,0.9449366276391375,0.16068432,"monoatomic ion transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9062455303621825,0.44310395,"monoatomic ion transport"),
                     c("GO:0051646","mitochondrion localization",0.04394385243028803,1,0.9298546266470985,0.24160098,"monoatomic ion transport"),
                     c("GO:0072583","clathrin-dependent endocytosis",0.06990989008540521,2,0.9308555064685979,0.26649209,"monoatomic ion transport"),
                     c("GO:0072699","protein localization to cortical microtubule cytoskeleton",0.00015774337066231612,1,0.9433731536634987,0.15867185,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,20,0.8876272230726193,0,"defense response"),
                     c("GO:0000165","MAPK cascade",0.14246691110973742,2,0.8588402128071748,0.45634702,"defense response"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8956770006647753,0.4896406,"defense response"),
                     c("GO:0006955","immune response",0.6433982378290884,2,0.9111344210371278,0.30172547,"defense response"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.12287469152544445,1,0.8637303629479638,0.44951817,"defense response"),
                     c("GO:0007300","ovarian nurse cell to oocyte transport",0.00031548674132463224,1,0.8738949051222475,0.46428331,"defense response"),
                     c("GO:0009611","response to wounding",0.16569462243976346,2,0.9035040983547749,0.45786068,"defense response"),
                     c("GO:0009641","shade avoidance",0.0003204162216578296,1,0.9195167266622747,0.49975625,"defense response"),
                     c("GO:0009744","response to sucrose",0.028226204387888188,3,0.8965158682779218,0.22702886,"defense response"),
                     c("GO:0010117","photoprotection",0.00019717921332789513,1,0.9222402904495552,0.48975475,"defense response"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9207574345423947,0.20406574,"defense response"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.252707344541197,2,0.8872073311603659,0.47637699,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,3,0.8183183765622872,0.37522238,"defense response"),
                     c("GO:0042742","defense response to bacterium",0.140633144425788,6,0.882390566826153,0.45104908,"defense response"),
                     c("GO:0046740","transport of virus in host, cell to cell",0.0009094891214749163,1,0.9813132397214319,0.48259105,"defense response"),
                     c("GO:0048573","photoperiodism, flowering",0.0014763793597926149,2,0.8198788496978948,0.46247226,"defense response"),
                     c("GO:0071284","cellular response to lead ion",9.366012633075019E-05,1,0.9183645836820575,0.26100449,"defense response"),
                     c("GO:0071456","cellular response to hypoxia",0.026404761404771757,4,0.859551573207622,0.39162542,"defense response"),
                     c("GO:0071497","cellular response to freezing",7.640694516455937E-05,1,0.9085482240172137,0.40112537,"defense response"),
                     c("GO:0007049","cell cycle",2.8073119376140743,2,0.9904146740988637,0.01433045,"cell cycle"),
                     c("GO:0007623","circadian rhythm",0.10049731555289494,1,0.9955431861403087,-0,"circadian rhythm"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,4,1,-0,"metabolic process"),
                     c("GO:0008219","cell death",0.4651211168388386,1,0.9918986929793846,0.01044746,"cell death"),
                     c("GO:0015979","photosynthesis",0.228607115192195,2,0.9607037950593874,0.06311855,"photosynthesis"),
                     c("GO:0016477","cell migration",0.49093927008395993,1,0.9918610395838999,0.01050372,"cell migration"),
                     c("GO:0016567","protein ubiquitination",1.2678648064385323,12,0.9027159603986793,-0,"protein ubiquitination"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,7,0.9465354534227138,0.15163831,"protein ubiquitination"),
                     c("GO:0005985","sucrose metabolic process",0.0829064649838801,1,0.928854251766739,0.44424972,"protein ubiquitination"),
                     c("GO:0006004","fucose metabolic process",0.06850745293061056,2,0.9071417508253294,0.43756433,"protein ubiquitination"),
                     c("GO:0006071","glycerol metabolic process",0.19450743498730216,1,0.8923138496616365,0.4768145,"protein ubiquitination"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,3,0.8025323564716771,0.41940643,"protein ubiquitination"),
                     c("GO:0006260","DNA replication",1.488685807444442,2,0.9033179992280823,0.21165279,"protein ubiquitination"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8896195004172376,0.37832162,"protein ubiquitination"),
                     c("GO:0006457","protein folding",1.174377211919444,4,0.8915054911061395,0.37292242,"protein ubiquitination"),
                     c("GO:0006508","proteolysis",5.2622572267907,5,0.9046114465127817,0.44941082,"protein ubiquitination"),
                     c("GO:0006520","amino acid metabolic process",5.369313215926915,1,0.9465016164309483,0.15178462,"protein ubiquitination"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,10,0.945311378466218,0.12094833,"protein ubiquitination"),
                     c("GO:0006694","steroid biosynthetic process",0.2949431320360321,2,0.8659733844576685,0.33431465,"protein ubiquitination"),
                     c("GO:0006807","nitrogen compound metabolic process",40.19370343674496,3,0.9575436664180312,0.12823919,"protein ubiquitination"),
                     c("GO:0008380","RNA splicing",0.9742871404547958,1,0.8911240835883005,0.39771984,"protein ubiquitination"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,2,0.9596258578399225,0.21768792,"protein ubiquitination"),
                     c("GO:0009073","aromatic amino acid family biosynthetic process",0.5154831526629496,1,0.8632746428917868,0.41955763,"protein ubiquitination"),
                     c("GO:0015995","chlorophyll biosynthetic process",0.07320031820781447,2,0.898647303381967,0.16277389,"protein ubiquitination"),
                     c("GO:0016310","phosphorylation",5.235381700014107,24,0.9201586061779947,0.46324041,"protein ubiquitination"),
                     c("GO:0018345","protein palmitoylation",0.0003918936864891916,1,0.927482938228154,0.37471442,"protein ubiquitination"),
                     c("GO:0030244","cellulose biosynthetic process",0.03957140337474196,4,0.8903472165832508,0.14129671,"protein ubiquitination"),
                     c("GO:0031146","SCF-dependent proteasomal ubiquitin-dependent protein catabolic process",0.12335531585793119,2,0.9054990986254112,0.46794754,"protein ubiquitination"),
                     c("GO:0033314","mitotic DNA replication checkpoint signaling",0.042962885843981745,1,0.8389238091593698,0.4495745,"protein ubiquitination"),
                     c("GO:0033511","luteolin biosynthetic process",2.711214183258558E-05,1,0.9474209773171544,0.45767193,"protein ubiquitination"),
                     c("GO:0042761","very long-chain fatty acid biosynthetic process",0.055542919654301456,1,0.8653016314022828,0.4941006,"protein ubiquitination"),
                     c("GO:0042853","L-alanine catabolic process",0.0402368832197236,1,0.8984918709082395,0.42825993,"protein ubiquitination"),
                     c("GO:0046488","phosphatidylinositol metabolic process",0.5598188987797269,2,0.8680597043051811,0.47629275,"protein ubiquitination"),
                     c("GO:0046777","protein autophosphorylation",0.001434478776960437,1,0.9247723766968365,0.40604291,"protein ubiquitination"),
                     c("GO:0046938","phytochelatin biosynthetic process",0.003766122974562797,1,0.9264348474927936,0.20437656,"protein ubiquitination"),
                     c("GO:0055047","generative cell mitosis",7.147746483136198E-05,1,0.8982388243261324,0.45916025,"protein ubiquitination"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.9161100139273048,0.45115494,"protein ubiquitination"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9444518331173675,0.45105203,"protein ubiquitination"),
                     c("GO:1901576","organic substance biosynthetic process",28.21764434528959,1,0.9094524352519981,0.21493739,"protein ubiquitination"),
                     c("GO:1990918","double-strand break repair involved in meiotic recombination",0.009479390680738558,1,0.7765808641048081,0.40319107,"protein ubiquitination"),
                     c("GO:0030010","establishment of cell polarity",0.11765190711242182,1,0.9927490803456864,0.00919456,"establishment of cell polarity"),
                     c("GO:0030154","cell differentiation",2.279586420543629,3,0.8990410185370589,0.01240155,"cell differentiation"),
                     c("GO:0009826","unidimensional cell growth",0.016627137163874758,2,0.9143084839834023,0.48841036,"cell differentiation"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,2,0.9219472873147521,0.48986076,"cell differentiation"),
                     c("GO:0055046","microgametogenesis",0.0014049018949612527,1,0.9243587561470956,0.40770388,"cell differentiation"),
                     c("GO:0032259","methylation",2.6278542060840238,1,0.97034765211817,0.05843136,"methylation"),
                     c("GO:0040011","locomotion",0.5148398554794674,1,1,-0,"locomotion"),
                     c("GO:0048544","recognition of pollen",0.04932191547380636,1,0.9931984117235034,0.00854634,"recognition of pollen"),
                     c("GO:0050896","response to stimulus",17.567785530535815,4,1,-0,"response to stimulus"),
                     c("GO:0060361","flight",9.612486649734888E-05,1,1,-0,"flight"),
                     c("GO:0070814","hydrogen sulfide biosynthetic process",0.03377679924306844,1,0.9402148354754484,0.03784896,"hydrogen sulfide biosynthetic process"),
                     c("GO:0000103","sulfate assimilation",0.10156701278519879,1,0.9563161258332863,0.48363934,"hydrogen sulfide biosynthetic process"),
                     c("GO:0042350","GDP-L-fucose biosynthetic process",0.02421607213683212,1,0.8835380103040614,0.34490714,"hydrogen sulfide biosynthetic process"),
                     c("GO:0046166","glyceraldehyde-3-phosphate biosynthetic process",0.03550458209985412,1,0.9116980138604488,0.11775102,"hydrogen sulfide biosynthetic process"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,8,0.9197627632454867,-0,"cell wall organization"),
                     c("GO:0006325","chromatin organization",1.3384303174082528,3,0.9293852555636286,0.4063955,"cell wall organization"),
                     c("GO:0006997","nucleus organization",0.12424015757774012,1,0.942052264292397,0.37383592,"cell wall organization"),
                     c("GO:0009662","etioplast organization",0.0006778035458146395,1,0.9517053599085655,0.22083297,"cell wall organization"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9401213080094034,0.32062999,"cell wall organization"),
                     c("GO:0051013","microtubule severing",0.012493767904488756,1,0.9408521597191466,0.26779102,"cell wall organization"),
                     c("GO:0070206","protein trimerization",0.00727591297179933,1,0.9495513901672242,0.4524177,"cell wall organization"),
                     c("GO:0140964","intracellular auxin homeostasis",0.00018239077232830298,1,0.9798816262649623,-0,"intracellular auxin homeostasis"),
                     c("GO:0098771","inorganic ion homeostasis",1.0570555799893464,1,0.9677054967424179,0.45565197,"intracellular auxin homeostasis"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,3,0.915753539696631,-0,"regulation of defense response to fungus"),
                     c("GO:0010105","negative regulation of ethylene-activated signaling pathway",0.005397780964851129,1,0.9119840385339306,0.43825323,"regulation of defense response to fungus"),
                     c("GO:1900457","regulation of brassinosteroid mediated signaling pathway",0.001054908791304239,1,0.9314577900279192,0.34087448,"regulation of defense response to fungus"),
                     c("GO:1900458","negative regulation of brassinosteroid mediated signaling pathway",0.0003844994659893955,1,0.9228448013437632,0.32535433,"regulation of defense response to fungus"),
                     c("GO:1901001","negative regulation of response to salt stress",0.00037710524548959944,1,0.9192906002379682,0.46916415,"regulation of defense response to fungus"),
                     c("GO:2000022","regulation of jasmonic acid mediated signaling pathway",0.011707515791343773,1,0.9222577951337201,0.42518681,"regulation of defense response to fungus"),
                     c("GO:2000028","regulation of photoperiodism, flowering",0.006859371883644151,1,0.9233401420289508,0.37952976,"regulation of defense response to fungus"),
                     c("GO:2000031","regulation of salicylic acid mediated signaling pathway",0.008145966250608667,2,0.9238018124888107,0.37732875,"regulation of defense response to fungus"),
                     c("GO:2000652","regulation of secondary cell wall biogenesis",0.0004904832931531391,1,0.9572395386714166,0.08761323,"regulation of secondary cell wall biogenesis"),
                     c("GO:1903527","positive regulation of membrane tubulation",0.0019298915504467734,1,0.9338975824675887,0.46246533,"regulation of secondary cell wall biogenesis"),
                     c("GO:1903553","positive regulation of extracellular exosome assembly",0.0003918936864891916,1,0.9364922075877135,0.40819627,"regulation of secondary cell wall biogenesis"),
                     c("GO:2001006","regulation of cellulose biosynthetic process",0.002356291599268347,2,0.947090576042956,0.0943689,"regulation of cellulose biosynthetic process"),
                     c("GO:0009789","positive regulation of abscisic acid-activated signaling pathway",0.00080843477464437,1,0.9214809266030887,0.38244242,"regulation of cellulose biosynthetic process"),
                     c("GO:0009934","regulation of meristem structural organization",0.0007295630893132119,1,0.9488957450726474,0.37175537,"regulation of cellulose biosynthetic process"),
                     c("GO:0031047","regulatory ncRNA-mediated gene silencing",0.26322685557224024,2,0.9070631389466772,0.33664393,"regulation of cellulose biosynthetic process"),
                     c("GO:0045893","positive regulation of DNA-templated transcription",1.588722216586183,2,0.899438868614809,0.23293868,"regulation of cellulose biosynthetic process"),
                     c("GO:0051302","regulation of cell division",0.16982799169914947,1,0.9414309705668448,0.18590747,"regulation of cellulose biosynthetic process"),
                     c("GO:0051510","regulation of unidimensional cell growth",0.00424921204721614,1,0.9265414734291049,0.45200802,"regulation of cellulose biosynthetic process"),
                     c("GO:0051726","regulation of cell cycle",0.9132256675674799,2,0.9331342057273261,0.13410528,"regulation of cellulose biosynthetic process"),
                     c("GO:0090333","regulation of stomatal closure",0.012828972567146176,1,0.9506863636529148,0.15074542,"regulation of cellulose biosynthetic process"),
                     c("GO:1901141","regulation of lignin biosynthetic process",0.00025140349699306625,1,0.9510187473965689,0.14312882,"regulation of cellulose biosynthetic process"),
                     c("GO:1905038","regulation of membrane lipid metabolic process",0.009506502822571143,1,0.9380563272407281,0.25241702,"regulation of cellulose biosynthetic process"),
                     c("GO:1905157","positive regulation of photosynthesis",0.00010598382716374363,1,0.9456428144310399,0.39735289,"regulation of cellulose biosynthetic process"),
                     c("GO:1905183","negative regulation of protein serine/threonine phosphatase activity",1.7253181166190823E-05,1,0.9463200482012567,0.1289891,"regulation of cellulose biosynthetic process"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9433671437873972,0.49983935,"regulation of cellulose biosynthetic process"),
                     c("GO:2000067","regulation of root morphogenesis",0.00304641884591598,1,0.9442925564043755,0.347328,"regulation of cellulose biosynthetic process"),
                     c("GO:2000306","positive regulation of photomorphogenesis",0.00015281389032911874,1,0.9191016589735107,0.33080078,"regulation of cellulose biosynthetic process"),
                     c("GO:2000369","regulation of clathrin-dependent endocytosis",0.010906475237199198,1,0.9296953924094261,0.14897448,"regulation of cellulose biosynthetic process"),
                     c("GO:2000377","regulation of reactive oxygen species metabolic process",0.02678186665026136,1,0.9396520077476884,0.26299316,"regulation of cellulose biosynthetic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_BP_Table.tsv", sep = "\t")
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

