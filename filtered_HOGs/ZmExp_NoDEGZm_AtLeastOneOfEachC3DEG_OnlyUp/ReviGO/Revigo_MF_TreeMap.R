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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,9,1,-0,"molecular_function"),
c("GO:0003682","chromatin binding",0.6287702097345026,6,0.9556687776402477,0.04661592,"chromatin binding"),
c("GO:0003729","mRNA binding",1.136762432364793,31,0.9300196903212681,0,"mRNA binding"),
c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,2,0.9476843341207317,0.42104634,"mRNA binding"),
c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,5,0.9193094385579699,0.31664803,"mRNA binding"),
c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,2,0.9517726307655306,0.19519762,"mRNA binding"),
c("GO:0003684","damaged DNA binding",0.33181460177613226,1,0.9403499025777665,0.40218543,"mRNA binding"),
c("GO:0003697","single-stranded DNA binding",0.40260363864918647,1,0.9393076541150771,0.41085789,"mRNA binding"),
c("GO:0003713","transcription coactivator activity",0.24395086691346182,3,0.9693393828269257,0.29461398,"mRNA binding"),
c("GO:0003723","RNA binding",6.099813894661886,48,0.927900773543416,0.40176997,"mRNA binding"),
c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9401925612517594,0.48968021,"mRNA binding"),
c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.939848294598848,0.47960958,"mRNA binding"),
c("GO:0003743","translation initiation factor activity",0.37315089336372126,8,0.9172452462340493,0.24940265,"mRNA binding"),
c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8881524585665113,0.3505363,"mRNA binding"),
c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8893166604050521,0.33108637,"mRNA binding"),
c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8885916556311917,0.3394079,"mRNA binding"),
c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9468722449079794,0.428381,"mRNA binding"),
c("GO:0008327","methyl-CpG binding",0.009309036644208929,1,0.9390418609120744,0.46750139,"mRNA binding"),
c("GO:0008878","glucose-1-phosphate adenylyltransferase activity",0.021398809818155354,1,0.8890101964432305,0.29083106,"mRNA binding"),
c("GO:0016780","phosphotransferase activity, for other substituted phosphate groups",0.32920461222629094,1,0.8716406785187323,0.35102513,"mRNA binding"),
c("GO:0020037","heme binding",1.453841791525211,2,0.9578412310689256,0.13060477,"mRNA binding"),
c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9488267691429012,0.41077037,"mRNA binding"),
c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8701342346445498,0.27195941,"mRNA binding"),
c("GO:0044183","protein folding chaperone",0.39473153762799984,1,0.9779221791394449,0.32443613,"mRNA binding"),
c("GO:0003735","structural constituent of ribosome",2.128498588986873,6,0.9957656723800591,-0,"structural constituent of ribosome"),
c("GO:0003774","cytoskeletal motor activity",0.3926027441124113,1,1,-0,"cytoskeletal motor activity"),
c("GO:0003824","catalytic activity",60.727864217389204,13,1,-0,"catalytic activity"),
c("GO:0005096","GTPase activator activity",0.43493469016718705,2,0.9847701438688098,-0,"GTPase activator activity"),
c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
c("GO:0005515","protein binding",8.610051728351934,64,0.9676666389088371,0.07766851,"protein binding"),
c("GO:0005543","phospholipid binding",0.6999096105403309,2,0.9593664947520258,0.04714527,"phospholipid binding"),
c("GO:0008270","zinc ion binding",3.7731179768939866,9,0.9407825034733968,0.05738816,"zinc ion binding"),
c("GO:0000166","nucleotide binding",18.476222463002568,48,0.9174794025435898,0.48160605,"zinc ion binding"),
c("GO:0003676","nucleic acid binding",20.580503808256374,36,0.9411278764904909,0.20275684,"zinc ion binding"),
c("GO:0005509","calcium ion binding",1.3422486614434648,4,0.9515337555035566,0.3693843,"zinc ion binding"),
c("GO:0005525","GTP binding",1.780363236924562,9,0.9294139947809361,0.19635374,"zinc ion binding"),
c("GO:0016208","AMP binding",0.08038413014592037,1,0.9433231958232332,0.29697106,"zinc ion binding"),
c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9588095165721036,0.26177282,"zinc ion binding"),
c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9751546291666777,0.2693746,"zinc ion binding"),
c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9706463507141544,0.11877911,"zinc ion binding"),
c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.9348882157731002,0.33990845,"zinc ion binding"),
c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,5,0.9495984971076281,0.18168089,"zinc ion binding"),
c("GO:0070403","NAD+ binding",0.16643173804060435,1,0.9466720577769359,0.29755636,"zinc ion binding"),
c("GO:0071949","FAD binding",0.630202710371034,1,0.9411153726895474,0.30272555,"zinc ion binding"),
c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9568488657554564,0.25647876,"zinc ion binding"),
c("GO:0008289","lipid binding",1.291468066123697,4,0.9740360747584395,0.05834827,"lipid binding"),
c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
c("GO:0015288","porin activity",0.12362746592455742,2,0.93884781391447,-0,"porin activity"),
c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9130671375083721,0.39044024,"porin activity"),
c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9417239126056423,0.47124298,"porin activity"),
c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,1,0.9411586434365969,0.31858441,"porin activity"),
c("GO:0032977","membrane insertase activity",0.038238453523758646,1,0.9623999851096069,0.18410598,"porin activity"),
c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9564962400079116,0.04318065,"CoA-ligase activity"),
c("GO:0016740","transferase activity",20.627439270288612,56,0.9472493489854454,0.0817378,"transferase activity"),
c("GO:0016491","oxidoreductase activity",11.23609814678324,19,0.9511880784179745,0.10406278,"transferase activity"),
c("GO:0016787","hydrolase activity",22.243200201100027,39,0.9467167859017094,0.12712323,"transferase activity"),
c("GO:0016757","glycosyltransferase activity",2.478461155483397,8,0.874378058566644,0.05667901,"glycosyltransferase activity"),
c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9146616131447959,0.20205602,"glycosyltransferase activity"),
c("GO:0008168","methyltransferase activity",2.5264987116585993,8,0.8279886835043281,0.33973871,"glycosyltransferase activity"),
c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,3,0.8786654680684459,0.2917808,"glycosyltransferase activity"),
c("GO:0008483","transaminase activity",0.6992399275802186,2,0.8785936626683389,0.28934569,"glycosyltransferase activity"),
c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8690028745887369,0.36085078,"glycosyltransferase activity"),
c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8844175700784569,0.30154932,"glycosyltransferase activity"),
c("GO:0016767","geranylgeranyl-diphosphate geranylgeranyltransferase activity",0.019117009268633918,1,0.9020886228737552,0.20439581,"glycosyltransferase activity"),
c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.8864348573497074,0.49223587,"glycosyltransferase activity"),
c("GO:0047216","inositol 3-alpha-galactosyltransferase activity",0.0009202596968429504,1,0.9192503600880888,0.43743078,"glycosyltransferase activity"),
c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.912186342655241,0.17480439,"glycosyltransferase activity"),
c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.8656573719390799,0.45431535,"glycosyltransferase activity"),
c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8670665535523243,0.39355689,"glycosyltransferase activity"),
c("GO:0016829","lyase activity",3.6221510367468346,4,0.9571217737769077,0.05997119,"lyase activity"),
c("GO:0016853","isomerase activity",2.413124934500793,6,0.958900996678351,0.05646078,"isomerase activity"),
c("GO:0016871","cycloartenol synthase activity",0.0006563780006397912,1,0.9622379275230969,0.02586138,"cycloartenol synthase activity"),
c("GO:0016874","ligase activity",3.248266153118885,3,0.9576146728115937,0.0589874,"ligase activity"),
c("GO:0016887","ATP hydrolysis activity",4.018519084389943,17,0.8653089897370425,-0,"ATP hydrolysis activity"),
c("GO:0000146","microfilament motor activity",0.092387421083296,1,0.9701164243091338,0.49829267,"ATP hydrolysis activity"),
c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,2,0.8780590553391917,0.42844795,"ATP hydrolysis activity"),
c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.915502130688966,0.3692388,"ATP hydrolysis activity"),
c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9398070340024738,0.19229254,"ATP hydrolysis activity"),
c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,2,0.9046890764331185,0.31839638,"ATP hydrolysis activity"),
c("GO:0004712","protein serine/threonine/tyrosine kinase activity",0.056286631048107494,1,0.8424023554143275,0.44384691,"ATP hydrolysis activity"),
c("GO:0004715","non-membrane spanning protein tyrosine kinase activity",0.060526478133321286,1,0.8395108893642167,0.44607451,"ATP hydrolysis activity"),
c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,3,0.8410457538068864,0.42141607,"ATP hydrolysis activity"),
c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,5,0.818388215743473,0.48775459,"ATP hydrolysis activity"),
c("GO:0008233","peptidase activity",4.167383840940452,14,0.8383411349272583,0.3656957,"ATP hydrolysis activity"),
c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,4,0.8511527842038522,0.45595338,"ATP hydrolysis activity"),
c("GO:0008568","microtubule severing ATPase activity",0.013493446398287598,1,0.8880141756101257,0.41923201,"ATP hydrolysis activity"),
c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9607133938713633,0.45552111,"ATP hydrolysis activity"),
c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8995729857130226,0.39048386,"ATP hydrolysis activity"),
c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.9103623179383683,0.32827501,"ATP hydrolysis activity"),
c("GO:0033919","glucan 1,3-alpha-glucosidase activity",0.000915824710352141,1,0.940637920264108,0.4867618,"ATP hydrolysis activity"),
c("GO:0052742","phosphatidylinositol kinase activity",0.07485148449863564,1,0.8728322439305342,0.40778693,"ATP hydrolysis activity"),
c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9359184516152639,0.17583997,"ATP hydrolysis activity"),
c("GO:0106310","protein serine kinase activity",0.08584138102286135,5,0.8377396328998755,0.37283533,"ATP hydrolysis activity"),
c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9351685479926555,0.12676211,"ATP hydrolysis activity"),
c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9371411775817016,0.39401957,"ATP hydrolysis activity"),
c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9746923776739819,0.04901613,"carbohydrate binding"),
c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.972489785549527,0.06283846,"protein-containing complex binding"),
c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9721699882542946,0.03325955,"ER retention sequence binding"),
c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.979940750299162,0.49249984,"ER retention sequence binding"),
c("GO:0047769","arogenate dehydratase activity",0.02556769711951619,1,0.9482814083541394,0.03410747,"arogenate dehydratase activity"),
c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.954541698643944,0.36924486,"arogenate dehydratase activity"),
c("GO:0051082","unfolded protein binding",0.5891679978648201,4,0.9154767838810479,0.0463004,"unfolded protein binding"),
c("GO:0000149","SNARE binding",0.26943873427614345,1,0.9206649941369229,0.40559978,"unfolded protein binding"),
c("GO:0003779","actin binding",0.7421262469463454,3,0.9093692418236538,0.44654025,"unfolded protein binding"),
c("GO:0005516","calmodulin binding",0.264039138223583,1,0.9207908819958428,0.40485813,"unfolded protein binding"),
c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9296308010906706,0.35368117,"unfolded protein binding"),
c("GO:0019900","kinase binding",0.4400460120978449,2,0.9080501720806812,0.42444044,"unfolded protein binding"),
c("GO:0030276","clathrin binding",0.11822565237875157,2,0.9254854729266847,0.37746303,"unfolded protein binding"),
c("GO:0030544","Hsp70 protein binding",0.07279586826014549,2,0.9280593508043852,0.36265291,"unfolded protein binding"),
c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9236883761991721,0.38789043,"unfolded protein binding"),
c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9284681538408958,0.36031388,"unfolded protein binding"),
c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,2,0.9477390078220179,0.25391209,"unfolded protein binding"),
c("GO:0032051","clathrin light chain binding",0.02013705616152008,1,0.9308049285888775,0.3284969,"unfolded protein binding"),
c("GO:0042393","histone binding",0.41579772345934446,3,0.9178683476412992,0.42217449,"unfolded protein binding"),
c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9166127381504338,0.42967902,"unfolded protein binding"),
c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9220851670854019,0.39725472,"unfolded protein binding"),
c("GO:0043621","protein self-association",0.05684543934594948,3,0.9293047565349698,0.3555383,"unfolded protein binding"),
c("GO:0046982","protein heterodimerization activity",0.29517274338906496,1,0.917948559934934,0.40897571,"unfolded protein binding"),
c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9103192038989392,0.47948026,"unfolded protein binding"),
c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9459807677703748,0.26331899,"unfolded protein binding"),
c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.936850692254367,0.31312228,"unfolded protein binding"),
c("GO:1990935","splicing factor binding",0.0007095978385295039,1,0.9458890993090746,0.26381104,"unfolded protein binding"),
c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9171112021374623,0.048252,"dioxygenase activity"),
c("GO:0004322","ferroxidase activity",0.06266857660838222,1,0.931427241325847,0.3160356,"dioxygenase activity"),
c("GO:0004601","peroxidase activity",0.4792135952914281,2,0.9129073152886723,0.37878622,"dioxygenase activity"),
c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9332768155056413,0.2837958,"dioxygenase activity"),
c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9043323631708091,0.44913707,"dioxygenase activity"),
c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9054338817098204,0.38257989,"dioxygenase activity"),
c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9150426384166358,0.41045836,"dioxygenase activity"),
c("GO:0016651","oxidoreductase activity, acting on NAD(P)H",0.7914588191768638,1,0.9166500814288117,0.39829104,"dioxygenase activity"),
c("GO:0016671","oxidoreductase activity, acting on a sulfur group of donors, disulfide as acceptor",0.10547063123118373,1,0.9237378380452589,0.33002637,"dioxygenase activity"));

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

