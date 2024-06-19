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
revigo.data <- rbind(c("GO:0000995","RNA polymerase III general transcription initiation factor activity",0.007712441507517546,1,0.9841871246288971,0.00722361,"RNA polymerase III general transcription initiation factor activity"),
                     c("GO:0003674","molecular_function",100,32,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,10,0.9689397995337248,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,42,0.9521522430374141,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,2,0.9641766955594843,0.42104634,"mRNA binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,6,0.9791567156099864,0.36500624,"mRNA binding"),
                     c("GO:0001006","RNA polymerase III type 3 promoter sequence-specific DNA binding",0.0190615719374988,1,0.9608605950710434,0.36290179,"mRNA binding"),
                     c("GO:0001217","DNA-binding transcription repressor activity",0.09621481442486451,1,0.9835529779383934,0.49577888,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,96,0.9454417456386672,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,3,0.966190227513994,0.19519762,"mRNA binding"),
                     c("GO:0003684","damaged DNA binding",0.33181460177613226,1,0.9582563633962112,0.4175852,"mRNA binding"),
                     c("GO:0003697","single-stranded DNA binding",0.40260363864918647,2,0.9575327891197729,0.42694223,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.959081297510557,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9591098450616944,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,11,0.9324790327408994,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,2,0.9636243827418178,0.428381,"mRNA binding"),
                     c("GO:0008327","methyl-CpG binding",0.009309036644208929,1,0.9540229178639519,0.34771157,"mRNA binding"),
                     c("GO:0019237","centromeric DNA binding",0.004625690909914204,1,0.9647495718932401,0.20059314,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,9,0.9680097931257933,0.13060477,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.964953719942075,0.41077037,"mRNA binding"),
                     c("GO:0035198","miRNA binding",0.029599099839661934,1,0.9644971613561167,0.41680219,"mRNA binding"),
                     c("GO:0042134","rRNA primary transcript binding",0.01629192287398833,1,0.9659286451165663,0.39794793,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9587398118533373,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,21,0.9458079172850538,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,6,0.9959828945185376,-0,"structural constituent of ribosome"),
                     c("GO:0003774","cytoskeletal motor activity",0.3926027441124113,1,1,-0,"cytoskeletal motor activity"),
                     c("GO:0003824","catalytic activity",60.727864217389204,34,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,62,0.7991312682332853,0,"protein kinase activity"),
                     c("GO:0000234","phosphoethanolamine N-methyltransferase activity",0.0010222643861315665,1,0.9005656057222542,0.45971198,"protein kinase activity"),
                     c("GO:0003755","peptidyl-prolyl cis-trans isomerase activity",0.32422412239711196,4,0.8861130940300226,0.4205441,"protein kinase activity"),
                     c("GO:0003975","UDP-N-acetylglucosamine-dolichyl-phosphate N-acetylglucosaminephosphotransferase activity",0.005718915079898721,1,0.9032388864135855,0.35348317,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8893158783180034,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8898915598540006,0.48555173,"protein kinase activity"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8899083254406018,0.46089842,"protein kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,2,0.8913250950727402,0.25688843,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8465994648490534,0.49507895,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,10,0.8545818422016799,0.41376722,"protein kinase activity"),
                     c("GO:0004801","transaldolase activity",0.029270910839342038,1,0.9113911394817553,0.21824965,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,22,0.8235425141772773,0.47753726,"protein kinase activity"),
                     c("GO:0008146","sulfotransferase activity",0.16979789278712867,1,0.8972810614263171,0.25770873,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,15,0.8467149998284565,0.35677507,"protein kinase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,7,0.8646788483075069,0.30425857,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,7,0.8647495345145602,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,3,0.871128081007062,0.35123485,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,6,0.8791448240759565,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8934390285663193,0.40329157,"protein kinase activity"),
                     c("GO:0008798","beta-aspartyl-peptidase activity",0.023041972313000234,1,0.8906679336776527,0.44886191,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.8506214055277963,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,11,0.8716753777217396,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,37,0.876628926563111,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,2,0.8859296476106466,0.31489566,"protein kinase activity"),
                     c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.04143767537877,2,0.8620848457185354,0.43304972,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.8988388139109963,0.49223587,"protein kinase activity"),
                     c("GO:0033855","nicotianamine aminotransferase activity",6.6524797362141E-06,1,0.9356878940872976,0.45711636,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8743065572322031,0.49253483,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9150530150150219,0.30679413,"protein kinase activity"),
                     c("GO:0047326","inositol tetrakisphosphate 5-kinase activity",0.0016542499610719059,1,0.9024036089977014,0.40502966,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,2,0.8849374424809161,0.48733868,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9222526746562715,0.29376839,"protein kinase activity"),
                     c("GO:0050734","hydroxycinnamoyltransferase activity",0.0034969868480032116,1,0.90309679849696,0.35899295,"protein kinase activity"),
                     c("GO:0051753","mannan synthase activity",0.005627997856837128,2,0.8976585684333968,0.49140545,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9039044082734732,0.47557181,"protein kinase activity"),
                     c("GO:0052667","phosphomethylethanolamine N-methyltransferase activity",0.00011530964876104439,1,0.9103636574993704,0.40189746,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9213773782071678,0.48656541,"protein kinase activity"),
                     c("GO:0070012","oligopeptidase activity",0.014826159838775825,1,0.8933720677284,0.31386666,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.9013348578283663,0.18579171,"protein kinase activity"),
                     c("GO:0102732","myo-inositol-1,2,3,4,6-heptakisphosphate 5-kinase activity",0.0006741179466030287,1,0.9086267095387811,0.38119148,"protein kinase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.8843734136147638,0.45431535,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8832947747905774,0.39355689,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,6,0.9746143375880698,-0,"GTPase activator activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,173,0.9765859638010934,0.07766851,"protein binding"),
                     c("GO:0097159","organic cyclic compound binding",39.22969739688024,1,0.9707619643887931,0.13134756,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,30,0.9548331501228443,0.05738816,"zinc ion binding"),
                     c("GO:0000035","acyl binding",0.03278342014006308,1,0.9779977426252987,0.11604884,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,148,0.9381184188920123,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,6,0.9622465842144501,0.3830934,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,2,0.9718546830497732,0.21260669,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,52,0.9580118674808056,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,12,0.9632929394898434,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,20,0.9471629003822944,0.19635374,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,5,0.9562520289102728,0.24645234,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9551793223880397,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.95672267465334,0.29697106,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,2,0.9651036010951621,0.26177282,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,3,0.9505238324807589,0.31800282,"zinc ion binding"),
                     c("GO:0030955","potassium ion binding",0.05969270067304912,1,0.9720836061172309,0.26193125,"zinc ion binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9821582854632487,0.2693746,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,2,0.956819596612351,0.2885161,"zinc ion binding"),
                     c("GO:0043325","phosphatidylinositol-3,4-bisphosphate binding",0.013358179310317913,1,0.9737742438971143,0.211085,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9775729510461595,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,7,0.9506049035463411,0.33990845,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,1,0.9564209130429503,0.34531063,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,5,0.955636599785464,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,12,0.9623105778002776,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9605050879749866,0.28523818,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,3,0.9590380608318773,0.29755636,"zinc ion binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,4,0.9552335844763332,0.30272555,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,6,0.9811618864945685,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,3,1,-0,"electron transfer activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,2,0.9870190666322757,0.02769872,"auxin binding"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,6,0.9143108305903176,-0,"antiporter activity"),
                     c("GO:0005338","nucleotide-sugar transmembrane transporter activity",0.0598301852542642,2,0.9275650340655242,0.3876045,"antiporter activity"),
                     c("GO:0015079","potassium ion transmembrane transporter activity",0.4915428577358782,3,0.9006833655355473,0.39988553,"antiporter activity"),
                     c("GO:0015267","channel activity",1.6825096949913436,3,0.9041584076811204,0.45334576,"antiporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9440732240270157,0.46922631,"antiporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9296392077666529,0.35182404,"antiporter activity"),
                     c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,2,0.9245893969014204,0.41895885,"antiporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,3,0.9555101850473171,0.20645071,"antiporter activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,2,0.9235065769550767,0.36857202,"antiporter activity"),
                     c("GO:0042910","xenobiotic transmembrane transporter activity",0.2653474592383718,1,0.9291784951893058,0.37757719,"antiporter activity"),
                     c("GO:0051119","sugar transmembrane transporter activity",0.13162596406073218,2,0.9178488318907265,0.35505294,"antiporter activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9343558497731906,0.4373877,"antiporter activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,7,0.8810902927317588,0.0531541,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004322","ferroxidase activity",0.06266857660838222,1,0.9179378657781734,0.34042691,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004324","ferredoxin-NADP+ reductase activity",0.018482806200448173,1,0.9239712537810244,0.30749203,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,5,0.8975107343242501,0.46480825,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,4,0.8877923428856705,0.41437051,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0010242","oxygen evolving activity",0.0031710153409287207,1,0.9312468395343216,0.26980743,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0015930","glutamate synthase activity",0.041380641452497105,2,0.9158182898465297,0.32846804,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9200321084537753,0.46927032,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,3,0.8964682192013551,0.41891471,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,2,0.8998981679879802,0.44913707,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016630","protochlorophyllide reductase activity",0.0022374506846133423,1,0.92786850794828,0.26342111,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,2,0.9090472749922742,0.39197077,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016651","oxidoreductase activity, acting on NAD(P)H",0.7914588191768638,1,0.90165282556963,0.4378256,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,6,0.8961932035919564,0.47359913,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9079736615097322,0.49506653,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.924129806024186,0.41198066,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0019139","cytokinin dehydrogenase activity",0.0057233500663895305,1,0.928972415421427,0.2813582,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0045543","gibberellin 2-beta-dioxygenase activity",0.006071496505918068,1,0.918140273601625,0.48062226,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0047037","salutaridine reductase (NADPH) activity",1.77399459632376E-05,1,0.9382493475303578,0.38914098,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0047501","(+)-neomenthol dehydrogenase activity",2.66099189448564E-05,1,0.9371855713015793,0.3973742,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0048529","magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase activity",0.002614424536332141,1,0.9256782569963189,0.33893217,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,6,0.9021566538274115,0.43460917,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102111","gibberellin A20,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102652","gibberellin A9,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102721","ubiquinol:oxygen oxidoreductase activity",0.0020112663735820627,1,0.9300274199055578,0.26152924,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102772","sphingolipid C4-monooxygenase activity",4.4349864908093996E-05,1,0.9387728340281034,0.20804942,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102924","gibberellin A44,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9355263830457752,0.35384313,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0106292","superoxide-generating NADPH oxidase activity",0.0027186467188661623,1,0.9270682711279755,0.26695083,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:1990136","linoleate 9S-lipoxygenase activity",0.00026609918944856393,1,0.9338910724870338,0.23015728,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,199,0.94795157278954,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,59,0.9514904785786901,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,96,0.9474747332767991,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,19,0.9568951737037287,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,6,0.9299383755313524,0.04651897,"carboxy-lyase activity"),
                     c("GO:0009978","allene oxide synthase activity",9.09172230615927E-05,1,0.9582979753288584,0.33685723,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.949980592318686,0.45277577,"carboxy-lyase activity"),
                     c("GO:0047769","arogenate dehydratase activity",0.02556769711951619,2,0.9441430404368263,0.47917335,"carboxy-lyase activity"),
                     c("GO:0106099","2-keto-3-deoxy-L-rhamnonate aldolase activity",0.00019513940559561356,1,0.9563425525954127,0.48677249,"carboxy-lyase activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9642486139462519,0.28175041,"carboxy-lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,18,0.9585373428327074,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,15,0.9573490421877254,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,28,0.8699916892191324,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9214819569250295,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004301","epoxide hydrolase activity",0.06371079843372243,1,0.9311110808384813,0.22109293,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9392140747595851,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,11,0.8900502416316907,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008568","microtubule severing ATPase activity",0.013493446398287598,1,0.9039060921890452,0.41923201,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.9014398579567025,0.38812549,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,11,0.9114266549752635,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9116529749120368,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9136811150660402,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,2,0.9316529135975979,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0033919","glucan 1,3-alpha-glucosidase activity",0.000915824710352141,1,0.9276848533170794,0.4867618,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,2,0.9322888130169236,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047668","amygdalin beta-glucosidase activity",2.4392425699451695E-05,1,0.9382288266081679,0.44127145,"ATP hydrolysis activity"),
                     c("GO:0047701","beta-L-arabinosidase activity",4.878485139890339E-05,1,0.9401401842035029,0.41027052,"ATP hydrolysis activity"),
                     c("GO:0047782","coniferin beta-glucosidase activity",2.66099189448564E-05,1,0.9380013667085084,0.44258599,"ATP hydrolysis activity"),
                     c("GO:0047840","dCTP diphosphatase activity",0.00782775115627859,1,0.9364963382573893,0.42129781,"ATP hydrolysis activity"),
                     c("GO:0047884","FAD diphosphatase activity",0.00013970207446049608,1,0.9456770136428667,0.31970647,"ATP hydrolysis activity"),
                     c("GO:0050224","prunasin beta-glucosidase activity",1.55224527178329E-05,1,0.9393837076440363,0.43456671,"ATP hydrolysis activity"),
                     c("GO:0080079","cellobiose glucosidase activity",7.982975683456919E-05,1,0.9349762316979057,0.45988357,"ATP hydrolysis activity"),
                     c("GO:0080083","beta-gentiobiose beta-glucosidase activity",7.76122635891645E-05,1,0.9350575361915731,0.42073911,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9402006720937442,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9387706512195393,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140326","ATPase-coupled intramembrane lipid transporter activity",0.056601515088954966,1,0.9035448513096077,0.47545686,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,2,0.9359460693424815,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,8,0.981632546991765,0.04901613,"carbohydrate binding"),
                     c("GO:0030674","protein-macromolecule adaptor activity",1.2206258094127533,3,0.9818006831662439,-0,"protein-macromolecule adaptor activity"),
                     c("GO:0032556","pyrimidine deoxyribonucleotide binding",9.09172230615927E-05,1,0.9763354793621403,0.07629645,"pyrimidine deoxyribonucleotide binding"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9477041521979613,0.03310522,"DNA demethylase activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,4,0.98155668749802,-0,"signaling receptor activity"),
                     c("GO:0009885","transmembrane histidine kinase cytokinin receptor activity",3.769738517187989E-05,1,0.8767415148035875,0.38178873,"signaling receptor activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,15,0.9437887856978878,0.04095538,"protein homodimerization activity"),
                     c("GO:0000149","SNARE binding",0.26943873427614345,2,0.942792002641888,0.40354822,"protein homodimerization activity"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,11,0.942883750636101,0.40281405,"protein homodimerization activity"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,13,0.9354328706465737,0.38709915,"protein homodimerization activity"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,2,0.9493128327061413,0.3218275,"protein homodimerization activity"),
                     c("GO:0019901","protein kinase binding",0.4024661540679714,5,0.9328854808754286,0.41867396,"protein homodimerization activity"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.9441555178342247,0.39266904,"protein homodimerization activity"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,3,0.9463011211239948,0.34139999,"protein homodimerization activity"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9481717059123352,0.32923903,"protein homodimerization activity"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,3,0.9449938446832636,0.34990764,"protein homodimerization activity"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9484686286287121,0.32731003,"protein homodimerization activity"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,2,0.9624225318942571,0.2370668,"protein homodimerization activity"),
                     c("GO:0032050","clathrin heavy chain binding",0.02621520514717436,2,0.9500809474712852,0.30626295,"protein homodimerization activity"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,2,0.9434110870812179,0.39860027,"protein homodimerization activity"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9521192316858261,0.48770765,"protein homodimerization activity"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,2,0.9455109213449056,0.33872441,"protein homodimerization activity"),
                     c("GO:0042393","histone binding",0.41579772345934446,6,0.9407521339487543,0.41995229,"protein homodimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,9,0.9398351958719575,0.42737733,"protein homodimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,4,0.9438266869026394,0.3952865,"protein homodimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,6,0.9490761271825151,0.32336445,"protein homodimerization activity"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,9,0.935226970665565,0.4679269,"protein homodimerization activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,11,0.9390050009843673,0.43412799,"protein homodimerization activity"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,2,0.9422377809164217,0.40798975,"protein homodimerization activity"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,3,0.9529726326090413,0.28789489,"protein homodimerization activity"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,2,0.9490154651867558,0.32375837,"protein homodimerization activity"),
                     c("GO:1990935","splicing factor binding",0.0007095978385295039,1,0.961086045252152,0.24567361,"protein homodimerization activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.980052909297772,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9969861288752431,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9808115239171371,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9862914642874989,0.49249984,"ER retention sequence binding"),
                     c("GO:1904408","melatonin binding",0.00021066185831344647,1,0.9793104873168712,0.43339648,"ER retention sequence binding"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,2,0.9623917017344634,0.02253254,"oxalate-CoA ligase activity"),
                     c("GO:0003972","RNA ligase (ATP) activity",0.012360307349885798,1,0.9369948658263713,0.4132865,"oxalate-CoA ligase activity"),
                     c("GO:0003989","acetyl-CoA carboxylase activity",0.08698117255099935,1,0.9486761697043934,0.40308719,"oxalate-CoA ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.955587030443094,0.28284902,"oxalate-CoA ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.963824697656016,0.2351635,"oxalate-CoA ligase activity"),
                     c("GO:0106286","(E)-caffeate-CoA ligase activity",2.4392425699451695E-05,1,0.9642654643218527,0.39917988,"oxalate-CoA ligase activity"),
                     c("GO:0050373","UDP-arabinose 4-epimerase activity",8.20472500799739E-05,1,0.9612338253353215,0.02256811,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0016871","cycloartenol synthase activity",0.0006563780006397912,1,0.9581354922059703,0.28013495,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:2001070","starch binding",0.024815966909324,1,0.9847875682274669,0.03482451,"starch binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file

pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No Zm DE, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_MF_Table.tsv", sep = "\t")
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

