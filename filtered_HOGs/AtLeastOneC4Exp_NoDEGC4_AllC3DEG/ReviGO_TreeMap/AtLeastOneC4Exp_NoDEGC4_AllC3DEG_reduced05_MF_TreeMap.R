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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,13,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,4,0.9611351788227029,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,23,0.9397627743394634,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9552139747970011,0.42104634,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,46,0.9304571034993059,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9582993692262955,0.19519762,"mRNA binding"),
                     c("GO:0003713","transcription coactivator activity",0.24395086691346182,2,0.9646157548972838,0.29461398,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9486811596841636,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.948376383859577,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,5,0.9273718332035783,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9545072611422076,0.428381,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,5,0.961702143109721,0.13060477,"mRNA binding"),
                     c("GO:0032934","sterol binding",0.11902838493358808,1,0.9602753423560685,0.10144903,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9500987454666109,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,11,0.9323975125272463,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,4,0.9963047271514212,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,20,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,28,0.8012889512105662,0,"protein kinase activity"),
                     c("GO:0003878","ATP citrate synthase activity",0.010080724293609766,1,0.90111202967131,0.32478201,"protein kinase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9196246748807231,0.20796201,"protein kinase activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,4,0.8692663224412006,0.41115081,"protein kinase activity"),
                     c("GO:0004659","prenyltransferase activity",0.25118876486646274,1,0.8956877844769857,0.26852283,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8507841213030327,0.49507895,"protein kinase activity"),
                     c("GO:0004721","phosphoprotein phosphatase activity",0.6445676316147657,3,0.827746251708972,0.4549882,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,8,0.8195039001097916,0.47753726,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,4,0.8652158530088412,0.35677507,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8946693927575179,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8985318456786945,0.40329157,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.854758657386247,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,5,0.877182289485356,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,9,0.8821684863712335,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,2,0.9012515924548455,0.21201958,"protein kinase activity"),
                     c("GO:0019786","protein-phosphatidylethanolamide deconjugating activity",0.010599617713034465,1,0.903575015657121,0.30543886,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.9146798450832966,0.49223587,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8797085259364104,0.49253483,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9147899061265138,0.30679413,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9229099827044737,0.29376839,"protein kinase activity"),
                     c("GO:0051741","2-methyl-6-phytyl-1,4-benzoquinone methyltransferase activity",0.0016742074002805483,1,0.9169995905573177,0.47516484,"protein kinase activity"),
                     c("GO:0051742","2-methyl-6-solanyl-1,4-benzoquinone methyltransferase activity",1.1087466227023499E-05,1,0.9344914229851815,0.12943404,"protein kinase activity"),
                     c("GO:0052630","UDP-N-acetylgalactosamine diphosphorylase activity",4.4349864908093996E-05,1,0.9201347359742593,0.25791543,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9138386017595478,0.3846796,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9282637979198199,0.49777323,"protein kinase activity"),
                     c("GO:0070204","2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylic-acid synthase activity",0.007958583257757468,1,0.9201019614889142,0.19601592,"protein kinase activity"),
                     c("GO:0070569","uridylyltransferase activity",0.1030402586342202,1,0.8813259139456104,0.45345646,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.8989056130651917,0.18579171,"protein kinase activity"),
                     c("GO:0102550","2-methyl-6-geranylgeranyl-1,4-benzoquinol methyltransferase activity",7.76122635891645E-05,1,0.9297380499560842,0.39293218,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.9000121762031924,0.39355689,"protein kinase activity"),
                     c("GO:1990714","hydroxyproline O-galactosyltransferase activity",0.0006475080276581724,1,0.9193389090433758,0.42279362,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,73,0.9709642464979593,0.07766851,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,13,0.9431261818443263,0.05738816,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,65,0.9198304274296102,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,3,0.951980980576073,0.3830934,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,26,0.9461399564400542,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,6,0.9533692929286738,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,7,0.9363475014305243,0.19635374,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9451857612930475,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9468571600713641,0.29697106,"zinc ion binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,2,0.9474421564492345,0.24057124,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,1,0.9450352440119275,0.2885161,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9722946711923405,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.9385909232649782,0.33990845,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,3,0.9443192341157636,0.34531063,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,2,0.9432731516283231,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,5,0.9495443006076708,0.18168089,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,2,0.9493950186232514,0.29755636,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,5,0.9768257194973491,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,1,0.9942543233776147,-0,"photoreceptor activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,6,0.8968261985078174,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,6,0.8982143348026331,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,2,0.8907659990398528,0.40806599,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9189858187898816,0.31075808,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,2,0.890425659736503,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.8961949187017461,0.41247223,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9007265043525489,0.44173969,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,1,0.9103083375357702,0.3863248,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016639","oxidoreductase activity, acting on the CH-NH2 group of donors, NAD or NADP as acceptor",0.09697097962154752,1,0.9115284848589319,0.34937139,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016711","flavonoid 3'-monooxygenase activity",8.20472500799739E-05,1,0.937270470200862,0.21343278,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9078332520701534,0.4215727,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0045486","naringenin 3-dioxygenase activity",0.00033705897330151436,1,0.9297661353959318,0.38276908,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0050664","oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor",0.059595130970251306,1,0.915726359435088,0.3347022,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,3,0.9030988232273004,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:1990136","linoleate 9S-lipoxygenase activity",0.00026609918944856393,1,0.9349402817373502,0.22819902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,74,0.9460377115566114,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,26,0.9499681698786101,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,38,0.9455067431342304,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,8,0.9559105526427042,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,2,0.9419042058348917,0.04651897,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9569626007807175,0.45277577,"carboxy-lyase activity"),
                     c("GO:0043748","O-succinylbenzoate synthase activity",0.006195676127660732,1,0.9555538612351807,0.43317182,"carboxy-lyase activity"),
                     c("GO:0070205","2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase activity",0.0028051289554369453,1,0.9579522593167777,0.41110877,"carboxy-lyase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9587081143100623,0.03072632,"magnesium chelatase activity"),
                     c("GO:0004088","carbamoyl-phosphate synthase (glutamine-hydrolyzing) activity",0.05676560958911491,1,0.9534244275930002,0.49740294,"magnesium chelatase activity"),
                     c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9467497923961009,0.43449047,"magnesium chelatase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,8,0.9576985621983991,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,7,0.9564055825885092,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,14,0.8623500997419943,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004427","inorganic diphosphate phosphatase activity",0.053847388478162325,1,0.9216085953602194,0.49693897,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,3,0.8986105214288743,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0004636","phosphoribosyl-ATP diphosphatase activity",0.02123693281124081,1,0.9263463719477826,0.45732449,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.903912906933785,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9041668247449802,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9019281482572177,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0017057","6-phosphogluconolactonase activity",0.03804774910465384,1,0.8969452217644882,0.49474271,"ATP hydrolysis activity"),
                     c("GO:0035529","NADH pyrophosphatase activity",0.024569825159084076,1,0.922203574775657,0.46310876,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,1,0.9261377070261394,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047714","galactolipase activity",0.0037120836928074678,1,0.9070072417778053,0.35781567,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9288473613191554,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9262021708941423,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9440469655202162,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,3,0.9774252496409356,0.04901613,"carbohydrate binding"),
                     c("GO:0030414","peptidase inhibitor activity",0.25608720744556174,2,0.9711923644749612,-0,"peptidase inhibitor activity"),
                     c("GO:0005094","Rho GDP-dissociation inhibitor activity",0.010056331867910313,1,0.9778961338048868,0.49485349,"peptidase inhibitor activity"),
                     c("GO:0043879","glycolate transmembrane transporter activity",0.00030379657462044384,1,0.9484034003314635,-0,"glycolate transmembrane transporter activity"),
                     c("GO:0008526","phosphatidylinositol transfer activity",0.022729305765398174,1,0.9456378334776795,0.2866085,"glycolate transmembrane transporter activity"),
                     c("GO:0009678","diphosphate hydrolysis-driven proton transmembrane transporter activity",0.01596151638042303,1,0.9316782497316297,0.42292952,"glycolate transmembrane transporter activity"),
                     c("GO:0015116","sulfate transmembrane transporter activity",0.08090967604508129,1,0.9275559067027012,0.28203837,"glycolate transmembrane transporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9161894325329503,0.35239237,"glycolate transmembrane transporter activity"),
                     c("GO:0015165","pyrimidine nucleotide-sugar transmembrane transporter activity",0.046425438585792796,1,0.9464063856271305,0.22592317,"glycolate transmembrane transporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9398362287280138,0.47124298,"glycolate transmembrane transporter activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9754100624782831,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9946516324784443,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8859491001696689,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9590108201390692,0.02641673,"carotenoid isomerase activity"),
                     c("GO:0004165","delta(3)-delta(2)-enoyl-CoA isomerase activity",0.021527424426388823,1,0.9516491001271186,0.49767449,"carotenoid isomerase activity"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9574994292018759,0.44355864,"carotenoid isomerase activity"),
                     c("GO:0009982","pseudouridine synthase activity",0.21083704027983347,1,0.9417551056910062,0.41191472,"carotenoid isomerase activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,1,0.9756286294269892,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9825458914650622,0.49249984,"ER retention sequence binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,5,0.925699410242818,0.0463004,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,3,0.9304151030336154,0.40485813,"unfolded protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,4,0.922196046603978,0.43412799,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9382433373927986,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,4,0.920971180283221,0.42444044,"unfolded protein binding"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.931962704354966,0.39461121,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,1,0.9345748502102412,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9338519918891887,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,2,0.9329831565563115,0.38789043,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9542223152176578,0.25391209,"unfolded protein binding"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,1,0.9310567427136419,0.4006017,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,2,0.9278226197158785,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,2,0.926708098819392,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,5,0.9328164464472001,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,2,0.9315624979994015,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,3,0.9379549488732737,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9211149610924724,0.4679269,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9505342250815753,0.26331899,"unfolded protein binding"),
                     c("GO:0070678","preprotein binding",3.326239868107049E-05,1,0.9592829924781565,0.22354599,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.9431309481722314,0.31312228,"unfolded protein binding"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,1,0.9378810432751246,0.35601456,"unfolded protein binding"),
                     c("GO:0140098","catalytic activity, acting on RNA",3.705433430564499,1,0.9134036898318094,0.06018061,"catalytic activity, acting on RNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9403686812373283,0.42852185,"catalytic activity, acting on RNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AtLeastOneC4Exp_NoDEGC4_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


hogs <- read.csv("..\\AtLeastOneC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("At Least One C4 Expanded, No DEG in C4, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
