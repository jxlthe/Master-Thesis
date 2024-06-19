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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,17,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,4,0.9649126240562929,0.04661592,"chromatin binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,28,0.9429737859643986,-0,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9573796320794878,0.42104634,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,46,0.9368707036252987,0.48153711,"mRNA binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9625789804256175,0.19519762,"mRNA binding"),
                     c("GO:0003714","transcription corepressor activity",0.25663049329068593,2,0.9768100542757125,0.29590577,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9512767475415395,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9511270322037229,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,8,0.9313035085598319,0.24940265,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9567183725746465,0.428381,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,6,0.9647944226121699,0.13060477,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9583097999647859,0.41077037,"mRNA binding"),
                     c("GO:0035198","miRNA binding",0.029599099839661934,1,0.9577632770366261,0.41680219,"mRNA binding"),
                     c("GO:0042134","rRNA primary transcript binding",0.01629192287398833,1,0.9594766474822041,0.39794793,"mRNA binding"),
                     c("GO:0043047","single-stranded telomeric DNA binding",0.03936272259917883,1,0.9556931252544846,0.33632722,"mRNA binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,12,0.9395673469359026,0.33073906,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,3,0.996953857194905,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,22,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,29,0.7966955598053842,0,"protein kinase activity"),
                     c("GO:0000234","phosphoethanolamine N-methyltransferase activity",0.0010222643861315665,1,0.9101961182379926,0.45971198,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,3,0.9041483348680777,0.34856872,"protein kinase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.917911386475061,0.20796201,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8943318146483252,0.45866285,"protein kinase activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,4,0.8580116616013308,0.41115081,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8945077899471767,0.48555173,"protein kinase activity"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.89424154920464,0.46089842,"protein kinase activity"),
                     c("GO:0004619","phosphoglycerate mutase activity",0.0367083831844294,1,0.9534616077201333,0.4684783,"protein kinase activity"),
                     c("GO:0004659","prenyltransferase activity",0.25118876486646274,1,0.8953678773240584,0.26852283,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,5,0.8356403557433328,0.41376722,"protein kinase activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,7,0.8225227726558936,0.47753726,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,7,0.860152910741148,0.35677507,"protein kinase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,2,0.86787225776421,0.30425857,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8606159888939169,0.35123485,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,4,0.8785653070669813,0.30161168,"protein kinase activity"),
                     c("GO:0008728","GTP diphosphokinase activity",0.028889502001132432,1,0.8974818261195343,0.40329157,"protein kinase activity"),
                     c("GO:0008798","beta-aspartyl-peptidase activity",0.023041972313000234,1,0.880325975233364,0.42550609,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,2,0.8540627345701525,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,7,0.8760345470755091,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,19,0.8809090870526365,0.35580262,"protein kinase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.908880760441058,0.49223587,"protein kinase activity"),
                     c("GO:0033855","nicotianamine aminotransferase activity",6.6524797362141E-06,1,0.9354078780944137,0.45711636,"protein kinase activity"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8810299774872884,0.49253483,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.8991384315701102,0.49478365,"protein kinase activity"),
                     c("GO:0047159","1-alkenylglycerophosphocholine O-acyltransferase activity",0.00010865716902483029,1,0.9171381390504977,0.42478125,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,2,0.8896139113429542,0.48733868,"protein kinase activity"),
                     c("GO:0050200","plasmalogen synthase activity",3.769738517187989E-05,1,0.9251136312087593,0.26985858,"protein kinase activity"),
                     c("GO:0052636","arabinosyltransferase activity",0.0034526369830951177,1,0.9069933912781002,0.47557181,"protein kinase activity"),
                     c("GO:0052667","phosphomethylethanolamine N-methyltransferase activity",0.00011530964876104439,1,0.9190552469379513,0.40189746,"protein kinase activity"),
                     c("GO:0052923","all-trans-nonaprenyl-diphosphate synthase (geranyl-diphosphate specific) activity",0.00042797619636310703,1,0.9277004896234339,0.49777323,"protein kinase activity"),
                     c("GO:0071618","lysophosphatidylethanolamine acyltransferase activity",0.003938268003838747,1,0.9042916666933052,0.35248487,"protein kinase activity"),
                     c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.9072941853896729,0.31861303,"protein kinase activity"),
                     c("GO:0106262","1-acylglycerophosphoethanolamine O-acyltransferase activity",0.0003703213719825849,1,0.9122456680807022,0.44573481,"protein kinase activity"),
                     c("GO:0106263","1-acylglycerophosphoserine O-acyltransferase activity",0.0004324111828539164,1,0.9115852939314706,0.15964368,"protein kinase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.892179336596866,0.45431535,"protein kinase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.891768907063418,0.39355689,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,3,0.9777415698833525,-0,"GTPase activator activity"),
                     c("GO:0005515","protein binding",8.610051728351934,88,0.9721879680375627,0.07766851,"protein binding"),
                     c("GO:0097159","organic cyclic compound binding",39.22969739688024,1,0.9650038045512186,0.13134756,"protein binding"),
                     c("GO:0005543","phospholipid binding",0.6999096105403309,3,0.9638905596726693,0.04714527,"phospholipid binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,15,0.9479648068525763,0.05738816,"zinc ion binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,88,0.9273531578387632,0.48160605,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,4,0.95568184563194,0.3830934,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,2,0.9688147268833677,0.21260669,"zinc ion binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,34,0.9506227862818133,0.20275684,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,6,0.9569447823796046,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,13,0.9382587711667798,0.19635374,"zinc ion binding"),
                     c("GO:0010181","FMN binding",0.4661059927178409,1,0.9480095235583085,0.34483719,"zinc ion binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9499164392995864,0.29697106,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,3,0.9439997884366288,0.31800282,"zinc ion binding"),
                     c("GO:0030955","potassium ion binding",0.05969270067304912,1,0.9674813988533001,0.26193125,"zinc ion binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9785687783374003,0.2693746,"zinc ion binding"),
                     c("GO:0043169","cation binding",18.365868911645002,2,0.9484806002806669,0.2885161,"zinc ion binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.9738119879027932,0.11877911,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,5,0.9428110545651036,0.33990845,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,3,0.9486566541871637,0.35257211,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,6,0.9538370946043542,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9546900860383924,0.28523818,"zinc ion binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,3,0.9529625669172973,0.29755636,"zinc ion binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,3,0.9483293052956454,0.30272555,"zinc ion binding"),
                     c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9653524532195201,0.25647876,"zinc ion binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,4,0.9777629687338549,0.05834827,"lipid binding"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,3,0.9931998046540986,-0,"photoreceptor activity"),
                     c("GO:0009882","blue light photoreceptor activity",0.009566265860675875,1,0.9807881196363898,0.45612612,"photoreceptor activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,2,0.9872072771190255,0.02769872,"auxin binding"),
                     c("GO:0015293","symporter activity",0.6379683717164414,3,0.923821166373216,-0,"symporter activity"),
                     c("GO:0005254","chloride channel activity",0.1887840699542837,2,0.9242613427716713,0.3683835,"symporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9161173207841773,0.46407837,"symporter activity"),
                     c("GO:0015165","pyrimidine nucleotide-sugar transmembrane transporter activity",0.046425438585792796,1,0.9438418196992429,0.32778933,"symporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.9445611913696884,0.47124298,"symporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9386913249698924,0.3537519,"symporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.96521610515402,0.20748322,"symporter activity"),
                     c("GO:0042910","xenobiotic transmembrane transporter activity",0.2653474592383718,1,0.9394495999257795,0.3797985,"symporter activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9246894291754755,0.43531693,"symporter activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,5,0.9032465887001682,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,4,0.9045234927405373,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004601","peroxidase activity",0.4792135952914281,1,0.9010181384485876,0.40806599,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0010242","oxygen evolving activity",0.0031710153409287207,1,0.9368825710449135,0.26712027,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9251361833992277,0.46927032,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,4,0.8900399704432891,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9055060857322698,0.41247223,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.906835546147857,0.44173969,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016630","protochlorophyllide reductase activity",0.0022374506846133423,1,0.9324827698857283,0.26085905,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors",0.2806015952735107,1,0.9156698087779124,0.3863248,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016639","oxidoreductase activity, acting on the CH-NH2 group of donors, NAD or NADP as acceptor",0.09697097962154752,1,0.9173384851662528,0.34937139,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9178622094095388,0.48602252,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0019139","cytokinin dehydrogenase activity",0.0057233500663895305,1,0.9347256405517566,0.27843728,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0019797","procollagen-proline 3-dioxygenase activity",0.004459378916508851,1,0.8891486089271369,0.3991339,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0048529","magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase activity",0.002614424536332141,1,0.9342795792957843,0.26361814,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9090203772696658,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,95,0.9463340990000323,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,33,0.9501098369352733,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,61,0.9458248194786081,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,12,0.9558474527626463,0.05997119,"lyase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,3,0.9312451570699782,0.04651897,"carboxy-lyase activity"),
                     c("GO:0009978","allene oxide synthase activity",9.09172230615927E-05,1,0.9586286303151453,0.33685723,"carboxy-lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9507275672268762,0.45277577,"carboxy-lyase activity"),
                     c("GO:0018812","3-hydroxyacyl-CoA dehydratase activity",0.012180690397008016,1,0.9470672300218994,0.4947816,"carboxy-lyase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9445868285842907,0.48814337,"carboxy-lyase activity"),
                     c("GO:0106099","2-keto-3-deoxy-L-rhamnonate aldolase activity",0.00019513940559561356,1,0.9569514944109843,0.48677249,"carboxy-lyase activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9644978415543609,0.28175041,"carboxy-lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,10,0.9575820482613279,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,11,0.9563272943543439,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,20,0.86154668418708,0.05972228,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.9066338370600061,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9303409957872354,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,4,0.872994834358962,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8863418508493948,0.38812549,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,4,0.8980043197460424,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.898268516476489,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.8995965306598724,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0033907","beta-D-fucosidase activity",0.0008093850345727155,1,0.9213380368262046,0.48296795,"ATP hydrolysis activity"),
                     c("GO:0047631","ADP-ribose diphosphatase activity",0.02218158493378321,2,0.9236516297613782,0.45903615,"ATP hydrolysis activity"),
                     c("GO:0047668","amygdalin beta-glucosidase activity",2.4392425699451695E-05,1,0.9297547974380621,0.44127145,"ATP hydrolysis activity"),
                     c("GO:0047701","beta-L-arabinosidase activity",4.878485139890339E-05,1,0.9311589858280174,0.41027052,"ATP hydrolysis activity"),
                     c("GO:0047884","FAD diphosphatase activity",0.00013970207446049608,1,0.9378284698963498,0.31970647,"ATP hydrolysis activity"),
                     c("GO:0050224","prunasin beta-glucosidase activity",1.55224527178329E-05,1,0.9310740794270776,0.43456671,"ATP hydrolysis activity"),
                     c("GO:0080079","cellobiose glucosidase activity",7.982975683456919E-05,1,0.9260382380344921,0.45988357,"ATP hydrolysis activity"),
                     c("GO:0080083","beta-gentiobiose beta-glucosidase activity",7.76122635891645E-05,1,0.9261311566013791,0.42073911,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9289246491016637,0.17583997,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9258252545029961,0.12676211,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9432055679184503,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,4,0.9783323314588592,0.04901613,"carbohydrate binding"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9477681450997801,0.03310522,"DNA demethylase activity"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9764184523015997,0.06283846,"protein-containing complex binding"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9955908316513356,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8901095611423163,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9775829143049821,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9839743343245027,0.49249984,"ER retention sequence binding"),
                     c("GO:0050373","UDP-arabinose 4-epimerase activity",8.20472500799739E-05,1,0.9657380793298711,0.02256811,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9591450397867443,0.44049754,"UDP-arabinose 4-epimerase activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,9,0.9291812586835051,0.0463004,"unfolded protein binding"),
                     c("GO:0000149","SNARE binding",0.26943873427614345,1,0.9335923777889343,0.40559978,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,4,0.93369926377063,0.40485813,"unfolded protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,5,0.9263343769046012,0.43412799,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9411896537168161,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,3,0.9240850272295623,0.42444044,"unfolded protein binding"),
                     c("GO:0019904","protein domain specific binding",0.1980687141727932,1,0.9351809171549877,0.39461121,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,2,0.9376807483820682,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,3,0.9398601683828287,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,2,0.9361576356371167,0.38789043,"unfolded protein binding"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9402061067010586,0.36031388,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9564547862130692,0.25391209,"unfolded protein binding"),
                     c("GO:0032050","clathrin heavy chain binding",0.02621520514717436,1,0.9439845079074578,0.33497261,"unfolded protein binding"),
                     c("GO:0032182","ubiquitin-like protein binding",0.2347371824788053,1,0.9343136186262746,0.4006017,"unfolded protein binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9439935311809833,0.31955157,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,3,0.931216122123806,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,6,0.9301481221369102,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,8,0.9359980831746438,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,3,0.934797809663528,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,5,0.9409138808360964,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,3,0.9247830160760373,0.4679269,"unfolded protein binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.9329467250005615,0.41008684,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.953210853019769,0.26331899,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,2,0.9449917264125887,0.48770765,"unfolded protein binding"),
                     c("GO:0097602","cullin family protein binding",0.05781226640094593,1,0.9408432064229358,0.35601456,"unfolded protein binding"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.963076983842748,0.02238092,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0003972","RNA ligase (ATP) activity",0.012360307349885798,1,0.9363054715674739,0.36034906,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9546792169747752,0.28284902,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,1,0.9633250098724788,0.2351635,"jasmonoyl-L-amino acid ligase activity"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:2001070","starch binding",0.024815966909324,1,0.9816736548671297,0.03482451,"starch binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3DEG_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No DEG in Zm, All C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG_reduced05_MF_Table.tsv", sep = "\t")
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