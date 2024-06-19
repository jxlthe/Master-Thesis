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
revigo.data <- rbind(c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,7,0.9419628176739984,0.066121,"transcription cis-regulatory region binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,5,0.9499660975376847,0.22593798,"transcription cis-regulatory region binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,24,0.94105616324916,0.45223656,"transcription cis-regulatory region binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9665228957965674,0.31120183,"transcription cis-regulatory region binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,16,0.9864138522601702,0.40628023,"transcription cis-regulatory region binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,3,0.9562940102124945,0.31664803,"transcription cis-regulatory region binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9430915716595631,0.27764537,"transcription cis-regulatory region binding"),
                     c("GO:0019237","centromeric DNA binding",0.004625690909914204,1,0.9667754581017608,0.49059442,"transcription cis-regulatory region binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9651443006653624,0.14801235,"transcription cis-regulatory region binding"),
                     c("GO:0003674","molecular_function",100,10,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,1,0.9738454051056479,0.0533165,"chromatin binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,11,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,23,0.7744996333649337,0,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,2,0.8845724263170435,0.34856872,"protein kinase activity"),
                     c("GO:0003975","UDP-N-acetylglucosamine-dolichyl-phosphate N-acetylglucosaminephosphotransferase activity",0.005718915079898721,1,0.8976980242425416,0.35348317,"protein kinase activity"),
                     c("GO:0004106","chorismate mutase activity",0.044957458057334886,1,0.9475028729130611,0.47467998,"protein kinase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.8991077512603162,0.25688843,"protein kinase activity"),
                     c("GO:0004709","MAP kinase kinase kinase activity",0.022518643907084728,1,0.8337664879965059,0.49507895,"protein kinase activity"),
                     c("GO:0004781","sulfate adenylyltransferase (ATP) activity",0.04230755362907627,1,0.8822900625187642,0.41713542,"protein kinase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,1,0.8725292313409576,0.35677507,"protein kinase activity"),
                     c("GO:0008236","serine-type peptidase activity",1.2663216927208079,2,0.8316623479874524,0.49481639,"protein kinase activity"),
                     c("GO:0008373","sialyltransferase activity",0.03638684666384572,1,0.8914279426635238,0.43994578,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,2,0.8750557656849572,0.30161168,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,1,0.8443049032880112,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,6,0.8672270504187616,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,12,0.8727657548687641,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,3,0.8851166610554807,0.21201958,"protein kinase activity"),
                     c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.04143767537877,1,0.8564528272494051,0.43304972,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.8910492855322953,0.37780808,"protein kinase activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9559696927130285,0.38075906,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,1,0.8730438925023081,0.48733868,"protein kinase activity"),
                     c("GO:0050734","hydroxycinnamoyltransferase activity",0.0034969868480032116,1,0.8981091164002102,0.18416899,"protein kinase activity"),
                     c("GO:0051753","mannan synthase activity",0.005627997856837128,1,0.8996393290260886,0.47286456,"protein kinase activity"),
                     c("GO:0061630","ubiquitin protein ligase activity",0.6665452071699717,5,0.8056777694326984,0.45681379,"protein kinase activity"),
                     c("GO:0005096","GTPase activator activity",0.43493469016718705,1,0.9952182751637113,-0,"GTPase activator activity"),
                     c("GO:0005515","protein binding",8.610051728351934,46,0.9716427288876132,0.07766851,"protein binding"),
                     c("GO:0008017","microtubule binding",0.5569810834077709,6,0.9166999090790998,0.05255315,"microtubule binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,5,0.9309748412630602,0.40281405,"microtubule binding"),
                     c("GO:0019901","protein kinase binding",0.4024661540679714,2,0.92603792544991,0.41867396,"microtubule binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9335958346372112,0.3860137,"microtubule binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9443981783814701,0.31827679,"microtubule binding"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,1,0.9357296086277153,0.37244819,"microtubule binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9283267326886423,0.41995229,"microtubule binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,3,0.9271877990120283,0.42737733,"microtubule binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,5,0.9307049517133249,0.38709915,"microtubule binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9321461934406067,0.3952865,"microtubule binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9386618294523559,0.35396093,"microtubule binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,5,0.9214704781073605,0.4651985,"microtubule binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9261568333118312,0.4679269,"microtubule binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.930172341394325,0.40798975,"microtubule binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,11,0.9405145586980309,-0,"zinc ion binding"),
                     c("GO:0000035","acyl binding",0.03278342014006308,1,0.9712133733070926,0.11604884,"zinc ion binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,1,0.9497634817826637,0.3830934,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,3,0.9512389854106326,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,2,0.9321986694039569,0.35257211,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,3,0.9451137395441458,0.13873458,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9527322215068251,0.20283352,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,1,0.9381654281197742,0.31639809,"zinc ion binding"),
                     c("GO:0043295","glutathione binding",0.01852715606535627,1,0.9635077978033577,0.17458707,"zinc ion binding"),
                     c("GO:0043325","phosphatidylinositol-3,4-bisphosphate binding",0.013358179310317913,1,0.964338608907905,0.17101403,"zinc ion binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,1,0.9341467725870477,0.31135007,"zinc ion binding"),
                     c("GO:0050661","NADP binding",0.7038257036117155,1,0.944874429694901,0.3462071,"zinc ion binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,2,0.9438052293524416,0.16293278,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9466715608305454,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.9479750821545041,0.28982977,"zinc ion binding"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,1,0.9821431196277406,0.02993406,"auxin binding"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,4,0.9090083182364842,-0,"antiporter activity"),
                     c("GO:0005338","nucleotide-sugar transmembrane transporter activity",0.0598301852542642,2,0.9193689970476834,0.33273016,"antiporter activity"),
                     c("GO:0008553","P-type proton-exporting transporter activity",0.01860920331543624,1,0.8966669279706063,0.47571857,"antiporter activity"),
                     c("GO:0015276","ligand-gated monoatomic ion channel activity",0.3826107195486177,2,0.8847836857049981,0.39050958,"antiporter activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9238202754911319,0.35182404,"antiporter activity"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.9525787792498923,0.20645071,"antiporter activity"),
                     c("GO:0035673","oligopeptide transmembrane transporter activity",0.2025591379947377,1,0.9109513286038611,0.3876045,"antiporter activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9055820684265667,0.43240123,"antiporter activity"),
                     c("GO:0051980","iron-nicotianamine transmembrane transporter activity",0.0019447415762199217,1,0.9263030778312501,0.28474685,"antiporter activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,3,0.8796902445130788,0.0531541,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004324","ferredoxin-NADP+ reductase activity",0.018482806200448173,1,0.9281374100439902,0.30749203,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9024394673204892,0.46480825,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9268905757653264,0.31440091,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,1,0.9058668394137597,0.41891471,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9011585081526731,0.47359913,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.906941954769149,0.49506653,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0045543","gibberellin 2-beta-dioxygenase activity",0.006071496505918068,1,0.9167721931620054,0.48062226,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,2,0.9069600712008533,0.43460917,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102111","gibberellin A20,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102652","gibberellin A9,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.39925386,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102772","sphingolipid C4-monooxygenase activity",4.4349864908093996E-05,1,0.9380050977762479,0.28012766,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0102924","gibberellin A44,2-oxoglutarate:oxygen oxidoreductase activity",4.2132371662689295E-05,1,0.9342040754960839,0.20747889,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,55,0.9437836522476242,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,18,0.9479842952457235,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,25,0.9432159328868454,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,7,0.9543080049624433,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,7,0.9562022310468056,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,5,0.9548328537234818,0.05784577,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,5,0.871929271156829,0.05997119,"ATP hydrolysis activity"),
                     c("GO:0004301","epoxide hydrolase activity",0.06371079843372243,1,0.924843077902784,0.22109293,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,3,0.8864327502617124,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8607554481442179,0.454689,"ATP hydrolysis activity"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,2,0.86706935078231,0.32020834,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,4,0.9048402323081485,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.9050971095270066,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.9058211895932066,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0140326","ATPase-coupled intramembrane lipid transporter activity",0.056601515088954966,1,0.9049973218123009,0.47545686,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,4,0.9781376477450502,0.05647972,"carbohydrate binding"),
                     c("GO:0030527","structural constituent of chromatin",0.19884040182219404,1,1,-0,"structural constituent of chromatin"),
                     c("GO:0030674","protein-macromolecule adaptor activity",1.2206258094127533,1,1,-0,"protein-macromolecule adaptor activity"),
                     c("GO:0038023","signaling receptor activity",2.0951718830016857,4,0.9740485483191771,-0,"signaling receptor activity"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9966774087712642,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8780311780345175,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0050203","oxalate-CoA ligase activity",7.982975683456919E-05,1,0.9719635630397723,0.02253254,"oxalate-CoA ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9645846351874537,0.44358297,"oxalate-CoA ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9539150239545147,0.33857975,"oxalate-CoA ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9712579065613889,0.4447539,"oxalate-CoA ligase activity"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9331203890479748,0.0567159,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9455587216231844,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9688129977521718,0.01868194,"pterocarpan synthase activity"),
                     c("GO:0004383","guanylate cyclase activity",0.041591303310810554,1,0.9425222402969877,0.49730148,"pterocarpan synthase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,1,0.9387243910696051,0.48814337,"pterocarpan synthase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9507451098580505,0.32867735,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches



topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, No DEG in Zm, At Least One Of Each C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_MF_Table.tsv", sep = "\t")
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

