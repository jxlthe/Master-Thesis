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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,4,1,-0,"molecular_function"),
                     c("GO:0003729","mRNA binding",1.136762432364793,22,0.919189989329702,0,"mRNA binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,32,0.9055575952249222,0.48160605,"mRNA binding"),
                     c("GO:0000340","RNA 7-methylguanosine cap binding",0.033608327627353635,1,0.9396185885442078,0.42104634,"mRNA binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,3,0.9162007849864184,0.31664803,"mRNA binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,27,0.9332895993562562,0.20275684,"mRNA binding"),
                     c("GO:0003713","transcription coactivator activity",0.24395086691346182,2,0.9753831979375028,0.29461398,"mRNA binding"),
                     c("GO:0003723","RNA binding",6.099813894661886,31,0.9194448330560689,0.40176997,"mRNA binding"),
                     c("GO:0003727","single-stranded RNA binding",0.19311705175580451,1,0.9309513931976233,0.48968021,"mRNA binding"),
                     c("GO:0003730","mRNA 3'-UTR binding",0.15417121788676175,1,0.9302157238216323,0.47960958,"mRNA binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,6,0.9123267158409292,0.24940265,"mRNA binding"),
                     c("GO:0003919","FMN adenylyltransferase activity",0.02766322823642363,1,0.9008570455896044,0.2955927,"mRNA binding"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8946958237802933,0.3505363,"mRNA binding"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8978500615159202,0.33108637,"mRNA binding"),
                     c("GO:0004349","glutamate 5-kinase activity",0.026552264120475875,1,0.8971382099233871,0.3394079,"mRNA binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,2,0.9460972935657948,0.3693843,"mRNA binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,7,0.9200294102727018,0.13371594,"mRNA binding"),
                     c("GO:0008143","poly(A) binding",0.04161126075001919,1,0.9386791003004495,0.428381,"mRNA binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,5,0.9345942704456097,0.19635374,"mRNA binding"),
                     c("GO:0016208","AMP binding",0.08038413014592037,1,0.9360695411920246,0.29697106,"mRNA binding"),
                     c("GO:0016780","phosphotransferase activity, for other substituted phosphate groups",0.32920461222629094,1,0.8811699811123633,0.35102513,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9526208862981007,0.13769902,"mRNA binding"),
                     c("GO:0030628","pre-mRNA 3'-splice site binding",0.024598652571274335,1,0.9409402074530071,0.41077037,"mRNA binding"),
                     c("GO:0032934","sterol binding",0.11902838493358808,1,0.9520409820555689,0.10144903,"mRNA binding"),
                     c("GO:0033862","UMP kinase activity",0.05586087234498979,1,0.8792496643222873,0.27195941,"mRNA binding"),
                     c("GO:0042834","peptidoglycan binding",0.04667158033603272,1,0.9713245237085395,0.2693746,"mRNA binding"),
                     c("GO:0044183","protein folding chaperone",0.39473153762799984,1,0.9817725532787365,0.32443613,"mRNA binding"),
                     c("GO:0050660","flavin adenine dinucleotide binding",1.7044163107626968,3,0.927582774792719,0.33990845,"mRNA binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,3,0.9412561574070847,0.16455763,"mRNA binding"),
                     c("GO:0070403","NAD+ binding",0.16643173804060435,1,0.9407046839109325,0.29755636,"mRNA binding"),
                     c("GO:0071949","FAD binding",0.630202710371034,1,0.9346852271099148,0.30272555,"mRNA binding"),
                     c("GO:1904047","S-adenosyl-L-methionine binding",0.05472329831009719,1,0.9550440586924542,0.25647876,"mRNA binding"),
                     c("GO:0003735","structural constituent of ribosome",2.128498588986873,3,0.9940744567948091,-0,"structural constituent of ribosome"),
                     c("GO:0003824","catalytic activity",60.727864217389204,10,1,-0,"catalytic activity"),
                     c("GO:0005515","protein binding",8.610051728351934,38,0.9644651205302437,0.0711246,"protein binding"),
                     c("GO:0008289","lipid binding",1.291468066123697,3,0.9716392748624806,0.05041733,"lipid binding"),
                     c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.0036233839629912796,1,0.9621416322696915,0.0291488,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0004165","delta(3)-delta(2)-enoyl-CoA isomerase activity",0.021527424426388823,1,0.958026799223716,0.44447646,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0004619","phosphoglycerate mutase activity",0.0367083831844294,1,0.9536718610095938,0.40080388,"3-hydroxybutyryl-CoA epimerase activity"),
                     c("GO:0009055","electron transfer activity",1.0648291689771099,2,1,-0,"electron transfer activity"),
                     c("GO:0015658","branched-chain amino acid transmembrane transporter activity",0.1658463198238175,1,0.9377674971012107,-0,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0015144","carbohydrate transmembrane transporter activity",0.45600974597151334,1,0.9003538444090412,0.39044024,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0015576","sorbitol transmembrane transporter activity",8.8699729816188E-06,1,0.930714653561947,0.47124298,"branched-chain amino acid transmembrane transporter activity"),
                     c("GO:0016405","CoA-ligase activity",0.2856463924068064,1,0.9532029424934207,0.04318065,"CoA-ligase activity"),
                     c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,2,0.9179776234311781,0.04634246,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016166","phytoene dehydrogenase activity",0.015531322690814519,1,0.9395229810745828,0.2770188,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.9191237651479918,0.44913707,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors",0.9894144411941413,1,0.9258026679233834,0.39379271,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0045174","glutathione dehydrogenase (ascorbate) activity",0.018870867518393994,1,0.9328564423684245,0.28092429,"oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor"),
                     c("GO:0016740","transferase activity",20.627439270288612,28,0.9473893227859949,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,12,0.9514633183740351,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,27,0.9468379296218662,0.12712323,"transferase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,4,0.8829813457657951,0.05667901,"glycosyltransferase activity"),
                     c("GO:0004076","biotin synthase activity",0.016586849475627156,1,0.9212959877218342,0.20205602,"glycosyltransferase activity"),
                     c("GO:0004364","glutathione transferase activity",0.16460674109963627,1,0.9074792207022762,0.24793639,"glycosyltransferase activity"),
                     c("GO:0008168","methyltransferase activity",2.5264987116585993,4,0.8474485826623853,0.33973871,"glycosyltransferase activity"),
                     c("GO:0008194","UDP-glycosyltransferase activity",0.7516503804353587,1,0.8922989155723076,0.2917808,"glycosyltransferase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8828842450274552,0.28934569,"glycosyltransferase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8778123282706567,0.36085078,"glycosyltransferase activity"),
                     c("GO:0030792","methylarsonite methyltransferase activity",0.002785171516228303,1,0.9005106240285965,0.49223587,"glycosyltransferase activity"),
                     c("GO:0106261","tRNA uridine(34) acetyltransferase activity",0.0023993276915278854,1,0.915281639990702,0.17480439,"glycosyltransferase activity"),
                     c("GO:0140903","histone H3R26 methyltransferase activity",3.9914878417284594E-05,1,0.879379968175537,0.45431535,"glycosyltransferase activity"),
                     c("GO:0140984","histone H4K12 methyltransferase activity",7.982975683456919E-05,1,0.8812924459963372,0.39355689,"glycosyltransferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,2,0.9575628370270132,0.05997119,"lyase activity"),
                     c("GO:0016844","strictosidine synthase activity",0.011741626734417886,1,0.9643519416613686,0.03194325,"strictosidine synthase activity"),
                     c("GO:0004300","enoyl-CoA hydratase activity",0.09504176049804544,1,0.9562995134925896,0.39952982,"strictosidine synthase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,3,0.9593802156632314,0.05646078,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,3,0.9580668665713258,0.0589874,"ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,13,0.859843567619985,-0,"ATP hydrolysis activity"),
                     c("GO:0004064","arylesterase activity",0.008042848001082846,1,0.899100324036126,0.3692388,"ATP hydrolysis activity"),
                     c("GO:0004190","aspartic-type endopeptidase activity",0.26350250485819504,3,0.8471534642533737,0.41870234,"ATP hydrolysis activity"),
                     c("GO:0004334","fumarylacetoacetase activity",0.013076557668151514,1,0.9334241023351099,0.19229254,"ATP hydrolysis activity"),
                     c("GO:0004430","1-phosphatidylinositol 4-kinase activity",0.022294677089298852,1,0.8963787319871025,0.37585892,"ATP hydrolysis activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8988180044165047,0.31839638,"ATP hydrolysis activity"),
                     c("GO:0004714","transmembrane receptor protein tyrosine kinase activity",0.12685170110337585,1,0.8401190602818799,0.47011379,"ATP hydrolysis activity"),
                     c("GO:0004842","ubiquitin-protein transferase activity",0.9578927747107135,2,0.8206292145514233,0.48775459,"ATP hydrolysis activity"),
                     c("GO:0008233","peptidase activity",4.167383840940452,9,0.8243494089336002,0.3656957,"ATP hydrolysis activity"),
                     c("GO:0016788","hydrolase activity, acting on ester bonds",6.029659060857018,1,0.8882722935976386,0.39048386,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.9004105219525277,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0070290","N-acylphosphatidylethanolamine-specific phospholipase D activity",0.030998338077512295,2,0.8774789740202524,0.20698894,"ATP hydrolysis activity"),
                     c("GO:0102549","1-18:1-2-16:0-monogalactosyldiacylglycerol lipase activity",3.769738517187989E-05,1,0.9233918951177618,0.23079392,"ATP hydrolysis activity"),
                     c("GO:0106310","protein serine kinase activity",0.08584138102286135,3,0.8445907312449777,0.37283533,"ATP hydrolysis activity"),
                     c("GO:0106375","deoxynucleoside triphosphate hydrolase activity",2.4392425699451695E-05,1,0.9221043593761462,0.3204136,"ATP hydrolysis activity"),
                     c("GO:0140567","membrane protein dislocase activity",0.006211198580378564,1,0.9372040574477077,0.39401957,"ATP hydrolysis activity"),
                     c("GO:0030295","protein kinase activator activity",0.14413927844455088,1,0.9873128641687623,-0,"protein kinase activator activity"),
                     c("GO:0043022","ribosome binding",0.5329368041478477,2,0.9513019234838447,0.04582217,"ribosome binding"),
                     c("GO:0044877","protein-containing complex binding",2.2272125182993086,1,0.9699060444846261,0.054576,"protein-containing complex binding"),
                     c("GO:0046923","ER retention sequence binding",0.013602103567312429,2,0.9666702280639078,0.03325955,"ER retention sequence binding"),
                     c("GO:0010209","vacuolar sorting signal binding",1.1087466227023499E-05,1,0.9757673014373125,0.49249984,"ER retention sequence binding"),
                     c("GO:0048039","ubiquinone binding",0.042507128021162695,1,0.96753560753899,0.03635569,"ubiquinone binding"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9183948627609377,0.0463004,"unfolded protein binding"),
                     c("GO:0003779","actin binding",0.7421262469463454,2,0.9105620046594465,0.44654025,"unfolded protein binding"),
                     c("GO:0017025","TBP-class protein binding",0.05320431543699496,1,0.9320630371182815,0.35368117,"unfolded protein binding"),
                     c("GO:0019900","kinase binding",0.4400460120978449,2,0.9155094327542793,0.42444044,"unfolded protein binding"),
                     c("GO:0030276","clathrin binding",0.11822565237875157,1,0.9280562478047915,0.37746303,"unfolded protein binding"),
                     c("GO:0030544","Hsp70 protein binding",0.07279586826014549,2,0.930543892534813,0.36265291,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9263198924398789,0.38789043,"unfolded protein binding"),
                     c("GO:0031369","translation initiation factor binding",0.06718339285602619,1,0.9309390668773161,0.36031388,"unfolded protein binding"),
                     c("GO:0031370","eukaryotic initiation factor 4G binding",0.0003658863854917755,1,0.9495734550393145,0.25391209,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9207010176079609,0.42217449,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.926138112437574,0.38898648,"unfolded protein binding"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9247713362959409,0.39725472,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9317478271669182,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9134297138568106,0.47948026,"unfolded protein binding"),
                     c("GO:0061649","ubiquitin modification-dependent histone binding",0.0006874229060754569,1,0.9478734934103231,0.26331899,"unfolded protein binding"),
                     c("GO:0070696","transmembrane receptor protein serine/threonine kinase binding",0.01030469111139564,1,0.9390444789707422,0.31312228,"unfolded protein binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Up_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Upregulated, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Up_reduced05_MF_Table.tsv", sep = "\t")
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
