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
revigo.data <- rbind(c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,3,0.9234219310168503,0.066121,"transcription cis-regulatory region binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.9356732319507387,0.22593798,"transcription cis-regulatory region binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,11,0.9208157986383229,0.45223656,"transcription cis-regulatory region binding"),
                     c("GO:0003680","minor groove of adenine-thymine-rich DNA binding",0.022786960589778697,1,0.9558081254739234,0.31120183,"transcription cis-regulatory region binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,7,0.9892278958137051,0.40628023,"transcription cis-regulatory region binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.9419095650686188,0.31664803,"transcription cis-regulatory region binding"),
                     c("GO:0003743","translation initiation factor activity",0.37315089336372126,1,0.9344428751834606,0.27764537,"transcription cis-regulatory region binding"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9559633361948003,0.14801235,"transcription cis-regulatory region binding"),
                     c("GO:0003674","molecular_function",100,5,1,-0,"molecular_function"),
                     c("GO:0003682","chromatin binding",0.6287702097345026,1,0.9750857260813484,0.0533165,"chromatin binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,4,1,-0,"catalytic activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,8,0.7712364958384572,0,"protein kinase activity"),
                     c("GO:0003756","protein disulfide isomerase activity",0.04975833093363606,2,0.8956709135511683,0.34856872,"protein kinase activity"),
                     c("GO:0004180","carboxypeptidase activity",0.4650748583587277,2,0.8278945384231705,0.43795071,"protein kinase activity"),
                     c("GO:0004781","sulfate adenylyltransferase (ATP) activity",0.04230755362907627,1,0.8917147765889218,0.41713542,"protein kinase activity"),
                     c("GO:0008483","transaminase activity",0.6992399275802186,1,0.8730524232571509,0.30161168,"protein kinase activity"),
                     c("GO:0015035","protein-disulfide reductase activity",0.26137371134260656,1,0.8518403620480547,0.4107919,"protein kinase activity"),
                     c("GO:0016746","acyltransferase activity",3.815886769118107,2,0.8739778927909299,0.37902661,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,5,0.879331752206891,0.35580262,"protein kinase activity"),
                     c("GO:0016760","cellulose synthase (UDP-forming) activity",0.02088878637171227,2,0.8896342995963697,0.21201958,"protein kinase activity"),
                     c("GO:0046027","phospholipid:diacylglycerol acyltransferase activity",0.010879021861955458,1,0.9061901827145281,0.20092889,"protein kinase activity"),
                     c("GO:0046608","carotenoid isomerase activity",0.0010289168658677808,1,0.9606524111978103,0.38075906,"protein kinase activity"),
                     c("GO:0047334","diphosphate-fructose-6-phosphate 1-phosphotransferase activity",0.018686815579025403,1,0.880676507859515,0.48733868,"protein kinase activity"),
                     c("GO:0061630","ubiquitin protein ligase activity",0.6665452071699717,3,0.8153959195366739,0.45681379,"protein kinase activity"),
                     c("GO:0005515","protein binding",8.610051728351934,20,0.9655182250160117,0.07766851,"protein binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,4,0.9259323685080195,-0,"zinc ion binding"),
                     c("GO:0000822","inositol hexakisphosphate binding",0.014655412858879661,1,0.9597188024449858,0.21260669,"zinc ion binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,2,0.94135628999035,0.3693843,"zinc ion binding"),
                     c("GO:0005525","GTP binding",1.780363236924562,2,0.9210495842930069,0.19635374,"zinc ion binding"),
                     c("GO:0005546","phosphatidylinositol-4,5-bisphosphate binding",0.08566176406998356,1,0.9537792124994369,0.24645234,"zinc ion binding"),
                     c("GO:0016597","amino acid binding",0.16394149312601486,1,0.9511244949970609,0.26177282,"zinc ion binding"),
                     c("GO:0030170","pyridoxal phosphate binding",1.0388268431814944,1,0.9278312378715666,0.31800282,"zinc ion binding"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9313622633952597,0.18168089,"zinc ion binding"),
                     c("GO:0070402","NADPH binding",0.10609152933989706,1,0.940630917527751,0.28523818,"zinc ion binding"),
                     c("GO:0009881","photoreceptor activity",0.06493041971869501,1,0.9932715597699995,-0,"photoreceptor activity"),
                     c("GO:0010011","auxin binding",0.0009269121765791645,1,0.9851425715980695,0.02993406,"auxin binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9091395406627808,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9103773871717047,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0008670","2,4-dienoyl-CoA reductase (NADPH) activity",0.02439020820620629,1,0.9339076030275826,0.31075808,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8921448713017424,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1,0.9225127918880266,0.35139584,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016740","transferase activity",20.627439270288612,20,0.9383201183959232,0.0817378,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,9,0.9431975421019412,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,14,0.9376609006114417,0.12712323,"transferase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,5,0.9504884412870329,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,4,0.9526540868602282,0.05541399,"isomerase activity"),
                     c("GO:0016874","ligase activity",3.248266153118885,3,0.9510894178819743,0.05784577,"ligase activity"),
                     c("GO:0016881","acid-amino acid ligase activity",0.3905693028063752,1,0.9491572006995631,0.04406353,"acid-amino acid ligase activity"),
                     c("GO:0016851","magnesium chelatase activity",0.008235769913433055,1,0.9612669306650197,0.44358297,"acid-amino acid ligase activity"),
                     c("GO:0080123","jasmonoyl-L-amino acid ligase activity",7.09597838529504E-05,1,0.9681016194434066,0.4447539,"acid-amino acid ligase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,3,0.8546466589755487,0.05997119,"ATP hydrolysis activity"),
                     c("GO:0004566","beta-glucuronidase activity",0.008790143224784231,1,0.9039648109878481,0.4776459,"ATP hydrolysis activity"),
                     c("GO:0016791","phosphatase activity",1.6807667453004556,2,0.8256388430999533,0.32020834,"ATP hydrolysis activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8803578753850313,0.32827501,"ATP hydrolysis activity"),
                     c("GO:0016810","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds",1.9557758226090554,1,0.8806899373247634,0.32699612,"ATP hydrolysis activity"),
                     c("GO:0016811","hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",0.9259763944296038,1,0.8794282273547864,0.29603328,"ATP hydrolysis activity"),
                     c("GO:0046556","alpha-L-arabinofuranosidase activity",0.046915504593027235,1,0.8929193777485857,0.2148739,"ATP hydrolysis activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9737960570931955,0.05647972,"carbohydrate binding"),
                     c("GO:0032977","membrane insertase activity",0.038238453523758646,2,0.9657882125052901,-0,"membrane insertase activity"),
                     c("GO:0015085","calcium ion transmembrane transporter activity",0.312553455446547,1,0.9213806366238194,0.40515784,"membrane insertase activity"),
                     c("GO:0015297","antiporter activity",0.5980379708464388,1,0.9248862352263971,0.43240123,"membrane insertase activity"),
                     c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,1,0.9433329921465046,0.40847859,"membrane insertase activity"),
                     c("GO:0046943","carboxylic acid transmembrane transporter activity",1.0773003509892658,1,0.9264479109891621,0.21624994,"membrane insertase activity"),
                     c("GO:0045158","electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity",0.003869525713231201,1,0.9943526289119874,-0,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0009496","plastoquinol--plastocyanin reductase activity",0.0010510917983218278,1,0.8916867778567433,0.41566947,"electron transporter, transferring electrons within cytochrome b6/f complex of photosystem II activity"),
                     c("GO:0051082","unfolded protein binding",0.5891679978648201,3,0.9279089282571554,0.05290415,"unfolded protein binding"),
                     c("GO:0005516","calmodulin binding",0.264039138223583,2,0.9325759945518426,0.40485813,"unfolded protein binding"),
                     c("GO:0031072","heat shock protein binding",0.16268417445587038,1,0.9351237887764527,0.38789043,"unfolded protein binding"),
                     c("GO:0033612","receptor serine/threonine kinase binding",0.01374402313501833,2,0.9456650303771515,0.31955157,"unfolded protein binding"),
                     c("GO:0035064","methylated histone binding",0.1065793778538861,1,0.9372015203400825,0.37419503,"unfolded protein binding"),
                     c("GO:0042393","histone binding",0.41579772345934446,1,0.9300078306032469,0.42217449,"unfolded protein binding"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,3,0.928905480268335,0.42967902,"unfolded protein binding"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.9349583140706972,0.38898648,"unfolded protein binding"),
                     c("GO:0043621","protein self-association",0.05684543934594948,2,0.9400609143876874,0.3555383,"unfolded protein binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9233978415106203,0.4679269,"unfolded protein binding"),
                     c("GO:0051087","protein-folding chaperone binding",0.3040693262896286,1,0.9317970257605913,0.41008684,"unfolded protein binding"),
                     c("GO:0140097","catalytic activity, acting on DNA",2.838298219401709,1,0.9332831307323307,0.0567159,"catalytic activity, acting on DNA"),
                     c("GO:0035514","DNA demethylase activity",0.020587207290337233,1,0.9528444094497709,0.41904136,"catalytic activity, acting on DNA"),
                     c("GO:0140311","protein sequestering activity",0.002454765022663003,1,1,-0,"protein sequestering activity"),
                     c("GO:0140859","pterocarpan synthase activity",2.2174932454047E-06,1,0.9653813869915613,0.01868194,"pterocarpan synthase activity"),
                     c("GO:0016831","carboxy-lyase activity",0.6244594028654359,1,0.9347660347962944,0.48814337,"pterocarpan synthase activity"),
                     c("GO:0030570","pectate lyase activity",0.03267476297103825,1,0.9458411415056922,0.32867735,"pterocarpan synthase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_NoDEGZm_AllC3Down_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, No DEG in Zm, All C3 Downregulated, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_NoDEGZm_AllC3Down_reduced05_MF_Table.tsv", sep = "\t")
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


