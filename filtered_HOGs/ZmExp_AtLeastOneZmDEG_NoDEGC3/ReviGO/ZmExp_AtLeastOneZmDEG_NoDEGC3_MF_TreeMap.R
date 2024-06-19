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
revigo.data <- rbind(c("GO:0003674","molecular_function",100,3,1,-0,"molecular_function"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9674763781434593,0.03485749,"phosphopyruvate hydratase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,5,0.7495939076464454,0,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8800192775174226,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8806390324819089,0.48555173,"protein kinase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9064507700838336,0.2197454,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.7911728437529135,0.41376722,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8246484213155818,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8330343102594516,0.35123485,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8843305413378043,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8940567007377429,0.31489566,"protein kinase activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9220076937755323,0.21182638,"protein kinase activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,10,0.9513977847096342,0.06960742,"protein binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.8906488280094385,0.05213064,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0003954","NADH dehydrogenase activity",0.34065796483880617,1,0.9039063723566827,0.39393019,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.8921535204635909,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8819424727684569,0.47359913,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.918923140274019,0.32214763,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9352520491503108,0.19910628,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8974800889555598,0.42767891,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,9,0.9359688855353118,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,4,0.9419777446730677,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,8,0.9366806644692963,0.12712323,"hydrolase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8877803497394471,0.05402108,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8790333722593702,0.29096067,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8793658005888485,0.32827501,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9011474878550547,0.17478095,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9498762419375432,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9522034682976344,0.05541399,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9551689980615783,0.03465455,"intramolecular lyase activity"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9367313112785355,0.05445618,"heme binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,2,0.9014601427576004,0.31664803,"heme binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9624150224777499,0.35731485,"heme binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,2,0.9061840449122336,0.19568843,"heme binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,8,0.8858787347320689,0.45223656,"heme binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.9170075268241012,0.13060477,"heme binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8945491782583571,0.24450519,"heme binding"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9636134301546603,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9505756254153493,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9168542342812899,0.05306855,"protein dimerization activity"),
                     c("GO:0008083","growth factor activity",0.15856185451266308,1,0.8944806339717133,0.41153522,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.922727368133231,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.929258472289716,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9279111023338894,0.42311414,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9348016294787308,0.37611117,"protein dimerization activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9131244255135076,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,7,0.8582856764593443,0.42621458,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9100641587982604,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.912898418532941,0.3693843,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8933372358306384,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9114834216265212,0.21206178,"iron-sulfur cluster binding"),
                     c("GO:0051287","NAD binding",0.8468961503119814,1,0.9046146732125713,0.15052371,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)



# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE, No C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_NoDEGC3_MF_Table.tsv", sep = "\t")
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


