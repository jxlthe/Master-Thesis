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
                     c("GO:0003729","mRNA binding",1.136762432364793,2,0.91463056486156,0.0528646,"mRNA binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,2,0.8983401523761327,0.31664803,"mRNA binding"),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9647688759912599,0.35731485,"mRNA binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,8,0.88239745029626,0.45223656,"mRNA binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8932722122445766,0.24450519,"mRNA binding"),
                     c("GO:0020037","heme binding",1.453841791525211,1,0.9360712416699458,0.13060477,"mRNA binding"),
                     c("GO:0003824","catalytic activity",60.727864217389204,2,1,-0,"catalytic activity"),
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9677288198492044,0.03485749,"phosphopyruvate hydratase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,5,0.7343532747717424,0,"protein kinase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8731727486406581,0.45866285,"protein kinase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8737872867029048,0.48555173,"protein kinase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9017348796980814,0.2197454,"protein kinase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.7781290636399422,0.41376722,"protein kinase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8142778824417033,0.44701273,"protein kinase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8227360115156692,0.35123485,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8790819384399416,0.35580262,"protein kinase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8892380938671747,0.31489566,"protein kinase activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9184520922391256,0.21182638,"protein kinase activity"),
                     c("GO:0005085","guanyl-nucleotide exchange factor activity",0.4464279576581196,1,1,-0,"guanyl-nucleotide exchange factor activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,9,0.9504356501827331,0.06960742,"protein binding"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,8,0.9360752783255815,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9421263523713012,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,8,0.9367913538542104,0.12712323,"hydrolase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,2,0.8827212565347011,0.05402108,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.8732006896593303,0.29096067,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8739323593148979,0.32827501,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.8957649409643302,0.17478095,"hydrolase activity, acting on glycosyl bonds"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9500792730745891,0.05879155,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.9524195269800306,0.05541399,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9543983189341552,0.03465455,"intramolecular lyase activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,2,0.9630106320612827,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9487311608463957,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,2,0.9222392509874762,0.05306855,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9277179683923644,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,2,0.933832199442618,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.9325696250840736,0.42311414,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9390304801000545,0.37611117,"protein dimerization activity"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.9246415584092138,0.0474854,"dioxygenase activity"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,1,0.9207737630804183,0.4568902,"dioxygenase activity"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,1,0.9196885870423279,0.42767891,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.9231628623785453,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9532195581978973,0.45929075,"dioxygenase activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9112094691370106,-0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,7,0.8592894703414051,0.28451694,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9078469500025347,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,2,0.9047787738628846,0.34784993,"iron-sulfur cluster binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.9107662553807204,0.3693843,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8902244543675969,0.18168089,"iron-sulfur cluster binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9107682460071076,0.11631471,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)

# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, Only Any Zm DE, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

if(length(topGO_data$GO.ID) > 0)
{
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
}


# also find the reduced form of topGO sig Terms and mark them too
not_plotted <- topGO_data[!topGO_data$Term %in% topGO_ReviGO_Intersect,]$GO.ID
ReviGO_table <- read.csv(".\\ZmExp_OnlyAnyZmDEG_NoDEGC3_reduced05_MF_Table.tsv", sep = "\t")
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



