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
                     c("GO:0004634","phosphopyruvate hydratase activity",0.0373204113201611,1,0.9675380968275328,0.03247589,"phosphopyruvate hydratase activity"),
                     c("GO:0004722","protein serine/threonine phosphatase activity",0.27943519382642784,2,0.775754863567482,0.03898172,"protein serine/threonine phosphatase activity"),
                     c("GO:0004103","choline kinase activity",0.008863320501882585,1,0.8871755471748216,0.45866285,"protein serine/threonine phosphatase activity"),
                     c("GO:0004305","ethanolamine kinase activity",0.017884083024188903,1,0.8899807862151725,0.48555173,"protein serine/threonine phosphatase activity"),
                     c("GO:0004559","alpha-mannosidase activity",0.07698914798720577,1,0.908677723345562,0.17913273,"protein serine/threonine phosphatase activity"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,2,0.7403563373645067,0.41376722,"protein serine/threonine phosphatase activity"),
                     c("GO:0008234","cysteine-type peptidase activity",0.5549587295679618,1,0.8171674534539154,0.44701273,"protein serine/threonine phosphatase activity"),
                     c("GO:0008474","palmitoyl-(protein) hydrolase activity",0.0540713552959482,1,0.8254607650217357,0.3938196,"protein serine/threonine phosphatase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.8781281494576612,0.32827501,"protein serine/threonine phosphatase activity"),
                     c("GO:0016887","ATP hydrolysis activity",4.018519084389943,1,0.8690471915223179,0.25703449,"protein serine/threonine phosphatase activity"),
                     c("GO:0005085","guanyl-nucleotide exchange factor activity",0.4464279576581196,1,1,-0,"guanyl-nucleotide exchange factor activity"),
                     c("GO:0005515","protein binding",8.610051728351934,4,0.9467882816804609,0.06960742,"protein binding"),
                     c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen",1.5459808533650217,2,0.9327176990849406,-0,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0004497","monooxygenase activity",1.3240297369392198,2,0.9335678771959217,0.4568902,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0018685","alkane 1-monooxygenase activity",0.003854003260513368,1,0.9562486155148083,0.27075548,"oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,6,0.9350402801182242,0.08125579,"hydrolase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9412670072024919,0.10619003,"hydrolase activity"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.9357744210752793,0.12712323,"hydrolase activity"),
                     c("GO:0016829","lyase activity",3.6221510367468346,1,0.9494660574333623,0.05562821,"lyase activity"),
                     c("GO:0016853","isomerase activity",2.413124934500793,1,0.951874095225543,0.04962823,"isomerase activity"),
                     c("GO:0016872","intramolecular lyase activity",0.034943258561087265,1,0.9510245110061127,0.03229967,"intramolecular lyase activity"),
                     c("GO:0020037","heme binding",1.453841791525211,2,0.9326437823545031,0.05445618,"heme binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,3,0.8792267248252142,0.48153711,"heme binding"),
                     c("GO:0003700","DNA-binding transcription factor activity",5.926133171202054,2,0.986595269489872,0.3967737,"heme binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.9119890358065028,0.33073906,"heme binding"),
                     c("GO:0003746","translation elongation factor activity",0.305002890945944,1,0.8995508020656475,0.28189022,"heme binding"),
                     c("GO:0043565","sequence-specific DNA binding",4.437827099656763,1,0.8897098336288308,0.15481958,"heme binding"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,1,0.9604459476809007,0.05209277,"carbohydrate binding"),
                     c("GO:0035091","phosphatidylinositol binding",0.4369681314732231,1,0.9423192183734204,0.04747271,"phosphatidylinositol binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9212899096591148,0.05306855,"protein dimerization activity"),
                     c("GO:0042802","identical protein binding",0.5005104004202948,1,0.9268058925071231,0.46009366,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9330046203362361,0.41374711,"protein dimerization activity"),
                     c("GO:0043130","ubiquitin binding",0.21361999430281634,1,0.931722119464749,0.42311414,"protein dimerization activity"),
                     c("GO:0050291","sphingosine N-acyltransferase activity",0.0206648195539264,1,0.9253715940546006,0.03095884,"sphingosine N-acyltransferase activity"),
                     c("GO:0004337","geranyltranstransferase activity",0.03165028109166128,1,0.9054205927040189,0.1977699,"sphingosine N-acyltransferase activity"),
                     c("GO:0016765","transferase activity, transferring alkyl or aryl (other than methyl) groups",0.9927517685284755,1,0.8987044811614465,0.19133236,"sphingosine N-acyltransferase activity"),
                     c("GO:0051536","iron-sulfur cluster binding",1.9432846831576909,2,0.9029540873696541,0,"iron-sulfur cluster binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,4,0.8524816754298336,0.26407481,"iron-sulfur cluster binding"),
                     c("GO:0000287","magnesium ion binding",1.7610111733719154,2,0.9047778837205973,0.3830934,"iron-sulfur cluster binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.8992681387269431,0.34784993,"iron-sulfur cluster binding"),
                     c("GO:0008270","zinc ion binding",3.7731179768939866,2,0.8847467633007552,0.18168089,"iron-sulfur cluster binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)




# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_MF_Table.tsv", sep = "\t")
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
