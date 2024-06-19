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
revigo.data <- rbind(c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.287885350986827,1,0.9423177302167721,-0,"DNA-binding transcription factor activity, RNA polymerase II-specific"),
                     c("GO:0004672","protein kinase activity",3.5248673905776853,3,0.7712757143447455,0,"protein kinase activity"),
                     c("GO:0016757","glycosyltransferase activity",2.478461155483397,1,0.8720970226076804,0.35580262,"protein kinase activity"),
                     c("GO:0005198","structural molecule activity",3.1132296844467198,1,1,-0,"structural molecule activity"),
                     c("GO:0005515","protein binding",8.610051728351934,6,0.939878503395764,0.06452231,"protein binding"),
                     c("GO:0016740","transferase activity",20.627439270288612,4,0.939315351942652,0.07956199,"transferase activity"),
                     c("GO:0016491","oxidoreductase activity",11.23609814678324,2,0.9445879769274047,0.10406278,"transferase activity"),
                     c("GO:0016787","hydrolase activity",22.243200201100027,3,0.9386239715185933,0.12712323,"transferase activity"),
                     c("GO:0030246","carbohydrate binding",1.0034689133835164,1,0.9556232196159946,0.04919141,"carbohydrate binding"),
                     c("GO:0046983","protein dimerization activity",1.1741382810160892,1,0.9013056627527537,0.04828594,"protein dimerization activity"),
                     c("GO:0008083","growth factor activity",0.15856185451266308,1,0.8796806122877222,0.41153522,"protein dimerization activity"),
                     c("GO:0042803","protein homodimerization activity",0.16806824805571302,1,0.9148935368882537,0.41374711,"protein dimerization activity"),
                     c("GO:0043621","protein self-association",0.05684543934594948,1,0.9212121829196258,0.37611117,"protein dimerization activity"),
                     c("GO:0051213","dioxygenase activity",0.7412059872495025,1,0.8874789418037062,0.0474854,"dioxygenase activity"),
                     c("GO:0003954","NADH dehydrogenase activity",0.34065796483880617,1,0.8939506935615665,0.36657586,"dioxygenase activity"),
                     c("GO:0016616","oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",1.7866720052077387,1,0.8648493115254148,0.43460917,"dioxygenase activity"),
                     c("GO:0016706","2-oxoglutarate-dependent dioxygenase activity",0.28018914152986546,1,0.8821277040160563,0.35993207,"dioxygenase activity"),
                     c("GO:0047998","hyoscyamine (6S)-dioxygenase activity",2.2174932454046998E-05,1,0.9255406726276421,0.45929075,"dioxygenase activity"),
                     c("GO:0051287","NAD binding",0.8468961503119814,1,0.8756421678359987,-0,"NAD binding"),
                     c("GO:0000166","nucleotide binding",18.476222463002568,3,0.8164320127475632,0.42621458,"NAD binding"),
                     c("GO:0000976","transcription cis-regulatory region binding",3.1655181751733625,1,0.8745510615030264,0.31664803,"NAD binding"),
                     c("GO:0003676","nucleic acid binding",20.580503808256374,1,0.872334675051111,0.18774095,"NAD binding"),
                     c("GO:0003677","DNA binding",12.252016067246458,5,0.8481322633103344,0.45223656,"NAD binding"),
                     c("GO:0003729","mRNA binding",1.136762432364793,1,0.8897315936146176,0.12297344,"NAD binding"),
                     c("GO:0005509","calcium ion binding",1.3422486614434648,1,0.9111581046147357,0.14439095,"NAD binding"),
                     c("GO:0031418","L-ascorbic acid binding",0.06531848103664084,1,0.9075733758048868,0.21206178,"NAD binding"),
                     c("GO:0046872","metal ion binding",18.074696526070646,4,0.8758810110223696,0.45949964,"NAD binding"),
                     c("GO:0052793","pectin acetylesterase activity",0.007730181453480783,1,0.9587271113193424,0.03057439,"pectin acetylesterase activity"),
                     c("GO:0004553","hydrolase activity, hydrolyzing O-glycosyl compounds",1.6123703836391927,1,0.937769329963123,0.17194065,"pectin acetylesterase activity"),
                     c("GO:0016798","hydrolase activity, acting on glycosyl bonds",2.011009144365596,1,0.9365366538460628,0.29096067,"pectin acetylesterase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_MF_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "MF" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("Zm Expanded, At Least One Zm DEG Only Downregulated, No DEG in C3, MF, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_MF_Table.tsv", sep = "\t")
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

