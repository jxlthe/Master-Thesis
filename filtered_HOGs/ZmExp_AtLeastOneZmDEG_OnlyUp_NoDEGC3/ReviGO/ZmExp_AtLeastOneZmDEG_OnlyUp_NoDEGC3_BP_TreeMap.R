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
revigo.data <- rbind(c("GO:0000911","cytokinesis by cell plate formation",0.008712856488926366,1,0.9943739120321408,0.00821801,"cytokinesis by cell plate formation"),
                     c("GO:0006879","intracellular iron ion homeostasis",0.421936404379863,1,1,-0,"intracellular iron ion homeostasis"),
                     c("GO:0007155","cell adhesion",1.1091577223710762,1,0.9917207501957832,0.01321071,"cell adhesion"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,1,1,-0,"metabolic process"),
                     c("GO:0008202","steroid metabolic process",0.582585703698599,2,0.8123060199978509,-0,"steroid metabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,1,0.9366026159428562,0.15163831,"steroid metabolic process"),
                     c("GO:0006096","glycolytic process",0.4557945400484291,1,0.8058495451850985,0.48590569,"steroid metabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,4,0.9352545421856103,0.10901633,"steroid metabolic process"),
                     c("GO:0008203","cholesterol metabolic process",0.2200421431132646,2,0.790056955917325,0.45443823,"steroid metabolic process"),
                     c("GO:0008299","isoprenoid biosynthetic process",0.527156162091961,1,0.7666090970461414,0.49264795,"steroid metabolic process"),
                     c("GO:0042759","long-chain fatty acid biosynthetic process",0.04469313344093403,1,0.7993405984293194,0.3981042,"steroid metabolic process"),
                     c("GO:0046856","phosphatidylinositol dephosphorylation",0.10928164950665267,1,0.7771081847394969,0.42785472,"steroid metabolic process"),
                     c("GO:0009793","embryo development ending in seed dormancy",0.0206816347379296,2,0.8483925358219894,-0,"embryo development ending in seed dormancy"),
                     c("GO:0007517","muscle organ development",0.07354784657130489,1,0.8854406117318956,0.40853039,"embryo development ending in seed dormancy"),
                     c("GO:0010015","root morphogenesis",0.013467340270295237,1,0.8560617638199274,0.43918043,"embryo development ending in seed dormancy"),
                     c("GO:0010090","trichome morphogenesis",0.006709022733481632,1,0.861821252260259,0.35446469,"embryo development ending in seed dormancy"),
                     c("GO:0015031","protein transport",3.093438694074183,3,0.927458842668181,0,"protein transport"),
                     c("GO:0016192","vesicle-mediated transport",2.6087919056355493,1,0.9427287900065358,0.38566904,"protein transport"),
                     c("GO:0051641","cellular localization",5.891744471119505,1,0.9405307129667297,0.41259033,"protein transport"),
                     c("GO:0016135","saponin biosynthetic process",9.858960666394757E-06,1,0.9554118934769257,0.03120971,"saponin biosynthetic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,3,0.8780480233264657,0.05779371,"phosphorylation"),
                     c("GO:0050896","response to stimulus",17.567785530535815,1,1,-0,"response to stimulus"),
                     c("GO:0070417","cellular response to cold",0.007012185773973271,1,0.9078628220050563,0.00808117,"cellular response to cold"),
                     c("GO:0006950","response to stress",6.919654498638822,1,0.8945346792326945,0.4896406,"cellular response to cold"),
                     c("GO:0006952","defense response",1.1604144588756624,1,0.8982880972287655,0.45786068,"cellular response to cold"),
                     c("GO:0009611","response to wounding",0.16569462243976346,1,0.9123593013600194,0.31137091,"cellular response to cold"),
                     c("GO:0009624","response to nematode",0.0011042035946362127,1,0.9274399796732652,0.34321484,"cellular response to cold"),
                     c("GO:0010274","hydrotropism",0.006815006560645376,1,0.9110428880140717,0.14965038,"cellular response to cold"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,1,0.8330004593099042,0.30465566,"cellular response to cold"),
                     c("GO:0043434","response to peptide hormone",0.19767216136121488,1,0.8792714644655645,0.18156262,"cellular response to cold"),
                     c("GO:0071555","cell wall organization",0.8943925879544993,1,0.9918913945769161,0.01286368,"cell wall organization"),
                     c("GO:0098734","macromolecule depalmitoylation",0.06083225205182225,1,0.9562316419110493,0.05291306,"macromolecule depalmitoylation"),
                     c("GO:0006013","mannose metabolic process",0.050640551462936674,1,0.9316211055970944,0.33438746,"macromolecule depalmitoylation"),
                     c("GO:0006412","translation",4.38869169324396,1,0.8436745944312076,0.25211911,"macromolecule depalmitoylation"),
                     c("GO:0045489","pectin biosynthetic process",0.012338489273993038,1,0.9108732274597405,0.10510175,"macromolecule depalmitoylation"),
                     c("GO:1905421","regulation of plant organ morphogenesis",0.005429822587016912,1,0.9545010570790557,0.09345563,"regulation of plant organ morphogenesis"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1,0.9470816734848149,0.1515582,"regulation of plant organ morphogenesis"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,3,0.9233290044708715,0.38614556,"regulation of plant organ morphogenesis"),
                     c("GO:0010030","positive regulation of seed germination",0.0013630013121290752,1,0.9522524048671861,0.39302598,"regulation of plant organ morphogenesis"),
                     c("GO:0032012","regulation of ARF protein signal transduction",0.05576721100946194,1,0.9542546721464553,0.14116446,"regulation of plant organ morphogenesis"),
                     c("GO:0040008","regulation of growth",0.2209540969749061,1,0.9512874129782445,0.12053119,"regulation of plant organ morphogenesis"),
                     c("GO:0045927","positive regulation of growth",0.06802189911779062,1,0.9501629219113467,0.14325584,"regulation of plant organ morphogenesis"),
                     c("GO:1902066","regulation of cell wall pectin metabolic process",0.00012816648866313183,1,0.9645507253381704,0.16810102,"regulation of plant organ morphogenesis"),
                     c("GO:2000012","regulation of auxin polar transport",0.002183759787606439,1,0.9424827065100794,-0,"regulation of auxin polar transport"),
                     c("GO:0010540","basipetal auxin transport",0.0002908393396586453,1,0.9285531757265586,0.45061093,"regulation of auxin polar transport"),
                     c("GO:0048209","regulation of vesicle targeting, to, from or within Golgi",2.4647401665986893E-06,1,0.9661947727929658,0.33066454,"regulation of auxin polar transport"),
                     c("GO:0099175","regulation of postsynapse organization",0.027792410118566816,1,0.9449895798698879,0.3459764,"regulation of auxin polar transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
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
  title = paste("Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_BP_Table.tsv", sep = "\t")
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

