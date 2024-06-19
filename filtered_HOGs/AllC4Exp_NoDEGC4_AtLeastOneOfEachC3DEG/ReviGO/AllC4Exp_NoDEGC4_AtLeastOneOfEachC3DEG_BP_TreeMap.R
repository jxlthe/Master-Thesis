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
revigo.data <- rbind(c("GO:0006811","monoatomic ion transport",4.7761710300947735,5,0.8693142036127602,0,"monoatomic ion transport"),
                     c("GO:0034220","monoatomic ion transmembrane transport",4.161151810543902,2,0.8166078163982737,0.44154323,"monoatomic ion transport"),
                     c("GO:0035672","oligopeptide transmembrane transport",0.23299928216907387,1,0.8627103740851154,0.43869347,"monoatomic ion transport"),
                     c("GO:0046907","intracellular transport",2.9683457363987995,1,0.8501483265362694,0.4188028,"monoatomic ion transport"),
                     c("GO:0051028","mRNA transport",0.25417879442065644,1,0.8648015434858708,0.30464429,"monoatomic ion transport"),
                     c("GO:0006952","defense response",1.1604144588756624,3,0.8381661356544009,-0,"defense response"),
                     c("GO:0009651","response to salt stress",0.07343446852364134,2,0.8545317476844199,0.4259419,"defense response"),
                     c("GO:0010039","response to iron ion",0.05924742412469929,1,0.8678249553799107,0.46879008,"defense response"),
                     c("GO:0019722","calcium-mediated signaling",0.15710007347883384,1,0.8452169740379526,0.46097507,"defense response"),
                     c("GO:0035556","intracellular signal transduction",4.1384737362710275,2,0.7888309661612708,0.37522238,"defense response"),
                     c("GO:0042221","response to chemical",4.8560360057130705,1,0.8500399690180803,0.46171648,"defense response"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1,0.8262373553023209,0.42151038,"defense response"),
                     c("GO:0008150","biological_process",100,2,1,-0,"biological_process"),
                     c("GO:0008152","metabolic process",57.597931274565454,2,1,-0,"metabolic process"),
                     c("GO:0009820","alkaloid metabolic process",0.01601588160255828,2,0.9436353225576943,0.05107695,"alkaloid metabolic process"),
                     c("GO:0006413","translational initiation",0.4742529791560868,1,0.8654486740322442,0.36848795,"alkaloid metabolic process"),
                     c("GO:0006457","protein folding",1.174377211919444,1,0.8677155725862653,0.27735821,"alkaloid metabolic process"),
                     c("GO:0006508","proteolysis",5.2622572267907,2,0.8680221978172313,0.44475214,"alkaloid metabolic process"),
                     c("GO:0030091","protein repair",0.06087415263465443,1,0.91015480010716,0.11683693,"alkaloid metabolic process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,2,0.8973857222610786,-0,"lipid catabolic process"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,2,0.9353800814722798,0.15163831,"lipid catabolic process"),
                     c("GO:0006629","lipid metabolic process",6.477701939366011,2,0.9338811362930068,0.12277613,"lipid catabolic process"),
                     c("GO:0006654","phosphatidic acid biosynthetic process",0.11421359458001666,1,0.9085864140814529,0.46381088,"lipid catabolic process"),
                     c("GO:0016137","glycoside metabolic process",0.12090536413233209,1,0.9517454394802891,0.06048027,"glycoside metabolic process"),
                     c("GO:0033491","coniferin metabolic process",1.7253181166190823E-05,1,0.9442566754145938,0.44674727,"glycoside metabolic process"),
                     c("GO:0016310","phosphorylation",5.235381700014107,1,0.9338745778886478,0.07741388,"phosphorylation"),
                     c("GO:0009058","biosynthetic process",29.004480601854056,1,0.9507663781277531,0.11172681,"phosphorylation"),
                     c("GO:0030261","chromosome condensation",0.1195448275603696,1,0.9798385633655547,-0,"chromosome condensation"),
                     c("GO:0033500","carbohydrate homeostasis",0.08834861127173001,1,0.9829765143020287,-0,"carbohydrate homeostasis"),
                     c("GO:0048364","root development",0.0376070054619628,1,0.9455995749748067,-0,"root development"),
                     c("GO:0010233","phloem transport",0.0002686566781592571,1,0.9004928829494658,0.29362419,"root development"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,1,0.9484621531515232,0.43353892,"root development"),
                     c("GO:0048316","seed development",0.03280076213709535,1,0.9264830084500535,0.46716202,"root development"),
                     c("GO:0080111","DNA demethylation",0.039697105123238485,1,0.8845783069732632,0.07105297,"DNA demethylation"),
                     c("GO:0006076","(1->3)-beta-D-glucan catabolic process",2.4647401665986893E-06,1,0.9449009301204617,0.22339313,"DNA demethylation"),
                     c("GO:0009699","phenylpropanoid biosynthetic process",0.0697373582737433,1,0.9011708594395732,0.15925058,"DNA demethylation"),
                     c("GO:0010731","protein glutathionylation",2.4647401665986892E-05,1,0.9307160410633457,0.23039539,"DNA demethylation"),
                     c("GO:0016926","protein desumoylation",0.02830261133305275,1,0.8950802362388776,0.33184604,"DNA demethylation"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1,0.893709612099068,0.15970956,"DNA demethylation"),
                     c("GO:0080151","positive regulation of salicylic acid mediated signaling pathway",7.887168533115806E-05,1,0.9670039793939942,0.08184898,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:0045739","positive regulation of DNA repair",0.0330817425160876,1,0.9461004114183711,0.41456906,"positive regulation of salicylic acid mediated signaling pathway"),
                     c("GO:1902456","regulation of stomatal opening",0.0015108857221249963,1,0.9714778756693939,-0,"regulation of stomatal opening"),
                     c("GO:0006355","regulation of DNA-templated transcription",11.048858347143273,1,0.9208134133543636,0.38658972,"regulation of stomatal opening"),
                     c("GO:0042325","regulation of phosphorylation",0.2008023813727952,1,0.94918833865427,0.11892868,"regulation of stomatal opening"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches

topGO_data <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("All C4 Expanded, No DEG in C4, At Least One Of Each C3 DE, BP, ReviGO 0.5, # HOGs", n_hogs, ", Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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
ReviGO_table <- read.csv(".\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_BP_Table.tsv", sep = "\t")
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

