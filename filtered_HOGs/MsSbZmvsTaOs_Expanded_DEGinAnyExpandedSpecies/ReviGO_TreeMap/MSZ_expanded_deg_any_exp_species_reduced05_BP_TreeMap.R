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
revigo.data <- rbind(c("GO:0006123","mitochondrial electron transport, cytochrome c to oxygen",0.047058804746273615,1,0.9708155666584158,0.05401416,"mitochondrial electron transport, cytochrome c to oxygen"),
c("GO:0006457","protein folding",1.0590625918025878,2,0.9915181497788998,0.01051422,"protein folding"),
c("GO:0006730","one-carbon metabolic process",0.3470171077448388,1,0.9520586567451994,0.02817366,"one-carbon metabolic process"),
c("GO:0006108","malate metabolic process",0.054951831581367426,1,0.9516425963768287,0.2431141,"one-carbon metabolic process"),
c("GO:0006865","amino acid transport",0.6097621009115332,2,0.9183573676658028,-0,"amino acid transport"),
c("GO:0006811","monoatomic ion transport",4.401907471506606,1,0.9333500838758619,0.33548775,"amino acid transport"),
c("GO:0042908","xenobiotic transport",0.24723501249579571,1,0.9503034467272828,0.24847794,"amino acid transport"),
c("GO:0098655","monoatomic cation transmembrane transport",2.923892461650467,1,0.9104411394803963,0.41450816,"amino acid transport"),
c("GO:0006952","defense response",1.0410213876080876,4,0.8625541555957023,-0,"defense response"),
c("GO:0009739","response to gibberellin",0.015047641551607423,1,0.8875402733962205,0.47319047,"defense response"),
c("GO:0010044","response to aluminum ion",0.00414442117004926,1,0.8964838355643986,0.4322357,"defense response"),
c("GO:0098869","cellular oxidant detoxification",0.8170663139318542,2,0.8466371269544272,0.34943531,"defense response"),
c("GO:0007017","microtubule-based process",0.8354933550892162,1,0.9917066152132382,0.01025273,"microtubule-based process"),
c("GO:0009653","anatomical structure morphogenesis",0.7824540774185857,2,0.969406932589886,-0,"anatomical structure morphogenesis"),
c("GO:0009877","nodulation",0.0006951717692939769,1,0.9758601487640157,0.45533059,"anatomical structure morphogenesis"),
c("GO:0055046","microgametogenesis",0.0015100860443036628,1,0.9658151639493533,0.47747485,"anatomical structure morphogenesis"),
c("GO:0009813","flavonoid biosynthetic process",0.015656332622328944,2,0.9411979066705443,0.03204325,"flavonoid biosynthetic process"),
c("GO:0006556","S-adenosylmethionine biosynthetic process",0.05360805457278003,1,0.9428016116093211,0.12361764,"flavonoid biosynthetic process"),
c("GO:0030244","cellulose biosynthetic process",0.03770891075830534,1,0.9329711041803669,0.11447671,"flavonoid biosynthetic process"),
c("GO:0010345","suberin biosynthetic process",0.0015433478514469152,1,0.9237423665884029,0.09640071,"suberin biosynthetic process"),
c("GO:0016310","phosphorylation",7.474579996508841,4,0.9484773228342563,0.08716535,"phosphorylation"),
c("GO:0016567","protein ubiquitination",0.8261966799926771,2,0.9297890717699825,-0,"protein ubiquitination"),
c("GO:0006260","DNA replication",1.444034747678592,1,0.9254773538525394,0.13081923,"protein ubiquitination"),
c("GO:0006412","translation",5.085673767131161,2,0.8731225150522767,0.35814401,"protein ubiquitination"),
c("GO:0006486","protein glycosylation",0.628352124923897,1,0.8976389749369134,0.45949342,"protein ubiquitination"),
c("GO:0006508","proteolysis",5.350747086797883,1,0.9245589314391058,0.46208387,"protein ubiquitination"),
c("GO:0042744","hydrogen peroxide catabolic process",0.1277353179722325,1,0.9346249079263118,0.47893235,"protein ubiquitination"),
c("GO:0046777","protein autophosphorylation",0.0549352006777958,1,0.926755760387361,0.36847714,"protein ubiquitination"),
c("GO:0051603","proteolysis involved in protein catabolic process",0.9636876860000256,1,0.8995835722780419,0.28918282,"protein ubiquitination"),
c("GO:0019375","galactolipid biosynthetic process",0.0015799358393044929,1,0.9301884269657771,0.09615063,"galactolipid biosynthetic process"),
c("GO:0016114","terpenoid biosynthetic process",0.30528351832219996,1,0.905325448217654,0.41553037,"galactolipid biosynthetic process"),
c("GO:0019953","sexual reproduction",0.35416507009992376,2,1,-0,"sexual reproduction"),
c("GO:0032259","methylation",3.1286355155207577,1,0.9807346094481942,0.03567345,"methylation"),
c("GO:0045893","positive regulation of DNA-templated transcription",0.723806859063603,3,0.8384389936430329,-0,"positive regulation of DNA-templated transcription"),
c("GO:0006357","regulation of transcription by RNA polymerase II",1.917576443615649,1,0.880174015900248,0.4674003,"positive regulation of DNA-templated transcription"),
c("GO:0010888","negative regulation of lipid storage",0.0023216741385990235,1,0.9219438450941123,0.41955325,"positive regulation of DNA-templated transcription"),
c("GO:0031047","RNA-mediated gene silencing",0.1393037744966557,1,0.8771826635315785,0.32456012,"positive regulation of DNA-templated transcription"),
c("GO:0040010","positive regulation of growth rate",0.00047897002286283576,1,0.9122786726752258,0.44241796,"positive regulation of DNA-templated transcription"),
c("GO:0042691","positive regulation of crystal cell differentiation",0.00015633049357328667,1,0.914839848767494,0.44286692,"positive regulation of DNA-templated transcription"),
c("GO:0050821","protein stabilization",0.036474897713290676,1,0.9366171809747011,0.17581023,"positive regulation of DNA-templated transcription"),
c("GO:0051726","regulation of cell cycle",0.5059719579017281,1,0.9184944712765444,0.24768465,"positive regulation of DNA-templated transcription"),
c("GO:0060966","regulation of gene silencing by RNA",0.009459657951541007,1,0.9222941851133702,0.26306049,"positive regulation of DNA-templated transcription"),
c("GO:1902290","positive regulation of defense response to oomycetes",0.00011974250571570894,1,0.8890214439396393,0.40937928,"positive regulation of DNA-templated transcription"),
c("GO:2000652","regulation of secondary cell wall biogenesis",0.0007616953835804819,1,0.9461555405419374,0.15137793,"positive regulation of DNA-templated transcription"),
c("GO:0050896","response to stimulus",14.674014998147983,1,1,-0,"response to stimulus"),
c("GO:0051085","chaperone cofactor-dependent protein refolding",0.016697427185912744,1,0.9817854419792842,0.00726942,"chaperone cofactor-dependent protein refolding"),
c("GO:0071555","cell wall organization",0.8657349901438613,7,0.8699903461137262,0,"cell wall organization"),
c("GO:0007005","mitochondrion organization",0.3670107800186479,1,0.9149178070077862,0.45467033,"cell wall organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="MSZ_expanded_deg_any_exp_species_reduced05_BP_TreeMap.pdf", width=16, height=9 ) # width and height are in inches


topGO_data <- read.csv("..\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")
topGO_data <- topGO_data[topGO_data$Ontology == "BP" & topGO_data$weight01 <= 0.05,]
topGO_ReviGO_Intersect <- intersect(topGO_data$Term, revigo.data[,2])

hogs <- read.csv("..\\MSZ_expanded_deg_any_exp_species_hogs.csv")
n_hogs <- length(hogs$X)


# check the tmPlot command documentation for all possible parameters - there are a lot more
t = treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = paste("C4 vs C3 Expanded, DEG in Any Expanded Species, BP, ReviGO 0.5, # HOGs", n_hogs, "Total Sig. GO-Terms", length(topGO_data$GO.ID)),
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

dev.off()

