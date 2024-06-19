setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\filtered_HOGs")

Zm_C4DEG <- read.csv(".\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_unique_GO.csv")
Zm_C3DEG <- read.csv(".\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_unique_GO.csv")
Zm_intersect_GOID <- intersect(Zm_C3DEG$GOID, Zm_C4DEG$GOID)

Zm_intersect <- Zm_C4DEG[Zm_C4DEG$GOID %in% Zm_intersect_GOID,]

Zm_C4DEG_unique <- Zm_C4DEG[!Zm_C4DEG$GOID %in% Zm_intersect_GOID,]
Zm_C3DEG_unique <- Zm_C3DEG[!Zm_C3DEG$GOID %in% Zm_intersect_GOID,]


All_C4DEG <- read.csv(".\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_unique_GO.csv")
All_C3DEG <- read.csv(".\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_unique_GO.csv")
All_intersect_GOID <- intersect(All_C3DEG$GOID, All_C4DEG$GOID)

All_intersect <- All_C4DEG[All_C4DEG$GOID %in% All_intersect_GOID,]

All_C4DEG_unique <- All_C4DEG[!All_C4DEG$GOID %in% All_intersect_GOID,]
All_C3DEG_unique <- All_C3DEG[!All_C3DEG$GOID %in% All_intersect_GOID,]
