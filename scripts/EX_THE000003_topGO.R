setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")

#51 HOGs
# MSZ_expanded_deg_any_exp_species_hogs ################################################
MSZ_expanded_deg_any_exp_species_hogs <- read.csv(".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MSZ_expanded_deg_any_exp_species_HOGs.csv")
MSZ_expanded_deg_any_exp_species_hogs <- MSZ_expanded_deg_any_exp_species_hogs$x
MSZ_expanded_deg_any_exp_species_enrich <- test_GO_enrich(MSZ_expanded_deg_any_exp_species_hogs)
plot_GO(MSZ_expanded_deg_any_exp_species_enrich, "Ms, Sb, Zm vs Ta, Os, Conserved, Expanded, DEG in Any Expanded Species\ntopGO Analysis",
        ".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.pdf")
write.csv(MSZ_expanded_deg_any_exp_species_enrich, file = ".\\filtered_HOGs\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies\\MsSbZmvsTaOs_Expanded_DEGinAnyExpandedSpecies_moreThan1Sig_topGO.csv")

#20 HOGs
# AllC4Exp_NoDEGC4_AllC3DEG_hogs ################################################
AllC4Exp_NoDEGC4_AllC3DEG_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs.csv")
AllC4Exp_NoDEGC4_AllC3DEG_hogs <- AllC4Exp_NoDEGC4_AllC3DEG_hogs$x
AllC4Exp_NoDEGC4_AllC3DEG_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3DEG_hogs)
plot_GO(AllC4Exp_NoDEGC4_AllC3DEG_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
#enrichemnt test with background conserved, AllC3DEG, NoDEGC4
AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs <- read.csv(".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs.csv")
AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs <- AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs$x
AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3DEG_hogs, background = AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs)
plot_GO(AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG vs Conserved, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\AllC4Exp_NoDEGC4_AllC3DEG\\AllC4Exp_NoDEGC4_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.csv")


# # AllC4Exp_NoDEGC4_AllC3Up_hogs ################################################
# AllC4Exp_NoDEGC4_AllC3Up_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs.csv")
# AllC4Exp_NoDEGC4_AllC3Up_hogs <- AllC4Exp_NoDEGC4_AllC3Up_hogs$x
# AllC4Exp_NoDEGC4_AllC3Up_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AllC3Up_hogs)
# plot_GO(AllC4Exp_NoDEGC4_AllC3Up_enrich, "All C4 Expanded, No DEG in C4, All C3 Upregulated\ntopGO Analysis",
#         ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs_moreThan1Sig_topGO.pdf")
# write.csv(AllC4Exp_NoDEGC4_AllC3Up_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AllC3Up\\AllC4Exp_NoDEGC4_AllC3Up_hogs_moreThan1Sig_topGO.csv")

# # AllC4Exp_DEGAnyC4_NoDEGC3 ################################################
# AllC4Exp_DEGAnyC4_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_hogs.csv")
# AllC4Exp_DEGAnyC4_NoDEGC3_hogs <- AllC4Exp_DEGAnyC4_NoDEGC3_hogs$x
# AllC4Exp_DEGAnyC4_NoDEGC3_enrich <- test_GO_enrich(AllC4Exp_DEGAnyC4_NoDEGC3_hogs)
# plot_GO(AllC4Exp_DEGAnyC4_NoDEGC3_enrich, "All C4 Expanded, DEG in Any C4, No DEG in C3\ntopGO Analysis",
#         ".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_moreThan1Sig_topGO.pdf")
# write.csv(AllC4Exp_DEGAnyC4_NoDEGC3_enrich, file = ".\\filtered_HOGs\\AllC4Exp_DEGAnyC4_NoDEGC3\\AllC4Exp_DEGAnyC4_NoDEGC3_moreThan1Sig_topGO.csv")

# # ZmExp_AllZmDEG_NoDEGC3_hogs ################################################
# ZmExp_AllZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs.csv")
# ZmExp_AllZmDEG_NoDEGC3_hogs <- ZmExp_AllZmDEG_NoDEGC3_hogs$x
# ZmExp_AllZmDEG_NoDEGC3_enrich <- test_GO_enrich(ZmExp_AllZmDEG_NoDEGC3_hogs)
# plot_GO(ZmExp_AllZmDEG_NoDEGC3_enrich, "Zm Expanded, All Zm DEG, No DEG in C3\ntopGO Analysis",
#         ".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
# write.csv(ZmExp_AllZmDEG_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_AllZmDEG_NoDEGC3\\ZmExp_AllZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")

#48 HOGs
# ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs ################################################
ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs.csv")
ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs <- ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs$x
ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich <- test_GO_enrich(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs)
plot_GO(ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich, "Zm Expanded, Only Any Zm DEG, No DEG in C3\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")
#enrichment test with background: Conserved Only Any Zm DEG, No DEG C3
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs <- read.csv(".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs.csv")
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs <- ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs$x
ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich <- test_GO_enrich(ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs, background=ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_hogs)
plot_GO(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich, "Zm Expanded, Only Any Zm DEG, No DEG in C3 vs Conserved, Only Any Zm DEG, No DEG in C3\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_bg_topGO.pdf")
write.csv(ZmExp_OnlyAnyZmDEG_NoDEGC3_bg_enrich, file = ".\\filtered_hogs\\ZmExp_OnlyAnyZmDEG_NoDEGC3\\ZmExp_OnlyAnyZmDEG_NoDEGC3_hogs_moreThan1Sig_bg_topGO.csv")


#513 HOGs
# ZmExp_NoDEGZm_AllC3DEG_hogs ################################################
ZmExp_NoDEGZm_AllC3DEG_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs.csv")
ZmExp_NoDEGZm_AllC3DEG_hogs <- ZmExp_NoDEGZm_AllC3DEG_hogs$x
ZmExp_NoDEGZm_AllC3DEG_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3DEG_enrich, "Zm Expanded, No Zm DEG, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3DEG_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG DEG
ZmExp_NoDEGZm_AllC3DEG_bg_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs.csv")
ZmExp_NoDEGZm_AllC3DEG_bg_hogs <- ZmExp_NoDEGZm_AllC3DEG_bg_hogs$x
ZmExp_NoDEGZm_AllC3DEG_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_hogs, background=ZmExp_NoDEGZm_AllC3DEG_bg_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3DEG_bg_enrich, "Zm Expanded, No Zm DEG, All C3 DEGregulated vs Conserved, No Zm DEG, All C3 DEG\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\ZmExp_NoDEGZm_AllC3DEG_bg_hogs_moreThan1Sig_bg_topGO.csv")

NoDEGZm_AllC3DEG_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3DEG_bg_hogs)
plot_GO(NoDEGZm_AllC3DEG_bg_enrich, "No Zm DEG, All C3 DEGregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(NoDEGZm_AllC3DEG_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3DEG\\NoDEGZm_AllC3DEG_hogs_moreThan1Sig_bg_topGO.csv")



#204 HOGs
# ZmExp_NoDEGZm_AllC3Up_hogs ################################################
ZmExp_NoDEGZm_AllC3Up_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs.csv")
ZmExp_NoDEGZm_AllC3Up_hogs <- ZmExp_NoDEGZm_AllC3Up_hogs$x
ZmExp_NoDEGZm_AllC3Up_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Up_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Up_enrich, "Zm Expanded, No Zm DEG, All C3 Upregulated\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Up_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG Up
NoDEGZm_AllC3Up_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\NoDEGZm_AllC3Up_hogs.csv")
NoDEGZm_AllC3Up_hogs <- NoDEGZm_AllC3Up_hogs$x
ZmExp_NoDEGZm_AllC3Up_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Up_hogs, background = NoDEGZm_AllC3Up_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Up_bg_enrich, "Zm Expanded, No Zm DEG, All C3 Upregulated vs Conserved, No Zm DEG, All C3 Upregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Up_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Up\\ZmExp_NoDEGZm_AllC3Up_bg_hogs_moreThan1Sig_bg_topGO.csv")

#105 HOGs
# ZmExp_NoDEGZm_AllC3Down_hogs ################################################
ZmExp_NoDEGZm_AllC3Down_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs.csv")
ZmExp_NoDEGZm_AllC3Down_hogs <- ZmExp_NoDEGZm_AllC3Down_hogs$x
ZmExp_NoDEGZm_AllC3Down_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Down_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Down_enrich, "Zm Expanded, No Zm DEG, All C3 Downregulated\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Down_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_hogs_moreThan1Sig_topGO.csv")
#enrichment with No DEG Zm All C3 DEG Down
NoDEGZm_AllC3Down_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\NoDEGZm_AllC3Down_hogs.csv")
NoDEGZm_AllC3Down_hogs <- NoDEGZm_AllC3Down_hogs$x
ZmExp_NoDEGZm_AllC3Down_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AllC3Down_hogs, background = NoDEGZm_AllC3Down_hogs)
plot_GO(ZmExp_NoDEGZm_AllC3Down_bg_enrich, "Zm Expanded, No Zm DEG, All C3 Downregulated vs Conserved, No Zm DEG, All C3 Downregulated\ntopGO Analysis",
        ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_bg_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AllC3Down_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AllC3Down\\ZmExp_NoDEGZm_AllC3Down_bg_hogs_moreThan1Sig_bg_topGO.csv")

