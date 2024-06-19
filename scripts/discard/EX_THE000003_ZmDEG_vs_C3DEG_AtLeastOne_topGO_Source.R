setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")



#52 HOGs
# ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs ################################################
ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs.csv")
ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs <- ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs$x
ZmExp_AtLeastOneZmDEG_NoDEGC3_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs)
plot_GO(ZmExp_AtLeastOneZmDEG_NoDEGC3_enrich, "Zm Expanded, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_AtLeastOneZmDEG_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG, NoDEGC3
# AtLeastOneZmDEG_NoDEGC3_hogs <- read.csv(".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\AtLeastOneZmDEG_NoDEGC3_hogs.csv")
# AtLeastOneZmDEG_NoDEGC3_hogs <- AtLeastOneZmDEG_NoDEGC3_hogs$x
# ZmExp_AtLeastOneZmDEG_NoDEGC3_bg_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_NoDEGC3_hogs, background = AtLeastOneZmDEG_NoDEGC3_hogs)
# plot_GO(ZmExp_AtLeastOneZmDEG_NoDEGC3_bg_enrich, "Zm Expanded, At Least One Zm DE, No C3 DE vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_AtLeastOneZmDEG_NoDEGC3_bg_enrich, file = ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_NoDEGC3\\ZmExp_AtLeastOneZmDEG_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.csv")
# 




#22 HOGs
# ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs ################################################
ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs$x
ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs)
plot_GO(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_enrich, "Zm Expanded, At Least One Zm DE Only Downregulated, No C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG_OnlyDown, NoDEGC3
# AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- read.csv(".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs.csv")
# AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs <- AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs$x
# ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_bg_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs, background = AtLeastOneZmDEG_OnlyDown_NoDEGC3_hogs)
# plot_GO(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_bg_enrich, "Zm Expanded, At Least One Zm DE Only Downregulated, No C3 DE vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_bg_enrich, file = ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyDown_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.csv")





#30 HOGs
# ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs ################################################
ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs$x
ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs)
plot_GO(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_enrich, "Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_enrich, file = ".\\filtered_HOGs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG_OnlyUp, NoDEGC3
# AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- read.csv(".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs.csv")
# AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs <- AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs$x
# ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_bg_enrich <- test_GO_enrich(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs, background = AtLeastOneZmDEG_OnlyUp_NoDEGC3_hogs)
# plot_GO(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_bg_enrich, "Zm Expanded, At Least One Zm DE Only Upregulated, No C3 DE vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_bg_enrich, file = ".\\filtered_hogs\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3\\ZmExp_AtLeastOneZmDEG_OnlyUp_NoDEGC3_bg_hogs_moreThan1Sig_bg_topGO.csv")



#974 HOGs
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs ################################################
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs$x
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs)
plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG, NoDEGC3
# NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\NoDEGZm_AtLeastOneOfEachC3DEG_hogs.csv")
# NoDEGZm_AtLeastOneOfEachC3DEG_hogs <- NoDEGZm_AtLeastOneOfEachC3DEG_hogs$x
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_hogs, background = NoDEGZm_AtLeastOneOfEachC3DEG_hogs)
# plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_bg_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 DE vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_bg_hogs_moreThan1Sig_bg_topGO.csv")


#236 HOGs
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs ################################################
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs$x
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs)
plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 DEG Only Downregualted\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG, NoDEGC3
# NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
# NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs$x
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs, background = NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_hogs)
# plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_bg_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 DE Only Downregulated vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyDown_bg_hogs_moreThan1Sig_bg_topGO.csv")



#345 HOGs
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_OnlyUp_hogs ################################################
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- read.csv(".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs$x
ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs)
plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 Only Upregulated DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs_moreThan1Sig_topGO.pdf")
write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_enrich, file = ".\\filtered_HOGs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs_moreThan1Sig_topGO.csv")

# #enrichemnt test with background conserved, ZmDEG, NoDEGC3
# NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- read.csv(".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
# NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs$x
# ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_bg_enrich <- test_GO_enrich(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs, background = NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_hogs)
# plot_GO(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_bg_enrich, "Zm Expanded, No Zm DE, At Least One Of Each C3 Only Upregulated DE vs Conserved, At Least One Zm DE, No C3 DEG\ntopGO Analysis",
#         ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_bg_hogs_moreThan1Sig_bg_topGO.pdf")
# write.csv(ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_bg_enrich, file = ".\\filtered_hogs\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp\\ZmExp_NoDEGZm_AtLeastOneOfEachC3DEG_OnlyUp_bg_hogs_moreThan1Sig_bg_topGO.csv")


