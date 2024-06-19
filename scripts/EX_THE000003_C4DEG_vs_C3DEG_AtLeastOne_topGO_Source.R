setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003\\scripts")

source(".\\EX_THE000003_topGO_Source.R")


setwd("C:\\Users\\jerem\\Desktop\\Studium\\Nutzpflanzen Master\\Masterarbeit\\EX_THE000003")







# 14 HOGs
# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp ################################################
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs.csv")
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs <- AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs$x
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_hogs)
plot_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyUp_moreThan1Sig_topGO.csv")



# 10 HOGs
# AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown ################################################
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- read.csv(".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs.csv")
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs <- AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs$x
AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich <- test_GO_enrich(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_hogs)
plot_GO(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich, "All C4 Expanded, No DEG in C4, All C3 DEG\ntopGO Analysis",
        ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.pdf")
write.csv(AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_enrich, file = ".\\filtered_HOGs\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown\\AllC4Exp_NoDEGC4_AtLeastOneOfEachC3DEG_OnlyDown_moreThan1Sig_topGO.csv")
