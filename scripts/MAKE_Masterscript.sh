{
echo "This it the Masterscript with every script used in EX_THE000003. This is made to represent the workflow order, and is not meant to be a working R script"
echo "These are the filtering scripts used to filter the A2TEA RData returned by the A2TEA Workflow in EX_THE000003"
cat ./filter/EX_THE000003_Filter_Source.R ./filter/*
echo "The script to start the A2TEA WebApp to parse throught EX_THE000003_A2TEA_finished.RData"
cat EX_THE000003_A2TEA_Webapp.R
echo "script to calculate go term enrichment with topGO"
cat EX_THE000003_topGO_Source.R EX_THE000003_topGO.R
echo "After filtering and calculating the topThe GO Terms and HOG counts were put into ReviGO for clustering and the significant terms were higlighted in a TreeMap. The following are all the Treemap scripts in EX_THE000003"
cat `find ../filtered_HOGs/ -name "*.R"`
echo "script for checking base mean"
cat EX_THE000003_BaseMean.R
echo "plotting the upset plots for HOGs containing DEG and HOGs species composition as well as HOGs where species are not DEG"
cat EX_THE000003_UpSet.R EX_THE000003_HOGs_Analysis.R EX_THE000003_UpSet_NoDEG.R
echo "Script to caclulate parameters of hog sizes"
cat EX_THE000003_HOG_sizes.R
echo "script to calculate set enrichments"
cat EX_THE000003_Enrichments.R
echo "script to compare the all vs the conserved set"
cat EX_THE000003_All_Conserved_Comparison.R
echo "script to dicern which hypotheses contain a certain HOG"
cat EX_THE000003_Which_Hypotheses.R
echo "Script to calculate GO-Term enrichment which are in the At Least One of each C3 DE sets"
cat EX_THE000003_C4DEG_vs_C3DEG_AtLeastOne_topGO_Source.R
} > EX_THE000003_Masterscript.R
