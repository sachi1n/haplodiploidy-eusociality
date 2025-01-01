# Goal: To understand how the evolution of traits Eusociality and Haplodiploidy are correlated across the insect phylogeny
# NOTE: This R script is a trimmed version of another script (hap_and_eus_20230209.R). This specific version is made to streamline the analysis 
# NOTE: This script is also optimized for running in SOL HPC at ASU. 
# NOTE: In this script, fitPagel is run with the following settings: method = "fitMk", model = "ARD", pi = "fitzjohn"
# NOTE: In this script, halictidae is pruned to "force" two independent origins of eusocilality
# NOTE: In this script, haplodiploidy is used in the strict sense (HD.arrhenotoky) and eusociality is used in strcit sense 

# Load necessary packages ------------------------------------------------------

library(phytools) 
library(geiger)
#-------------------------------------------------------------------------------

# Specify relevant directories for different files------------------------------

# Data directory - contains .csv files of phenotypes and tree files
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for ASU desktop
data = ("/home/ssures53/haplodiploidy_and_eusociality/data/") # Path for running in SOL

# Output directory for the FitPagel analyses
fitpagel_halictidae_pruned = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel") #directory for ASU desktop
fitpagel_halictidae_pruned = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel_halictidae_pruned/") # Path for running in SOL

# Output directory ofr figures
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for ASU desktop
#-------------------------------------------------------------------------------


# Set working directory---------------------------------------------------------

setwd(data)
#setwd("C:/Users/tlink/Dropbox/PC/Documents/manuscripts/haplodiploidy eusociality")
#-------------------------------------------------------------------------------


# Load relevant files-----------------------------------------------------------
insect.tree = read.tree("halictidae_pruned_insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset
pheno.df = read.csv("halictidae_pruned_processed_pheno_df.csv") # This CSV file is a data-wrangled version of pheno_df_16102023.csv
row.names(pheno.df)=pheno.df$name
#-------------------------------------------------------------------------------

# Change the working directory
setwd(fitpagel_halictidae_pruned)

# Run FitPagel analysis---------------------------------------------------------
# Refer to this link to learn more about this analysis on: http://phytools.org/mexico2018/ex/8/Pagel94-method.html 
# Refer this too: http://blog.phytools.org/2016/05/fix-to-peculiar-quirk-of-fitpagel.html

x<-setNames(pheno.df$Eusocial.strict,rownames(pheno.df))
y<-setNames(pheno.df$HD.arrhenotoky,rownames(pheno.df))

fit_dep_x_calib<-fitPagel(insect.tree,x,y,dep.var="x", method = "fitMk", model = "ARD", pi = "fitzjohn") # x is dependent on y

# Save the result
saveRDS(fit_dep_x_calib,"halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib.RDS")
