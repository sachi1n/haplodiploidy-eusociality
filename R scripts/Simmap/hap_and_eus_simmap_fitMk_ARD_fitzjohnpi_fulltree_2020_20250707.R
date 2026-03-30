# Goal: To understand how the evolution of traits Eusociality and Haplodiploidy are correlated across the insect phylogeny
# NOTE: In this script, a simmap is run using the following settings: method = "fitMk", model = "ARD", pi = "fitzjohn"
# NOTE: In this script, haplodiploidy is used in the broad sense (HD.all). Eusociality is also used in a broad sense (Eusocial)
# NOTE: This simmap analysis is for the Chester's recent phylogeny (2020)


# Load necessary packages ------------------------------------------------------

library(phytools) 
library(dplyr)
#-------------------------------------------------------------------------------

# Set directory shortcuts-------------------------------------------------------
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for laptop

# Data directory for the phenotypic dataset and tree files
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for laptop

# Data directory for simmap results
simmap = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Simmap") #directory for ASU desktop
simmap = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") #directory for laptop
#-------------------------------------------------------------------------------


# Set working directory---------------------------------------------------------

setwd(data)
#setwd("C:/Users/tlink/Dropbox/PC/Documents/manuscripts/haplodiploidy eusociality")
#-------------------------------------------------------------------------------


# Load relevant files-----------------------------------------------------------
insect.tree = read.tree("insect_phylogeny_2020_rooted_calibrated.tre")

# Load the phenotype dataset
pheno.df = read.csv("pheno_processed_df_2020_20250626.csv") # This CSV file is a data-wrangled version of pheno_df_16102023.csv
row.names(pheno.df)=pheno.df$name
#-------------------------------------------------------------------------------

# Change the working directory
setwd(simmap)
list.files()

# Run Stochastic character mapping analysis-------------------------------------
# More about this analysis: http://blog.phytools.org/2019/07/stochastic-character-mapping-with.html

pheno.df$ploidy.sociality=interaction(pheno.df$HD.all,pheno.df$Eusocial)
ps.ss =setNames(pheno.df$ploidy.sociality,rownames(pheno.df))
ps.ss.trees=make.simmap(insect.tree,ps.ss,model="ARD",nsim=100,pi="fitzjohn")


# Save the result
saveRDS(ps.ss.trees,"hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_fulltree_2020.RDS")
#-------------------------------------------------------------------------------