#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: This analysis is performed on the full tree
# NOTE: HD.arrhenotoky is used, and the definition of Eusociality used for this specific analysis is from Boomsma and Gawne:https://onlinelibrary.wiley.com/doi/10.1111/brv.12330
# NOTE: This analysis is for the Chesters recent phylogeny (2020)

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# Specify relevant directories for different files------------------------------

# Directory for the FitPagel analyses
fitpagel = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel/") # Path for running in SOL
fitpagel = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final data and codes/R scripts/FitPagel/")
fitpagel = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")

setwd(fitpagel)

# Directory for figures
#figures = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/")
#-------------------------------------------------------------------------------

# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdarrhenotoky_EusKoos_calib_2020.RDS")
full_dep_xy

#Independent model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000599     8.0e-06    0.000591    0.000000
#eusocial|HD    0.000019    -6.1e-04    0.000000    0.000591
#solitary|DD    0.000006     0.0e+00   -0.000013    0.000008
#solitary|HD    0.000000     6.0e-06    0.000019   -0.000024

#Dependent (x & y) model rate matrix:
#               eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD       0e+00    0.000000     0.0e+00    0.000000
#eusocial|HD       0e+00   -0.000174     0.0e+00    0.000174
#solitary|DD       4e-06    0.000000    -1.2e-05    0.000008
#solitary|HD       0e+00    0.000044     1.9e-05   -0.000063

#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdarrhenotoky_EusKoos_calib_2020.RDS")
full_dep_x

#Independent model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000599     8.0e-06    0.000591    0.000000
#eusocial|HD    0.000019    -6.1e-04    0.000000    0.000591
#solitary|DD    0.000006     0.0e+00   -0.000013    0.000008
#solitary|HD    0.000000     6.0e-06    0.000019   -0.000024

#Dependent (x only) model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -8.0e-06    0.000008     0.0e+00    0.000000
#eusocial|HD     1.8e-05   -0.000246     0.0e+00    0.000228
#solitary|DD     4.0e-06    0.000000    -1.1e-05    0.000008
#solitary|HD     0.0e+00    0.000041     1.8e-05   -0.000059

#Model fit:
#  log-likelihood      AIC
#independent      -108.0578 224.1155
#dependent        -106.8223 225.6446

#Hypothesis test result:
#  likelihood-ratio:  2.47091 
#  p-value:  0.290703 

#Model fitting method used was fitMk 

#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdarrhenotoky_EusKoos_calib_2020.RDS")
full_dep_y

#Independent model rate matrix:
#  eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000599     8.0e-06    0.000591    0.000000
#eusocial|HD    0.000019    -6.1e-04    0.000000    0.000591
#solitary|DD    0.000006     0.0e+00   -0.000013    0.000008
#solitary|HD    0.000000     6.0e-06    0.000019   -0.000024

#Dependent (y only) model rate matrix:
#  eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000489    0.000000    0.000489    0.000000
#eusocial|HD    0.000000   -0.000489    0.000000    0.000489
#solitary|DD    0.000006    0.000000   -0.000014    0.000008
#solitary|HD    0.000000    0.000006    0.000020   -0.000026

#Model fit:
# log-likelihood      AIC
#independent      -108.0578 224.1155
#dependent        -107.9768 227.9537

#Hypothesis test result:
#  likelihood-ratio:  0.161849 
#p-value:  0.922263 

#Model fitting method used was fitMk 

#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#               log(L)     d.f.    AIC       weight
#independent  -108.0578    4    224.1155   0.59790851
#full_dep_xy  -106.8677    8    229.7354   0.03599969
#full_dep_x   -106.8223    6    225.6446   0.27835319
#full_dep_y   -107.9768    6    227.9537   0.08773862


# Plot the best supported model-------------------------------------------------
plot(full_dep_x, cex.main=1,cex.sub=0.8,cex.traits=0.7,cex.rates=0.7,lwd.by.rate=TRUE,max.lwd=6)

# Save as pdf after this
#-------------------------------------------------------------------------------


# Transition rate from solitary to eusocial in haplodiploids based on the dep_x model  - 0.00004
# Transition rate from solitary to eusocial in diploids based on the dep_x model  - 0.000004