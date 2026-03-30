#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: Halictidae is pruned to "force" two origins in this analysis
# NOTE: Broad definition of both eusociality and haplodiploidy is used
# NOTE: This analysis is for the Chesters recent phylogeny (2020)

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# Specify relevant directories for different files------------------------------

# Directory for the FitPagel analyses
fitpagel = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel/") # Path for running in SOL
fitpagel = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final data and codes/R scripts/FitPagel/")
fitpagel = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel_halictidae_pruned/")

setwd(fitpagel)

# Directory for figures
#figures = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/")
#-------------------------------------------------------------------------------


# ----- Look at the models -----------------------------------------------------

# Read model where x and y are interdependent over the other
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_calib_2020.RDS")
full_dep_xy

# Dependent (x & y) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD     0.0e+00     0.000000      0.0e+00     0.000000
#eusocial|HD     0.0e+00    -0.000248      0.0e+00     0.000248
#solitary|DD     3.3e-05     0.000000     -6.0e-05     0.000027
#solitary|HD     0.0e+00     0.000078      9.5e-05    -0.000172

#Model fit:
#  log-likelihood      AIC
#independent      -336.2405 680.4811
#dependent        -333.4023 682.8047

#Hypothesis test result:
#  likelihood-ratio:  5.67645 
# p-value:  0.22465 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_calib_2020.RDS")
full_dep_x

#Dependent (x only) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD    -2.6e-05     0.000026     0.0e+00     0.000000
#eusocial|HD     8.5e-05    -0.000338     0.0e+00     0.000252
#solitary|DD     3.3e-05     0.000000    -6.0e-05     0.000026
#solitary|HD     0.0e+00     0.000078     8.5e-05    -0.000163

#Model fit:
#  log-likelihood      AIC
#independent      -336.2405 680.4811
#dependent        -334.4241 680.8482

#Hypothesis test result:
#  likelihood-ratio:  3.63294 
#p-value:  0.162599 

#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_calib_2020.RDS")
full_dep_y

#Dependent (y only) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD   -0.000223    0.000000     0.000223     0.000000
#eusocial|HD    0.000000   -0.000223     0.000000     0.000223
#solitary|DD    0.000045    0.000000    -0.000071     0.000027
#solitary|HD    0.000000    0.000045     0.000095    -0.000139

#Model fit:
#  log-likelihood      AIC
#independent      -336.2405 680.4811
#dependent        -335.2179 682.4358

#Hypothesis test result:
#  likelihood-ratio:  2.0453 
#p-value:  0.35964 
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#              log(L)   d.f.      AIC     weight
#independent -336.2405    4    680.4811  0.3965795
#full_dep_xy -333.4023    8    682.8047  0.1241014
#full_dep_x  -334.4241    6    680.8482  0.3300834
#full_dep_y  -335.2179    6    682.4358  0.1492357


