#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: Halictidae is pruned to "force" two origins in this analysis
# NOTE: Strict definition of haplodiploidy and eusociality is used
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
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib_2020.RDS")
full_dep_xy

#Dependent (x & y) model rate matrix:
#               eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD       0e+00      0.0e+00      0.0e+00      0.0e+00
#eusocial|HD       0e+00      0.0e+00      0.0e+00      0.0e+00
#solitary|DD       4e-06      0.0e+00     -1.1e-05      8.0e-06
#solitary|HD       0e+00      7.4e-05      2.0e-05     -9.4e-05

#Model fit:
#                log-likelihood      AIC
#independent      -108.6137       225.2273
#dependent        -105.9959       227.9919

#Hypothesis test result:
#  likelihood-ratio:  5.23541 
# p-value:  0.263984 


#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib_2020.RDS")
full_dep_x

#Dependent (x only) model rate matrix:
#                eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -3.0e-05     8.0e-06     2.2e-05     0.0e+00
#eusocial|HD     1.9e-05    -3.4e-05     0.0e+00     1.6e-05
#solitary|DD     4.0e-06     0.0e+00    -1.1e-05     8.0e-06
#solitary|HD     0.0e+00     7.2e-05     1.9e-05    -9.1e-05

#Model fit:
#  log-likelihood      AIC
#independent      -108.6137 225.2273
#dependent        -106.2875 224.5749

#Hypothesis test result:
#  likelihood-ratio:  4.65238 
#p-value:  0.0976671 


#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib_2020.RDS")
full_dep_y

#Dependent (y only) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD   -0.000535    0.000000     0.000535     0.000000
#eusocial|HD    0.000000   -0.000535     0.000000     0.000535
#solitary|DD    0.000006    0.000000    -0.000014     0.000008
#solitary|HD    0.000000    0.000006     0.000020    -0.000026

#Model fit:
#  log-likelihood      AIC
#independent      -108.6137 225.2273
#dependent        -108.4518 228.9037

#Hypothesis test result:
#  likelihood-ratio:  0.323619 
#p-value:  0.850603 
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#log(L) d.f.      AIC     weight
#independent  -108.6137    4  225.2273 0.35768097
#full_dep_xy  -105.9959    8  227.9919 0.08977878
#full_dep_x   -106.2875    6  224.5749 0.49563137
#full_dep_y   -108.4518    6  228.9037 0.05690887