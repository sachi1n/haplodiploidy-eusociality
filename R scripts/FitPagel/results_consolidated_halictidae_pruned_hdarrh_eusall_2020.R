#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: Halictidae is pruned to "force" two origins in this analysis
# NOTE: Broad definition of eusociality and strict definition of haplodiploidy is used
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
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib_2020.RDS")
full_dep_xy


#Dependent (x & y) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD     0.0e+00    0.000000     0.0e+00      0.000000
#eusocial|HD     0.0e+00   -0.000205     0.0e+00      0.000205
#solitary|DD     3.3e-05    0.000000    -4.1e-05      0.000008
#solitary|HD     0.0e+00    0.000078     2.1e-05     -0.000099

#Model fit:
#  log-likelihood      AIC
#independent      -219.6098 447.2197
#dependent        -217.5815 451.1629

#Hypothesis test result:
#  likelihood-ratio:  4.05673 
#p-value:  0.398382 


#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib_2020.RDS")
full_dep_x

#Dependent (x only) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD    -8.0e-06     0.000008     0.0e+00      0.000000
#eusocial|HD     1.9e-05    -0.000266     0.0e+00      0.000248
#solitary|DD     3.3e-05     0.000000    -4.1e-05      0.000008
#solitary|HD     0.0e+00     0.000078     1.9e-05     -0.000097

#Model fit:
#              log-likelihood      AIC
#independent      -219.6098      447.2197
#dependent        -217.7606      447.5211

#Hypothesis test result:
#  likelihood-ratio:  3.69852 
# p-value:  0.157354 


#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib_2020.RDS")
full_dep_y


#Dependent (y only) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD   -0.000224     0.000000     0.000224     0.000000
#eusocial|HD    0.000000    -0.000224     0.000000     0.000224
#solitary|DD    0.000045     0.000000    -0.000052     0.000008
#solitary|HD    0.000000     0.000045     0.000021    -0.000065

#Model fit:
#               log-likelihood      AIC
#independent      -219.6098      447.2197
#dependent        -219.3814      450.7628

#Hypothesis test result:
#likelihood-ratio:  0.45687 
#p-value:  0.795778 
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#               log(L)    d.f.      AIC      weight
#independent  -219.6098    4     447.2197  0.46096393
#full_dep_xy  -217.5815    8     451.1629  0.06417968
#full_dep_x   -217.7606    6     447.5211  0.39646181
#full_dep_y   -219.3814    6     450.7628  0.07839458