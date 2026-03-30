# This is an R script for analysing the results of FitPagel which was run on the full tree after tree calibration
# NOTE: fitPagel was run with the following settings: method = "fitMk", model = "ARD", pi = "fitzjohn"
# NOTE: This R script is for analyzing the results of FitPagel for the recent Chester phylogeny (2020)


# Load necessary packages ------------------------------------------------------

library(phytools) 
#-------------------------------------------------------------------------------


# Specify relevant directories for different files------------------------------

# Directory for the FitPagel analyses
fitpagel = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel/") # Path for running in SOL
fitpagel = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")

# Directory for figures
#figures = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/")
#-------------------------------------------------------------------------------


# Set working directory---------------------------------------------------------

setwd(fitpagel)
#setwd("C:/Users/tlink/Dropbox/PC/Documents/manuscripts/haplodiploidy eusociality")
list.files()
#-------------------------------------------------------------------------------


# Load relevant files-----------------------------------------------------------

#-------------------------------------------------------------------------------
dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_calib_2020.RDS") # Both x and y are interdependent on each other
dep_xy

#Independent model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000829    0.000026    0.000803    0.000000
#eusocial|HD    0.000084   -0.000887    0.000000    0.000803
#solitary|DD    0.000067    0.000000   -0.000093    0.000026
#solitary|HD    0.000000    0.000067    0.000084   -0.000151

#Dependent (x & y) model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD     0.0e+00    0.000000     0.0e+00    0.000000
#eusocial|HD     0.0e+00   -0.000670     0.0e+00    0.000670
#solitary|DD     3.3e-05    0.000000    -6.0e-05    0.000027
#solitary|HD     0.0e+00    0.000179     9.3e-05   -0.000272

#Model fit:
#  log-likelihood      AIC
#independent      -460.7482 929.4963
#dependent        -450.2870 916.5740

#Hypothesis test result:
#  likelihood-ratio:  20.9224 
#p-value:  0.000328092 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_calib_2020.RDS") # x is dependent on y
dep_x

#Independent model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000829    0.000026    0.000803    0.000000
#eusocial|HD    0.000084   -0.000887    0.000000    0.000803
#solitary|DD    0.000067    0.000000   -0.000093    0.000026
#solitary|HD    0.000000    0.000067    0.000084   -0.000151

#Dependent (x only) model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -2.6e-05    0.000026     0.0e+00    0.000000
#eusocial|HD     8.5e-05   -0.000664     0.0e+00    0.000579
#solitary|DD     3.4e-05    0.000000    -6.0e-05    0.000026
#solitary|HD     0.0e+00    0.000179     8.5e-05   -0.000263

# Model fit:
# log-likelihood      AIC
# independent      -460.7482 929.4963
# dependent        -451.4739 914.9477

# Hypothesis test result:
# likelihood-ratio:  18.5486 
# p-value:  9.38041e-05 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_calib_2020.RDS") # y is dependent on x
dep_y

# Independent model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000829    0.000026    0.000803    0.000000
#eusocial|HD    0.000084   -0.000887    0.000000    0.000803
#solitary|DD    0.000067    0.000000   -0.000093    0.000026
#solitary|HD    0.000000    0.000067    0.000084   -0.000151

#Dependent (y only) model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000804    0.000000    0.000804    0.000000
#eusocial|HD    0.000000   -0.000804    0.000000    0.000804
#solitary|DD    0.000067    0.000000   -0.000093    0.000027
#solitary|HD    0.000000    0.000067    0.000093   -0.000160

# Model fit:
# log-likelihood      AIC
# independent      -460.7482 929.4963
# dependent        -459.7428 931.4856

# Hypothesis test result:
# likelihood-ratio:  2.0107 
# p-value:  0.365916 

#-------------------------------------------------------------------------------


# Do ANOVA to compare different models------------------------------------------
anova(dep_x,dep_y,dep_xy)
#output
#              log(L)    d.f.      AIC       weight
#independent -460.7482    4   929.4963    0.0004798624
#dep_x       -451.4739    6   914.9477    0.6923183813
#dep_y       -459.7428    6   931.4856    0.0001774788
#dep_xy      -450.2870    8   916.5740    0.3070242775
#-------------------------------------------------------------------------------

# Plot the best supported model ------------------------------------------------
plot(dep_x,lwd.by.rate=TRUE)
