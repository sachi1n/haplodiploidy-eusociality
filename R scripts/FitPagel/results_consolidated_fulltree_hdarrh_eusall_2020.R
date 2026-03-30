# This is an R script for analysing the results of FitPagel which was run on the full tree after tree calibration
# NOTE: fitPagel was run with the following settings: method = "fitMk", model = "ARD", pi = "fitzjohn"
# NOTE: This R script is for analyzing the results of FitPagel from Chesters recent phylogeny (2020)


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
dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib_2020.RDS") # Both x and y are interdependent on each other
dep_xy

#Independent model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000644    0.000008    0.000636    0.000000
#eusocial|HD    0.000018   -0.000655    0.000000    0.000636
#solitary|DD    0.000068    0.000000   -0.000076    0.000008
#solitary|HD    0.000000    0.000068    0.000018   -0.000086

#Dependent (x & y) model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD     0.0e+00    0.000000     0.0e+00    0.000000
#eusocial|HD     0.0e+00   -0.000489     0.0e+00    0.000489
#solitary|DD     3.3e-05    0.000000    -4.1e-05    0.000008
#solitary|HD     0.0e+00    0.000186     2.0e-05   -0.000206

#Model fit:
#  log-likelihood      AIC
#independent      -344.1890 696.3781
#dependent        -334.7321 685.4641

#Hypothesis test result:
#  likelihood-ratio:  18.914 
# p-value:  0.000817126 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib_2020.RDS") # x is dependent on y
dep_x

#Independent model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000644    0.000008    0.000636    0.000000
#eusocial|HD    0.000018   -0.000655    0.000000    0.000636
#solitary|DD    0.000068    0.000000   -0.000076    0.000008
#solitary|HD    0.000000    0.000068    0.000018   -0.000086

#Dependent (x only) model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -8.0e-06    0.000008     0.0e+00    0.000000
#eusocial|HD     1.8e-05   -0.000507     0.0e+00    0.000488
#solitary|DD     3.4e-05    0.000000    -4.2e-05    0.000008
#solitary|HD     0.0e+00    0.000182     1.8e-05   -0.000200

#Model fit:
#  log-likelihood      AIC
#independent      -344.1890 696.3781
#dependent        -334.9653 681.9306

#Hypothesis test result:
# likelihood-ratio:  18.4475 
# p-value:  9.86693e-05 

# Model fitting method used was fitMk
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib_2020.RDS") # y is dependent on x
dep_y

#Independent model rate matrix:
#              eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000644    0.000008    0.000636    0.000000
#eusocial|HD    0.000018   -0.000655    0.000000    0.000636
#solitary|DD    0.000068    0.000000   -0.000076    0.000008
#solitary|HD    0.000000    0.000068    0.000018   -0.000086

#Dependent (y only) model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000803    0.000000    0.000803    0.000000
#eusocial|HD    0.000000   -0.000803    0.000000    0.000803
#solitary|DD    0.000067    0.000000   -0.000074    0.000008
#solitary|HD    0.000000    0.000067    0.000020   -0.000087

#Model fit:
#  log-likelihood      AIC
#independent      -344.1890 696.3781
#dependent        -343.7912 699.5823

#Hypothesis test result:
#  likelihood-ratio:  0.79576 
#p-value:  0.671743 
#-------------------------------------------------------------------------------


# Do ANOVA to compare different models------------------------------------------
anova(dep_x,dep_y,dep_xy)
#output
#               log(L)   d.f.      AIC       weight
#independent -344.1890    4    696.3781   0.0006222015
#dep_x       -334.9653    6    681.9306   0.8534143584
#dep_y       -343.7912    6    699.5823   0.0001253543
#dep_xy      -334.7321    8    685.4641   0.1458380858
#-------------------------------------------------------------------------------

# Transition rate from solitary to eusocial in haplodiploids based on the dep_x model (least AIC value) - 0.0002
# Transition rate from solitary to eusocial in diploids based on the dep_x model (least AIC value) - 0.00003
