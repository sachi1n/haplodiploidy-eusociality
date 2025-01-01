# This is an R script for analysing the results of FitPagel which was run on the full tree after tree calibration
# NOTE: fitPagel was run with the following settings: method = "fitMk", model = "ARD", pi = "fitzjohn"


# Load necessary packages ------------------------------------------------------

library(phytools) 
#-------------------------------------------------------------------------------


# Specify relevant directories for different files------------------------------

# Data directory - contains .csv files of phenotypes and tree files
data = ("/home/ssures53/haplodiploidy_and_eusociality/data/") # Path for running in SOL
data = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data/")

# Output directory for the FitPagel analyses
fitpagel = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel/") # Path for running in SOL
fitpagel = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")

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
dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_calib.RDS") # Both x and y are interdependent on each other
dep_xy
#output:
#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000491    0.000034    0.000457    0.000000
#eusocial|HD    0.000350   -0.000807    0.000000    0.000457
#solitary|DD    0.000058    0.000000   -0.000092    0.000034
#solitary|HD    0.000000    0.000058    0.000350   -0.000408

#Dependent (x & y) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000776    0.000000    0.000776    0.000000
#eusocial|HD    0.000000    0.000000    0.000000    0.000000
#solitary|DD    0.000020    0.000000   -0.000054    0.000034
#solitary|HD    0.000000    0.000643    0.000383   -0.001025

#Model fit:
#            log-likelihood      AIC
#independent      -530.5190 1069.038
#dependent        -494.9995 1005.999

#Hypothesis test result:
#  likelihood-ratio:  71.0391 
#  p-value:  1.36957e-14 

#Model fitting method used was fitMk
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_calib.RDS") # x is dependent on y
dep_x
#output:
#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000491    0.000034    0.000457    0.000000
#eusocial|HD    0.000350   -0.000807    0.000000    0.000457
#solitary|DD    0.000058    0.000000   -0.000092    0.000034
#solitary|HD    0.000000    0.000058    0.000350   -0.000408

#Dependent (x only) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000810    0.000034    0.000776    0.000000
#eusocial|HD    0.000349   -0.000349    0.000000    0.000000
#solitary|DD    0.000020    0.000000   -0.000054    0.000034
#solitary|HD    0.000000    0.000643    0.000349   -0.000993

#Model fit:
#            log-likelihood      AIC
#independent      -530.5190 1069.038
#dependent        -496.1077 1004.215

#Hypothesis test result:
#  likelihood-ratio:  68.8227 
#  p-value:  1.13593e-15 

#Model fitting method used was fitMk 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_calib.RDS") # y is dependent on x
dep_y
#output:
#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000491    0.000034    0.000457    0.000000
#eusocial|HD    0.000350   -0.000807    0.000000    0.000457
#solitary|DD    0.000058    0.000000   -0.000092    0.000034
#solitary|HD    0.000000    0.000058    0.000350   -0.000408

#Dependent (y only) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000457    0.000000    0.000457    0.000000
#eusocial|HD    0.000000   -0.000457    0.000000    0.000457
#solitary|DD    0.000058    0.000000   -0.000092    0.000034
#solitary|HD    0.000000    0.000058    0.000384   -0.000442

#Model fit:
#            log-likelihood      AIC
#independent      -530.5190 1069.038
#dependent        -529.3804 1070.761

#Hypothesis test result:
#  likelihood-ratio:  2.27714 
#  p-value:  0.320276 

#Model fitting method used was fitMk 
 
#-------------------------------------------------------------------------------


# Do ANOVA to compare different models------------------------------------------
anova(dep_x,dep_y,dep_xy)
#output
#              log(L)     d.f.    AIC       weight
#independent  -530.5190    4   1069.038   5.953116e-15
#dep_x        -496.1077    6   1004.215   7.092585e-01
#dep_y        -529.3804    6   1070.761   2.515536e-15
#dep_xy       -494.9995    8   1005.999   2.907415e-01
#-------------------------------------------------------------------------------

# Plot the best supported model ------------------------------------------------
plot(dep_x,lwd.by.rate=TRUE)
