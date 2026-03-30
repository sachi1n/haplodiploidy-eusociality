#---- Consolidated results from the FitPagel analysis --------------------------
# Haplodiploidy is used in a strict broad sense here (HD.all) and eusociality in a strict sense (Eus.strict)
# This analysis is performed for the recent Chesters phylogeny (2020)

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# Specify relevant directories for different files------------------------------

# Directory for the FitPagel analyses
fitpagel = ("/home/ssures53/haplodiploidy_and_eusociality/fitpagel/") # Path for running in SOL
fitpagel = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")
fitpagel = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/FitPagel/")

# Directory for figures
#figures = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/")
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
setwd(fitpagel)

# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib_2020.RDS")
full_dep_xy

#Independent model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000553    0.000026    0.000527    0.000000
#eusocial|HD    0.000087   -0.000614    0.000000    0.000527
#solitary|DD    0.000006    0.000000   -0.000032    0.000026
#solitary|HD    0.000000    0.000006    0.000087   -0.000093

#Dependent (x & y) model rate matrix:
#               eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD       0e+00    0.000000     0.0e+00    0.000000
#eusocial|HD       0e+00   -0.000252     0.0e+00    0.000252
#solitary|DD       4e-06    0.000000    -3.0e-05    0.000027
#solitary|HD       0e+00    0.000039     9.1e-05   -0.000130

#Model fit:
#              log-likelihood      AIC
#independent      -224.8479     457.6958
#dependent        -222.5941     461.1883

# Hypothesis test result:
# likelihood-ratio:  4.50757 
# p-value:  0.341651 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
#full_dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib_2020.RDS")
#full_dep_x

#Independent model rate matrix:
#               eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000553    0.000026    0.000527    0.000000
#eusocial|HD    0.000087   -0.000614    0.000000    0.000527
#solitary|DD    0.000006    0.000000   -0.000032    0.000026
#solitary|HD    0.000000    0.000006    0.000087   -0.000093

#Dependent (x only) model rate matrix:
#             eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -2.6e-05    0.000026     0.0e+00    0.000000
#eusocial|HD     8.4e-05   -0.000412     0.0e+00    0.000328
#solitary|DD     4.0e-06    0.000000    -3.0e-05    0.000026
#solitary|HD     0.0e+00    0.000038     8.4e-05   -0.000121

#Model fit:
#  log-likelihood      AIC
#independent      -224.8479 457.6958
#dependent        -223.4092 458.8185

#Hypothesis test result:
#  likelihood-ratio:  2.87735 
#p-value:  0.237242 
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
#full_dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib_2020.RDS")
#full_dep_y

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000553    0.000026    0.000527    0.000000
#eusocial|HD    0.000087   -0.000614    0.000000    0.000527
#solitary|DD    0.000006    0.000000   -0.000032    0.000026
#solitary|HD    0.000000    0.000006    0.000087   -0.000093

#Dependent (y only) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000572    0.000000    0.000572    0.000000
#eusocial|HD    0.000000   -0.000572    0.000000    0.000572
#solitary|DD    0.000006    0.000000   -0.000032    0.000027
#solitary|HD    0.000000    0.000006    0.000094   -0.000100

#Model fit:
#            log-likelihood      AIC
#independent      -224.8479 457.6958
#dependent        -223.8708 459.7417

#Hypothesis test result:
#  likelihood-ratio:  1.95417 
#  p-value:  0.376406 

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)
#              log(L)    d.f.      AIC       weight
#independent -224.8479    4     457.6958   0.47518766
#full_dep_xy -222.5941    8     461.1883   0.08288818
#full_dep_x  -223.4092    6     458.8185   0.27107253
#full_dep_y  -223.8708    6     459.7417   0.17085164

# Plot the best supported model ------------------------------------------------
plot(full_dep_x,lwd.by.rate=TRUE)
