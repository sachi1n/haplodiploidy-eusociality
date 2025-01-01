#---- Consolidated results from the FitPagel analysis --------------------------
# Halictidae is pruned to "force" two origins in this analysis
# Haplodiploidy is used in a strict broad sense here (HD.all) and eusociality in a strict sense (Eus.strict)

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# ----- Set working directory --------------------------------------------------
setwd("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel/")
list.files()

# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib.RDS")
full_dep_xy

#Model fit:
#              log-likelihood      AIC
#independent      -272.5543      553.1086
#dependent        -262.7932      541.5864

# Hypothesis test result:
#  likelihood-ratio:  19.52225 
#  p-value:  0.000620 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib.RDS")
full_dep_x

#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000587    0.000034    0.000553    0.000000
#eusocial|HD    0.000350   -0.000903    0.000000    0.000553
#solitary|DD    0.000008    0.000000   -0.000042    0.000034
#solitary|HD    0.000000    0.000008    0.000350   -0.000358

#Dependent (x only) model rate matrix:
#               eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD    -3.4e-05     0.000034      0.0e+00     0.000000
#eusocial|HD     3.5e-04    -0.000350      0.0e+00     0.000000
#solitary|DD     2.0e-06     0.000000     -3.6e-05     0.000034
#solitary|HD     0.0e+00     0.000169      3.5e-04    -0.000518

#Model fit:
#               log-likelihood      AIC
#independent      -272.5543       553.1086
#dependent        -263.1942       539.4545

#Hypothesis test result:
#  likelihood-ratio:  17.65413 
#  p-value:  0.0001467083 
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdall_EusKoos_calib.RDS")
full_dep_y

#Model fit:
#               log-likelihood      AIC
#independent      -272.5543      553.1086
#dependent        -271.5884      555.1767

#Hypothesis test result:  
#  likelihood-ratio:  1.931897 
#  p-value:  0.3806221 
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#               log(L)    d.f.      AIC        weight
#independent  -272.5543    4     553.1086   0.0008054511
#full_dep_xy  -262.7932    8     541.5864   0.2558964745
#full_dep_x   -263.7272    6     539.4545   0.7430116854
#full_dep_y   -271.5884    6     555.1767   0.0002863889

# Plot the best supported model ------------------------------------------------
plot(full_dep_x,lwd.by.rate=TRUE)
