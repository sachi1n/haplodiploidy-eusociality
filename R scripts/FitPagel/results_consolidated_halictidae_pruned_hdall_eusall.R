#---- Consolidated results from the FitPagel analysis --------------------------
# Halictidae is pruned to "force" two origins in this analysis
# Haplodiploidy is used in broad sense (HD.all) and eusociality in a broad sense

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# ----- Set working directory --------------------------------------------------
setwd("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel/")
list.files()

# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_calib.RDS")
full_dep_xy

#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000549    0.000034    0.000515    0.000000
#eusocial|HD    0.000350   -0.000865    0.000000    0.000515
#solitary|DD    0.000037    0.000000   -0.000070    0.000034
#solitary|HD    0.000000    0.000037    0.000350   -0.000387

#Dependent (x & y) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000776    0.000000    0.000776    0.000000
#eusocial|HD    0.000000   -0.000230    0.000000    0.000230
#solitary|DD    0.000020    0.000000   -0.000054    0.000034
#solitary|HD    0.000000    0.000289    0.000384   -0.000674

#Model fit:
#            log-likelihood      AIC
#independent       -403.161   814.3220
#dependent         -390.367   796.7341

#Hypothesis test result:
#  likelihood-ratio:  25.5879 
#  p-value:  3.83129e-05 

#Model fitting method used was fitMk 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_calib.RDS")
full_dep_x

#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000549    0.000034    0.000515    0.000000
#eusocial|HD    0.000350   -0.000865    0.000000    0.000515
#solitary|DD    0.000037    0.000000   -0.000070    0.000034
#solitary|HD    0.000000    0.000037    0.000350   -0.000387

#Dependent (x only) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD    -0.00081    0.000034    0.000776    0.000000
#eusocial|HD     0.00035   -0.000578    0.000000    0.000228
#solitary|DD     0.00002    0.000000   -0.000054    0.000034
#solitary|HD     0.00000    0.000290    0.000350   -0.000640

#Model fit:
#            log-likelihood      AIC
#independent      -403.1610   814.3220
#dependent        -391.4956   794.9912

#Hypothesis test result:
#  likelihood-ratio:  23.3308 
#  p-value:  8.58587e-06 

# Model fitting method used was fitMk 
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_calib.RDS")
full_dep_y

#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000549    0.000034    0.000515    0.000000
#eusocial|HD    0.000350   -0.000865    0.000000    0.000515
#solitary|DD    0.000037    0.000000   -0.000070    0.000034
#solitary|HD    0.000000    0.000037    0.000350   -0.000387

#Dependent (y only) model rate matrix:
#            eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000515    0.000000    0.000515    0.000000
#eusocial|HD    0.000000   -0.000515    0.000000    0.000515
#solitary|DD    0.000037    0.000000   -0.000070    0.000034
#solitary|HD    0.000000    0.000037    0.000385   -0.000422

#Model fit:
#            log-likelihood      AIC
#independent      -403.1610 814.3220
#dependent        -402.0159 816.0318

#Hypothesis test result:
#  likelihood-ratio:  2.29019 
#  p-value:  0.318193 

# Model fitting method used was fitMk
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#               log(L)    d.f.      AIC        weight
#independent  -403.1610    4    814.3220    4.472624e-05
#full_dep_xy  -390.3670    8    796.7341    2.949362e-01
#full_dep_x   -391.4956    6    794.9912    7.050000e-01
#full_dep_y   -402.0159    6    816.0318    1.902314e-05

# Plot the best supported model ------------------------------------------------
plot(full_dep_x,lwd.by.rate=TRUE)
