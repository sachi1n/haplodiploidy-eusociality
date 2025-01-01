#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: This analysis is performed on the full tree
# NOTE: HD.arrhenotoky is used, and the definition of Eusociality used for this specific analysis is from Boomsma and Gawne:https://onlinelibrary.wiley.com/doi/10.1111/brv.12330


# ---- Load libraries ----------------------------------------------------------
library(phytools)

# ----- Set working directory --------------------------------------------------
setwd("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel/")
list.files()


# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdarrhenotoky_euskoos_calib.RDS")
full_dep_xy

# Model fit:
#                 log-likelihood      AIC
#  independent      -98.49477      204.9895
#  dependent        -89.49842      194.9968

# Hypothesis test result:
# likelihood-ratio:  17.9927 
# p-value:  0.00123816 

# Model fitting method used was fitMk

plot(full_dep_xy, lwd.by.rate=TRUE)
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdarrhenotoky_euskoos_calib.RDS")
full_dep_x

#Independent model rate matrix:
#               eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD   -0.000562    0.000009    0.000553    0.000000
#eusocial|HD    0.000000   -0.000553    0.000000    0.000553
#solitary|DD    0.000008    0.000000   -0.000017    0.000009
#solitary|HD    0.000000    0.000008    0.000000   -0.000008

#Dependent (x only) model rate matrix:
#                eusocial|DD eusocial|HD solitary|DD solitary|HD
#eusocial|DD      -9e-06    0.000009     0.0e+00    0.000000
#eusocial|HD       0e+00    0.000000     0.0e+00    0.000000
#solitary|DD       2e-06    0.000000    -1.1e-05    0.000009
#solitary|HD       0e+00    0.000174     0.0e+00   -0.000174

#Model fit:
#                log-likelihood      AIC
#independent      -98.49477       204.9895
#dependent        -89.50875       191.0175

#Hypothesis test result:
#  likelihood-ratio:  17.972 
#  p-value:  0.000125147 

# Model fitting method used was fitMk  

plot(full_dep_x)
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdarrhenotoky_euskoos_calib.RDS")
full_dep_y

#Model fit:
#                log-likelihood      AIC
#independent      -98.49477       204.9895
#dependent        -98.48441       208.9688

#Hypothesis test result:
#  likelihood-ratio:  0.0207085 
# p-value:  0.989699 

# Model fitting method used was fitMk 

plot(full_dep_y)
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#                log(L)     d.f.      AIC        weight
# independent  -98.49477     4     204.9895   0.0008127335
# full_dep_xy  -89.49842     8     194.9968   0.1201802647
# full_dep_x   -89.50875     6     191.0175   0.8788958655
# full_dep_y   -98.48441     6     208.9688   0.0001111363


# Plot the best supported model-------------------------------------------------
plot(full_dep_x, cex.main=1,cex.sub=0.8,cex.traits=0.7,cex.rates=0.7,lwd.by.rate=TRUE,max.lwd=6)

# Save as pdf after this
#-------------------------------------------------------------------------------
