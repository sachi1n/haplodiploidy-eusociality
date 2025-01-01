#---- Consolidated results from the FitPagel analysis --------------------------
# NOTE: This analysis is performed on the full tree

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# ----- Set working directory --------------------------------------------------
setwd("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel/")
list.files()


# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("fulltree_dep_xy_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib.RDS")
full_dep_xy

#Model fit:
#                 log-likelihood      AIC
#  independent      -356.4595      720.9189
#  dependent        -321.4300      658.8599

#Hypothesis test result:
#  likelihood-ratio:  70.059 
#  p-value:  2.20562e-14  

plot(full_dep_xy, lwd.by.rate=TRUE)
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("fulltree_dep_x_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib.RDS")
full_dep_x

#Model fit:
#                 log-likelihood      AIC
# independent      -356.4595        720.9189
# dependent        -321.4546        654.9091

#Hypothesis test result:
#  likelihood-ratio:  70.0098 
#  p-value:  6.27425e-16 

plot(full_dep_x)
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("fulltree_dep_y_fitMk_ARD_fitzjohnpi_hdarrhenotoky_calib.RDS")
full_dep_y

#Model fit:
#             log-likelihood      AIC
#independent      -356.4595     720.9189
#dependent        -356.4349     724.8698

#Hypothesis test result:
#  likelihood-ratio:  0.0491746 
#  p-value:  0.975713

plot(full_dep_y)
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#              log(L)    d.f.      AIC        weight
#independent -356.4595    4     720.9189   4.071350e-15
#full_dep_xy -321.4300    8     658.8599   1.218117e-01
#full_dep_x  -321.4546    6     654.9091   8.781883e-01
#full_dep_y  -356.4349    6     724.8698   5.647127e-16
