#---- Consolidated results from the FitPagel analysis --------------------------
# Halictidae is pruned to "force" two origins in this analysis
# Haplodiploidy is used in a strict sense here (HD.arrhenotoky) and eusociality is also used in a strict sense (Eusocial.strict)

# ---- Load libraries ----------------------------------------------------------
library(phytools)

# Specify relevant directories for different files------------------------------

data =("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for ASU desktop
data = ("/home/ssures53/haplodiploidy_and_eusociality/data/") # Path for running in SOL

# Output directory for the FitPagel analyses
fitpagel = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/FitPagel") #directory for ASU desktop

fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for ASU desktop
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
setwd(fitpagel)
list.files()
#-------------------------------------------------------------------------------


# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib.RDS")
full_dep_xy


#Dependent (x & y) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD       0e+00     0.000000      0.0e+00     0.000000
#eusocial|HD       0e+00     0.000000      0.0e+00     0.000000
#solitary|DD       2e-06     0.000000     -1.1e-05     0.000009
#solitary|HD       0e+00     0.000174      0.0e+00    -0.000174

#Model fit:
#  log-likelihood      AIC
#independent     -106.08090 220.1618
#dependent        -94.54137 205.0827

#Hypothesis test result:
#  likelihood-ratio:  23.0791 
#p-value:  0.000122103 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib.RDS")
full_dep_x


#Dependent (x only) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD      -9e-06     0.000009      0.0e+00     0.000000
#eusocial|HD       0e+00     0.000000      0.0e+00     0.000000
#solitary|DD       2e-06     0.000000     -1.1e-05     0.000009
#solitary|HD       0e+00     0.000168      0.0e+00    -0.000168

#Model fit:
#  log-likelihood      AIC
#independent     -106.08090 220.1618
#dependent        -94.56735 201.1347

#Hypothesis test result:
#  likelihood-ratio:  23.0271 
#p-value:  9.99372e-06 
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_calib.RDS")
full_dep_y

#Dependent (y only) model rate matrix:
#              eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD     0.0e+00      0.0e+00      0.0e+00      0.0e+00
#eusocial|HD     0.0e+00      0.0e+00      0.0e+00      0.0e+00
#solitary|DD     1.2e-05      0.0e+00     -2.1e-05      9.0e-06
#solitary|HD     0.0e+00      1.2e-05      0.0e+00     -1.2e-05

#Model fit:
#  log-likelihood      AIC
#independent      -106.0809 220.1618
#dependent        -106.0574 224.1148

#Hypothesis test result:
#  likelihood-ratio:  0.0469695 
#p-value:  0.976789 
#-------------------------------------------------------------------------------


# Doing ANOVA to compare the results
anova(full_dep_xy, full_dep_x, full_dep_y)

#              log(L)     d.f.    AIC         weight
#independent -106.08090    4    220.1618    6.483349e-05
#full_dep_xy  -94.54137    8    205.0827    1.219487e-01
#full_dep_x   -94.56735    6    201.1347    8.779774e-01
#full_dep_y  -106.05742    6    224.1148    8.982759e-06

