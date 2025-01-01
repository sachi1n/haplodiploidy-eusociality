#---- Consolidated results from the FitPagel analysis --------------------------
# Halictidae is pruned to "force" two origins in this analysis
# Haplodiploidy is used in a strict sense here (HD.arrhenotoky) and eusociality in a broad sense

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


# ----- Look at the models -----------------------------------------------------

#Read model where x and y are interdependent over the other
full_dep_xy = readRDS("halictidae_pruned_dep_xy_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib.RDS")
full_dep_xy

#Dependent (x & y) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD solitary|HD
#eusocial|DD   -0.000776    0.000000    0.000776    0.000000
#eusocial|HD    0.000000   -0.000221    0.000000    0.000221
#solitary|DD    0.000020    0.000000   -0.000029    0.000009
#solitary|HD    0.000000    0.000300    0.000000   -0.000300

#Model fit:
#              log-likelihood      AIC
#independent      -229.1166      466.2331
#dependent        -217.1696      450.3392

# Hypothesis test result:
#  likelihood-ratio:  23.894 
#  p-value:  8.38808e-05 
#-------------------------------------------------------------------------------

#Read model where x (eusociality) is the dependent trait
full_dep_x = readRDS("halictidae_pruned_dep_x_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib.RDS")
full_dep_x

#Dependent (x only) model rate matrix:
#             eusocial|DD  eusocial|HD  solitary|DD  solitary|HD
#eusocial|DD   -0.000785    0.000009     0.000776     0.000000
#eusocial|HD    0.000000   -0.000221     0.000000     0.000221
#solitary|DD    0.000020    0.000000    -0.000029     0.000009
#solitary|HD    0.000000    0.000300     0.000000    -0.000300

#Model fit:
#               log-likelihood      AIC
#independent      -229.1166       466.2331
#dependent        -217.1942       446.3884

#Hypothesis test result:
#  likelihood-ratio:  23.8447 
#  p-value:  6.64032e-06 
#-------------------------------------------------------------------------------

#Read model where y (haplodiploidy) is the dependent trait
full_dep_y = readRDS("halictidae_pruned_dep_y_fitMk_ARD_fitzjohnpi_haparrhenotoky_calib.RDS")
full_dep_y

#Model fit:
#               log-likelihood      AIC
#independent      -229.1166      466.2331
#dependent        -229.0920      470.1839

#Hypothesis test result:  
#  likelihood-ratio:  0.0492036 
#  p-value:  0.975698 
#-------------------------------------------------------------------------------

# Doing ANOVA to compare the results

anova(full_dep_xy, full_dep_x, full_dep_y)

#              log(L)    d.f.      AIC       weight
#independent -229.1166    4    466.2331   4.308676e-05
#full_dep_xy -217.1696    8    450.3392   1.218068e-01
#full_dep_x  -217.1942    6    446.3884   8.781441e-01
#full_dep_y  -229.0920    6    470.1839   5.976395e-06


plot(full_dep_x,lwd.by.rate=TRUE)
