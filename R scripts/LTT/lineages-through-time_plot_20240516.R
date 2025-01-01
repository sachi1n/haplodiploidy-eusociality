################################################################################
# GOAL: To create a lineage-through-time (LTT) plot for an insect phylogeny 
################################################################################

# Load necessary packages ------------------------------------------------------
library(phytools)   # For phylogenetic tree manipulation and LTT plotting
library(tidyverse)  # For data manipulation and plotting utilities
#-------------------------------------------------------------------------------


# Specify relevant directories for different files------------------------------
# NOTE: You can comment/uncomment the appropriate lines depending on whether
# you are using a desktop or a laptop. The same applies for all directory paths.

# Data directory containing the insect phylogeny
tree = ("C:/Users/ssures53/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data/Old phylogenies and datasets/") # directory in desktop
tree = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data/") # directory in laptop

# Data directory - contains RDS files of simmap
data = ("C:/Users/ssures53/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap") # directory in desktop
data = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap")  # directory in laptop

# Directory for the output LTT files
output = ("C:/Users/ssures53/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/LTT") # directory in desktop
output = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/LTT")  # directory in laptop

# Directory for figures
fig = ("C:/Users/ssures53/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/") # directory in desktop
fig = ("C:/Users/sachi/OneDrive - Texas Tech University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/")  # directory in laptop
#-------------------------------------------------------------------------------


# Set working directory---------------------------------------------------------
# The following lines are commented out for reference; uncomment and adapt
# if you need to change your working directory on the fly.

# setwd(data)
# list.files()
# setwd("C:/Users/tlink/Dropbox/PC/Documents/manuscripts/haplodiploidy eusociality")
#-------------------------------------------------------------------------------


# Load relevant files-----------------------------------------------------------
# Here, we would typically load the RDS file(s) from the previous simmap runs
# or from the final LTT computations if they have been done elsewhere.

# insect.maps = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree.RDS")
#-------------------------------------------------------------------------------


# Compute LTTs------------------------------------------------------------------
# Typically, we'd compute LTT from the simmap object:
# insect.ltts <- ltt(insect.maps)
# However, if this step was already done previously, we can just load the results.

# Save the LTT file (after computing above)
# setwd(output)
# saveRDS(insect.ltts, file = "insect_ltt_updated.rds")
#-------------------------------------------------------------------------------


# Load the LTT file-------------------------------------------------------------
setwd(output)  
# Switch to the 'output' directory to read the previously saved LTT results
insect.ltts <- readRDS("insect_ltt_updated.rds")

# Switch to 'data' directory to read the simmap object (alternative approach)
setwd(data)
insects.maps = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree.RDS")

# Switch to 'tree' directory to read the phylogeny
setwd(tree)
insect.tree = read.tree("insect_phylogeny_2017_collapsed_calibrated.tre")
#-------------------------------------------------------------------------------


# Plot the LTT for all mapped states combined------------------------------------
# The built-in plot method for LTT objects from phytools
plot(insect.ltts,
     show.total=FALSE,    # Do not plot the total lineage count as a separate line
     bty="n",             # Remove the box around the plot
     cex.axis=0.8,        # Scale axis text
     las=1,               # Rotate y-axis labels for readability
     xlab="millions of years (from the root)",  # Label for x-axis
     axes=FALSE)          # Suppress default axes for custom axis usage

# Manually adding the x-axis at specified intervals
axis(1, 
     at=round(seq(0, max(nodeHeights(insect.tree)), length.out=4), 1), 
     cex.axis=0.8)

# Manually adding the y-axis 
axis(2, las=1, cex.axis=0.8)

# Restrict any drawing to the region bounded by data range on x-axis and 
# the number of tips (Ntip) on y-axis
clip(0, max(nodeHeights(insect.tree)), 0, Ntip(insect.tree))

# Add a grid to the plot
grid()
#-------------------------------------------------------------------------------


# Example of plotting state-specific LTT as a cumulative area chart-------------
# Extract the LTT data from the first element of insect.ltts
obj <- insect.ltts[[1]]

# Calculate the cumulative lineage-through-time for each state
# The function 'apply' loops over each row, summing across states 
# in order to get the cumulative count at each time step.
# 'c(0, x[1:4])' is a quick way to prepend a zero before the partial sums.
cumul_ltt <- t(
  apply(obj$ltt, 1, function(x) cumsum(c(0, x[1:4])))
)

# Plot an empty chart, then we’ll fill it with polygons for the cumulative areas
plot(NA,
     xlim = range(obj$times),          # x-axis range from earliest to latest time
     ylim = range(cumul_ltt),          # y-axis range covers all cumulative counts
     xlab = "time",                    # Label for the x-axis
     ylab = "cumulative lineages",     # Label for the y-axis
     bty = "n",                        # No box around the plot
     las = 1)                          # Make y-axis labels horizontal

# Colors for different states to be stacked cumulatively. Adjust to taste.
cols <- c("black","#b1d9e1", "#005c77", "#f4c2c2", "#970000")

# Add each state’s area incrementally as a polygon
for(i in 2:ncol(cumul_ltt)) {
  polygon(
    x = c(obj$times, obj$times[length(obj$times):1]), 
    # This forms the 'top' and 'bottom' x-coordinates for the polygon
    y = c(cumul_ltt[, i-1], cumul_ltt[nrow(cumul_ltt):1, i]), 
    # This forms the 'top' and 'bottom' y-coordinates 
    col = cols[i-1],       # Fill color
    border = "transparent" # No border lines for the polygons
  )
}
#-------------------------------------------------------------------------------


# Plot the proportion of lineages in each state over time ----------------------
# Convert LTT matrix to data frame for easier manipulation
pdat = data.frame(obj$ltt)

# Reorder columns if necessary (this line picks columns and reorders them)
# colnames(pdat)[c(2,1,4,3,5)] is an example — adjust to your actual data structure
pdat = pdat[colnames(pdat)[c(2,1,4,3,5)]]

# Compute the cumulative proportion of lineages across states at each time
prop_ltt <- t(
  apply(pdat, 1, function(x) cumsum(c(0, x[1:4])) / x[5])
)
# The [5]-th column (x[5]) presumably contains the total number of lineages
# at each time point, so we divide the cumulative sums by that total
# to get proportions.

# Define colors for each state (matching the columns). The names (e.g. colnames(prop_ltt)[1:5])
# ensure each column color is properly matched.
cols <- setNames(
  c("black","#b1d9e1", "#005c77", "#f4c2c2", "#970000"),
  colnames(prop_ltt)[1:5]
)

# Initialize the plot with no data (type=NA) to set up axes
plot(NA,
     xlim = range(obj$times),
     ylim = range(prop_ltt),
     xlab = "time",
     ylab = "fraction of lineages",
     bty  = "n",
     las  = 1)

# Add stacked polygons to represent the fraction of each state over time
for(i in 2:ncol(prop_ltt)) {
  polygon(
    x = c(obj$times, obj$times[length(obj$times):1]),
    y = c(prop_ltt[, i-1], prop_ltt[nrow(prop_ltt):1, i]),
    col = cols[i],       # Use the matching color
    border = "transparent"
  )
}
#-------------------------------------------------------------------------------
