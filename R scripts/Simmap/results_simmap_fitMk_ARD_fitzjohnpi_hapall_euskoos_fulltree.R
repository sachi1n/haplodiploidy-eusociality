################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis
#       using a full insect phylogeny, focusing on transitions to eusociality.
#
# NOTE 1: "Haplodiploidy" here is broadly defined; "Eusociality" is defined strictly.
# NOTE 2: Adjust any directory paths and file names to match your local setup.
################################################################################


# ---- Load libraries ----------------------------------------------------------
library(phytools)    # For phylogenetic tree manipulation and simmap objects
library(ggplot2)     # For plotting
library(RColorBrewer)# For color palettes
library(viridis)     # For color palettes
library(patchwork)   # For combining multiple ggplot2 plots
library(gridExtra)   # For arranging multiple grid-based plots
library(ggridges)    # For ridgeline plots
library(svglite)     # For saving plots in SVG format
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Choose the directory paths that apply to your environment (comment/uncomment 
# based on whether you're on a desktop or a laptop).

# Directory for figure output
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") # ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures")   # Laptop

# Directory for data files (phylogeny, phenotype, etc.)
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") # ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data")   # Laptop

# Directory for Simmap results
res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") # ASU desktop
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results")   # Laptop
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
# Set working directory to where Simmap results are stored
setwd(res)
list.files()  # Check what files are available in this directory

# Optionally load or modify the simmap object (examples shown below).
# simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hdall_euskoos_fulltree.RDS")
# simmap.trees=drop.tip(simmap.trees,"Paraponera_clavata")
# saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hdall_euskoos_fulltree_paraponera_removed.RDS")

# Load the updated simmap object (with a tip "Paraponera_clavata" removed).
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hdall_euskoos_fulltree_paraponera_removed.RDS")
#-------------------------------------------------------------------------------


# ----- Set working directory for data files -----------------------------------
setwd(data)
list.files()  # Inspect the files in the data directory

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset (pre-processed CSV)
pheno.df = read.csv("pheno_processed_df_20231023.csv")
row.names(pheno.df)=pheno.df$name  # Use the 'name' column as row names for easier subsetting
#-------------------------------------------------------------------------------


# ------ Identify clade members from different orders/families -----------------
# We’ll use functions from caper and phangorn to detect common ancestors
# and identify all tips (taxa) within certain clades.

# Halictidae (some are eusocial bees)
halictids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Halictidae" & Eusocial == "eusocial")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Thysanoptera (thrips)
thrips = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, order == "Thysanoptera")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Aphididae (aphids)
aphids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Aphididae")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Blattodea (includes termites)
termites = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, order == "Blattodea")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Formicidae (ants), removing one tip manually if needed
ants = caper::clade.members(
  phangorn::Ancestors(insect.tree, as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Formicidae")))), type = "all")[2], 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = FALSE
)
ants = ants[ants != "Paraponera_clavata"]  # Remove a specific tip if desired

# Check which ant tips are not labeled as "eusocial"
subset(pheno.df, row.names(pheno.df) %in% ants & Eusocial != "eusocial")

# Platypodine ambrosia beetles
platypodine = c("Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus",
                "Cylindropalpus","Dendroplatypus","Dinoplatypus","Doliopygus","Epiplatypus",
                "Euplatypus","Megaplatypus","Mesoplatypus","Myoplatypus","Neotrachyostus",
                "Oxoplatypus","Pereioplatypus","Peroplatypus","Platyphysus","Platypus",
                "Teloplatypus","Trachyostus","Treptoplatypus","Triozastus")
platypodines = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, genus %in% platypodine)))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Vespidae (vespid wasps, e.g. paper wasps)
vespids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Vespidae")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Apidae (honey bees, bumble bees, stingless bees, etc.)
apids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Apidae")))), 
  insect.tree, 
  tip.labels = TRUE, 
  include.nodes = TRUE
)$tips

# Identify the haplodiploid/diploid sets (HD/DD)
HD.clade = row.names(subset(pheno.df, HD.all == "HD"))
DD.clade = row.names(subset(pheno.df, HD.all == "DD"))

# Create a list of clades to compare
clades = list(
  HD.clade, 
  DD.clade, 
  vespids, 
  platypodines, 
  ants, 
  termites, 
  aphids, 
  thrips, 
  halictids, 
  apids
)
names(clades) = c(
  "all haplodiploids","all diploids",
  "vespids","platypodines","ants",
  "termites","aphids","thrips","halictids","apids"
)

# Optionally reorder or subset the list to specific clades
clades = clades[c(1,2,3,5,6,10)]  # picking certain indices for analysis

# Initialize an empty data frame to collect results
clades.df = data.frame(matrix(nrow=0, ncol=3))
colnames(clades.df) = c("eusocial","clade","ploidy")

# Loop over each clade, subset the simmap trees, and describe the transitions
# (simulate or describe states) to gather a data frame of counts
for(c in 1:length(clades)){
  # Keep only the tips that are in this clade from the simmap object
  focal.trees = keep.tip(simmap.trees, clades[[c]])
  
  # Summarize the simmap for these focal trees
  obj = describe.simmap(focal.trees)
  
  # Extract a data frame of counts by state
  df = data.frame(obj$count)
  
  # Determine whether the clade is haplodiploid or diploid from the phenotype
  ploidy = unique(pheno.df[clades[[c]], ]$HD.all)
  
  # Subset to relevant columns (e.g. "haplodiploid.solitary", "haplodiploid.eusocial", etc.)
  df = df[c(paste(ploidy,"solitary", ploidy,"eusocial", sep="."))]
  
  # Add clade and ploidy info
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Name the columns to be consistent across all appended data
  colnames(df) = c("eusocial","clade","ploidy")
  
  # Append the newly created data frame to the master clades.df
  clades.df = rbind(clades.df, df)
}


# ---- Custom theme for ggplot2 ------------------------------------------------
custom_theme <- theme(
  strip.background = element_blank(),
  strip.placement = "outside",
  strip.text = element_text(family = "Arial", colour = "black", face = "italic", size = 14),
  panel.background = element_rect(fill = NA),
  # panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "bottom",
  legend.title = element_text(family = "Arial", colour = "black", size = 14, face = "italic"),
  legend.text = element_text(family = "Arial", colour = "black", size = 16),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(color = NA, fill = NA, linewidth = NA),
  axis.title = element_text(family = "Arial", colour = "black", face = "bold", size = 22),
  axis.title.y = element_text(family = "Arial", colour = "black", face = "bold", size = 22, vjust = 2),
  axis.text = element_text(family = "Arial", colour = "black", size = 20),
  axis.ticks = element_line(colour = "black", linewidth = 1),
  plot.background = element_rect(fill = "white"),
  plot.title = element_text(
    family = "Arial", colour = "black", face = "bold", size = 24, 
    hjust = 0.5, margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit = "cm")
  )
)


# ---- Simple demonstration ridgeline plot -------------------------------------
# This example shows the original ridgeline approach on the aggregated data
ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) + 
  geom_density_ridges(
    stat = "binline", 
    bins = 20, 
    scale = 0.95, 
    draw_baseline = FALSE
  ) +
  theme_ridges()


# ---- Adjust factor levels for consistent ordering ----------------------------
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c("all diploids","termites","all haplodiploids","apids","vespids","ants"))
)
clades.df$ploidy[clades.df$clade == "all haplodiploids"] = "HD"


# ---- Add more clades to be consistent across definitions ---------------------
# Here we artificially add some new rows (e.g., aphids, platypodines, etc.)
# for demonstration or uniform coverage in the final plot.
new_clades <- data.frame(
  clade   = rep(c("aphids", "platypodines", "halictids", "thrips"), each = 100),
  eusocial= rep(0, 4 * 100),   # e.g., zero transitions or dummy data
  ploidy  = rep(c("DD", "DD", "HD", "HD"), each = 100)
)
# Combine the new data with the existing clades.df
clades.df <- rbind(clades.df, new_clades)

# Update factor levels again, to accommodate the new clades
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c("all diploids","aphids","platypodines","termites",
                 "all haplodiploids","halictids","apids","vespids",
                 "ants","thrips"))
)

# Determine number of bins for binning the 'eusocial' variable in ridgelines
bins = max(clades.df$eusocial) - min(clades.df$eusocial)

# Remove any rows with missing values
clades.df = na.omit(clades.df)


# ---- Final ridgeline plot of eusocial origins across multiple clades --------
p1 = ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density),
                           fill = ploidy)) +
  # Use a binning ridgeline
  geom_density_ridges(
    stat = "binline", 
    binwidth = 1,  
    bins = bins, 
    scale = 0.7, 
    binwidth = 0.5,
    draw_baseline = TRUE, 
    color = "grey70"
  ) +
  # Adjust x-axis breaks and limits
  scale_x_continuous(
    breaks = seq(0, 20, by = 2),
    name = "Origins of eusociality",
    expand = c(0, 0)
  ) +
  # Adjust y-axis to have no extra expansion
  scale_y_discrete(
    expand = c(0, 0),
    name = "Insect clade"
  ) +
  # Constrain x-limits 
  coord_cartesian(xlim = c(-1, 20)) +
  # Apply our custom theme
  custom_theme +
  # Fine-tune some minimal theme elements
  theme(
    axis.text.y = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # Title
  ggtitle("Superorganismality in full phylogeny") +
  # Manual fill scale for ploidy
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid", "HD" = "haplodiploid"), 
    guide = "none"
  ) +
  # Optional alpha scale if you wish to vary transparency by clade
  scale_alpha_manual(values = rev(c(1, rep(0.4,3), 1, rep(0.4,5))), guide = FALSE)

# Display the final plot
p1

# ---- Save the plot as an SVG -----------------------------------------------
setwd(fig)  # Change to figure output directory
ggsave(
  filename = "eusocial_origins_clades_hapall_euskoos_fultree.svg",
  plot = p1,
  device = svglite,
  path = fig,
  width = 8,
  height = 6
)
#-------------------------------------------------------------------------------