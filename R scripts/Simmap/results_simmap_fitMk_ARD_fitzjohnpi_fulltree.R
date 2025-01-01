################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis using an 
#       insect phylogeny, where both haplodiploidy and eusociality are defined 
#       in a broad sense.
################################################################################

# ---- Load libraries ----------------------------------------------------------
library(phytools)     # For phylogenetic tree manipulation, simmap objects, etc.
library(ggplot2)      # For visualization
library(RColorBrewer) # Additional color palettes
library(viridis)      # Additional color palettes
library(patchwork)    # Combine multiple ggplot-based plots
library(gridExtra)    # Arrange multiple grid-based plots
library(ggridges)     # Ridgeline plots
library(svglite)      # For saving plots as SVG
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Modify these for your environment (e.g., desktop vs. laptop).

fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") 
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures")  

data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") 
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data")  

res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") 
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results")
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
# Switch to the data directory to load the phylogeny and phenotype file
setwd(data)
list.files()

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset (already wrangled into a CSV)
pheno.df = read.csv("pheno_processed_df_20231023.csv") 
row.names(pheno.df) = pheno.df$name  # Set 'name' as the row names
#-------------------------------------------------------------------------------


# ------ Identify clade members from different orders/families -----------------
# We'll use caper::clade.members and findMRCA to locate tips within major clades.

halictids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, family == "Halictidae" & Eusocial == "eusocial")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

thrips = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, order == "Thysanoptera")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

aphids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, family == "Aphididae")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

termites = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, order == "Blattodea")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

ants = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, family == "Formicidae")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Optionally remove a specific tip (Paraponera_clavata) if it's known to cause issues
ants = ants[ants != "Paraponera_clavata"]

# Platypodine ambrosia beetles
platypodine = c(
  "Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus",
  "Cylindropalpus","Dendroplatypus","Dinoplatypus","Doliopygus","Epiplatypus",
  "Euplatypus","Megaplatypus","Mesoplatypus","Myoplatypus","Neotrachyostus",
  "Oxoplatypus","Pereioplatypus","Peroplatypus","Platyphysus","Platypus",
  "Teloplatypus","Trachyostus","Treptoplatypus","Triozastus"
)

platypodines = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, genus %in% platypodine)))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

vespids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, family == "Vespidae")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

apids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, 
                      tips = row.names(subset(pheno.df, family == "Apidae")))
  ), 
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Collect these major clades in a list if needed (for reference)
clades = list(
  vespids, platypodines, ants, termites, aphids, 
  thrips, halictids, apids
)
names(clades) = c(
  "vespids","platypodines","ants","termites",
  "aphids","thrips","halictids","apids"
)
#-------------------------------------------------------------------------------


# ----- Load the simmap results ------------------------------------------------
setwd(res)
list.files()  # Inspect available files

simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_fulltree.RDS")
# Example of removing a tip if necessary
simmap.trees = drop.tip(simmap.trees, "Paraponera_clavata")
saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_fulltree_paraponera_removed.RDS")

# Re-load the updated simmap object
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_fulltree_paraponera_removed.RDS")

# Identify haplodiploid/diploid sets in the phenotype data
HD.clade = row.names(subset(pheno.df, HD.all == "HD"))
DD.clade = row.names(subset(pheno.df, HD.all == "DD"))

# Combine all clade vectors (HD, DD, plus the insect families/orders identified above)
clades = list(
  HD.clade, DD.clade, vespids, platypodines, ants, 
  termites, aphids, thrips, halictids, apids
)
names(clades) = c(
  "all haplodiploids","all diploids","vespids","platypodines","ants",
  "termites","aphids","thrips","halictids","apids"
)
#-------------------------------------------------------------------------------


# ----- Summarize simmap results for each clade --------------------------------
# We'll create a data frame (clades.df) to hold the counts of lineage states (e.g., solitary/eusocial).

clades.df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(clades.df) = c("eusocial","clade","ploidy")

# Loop over each clade
for (c in 1:length(clades)) {
  # Subset the simmap trees to only the tips in this specific clade
  focal.trees = keep.tip(simmap.trees, clades[[c]])
  
  # Describe the simmap object to get counts of mapped states
  obj = describe.simmap(focal.trees)
  df  = data.frame(obj$count)
  
  # Determine whether the clade is haplodiploid or diploid based on phenotype data
  ploidy = unique(pheno.df[clades[[c]], ]$HD.all)
  
  # Subset columns relevant to your states of interest (here: "solitary" and "eusocial")
  df = df[c(paste(ploidy, "solitary", ploidy, "eusocial", sep = "."))]
  
  # Add columns for the clade name and ploidy
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Rename columns consistently
  colnames(df) = c("eusocial","clade","ploidy")
  
  # Append to the master data frame
  clades.df = rbind(clades.df, df)
}
#-------------------------------------------------------------------------------


# ---- Define a custom theme for ggplot2 ---------------------------------------
custom_theme <- theme(
  strip.background = element_blank(),
  strip.placement  = "outside",
  strip.text       = element_text(family = "Arial", colour = "black", 
                                  face = "italic", size = 14),
  panel.background = element_rect(fill = NA),
  panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position  = "bottom",
  legend.title     = element_text(family = "Arial", colour = "black", 
                                  size = 14, face = "italic"),
  legend.text      = element_text(family = "Arial", colour = "black", size = 16),
  legend.background= element_rect(fill = NA),
  legend.key       = element_rect(color = NA, fill = NA, linewidth = NA),
  axis.title       = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 18),
  axis.title.y     = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 18, vjust = 2),
  axis.text        = element_text(family = "Arial", colour = "black", size = 14),
  axis.ticks       = element_line(colour = "black", linewidth = 1),
  plot.background  = element_rect(fill = "white"),
  plot.title       = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 20, hjust = 0.5,
                                  margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, 
                                                  unit = "cm"))
)


# ---- Quick check with a ridgeline plot ---------------------------------------
# This initial plot uses stat = "binline" to bin the eusocial values for each clade.
ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, 
                      draw_baseline = FALSE) +
  theme_ridges()


# ---- Prepare factor levels & final ridgeline plot ----------------------------
# Reorder the factor levels to control the vertical ordering of clades
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c(
    "all diploids","aphids","platypodines","termites","all haplodiploids",
    "halictids","apids","vespids","ants","thrips"
  ))
)

# Ensure the "all haplodiploids" row is labeled as HD
clades.df$ploidy[clades.df$clade == "all haplodiploids"] = "HD"

# Calculate the bin range for the ridgeline plot
bins = max(clades.df$eusocial) - min(clades.df$eusocial)

# Build the final ridgeline plot
p1 = ggplot(clades.df, aes(x = eusocial, y = clade, height = stat(density),
                           fill = ploidy)) +
  geom_density_ridges(
    stat = "binline", 
    bins = bins,
    from = 0, 
    to   = 20,   # set range if needed
    scale        = 0.7, 
    draw_baseline= TRUE,
    color        = "grey70"
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 2),
    name   = "Origins of eusociality",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0),
    name   = "Insect clade"
  ) +
  coord_cartesian(xlim = c(-0.6, 20)) +
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Eusociality in pruned phylogeny") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid",  "HD" = "haplodiploid"),
    guide  = "none"
  ) +
  scale_alpha_manual(values = rev(c(1, rep(0.4,3), 1, rep(0.4,5))), guide = FALSE)

# Display the final plot
p1

# ---- Save the plot -----------------------------------------------------------
setwd(fig)
ggsave(
  filename = "eusocial_origins_clades_hapall_eusall.svg",
  plot     = p1,
  device   = svglite,
  path     = fig,
  width    = 8,
  height   = 6
)
#-------------------------------------------------------------------------------