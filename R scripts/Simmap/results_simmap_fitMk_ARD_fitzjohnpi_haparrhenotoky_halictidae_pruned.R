################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis on a pruned 
#       insect phylogeny (where Halictidae has been pruned to force only 2 
#       eusociality origins).
#
# NOTE 1: Haplodiploidy is used in the strict sense (HD.arrhenotoky).
# NOTE 2: Eusociality is used in the broad sense.
################################################################################

# ---- Load libraries ----------------------------------------------------------
library(phytools)    # For phylogenetic tree manipulation, simmap objects, etc.
library(ggplot2)     # For visualization
library(RColorBrewer)# Additional color palettes
library(viridis)     # Additional color palettes
library(patchwork)   # Combine multiple ggplot-based plots
library(gridExtra)   # Arrange multiple grid-based plots
library(ggridges)    # Ridgeline plots
library(svglite)     # Exporting plots as SVG
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Adjust these for your desktop or laptop environment. Comment/uncomment as needed.

# Directory for figure output
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") # ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures")   # laptop

# Directory for data files (phylogeny, phenotype, etc.)
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") # ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data")   # laptop

# Directory for Simmap results
res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") # ASU desktop
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results")   # laptop
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
setwd(res)
list.files()  # Inspect the files available in this directory

# Example lines (commented out) if you needed to remove a tip and re-save:
# simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_halictidae_pruned.RDS")
# simmap.trees = drop.tip(simmap.trees, "Paraponera_clavata")
# saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_halictidae_pruned_paraponera_removed.RDS")

# Load the updated Simmap object (with "Paraponera_clavata" removed)
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_halictidae_pruned_paraponera_removed.RDS")
#-------------------------------------------------------------------------------


# ----- Set working directory for data files -----------------------------------
setwd(data)
list.files()

# Load the pruned insect phylogeny (halictidae pruned to enforce certain eusociality assumptions)
insect.tree = read.tree("halictidae_pruned_insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset (already tailored for the pruned tree)
pheno.df = read.csv("halictidae_pruned_processed_pheno_df.csv")
row.names(pheno.df) = pheno.df$name
#-------------------------------------------------------------------------------


# ---- Identify clade members from different orders/families -------------------
# We'll use 'caper::clade.members' and 'findMRCA' to locate the tips in major clades.

halictids = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, family == "Halictidae" & Eusocial == "eusocial"))
  )),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

thrips = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, order == "Thysanoptera"))
  )),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

aphids = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, family == "Aphididae"))
  )),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

termites = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, order == "Blattodea"))
  )),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

ants = caper::clade.members(
  phangorn::Ancestors(
    insect.tree, 
    as.numeric(findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, family == "Formicidae"))
    )),
    type="all"
  )[2],
  insect.tree, tip.labels=TRUE, include.nodes=FALSE
)
ants = ants[ants != "Paraponera_clavata"]  # Remove a specific tip if desired

# Check ants in the phenotype df that are not labeled eusocial (sanity check)
subset(pheno.df, row.names(pheno.df) %in% ants & Eusocial != "eusocial")

# Platypodine ambrosia beetles
platypodine = c(
  "Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus","Cylindropalpus",
  "Dendroplatypus","Dinoplatypus","Doliopygus","Epiplatypus","Euplatypus","Megaplatypus",
  "Mesoplatypus","Myoplatypus","Neotrachyostus","Oxoplatypus","Pereioplatypus","Peroplatypus",
  "Platyphysus","Platypus","Teloplatypus","Trachyostus","Treptoplatypus","Triozastus"
)

platypodines = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, genus %in% platypodine))
  )),
  insect.tree, tip.labels=TRUE, include.nodes=TRUE
)$tips

vespids = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, family == "Vespidae"))
  )),
  insect.tree, tip.labels=TRUE, include.nodes=TRUE
)$tips

apids = caper::clade.members(
  as.numeric(findMRCA(
    insect.tree, 
    tips = row.names(subset(pheno.df, family == "Apidae"))
  )),
  insect.tree, tip.labels=TRUE, include.nodes=TRUE
)$tips

# Identify haplodiploid vs. diploid sets (HD.arrhenotoky)
HD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "HD"))
DD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "DD"))

# Collect all clades in a list, naming each element
clades = list(
  HD.clade, DD.clade, vespids, platypodines, ants, termites, 
  aphids, thrips, halictids, apids
)
names(clades) = c(
  "all haplodiploids","all diploids","vespids","platypodines","ants",
  "termites","aphids","thrips","halictids","apids"
)

# Create an empty dataframe to store summarized simmap counts
clades.df = data.frame(matrix(nrow=0, ncol=3))
colnames(clades.df) = c("eusocial","clade","ploidy")

# Loop over each clade, keep only the relevant tips from the Simmap object,
# then describe and extract the count data for each social state
for(c in 1:length(clades)){
  focal.trees = keep.tip(simmap.trees, clades[[c]])
  obj = describe.simmap(focal.trees)
  
  df = data.frame(obj$count)
  
  # Determine if the clade is haplodiploid (HD) or diploid (DD)
  ploidy = unique(pheno.df[clades[[c]], ]$HD.arrhenotoky)
  
  # Subset columns relevant to your states, e.g. "HD.solitary" and "HD.eusocial"
  df = df[c(paste(ploidy,"solitary",ploidy,"eusocial",sep="."))]
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Rename columns for consistency
  colnames(df) = c("eusocial","clade","ploidy")
  
  # Bind the row(s) to the master dataframe
  clades.df = rbind(clades.df, df)
}


# ---- Custom ggplot2 theme ----------------------------------------------------
custom_theme <- theme(
  strip.background = element_blank(),
  strip.placement  = "outside",
  strip.text       = element_text(family = "Arial", colour = "black", 
                                  face = "italic", size = 14),
  panel.background = element_rect(fill = NA),
  # panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position  = "bottom",
  legend.title     = element_text(family = "Arial", colour = "black", 
                                  size = 14, face = "italic"),
  legend.text      = element_text(family = "Arial", colour = "black", size = 16),
  legend.background= element_rect(fill = NA),
  legend.key       = element_rect(color = NA, fill = NA, linewidth = NA),
  axis.title       = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 22),
  axis.title.y     = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 22, vjust = 2),
  axis.text        = element_text(family = "Arial", colour = "black", size = 18),
  axis.ticks       = element_line(colour = "black", linewidth = 1),
  plot.background  = element_rect(fill = "white"),
  plot.title       = element_text(family = "Arial", colour = "black", 
                                  face = "bold", size = 24, hjust = 0.5,
                                  margin = margin(t=0.5, r=0.5, b=0.5, 
                                                  l=0.5, unit = "cm"))
)


# ---- Quick ridgeline demonstration ------------------------------------------
ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, 
                      draw_baseline = FALSE) +
  theme_ridges()


# ---- Adjust factor levels & finalize the ridgeline plot ----------------------
# Reorder the factor levels to ensure a logical vertical plot order
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c("all diploids","aphids","platypodines","termites",
                 "all haplodiploids","halictids","apids","vespids","ants","thrips"))
)
clades.df$ploidy[clades.df$clade == "all haplodiploids"] = "HD"

# Calculate the number of bins for binned ridgeline plot
bins = max(clades.df$eusocial) - min(clades.df$eusocial)

# Build the final plot
p1 = ggplot(clades.df, aes(x = eusocial, y = clade, height = stat(density),
                           fill = ploidy)) +
  geom_density_ridges(
    stat = "binline", 
    bins = bins, 
    scale = 0.7, 
    draw_baseline = TRUE
  ) +
  scale_x_continuous(
    expand = c(0,0),
    breaks = seq(0,20,by=2),
    name   = "origins of eusociality"
  ) +
  scale_y_discrete(
    expand = c(0,0),
    name   = "insect clade"
  ) +
  coord_cartesian(clip = "off") + 
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Eusociality in pruned phylogeny") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid", "HD" = "haplodiploid"),
    guide  = "none"
  )
# scale_alpha_manual could be used to adjust transparency if desired

# Display the final plot
p1

# Save the plot as an SVG
setwd(fig)
ggsave(
  filename = "eusocial_origins_clades_haparrh_eusall_halpruned.svg",
  plot     = p1,
  device   = svglite,
  path     = fig,
  width    = 8,
  height   = 6
)
#-------------------------------------------------------------------------------

