################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis on a full 
#       insect phylogeny, focusing on strict haplodiploidy (HD.arrhenotoky)
#       and a broad definition of eusociality.
################################################################################

# ---- Load libraries ----------------------------------------------------------
library(phytools)    # For phylogenetic tree manipulation, simmap objects, etc.
library(ggplot2)     # For data visualization
library(RColorBrewer)# Color palettes
library(viridis)     # Additional color palettes
library(patchwork)   # Combining multiple ggplot2 plots
library(gridExtra)   # Arranging multiple grid-based plots
library(ggridges)    # Ridgeline plots
library(svglite)     # Saving plots in SVG format
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Modify according to your local environment (desktop vs. laptop).

# Directory for figures
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") # ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures")   # laptop

# Directory for data (phylogeny, phenotype data)
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") # ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data")   # laptop

# Directory for simmap results
res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") # ASU desktop
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results")   # laptop
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
setwd(res)  
list.files()  # Check available files in the directory

# Optionally load/edit the original trees, then save a version with certain tips dropped:
# simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree.RDS")
# simmap.trees=drop.tip(simmap.trees,"Paraponera_clavata")
# saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree_paraponera_removed.RDS")

# Load the updated simmap object (with "Paraponera_clavata" removed)
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree_paraponera_removed.RDS")
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
setwd(data)
list.files()

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset (CSV with pre-wrangled data)
pheno.df = read.csv("pheno_processed_df_20231023.csv")
row.names(pheno.df) = pheno.df$name
#-------------------------------------------------------------------------------


# ------ Identify clade members from different orders/families -----------------
# Using 'caper::clade.members' along with findMRCA to locate tips in major clades.

halictids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Halictidae" & Eusocial == "eusocial")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

thrips = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, order == "Thysanoptera")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

aphids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Aphididae")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

termites = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, order == "Blattodea")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

ants = caper::clade.members(
  phangorn::Ancestors(
    insect.tree, 
    as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Formicidae")))), 
    type = "all"
  )[2],
  insect.tree, tip.labels = TRUE, include.nodes = FALSE
)
# Remove a specific tip if needed
ants = ants[ants != "Paraponera_clavata"]

# Check for ants in the phenotype dataframe that are not labeled as eusocial
subset(pheno.df, row.names(pheno.df) %in% ants & Eusocial != "eusocial")

# Platypodine ambrosia beetles
platypodine = c("Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus",
                "Cylindropalpus","Dendroplatypus","Dinoplatypus","Doliopygus","Epiplatypus",
                "Euplatypus","Megaplatypus","Mesoplatypus","Myoplatypus","Neotrachyostus",
                "Oxoplatypus","Pereioplatypus","Peroplatypus","Platyphysus","Platypus",
                "Teloplatypus","Trachyostus","Treptoplatypus","Triozastus")

platypodines = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, genus %in% platypodine)))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Vespids (paper wasps, etc.)
vespids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Vespidae")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Apids (honey bees, stingless bees, bumble bees, etc.)
apids = caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Apidae")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Haplodiploid/diploid sets (HD.arrhenotoky == "HD" or "DD")
HD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "HD"))
DD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "DD"))

# Create a list of clades
clades = list(
  HD.clade, DD.clade, vespids, platypodines, ants, termites, 
  aphids, thrips, halictids, apids
)
names(clades) = c(
  "all haplodiploids","all diploids","vespids","platypodines","ants",
  "termites","aphids","thrips","halictids","apids"
)

# Initialize a dataframe to store counts from simmap
clades.df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(clades.df) = c("eusocial","clade","ploidy")

# Loop over each clade, describe the simmap, and extract the count data
for(c in 1:length(clades)){
  focal.trees = keep.tip(simmap.trees, clades[[c]])
  obj         = describe.simmap(focal.trees)
  
  # Convert the resulting count matrix to a data frame
  df = data.frame(obj$count)
  
  # Identify ploidy (HD or DD) by looking at the phenotype info for this clade
  ploidy = unique(pheno.df[clades[[c]], ]$HD.arrhenotoky)
  
  # Subset to the relevant columns: e.g. "HD.solitary" and "HD.eusocial"
  df = df[c(paste(ploidy, "solitary", ploidy, "eusocial", sep = "."))]
  
  # Store extra info on clade and ploidy
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Rename columns in a consistent manner
  colnames(df) = c("eusocial","clade","ploidy")
  
  # Append to clades.df
  clades.df = rbind(clades.df, df)
}


# ---- Custom ggplot2 theme ----------------------------------------------------
# This theme is applied to the ridgeline plots below.
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
  axis.text        = element_text(family = "Arial", colour = "black", size = 20),
  axis.ticks       = element_line(colour = "black", linewidth = 1),
  plot.background  = element_rect(fill = "white"),
  plot.title       = element_text(
    family = "Arial", colour = "black", face = "bold", size = 24,
    hjust = 0.5, margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit = "cm")
  )
)


# ---- Quick demonstration of ridgeline plot -----------------------------------
ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, 
                      draw_baseline = FALSE) +
  theme_ridges()


# ---- Adjust factor levels for consistent ordering ----------------------------
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c("all diploids","aphids","platypodines","termites",
                 "all haplodiploids","halictids","apids","vespids","ants","thrips"))
)
clades.df$ploidy[clades.df$clade == "all haplodiploids"] = "HD"


# Determine the number of bins for ridgeline binning
bins = max(clades.df$eusocial) - min(clades.df$eusocial)

# ---- Final ridgeline plot showing eusocial origins across clades ------------
p1 = ggplot(clades.df, aes(x = eusocial, y = clade, height = stat(density),
                           fill = ploidy)) +
  geom_density_ridges(
    stat          = "binline",
    bins          = bins,
    scale         = 0.7,
    draw_baseline = TRUE,
    color         = "grey70"
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 2),
    name   = "Origins of eusociality",
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0), name = "Insect clade") +
  coord_cartesian(xlim = c(-0.5, 20)) +
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Eusociality in full phylogeny") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid", "HD" = "haplodiploid"),
    guide  = "none"
  ) +
  scale_alpha_manual(values = rev(c(1, rep(0.4,3), 1, rep(0.4,5))), guide = FALSE)

# Display the final plot
p1

# ---- Save the plot -----------------------------------------------------------
setwd(fig)
ggsave(
  filename = "eusocial_origins_clades_haparrh_eusall_fulltree.svg",
  plot     = p1,
  device   = svglite,
  path     = fig,
  width    = 8,
  height   = 6
)
#-------------------------------------------------------------------------------

