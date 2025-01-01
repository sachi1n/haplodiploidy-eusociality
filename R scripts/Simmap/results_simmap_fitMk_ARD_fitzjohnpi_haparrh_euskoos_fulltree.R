################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis
#       on a full insect phylogeny, focusing on haplodiploidy and eusociality
#       (both used in a strict sense).
################################################################################

# ---- Load libraries ----------------------------------------------------------
library(phytools)    # For phylogenetic tree manipulation, Simmap objects, etc.
library(ggplot2)     # For data visualization
library(RColorBrewer)# For color palettes
library(viridis)     # Additional color scales
library(patchwork)   # For combining multiple ggplot plots
library(gridExtra)   # For arranging multiple grid-based plots
library(ggridges)    # For ridgeline plots
library(svglite)     # For saving plots in SVG format
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Modify these paths for your local environment. You can comment/uncomment
# depending on whether you're on a desktop or laptop.

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
list.files()  # Inspect the files available in the 'res' directory

# Optional loading/editing steps:
# simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_fulltree.RDS")
# simmap.trees = drop.tip(simmap.trees, "Paraponera_clavata")
# saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_fulltree_paraponera_removed.RDS")

# Load the updated Simmap object with "Paraponera_clavata" removed
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_euskoos_fulltree_paraponera_removed.RDS")
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
setwd(data)
list.files()

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset (already wrangled)
pheno.df = read.csv("pheno_processed_df_20231023.csv")
row.names(pheno.df) = pheno.df$name
#-------------------------------------------------------------------------------


# ------ Find clade members from different orders/families ---------------------
# We’ll identify major groups (Halictidae, Thysanoptera, Aphididae, etc.)
# using caper::clade.members and findMRCA.

halictids = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, family == "Halictidae" & Eusocial == "eusocial"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

thrips = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, order == "Thysanoptera"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

aphids = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, family == "Aphididae"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

termites = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, order == "Blattodea"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

ants = caper::clade.members(
  phangorn::Ancestors(
    insect.tree, 
    as.numeric(
      findMRCA(
        insect.tree, 
        tips = row.names(subset(pheno.df, family == "Formicidae"))
      )
    ), 
    type = "all"
  )[2],
  insect.tree, tip.labels = TRUE, include.nodes = FALSE
)
ants = ants[ants != "Paraponera_clavata"]  # remove a specific tip if needed

# Check for ants that are not labeled as eusocial in the dataset
subset(pheno.df, row.names(pheno.df) %in% ants & Eusocial != "eusocial")

# Platypodine ambrosia beetles
platypodine = c(
  "Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus",
  "Cylindropalpus","Dendroplatypus","Dinoplatypus","Doliopygus","Epiplatypus",
  "Euplatypus","Megaplatypus","Mesoplatypus","Myoplatypus","Neotrachyostus",
  "Oxoplatypus","Pereioplatypus","Peroplatypus","Platyphysus","Platypus",
  "Teloplatypus","Trachyostus","Treptoplatypus","Triozastus"
)
platypodines = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, genus %in% platypodine))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Vespids (e.g., wasps)
vespids = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, family == "Vespidae"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Apids (e.g., honey bees, bumble bees)
apids = caper::clade.members(
  as.numeric(
    findMRCA(
      insect.tree, 
      tips = row.names(subset(pheno.df, family == "Apidae"))
    )
  ),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

# Identify the haplodiploid/diploid sets (HD, DD), here specifically for arrhenotoky
HD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "HD"))
DD.clade = row.names(subset(pheno.df, HD.arrhenotoky == "DD"))

# Create a list of the clades to explore
clades = list(
  HD.clade, DD.clade, vespids, platypodines, ants, termites, 
  aphids, thrips, halictids, apids
)
names(clades) = c(
  "all haplodiploids","all diploids","vespids","platypodines","ants",
  "termites","aphids","thrips","halictids","apids"
)

# Subset that list to the clades of interest
clades = clades[c(1,2,3,5,6,10)]

# Initialize a data frame to store the counts of eusocial vs. solitary across states
clades.df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(clades.df) = c("eusocial","clade","ploidy")

# For each clade, subset the simmap trees, describe them, and extract counts
for(c in 1:length(clades)){
  focal.trees = keep.tip(simmap.trees, clades[[c]])
  
  # Summarize the simmap for these focal trees
  obj = describe.simmap(focal.trees)
  
  # Convert the 'count' matrix to a data frame
  df = data.frame(obj$count)
  
  # Identify whether the clade is haplodiploid (HD) or diploid (DD)
  ploidy = unique(pheno.df[clades[[c]], ]$HD.arrhenotoky)
  
  # Subset relevant columns: e.g. "HD.solitary", "HD.eusocial"
  df = df[c(paste(ploidy,"solitary", ploidy,"eusocial", sep = "."))]
  
  # Add clade and ploidy as columns
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Rename columns to ensure consistency
  colnames(df) = c("eusocial","clade","ploidy")
  
  # Append the newly created data frame to clades.df
  clades.df = rbind(clades.df, df)
}


# ---- Custom theme for ggplot2 ------------------------------------------------
custom_theme <- theme(
  strip.background = element_blank(),
  strip.placement  = "outside",
  strip.text = element_text(family = "Arial", colour = "black", 
                            face = "italic", size = 14),
  panel.background  = element_rect(fill = NA),
  # panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  legend.position   = "bottom",
  legend.title      = element_text(family = "Arial", colour = "black", 
                                   size = 14, face = "italic"),
  legend.text       = element_text(family = "Arial", colour = "black", size = 16),
  legend.background = element_rect(fill = NA),
  legend.key        = element_rect(color = NA, fill = NA, linewidth = NA),
  axis.title        = element_text(family = "Arial", colour = "black", 
                                   face = "bold", size = 22),
  axis.title.y      = element_text(family = "Arial", colour = "black", 
                                   face = "bold", size = 22, vjust = 2),
  axis.text         = element_text(family = "Arial", colour = "black", size = 20),
  axis.ticks        = element_line(colour = "black", linewidth = 1),
  plot.background   = element_rect(fill = "white"),
  plot.title        = element_text(
    family = "Arial", colour = "black", face = "bold", size = 24,
    hjust = 0.5, 
    margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit = "cm")
  )
)


# ---- Quick ridgeline demonstration ------------------------------------------
ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, 
                      draw_baseline = FALSE) +
  theme_ridges()


# ---- Adjust factor levels for consistent ordering ----------------------------
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c(
    "all diploids","aphids","platypodines","termites",
    "all haplodiploids","halictids","apids","vespids","ants","thrips"
  ))
)
clades.df$ploidy[clades.df$clade == "all haplodiploids"] = "HD"


# ---- Add more clades to ensure consistency across definitions (dummy data) ---
new_clades <- data.frame(
  clade    = rep(c("aphids", "platypodines", "halictids", "thrips"), each = 100),
  eusocial = rep(0, 4 * 100),  # e.g. zero or dummy data
  ploidy   = rep(c("DD", "DD", "HD", "HD"), each = 100)
)
clades.df <- rbind(clades.df, new_clades)

# Re-factor to ensure the newly added clades fit into the desired order
clades.df$clade = factor(
  clades.df$clade,
  levels = rev(c(
    "all diploids","aphids","platypodines","termites","all haplodiploids",
    "halictids","apids","vespids","ants","thrips"
  ))
)

# Determine bin range for the ridgeline binning
bins = max(clades.df$eusocial) - min(clades.df$eusocial)

# ---- Final ridgeline plot ----------------------------------------------------
p1 = ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density),
                           fill = ploidy)) +
  geom_density_ridges(
    stat         = "binline",
    binwidth     = 1,
    bins         = bins,
    scale        = 0.7,
    binwidth     = 0.5, 
    draw_baseline= TRUE,
    color        = "grey70"
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 2),
    name   = "Origins of eusociality",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    expand = c(0,0),
    name   = "Insect clade"
  ) +
  coord_cartesian(xlim = c(-1, 20)) +
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Superorganismality in full phylogeny") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid", "HD" = "haplodiploid"),
    guide  = "none"
  ) +
  scale_alpha_manual(
    values = rev(c(1, rep(0.4,3), 1, rep(0.4,5))),
    guide  = FALSE
  )

# Display the final plot
p1

# ---- Save the plot -----------------------------------------------------------
setwd(fig)  # Switch to the directory for figure outputs
ggsave(
  filename = "eusocial_origins_clades_haparrh_euskoos_fulltree.svg",
  plot     = p1,
  device   = svglite,
  path     = fig,
  width    = 8,
  height   = 6
)
#-------------------------------------------------------------------------------
