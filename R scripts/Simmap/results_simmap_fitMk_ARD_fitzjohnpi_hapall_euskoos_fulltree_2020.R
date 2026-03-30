################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis
#       using a full insect phylogeny, focusing on transitions to eusociality.
#
# NOTE 1: "Haplodiploidy" here is broadly defined; "Eusociality" is defined strictly.
# NOTE 2: Adjust any directory paths and file names to match your local setup.
################################################################################


# ---- Load libraries ----------------------------------------------------------
library(phytools)     # For phylogenetic tree manipulation and simmap objects
library(ggplot2)      # For plotting
library(RColorBrewer) # For color palettes
library(viridis)      # For color palettes
library(patchwork)    # For combining multiple ggplot2 plots
library(gridExtra)    # For arranging multiple grid-based plots
library(ggridges)     # For ridgeline plots
library(svglite)      # For saving plots in SVG format
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


# ----- Set working directory for data files -----------------------------------
setwd(data)
list.files()  # Inspect the files in the data directory

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2020_rooted_calibrated.tre")

# Load the phenotype dataset (already wrangled into a CSV)
pheno.df = read.csv("pheno_processed_df_2020_20250626.csv") 
row.names(pheno.df) = pheno.df$name  # Set 'name' as the row names
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
# Set working directory to where Simmap results are stored
setwd(res)
list.files()  # Check what files are available in this directory

# Load or modify the simmap object 
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hdall_euskoos_fulltree_2020.RDS")
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

ants=caper::clade.members(
  as.numeric(findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Formicidae")))),
  insect.tree, tip.labels = TRUE, include.nodes = TRUE
)$tips

#add a solitary outgroup to ant clade
ant.out=subset(pheno.df,family=="Ampulicidae")$name[1]
#ant.out=caper::clade.members(getParent(insect.tree,findMRCA(insect.tree, tips = row.names(subset(pheno.df, family == "Formicidae")))),insect.tree,tip.labels=TRUE)[1]
ants=c(ants,ant.out)

# Platypodine ambrosia beetles
platypodine = c("Austroplatypus","Doliopygus","Mesoplatypus","Platypus")

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
HD.clade = row.names(subset(pheno.df, HD.all == "HD"))
DD.clade = row.names(subset(pheno.df, HD.all == "DD"))

# Create a list of clades
clades = list(
  HD.clade, DD.clade, vespids, ants, termites,apids)

names(clades) = c(
  "all haplodiploids","all diploids","vespids","ants",
  "termites","apids"
)

#for trial plotting, exclude clades with few spp
#clades = list(
#  HD.clade, DD.clade, vespids, ants, termites, 
#  aphids, halictids, apids
#)
#names(clades) = c(
#  "all haplodiploids","all diploids","vespids","ants",
#  "termites","aphids","halictids","apids"
#)

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
  ploidy = unique(pheno.df[clades[[c]], ]$HD.all)
  
  #subset solitary->eusocial and eusocial->solitary
  df = df[c(paste(ploidy, "solitary", ploidy, "eusocial", sep = "."),
            paste(ploidy, "eusocial", ploidy, "solitary", sep = "."))]
  
  # Store extra info on clade and ploidy
  df$clade  = names(clades)[c]
  df$ploidy = ploidy
  
  # Rename columns in a consistent manner
  colnames(df) = c("eusocial","solitary","clade","ploidy")
  
  # Append to clades.df
  clades.df = rbind(clades.df, df)
}

write.csv(clades.df, file = "clades_df_simmap_hapall_euskoos_fulltree_2020.csv", row.names = FALSE)
clades.df = read.csv("clades_df_simmap_hapall_euskoos_fulltree_2020.csv")

# ---- Custom theme for ggplot2 ------------------------------------------------
# This theme is applied to the ridgeline plots below.
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


# ----- PLOT 2 -----------------------------------------------------------------

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
    "all haplodiploids","halictids","apids","vespids","ants","thrips",
    "all diploids","aphids","platypodines","termites"
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
    color        = "grey"
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 1),
    name   = "Origins of eusociality",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.05, 0.1)),
    name   = "Insect clade"
  ) +
  coord_cartesian(xlim = c(-1.5, 20)) +
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin       = margin(t = 14, r = 10, b = 10, l = 12) 
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

setwd(fig)

ggsave(filename = "eusocial_origins_clades_hapall_euskoos_2020.svg", plot = p1, device = svglite, path = fig, width = 8, height = 6)
