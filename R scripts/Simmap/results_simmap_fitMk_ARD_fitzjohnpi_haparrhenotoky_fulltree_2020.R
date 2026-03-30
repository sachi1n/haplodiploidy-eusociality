################################################################################
# GOAL: Consolidate and visualize results from a Simmap analysis on a full 
#       insect phylogeny, focusing on strict haplodiploidy (HD.arrhenotoky)
#       and a broad definition of eusociality.
################################################################################

# ---- Load libraries ----------------------------------------------------------
library(phytools)     # For phylogenetic tree manipulation, simmap objects, etc.
library(ggplot2)      # For data visualization
library(RColorBrewer) # Color palettes
library(viridis)      # Additional color palettes
library(patchwork)    # Combining multiple ggplot2 plots
library(gridExtra)    # Arranging multiple grid-based plots
library(ggridges)     # Ridgeline plots
library(svglite)      # Saving plots in SVG format
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Modify according to your local environment

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


# ----- Set working directory --------------------------------------------------
setwd(data)
list.files()

# Load the insect phylogeny
insect.tree = read.tree("insect_phylogeny_2020_rooted_calibrated.tre")

# Load the phenotype dataset (already wrangled into a CSV)
pheno.df = read.csv("pheno_processed_df_2020_20250626.csv") 
row.names(pheno.df) = pheno.df$name  # Set 'name' as the row names
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
setwd(res)  
list.files()  # Check available files in the directory

# Load/edit the simmap:
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_haparrhenotoky_fulltree_2020.RDS")
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
  ploidy = unique(pheno.df[clades[[c]], ]$HD.arrhenotoky)
  
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

write.csv(clades.df, file = "clades_df_simmap_haparrhenotoky_fulltree_2020.csv", row.names = FALSE)
clades.df = read.csv("clades_df_simmap_haparrhenotoky_fulltree_2020.csv")

# ---- Custom ggplot2 theme ----------------------------------------------------
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


#show both gains and losses of eusociality
clades.df.eusocial=clades.df
clades.df.eusocial$transition="eusocial"
clades.df.eusocial=clades.df.eusocial[c("eusocial","clade","ploidy","transition")]
colnames(clades.df.eusocial)=c("count","clade","ploidy","transition")

clades.df.solitary=clades.df
clades.df.solitary$transition="solitary"
clades.df.solitary=clades.df.solitary[c("solitary","clade","ploidy","transition")]
colnames(clades.df.solitary)=c("count","clade","ploidy","transition")

clades.df.combined=rbind(clades.df.eusocial,clades.df.solitary)

# ---- Factor order (same as before) -------------------------------------------
clades.df.combined$clade <- factor(
  clades.df.combined$clade,
  levels = rev(c(
    "all haplodiploids","halictids","apids","vespids","ants","thrips",
    "all diploids","aphids","platypodines","termites"
  ))
)


# ---- Ensure "all haplodiploids" labeled as HD --------------------------------
clades.df.combined$ploidy[clades.df.combined$clade == "all haplodiploids"] <- "HD"

# ---- Axis & bin settings (as before, data-aware) ------------------------------
x_max   <- max(clades.df.combined$count, na.rm = TRUE)
x_upper <- max(20, ceiling(x_max / 2) * 2)      # at least 0–20, extend if needed
brks    <- seq(0, x_upper, by = 2)
bins    <- max(1, diff(range(clades.df.combined$count, na.rm = TRUE)))

# ---- One-panel overlay of gains & losses -------------------------------------
p1 <-
  ggplot(
    clades.df.combined,
    aes(
      x = count,
      y = clade,
      height = after_stat(density),
      fill = ploidy,
      linetype = transition,
      color = transition,
      alpha = transition,
      group = interaction(clade, transition)  # ensure 2 ridges per clade
    )
  ) +
  geom_density_ridges(
    stat          = "binline",
    bins          = bins,
    from          = 0,
    to            = x_upper,
    scale         = 0.7,
    draw_baseline = TRUE,
    color         = NA # we'll re-add outlines in a second layer
  ) +
  # Thin outline pass so the two transitions are visually separable
  geom_density_ridges(
    stat          = "binline",
    bins          = bins,
    from          = 0,
    to            = x_upper,
    scale         = 0.7,
    draw_baseline = FALSE,
    fill          = NA
  ) +
  scale_x_continuous(
    breaks = brks,
    name   = "Number of transitions (gains & losses)",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0),
    name   = "Insect clade"
  ) +
  coord_cartesian(xlim = c(-0.6, x_upper)) +
  custom_theme +
  theme(
    axis.text.y      = element_text(vjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Hap_arrh Eus_all") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid",  "HD" = "haplodiploid"),
    guide  = "none"   # keep ploidy legend off, as in your original
  ) +
  scale_linetype_manual(
    values = c(eusocial = "solid", solitary = "22"),
    labels = c(eusocial = "gains", solitary = "losses"),
    name   = "Transition"
  ) +
  scale_color_manual(
    values = c(eusocial = "black", solitary = "grey35"),
    labels = c(eusocial = "gains", solitary = "losses"),
    name   = "Transition"
  ) +
  scale_alpha_manual(
    values = c(eusocial = 0.9, solitary = 0.5),
    labels = c(eusocial = "gains", solitary = "losses"),
    name   = "Transition"
  )

p1






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
  filename = "eusocial_origins_clades_haparrh_eusall_fulltree_no_outgroup.svg",
  plot     = p1,
  device   = svglite,
  path     = fig,
  width    = 8,
  height   = 6
)
#-------------------------------------------------------------------------------

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
  coord_cartesian(xlim = c(-0.6, 20)) +
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

ggsave(filename = "eusocial_origins_clades_haparrh_eusall_2020.svg", plot = p1, device = svglite, path = fig, width = 8, height = 6)
