#---- Consolidated results from the Simmap--------------------------------------
# NOTE: This analysis is performed on the tree where halictidae is pruned to force only 2 eusociality origins
# NOTE: This is for the Chester's 2020 phylogeny
# NOTE: Haplodiploidy and Eusociality are used in broad sense

# ---- Load libraries ----------------------------------------------------------
library(phytools)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(gridExtra)
library(ggridges)
library(svglite)
#-------------------------------------------------------------------------------


# Set directory shortcuts-------------------------------------------------------

# Output directory for the figures
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures") #directory for laptop

# Data directory for the phenotypic dataset and tree files
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") #directory for laptop

# Data directory for simmap results
res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") #directory for ASU desktop
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/Simmap/updated final results") #directory for laptop
#-------------------------------------------------------------------------------


# ----- Set working directory --------------------------------------------------
setwd(data)
list.files()

# Load the insect phylogeny
insect.tree = read.tree("halictidae_pruned_insect_phylogeny_2020.tre")

# Load the phenotype dataset
pheno.df = read.csv("halictidae_pruned_processed_pheno_df_2020.csv") # This CSV file is a data-wrangled version of pheno_df_16102023.csv
row.names(pheno.df)=pheno.df$name
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
setwd(res)
list.files()

simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hapall_halictidae_pruned_2020.RDS")
#-------------------------------------------------------------------------------


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

# Haplodiploid/diploid sets (HD.all == "HD" or "DD")
HD.clade = row.names(subset(pheno.df, HD.all == "HD"))
DD.clade = row.names(subset(pheno.df, HD.all == "DD"))

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

write.csv(clades.df, file = "clades_df_simmap_hapall_eusall_halpruned_2020.csv", row.names = FALSE)
clades.df = read.csv("clades_df_simmap_hapall_eusall_halpruned_2020.csv")


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


# ----- PLOT 1 -----------------------------------------------------------------

# --- carry over factor order ---------------------------------------------------
clades.df.combined$clade <- factor(
  clades.df.combined$clade,
  levels = rev(c(
    "all haplodiploids","halictids","apids","vespids","ants","thrips",
    "all diploids","aphids","platypodines","termites"
  ))
)

# --- ensure "all haplodiploids" labeled as HD ---------------------------------
clades.df.combined$ploidy[clades.df.combined$clade == "all haplodiploids"] <- "HD"

# --- unit-width bins from 0..x_upper ------------------------------------------
x_max   <- max(clades.df.combined$count, na.rm = TRUE)
x_upper <- max(20, ceiling(x_max))                 # at least 0–20, extend if needed
bins    <- x_upper                                 # 1 bin per integer
brks    <- seq(0, x_upper, by = 2)

# --- base layer: fills only (no borders), alpha by transition ------------------
p2 <-
  ggplot(
    clades.df.combined,
    aes(
      x      = count,
      y      = clade,
      height = after_stat(density),
      fill   = ploidy,
      alpha  = transition
    )
  ) +
  geom_density_ridges(
    stat          = "binline",
    bins          = bins,
    from          = 0,
    to            = x_upper,
    scale         = 0.7,
    draw_baseline = TRUE,
    color         = NA                  # no border on this layer (affects losses)
  ) +
  # gains-only outline: thicker border, losses still have none
  geom_density_ridges(
    data = subset(clades.df.combined, transition == "eusocial"),
    aes(x = count, y = clade, height = after_stat(density)),
    stat          = "binline",
    bins          = bins,
    from          = 0,
    to            = x_upper,
    scale         = 0.7,
    draw_baseline = FALSE,
    fill          = NA,
    color         = "grey30",
    size          = 0.6
  ) +
  # axes, padding, labels -------------------------------------------------------
scale_x_continuous(
  breaks = brks,
  name   = "Number of transitions (gains & losses)",
  expand = expansion(mult = 0, add = c(0.8, 0))   # add left space so bars don't touch y-axis
) +
  scale_y_discrete(
    name   = "Insect clade",
    expand = expansion(mult = c(0.06, 0.10))        # add some headroom/top padding
  ) +
  coord_cartesian(xlim = c(0, x_upper)) +           # keep domain >= 0
  custom_theme +
  theme(
    axis.text.y       = element_text(vjust = 0),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(t = 14, r = 10, b = 10, l = 12)  # extra top margin
  ) +
  ggtitle("Hap_all Eus_all") +
  scale_fill_manual(
    values = c("DD" = "#4c96ac", "HD" = "#b05656"),
    labels = c("DD" = "diploid",  "HD" = "haplodiploid"),
    guide  = "none"
  ) +
  scale_alpha_manual(
    values = c(eusocial = 0.95, solitary = 0.5),
    labels = c(eusocial = "gains", solitary = "losses"),
    name   = "Transition"
  ) +
  guides(alpha = guide_legend(override.aes = list(color = c("grey30", NA), size = c(0.6, 0))))

p2


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

ggsave(filename = "eusocial_origins_clades_hapall_eusall_halpruned_2020.svg", plot = p1, device = svglite, path = fig, width = 8, height = 6)





































