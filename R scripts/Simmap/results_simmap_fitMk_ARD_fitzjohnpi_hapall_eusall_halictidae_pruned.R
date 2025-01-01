#---- Consolidated results from the Simmap--------------------------------------
# NOTE: This analysis is performed on the tree where halictidae is pruned to force only 2 eusociality origins
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
insect.tree = read.tree("halictidae_pruned_insect_phylogeny_2017_collapsed_calibrated.tre")

# Load the phenotype dataset
pheno.df = read.csv("halictidae_pruned_processed_pheno_df.csv") # This CSV file is a data-wrangled version of pheno_df_16102023.csv
row.names(pheno.df)=pheno.df$name
#-------------------------------------------------------------------------------


# ----- Read the .rds dataset and plot density ---------------------------------
setwd(res)
list.files()

#simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hapall_halictidae_pruned.RDS")

#simmap.trees=drop.tip(simmap.trees,"Paraponera_clavata")
#saveRDS(simmap.trees, file = "hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hapall_halictidae_pruned_paraponera_removed.RDS")
simmap.trees = readRDS("hap_and_eus_simmap_fitMk_ARD_fitzjohnpi_hapall_halictidae_pruned_paraponera_removed.RDS")
#-------------------------------------------------------------------------------



halictids=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,family=="Halictidae" & Eusocial=="eusocial")))),
                               insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

thrips=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,order=="Thysanoptera")))),
                            insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

aphids=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,family=="Aphididae")))),
                            insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

termites=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,order=="Blattodea")))),
                              insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips

ants=caper::clade.members(phangorn::Ancestors(insect.tree,as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,family=="Formicidae")))),type="all")[2]
                          ,insect.tree,tip.labels=TRUE,include.nodes=FALSE)
ants=ants[ants!="Paraponera_clavata"]

subset(pheno.df,row.names(pheno.df) %in% ants & Eusocial!="eusocial")

platypodine=c("Austroplatypus","Baiocis","Carchesiopygus","Costaroplatus","Crossotarsus","Cylindropalpus","Dendroplatypus","Dinoplatypus","Doliopygus",
              "Epiplatypus","Euplatypus","Megaplatypus","Mesoplatypus","Myoplatypus","Neotrachyostus","Oxoplatypus","Pereioplatypus","Peroplatypus","Platyphysus","Platypus","Teloplatypus","Trachyostus","Treptoplatypus","Triozastus")

platypodines=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,genus %in% platypodine)))),
                                  insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

vespids=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,family=="Vespidae")))),
                             insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

apids=caper::clade.members(as.numeric(findMRCA(insect.tree,tips=row.names(subset(pheno.df,family=="Apidae")))),
                           insect.tree,tip.labels=TRUE,include.nodes=TRUE)$tips # 238

HD.clade=row.names(subset(pheno.df,HD.all=="HD"))
DD.clade=row.names(subset(pheno.df,HD.all=="DD"))

clades=list(HD.clade,DD.clade,vespids,platypodines,ants,termites,aphids,thrips,halictids,apids)
names(clades)=c("all haplodiploids","all diploids","vespids","platypodines","ants","termites","aphids","thrips","halictids","apids")

clades.df=data.frame(matrix(nrow=0,ncol=3))
colnames(clades.df)=c("eusocial","clade","ploidy")

for(c in 1:length(clades)){
  focal.trees=keep.tip(simmap.trees,clades[[c]])
  obj=describe.simmap(focal.trees)
  df=data.frame(obj$count)
  ploidy=unique(pheno.df[clades[[c]],]$HD.all)
  df=df[c(paste(ploidy,"solitary",ploidy,"eusocial",sep="."))]
  df$clade=names(clades)[c]
  df$ploidy=ploidy
  colnames(df)=c("eusocial","clade","ploidy")
  clades.df=rbind(clades.df,df)
}


# Use a custom theme for ggplot2
custom_theme <- theme(
  strip.background = element_blank(),
  strip.placement = "outside",
  strip.text = element_text(family = "Arial", colour = "black", face = "italic", size = 14),
  panel.background = element_rect(fill = NA),
  #panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.), # Updated from size to linewidth
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "bottom",
  legend.title = element_text(family = "Arial", colour = "black", size = 14, face = "italic"),
  legend.text = element_text(family = "Arial", colour = "black", size = 16),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(color = NA, fill = NA, linewidth = NA), # Updated for consistency, though likely unnecessary for legend.key
  axis.title = element_text(family = "Arial", colour = "black", face = "bold", size = 22),
  axis.title.y = element_text(family = "Arial", colour = "black", face = "bold", size = 22, vjust = 2),
  axis.text = element_text(family = "Arial", colour = "black", size = 20),
  axis.ticks = element_line(colour = "black", linewidth = 1), # Updated from size to linewidth
  plot.background = element_rect(fill = "white"),
  plot.title = element_text(family = "Arial", colour = "black", face = "bold", size = 24, hjust = 0.5, margin = margin(t=0.5, r=0.5, b=0.5, l=0.5, unit = "cm"))
)



ggplot(clades.df, aes(x = eusocial, y = clade, height = after_stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE)+
  theme_ridges()


clades.df$clade=factor(clades.df$clade,levels=rev(c("all diploids","aphids","platypodines","termites",
                                                    "all haplodiploids","halictids","apids","vespids","ants","thrips")))
clades.df$ploidy[clades.df$clade=="all haplodiploids"]="HD"


bins=max(clades.df$eusocial)-min(clades.df$eusocial)

p1=ggplot(clades.df, aes(x = eusocial, y = clade, height = stat(density),
                         fill=ploidy))+
  #alpha=setNames(c(1,rep(0.5,3),1,rep(0.5,5)),clades.df$clade))) + 
  geom_density_ridges(stat = "binline", bins=bins,from = 0, to = 20, scale = 0.7, draw_baseline = TRUE, color = "grey70")+
  #  theme_ridges()+
  #  scale_x_continuous(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 2), 
    name = "Origins of eusociality", 
    expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0),name="Insect clade") +
  coord_cartesian(xlim = c(-0.6, 20)) + 
  #theme_ridges(grid = TRUE, center_axis_labels = TRUE)
  custom_theme +
  #theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_text(vjust = 0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  ggtitle("Eusociality in pruned phylogeny")+
  scale_fill_manual(values = c("DD" = "#4c96ac", 
                               "HD" = "#b05656"),
                    labels = c("DD" = "diploid",
                               "HD" = "haplodiploid"), guide="none")+ 
  scale_alpha_manual(values=rev(c(1,rep(0.4,3),1,rep(0.4,5))),guide=FALSE)

p1
setwd(fig)

ggsave(filename = "eusocial_origins_clades_hapall_eusall_halpruned.svg", plot = p1, device = svglite, path = fig, width = 8, height = 6)





































