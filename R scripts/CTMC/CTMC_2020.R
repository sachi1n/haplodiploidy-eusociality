#-------------------------------------------------------------------------------
# Regime-dependent eusociality transition rates (non-corHMM)
# via:
#   (1) stochastic mapping of ploidy (DD<->HD) on the tree
#   (2) CTMC model for eusociality (S<->E) whose q01 depends on regime
#
# CTMC = Continuous-Time Markov Chain.
# In this script the CTMC is for SOCIALITY (solitary vs eusocial).
#
# The CTMC is defined by a 2x2 rate matrix Q:
#     state 0 = solitary (S)
#     state 1 = eusocial (E)
#
#        [ -q01   q01 ]
#   Q  = [  q10  -q10 ]
#
# where:
#   q01 = instantaneous rate of S -> E (origins of eusociality)
#   q10 = instantaneous rate of E -> S (losses/reversals)
#
# For a branch segment of length t, transition probabilities are:
#   P(t) = expm(Q * t)
# which is exactly what matexpo(Q*t) computes.
#
# Here q01 is allowed to vary by "regime":
#   - DD        : diploid
#   - HDother   : haplodiploid outside Aculeata
#   - Acu       : haplodiploid inside Aculeata
#
# IMPORTANT: HD is NOT treated as a single clade regime.
# Instead, we allow DD<->HD transitions via make.simmap on x,
# then define regimes segment-by-segment on each mapped edge.
#-------------------------------------------------------------------------------


# ---- Load libraries ----------------------------------------------------------
library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(phytools)
library(caper)
#-------------------------------------------------------------------------------


# Set directory shortcuts ------------------------------------------------------
# Modify according to your local environment

# Directory for figures
fig = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures/2020 phylogeny") # ASU desktop
fig = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Figures/updated figures/2020 phylogeny")   # laptop

# Directory for data (phylogeny, phenotype data)
data = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data") # ASU desktop
data = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Data")   # laptop

# Directory for simmap results
res = ("C:/Users/ssures53/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/CTMC") # ASU desktop
res = ("C:/Users/sachi/OneDrive - Arizona State University/Linksvayer Lab/Projects/Haplodiploidy and Eusociality/Final results/CTMC")   # laptop
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


set.seed(1)

# =============================================================
# 1) Recoding helpers
# =============================================================
# Convert pheno.df$Eusocial: "solitary"/"eusocial" -> 0/1
recode_eusocial01 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x == "eusocial"] <- 1L
  out[x == "solitary"] <- 0L
  out
}


# Convert pheno.df$HD.arrhenotoky: "DD"/"HD" -> 0/1
recode_hd01 <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x == "HD"] <- 1L
  out[x == "DD"] <- 0L
  out
}

# =============================================================
# 2) Align tree + data (keeps tips with non-missing x and y)
# =============================================================
align_pair <- function(tree, df,
                       x_col = "HD.arrhenotoky",
                       y_col = "Eusocial",
                       name_col = "name",
                       min_tips = 50) {
  stopifnot(inherits(tree, "phylo"))
  stopifnot(all(c(name_col, x_col, y_col) %in% names(df)))
  
  dat <- df %>%
    transmute(
      .taxon = as.character(.data[[name_col]]),
      x = recode_hd01(.data[[x_col]]),
      y = recode_eusocial01(.data[[y_col]])
    ) %>%
    filter(!is.na(.taxon)) %>%
    filter(.taxon %in% tree$tip.label) %>%
    distinct(.taxon, .keep_all = TRUE) %>%
    filter(!is.na(x), !is.na(y))
  
  common <- intersect(tree$tip.label, dat$.taxon)
  if (length(common) < min_tips) stop("Too few overlapping tips after filtering.")
  
  tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))
  if (is.null(tree2) || Ntip(tree2) < min_tips) stop("Tree became too small after pruning.")
  
  dat2 <- dat %>%
    filter(.taxon %in% tree2$tip.label) %>%
    mutate(.taxon = factor(.taxon, levels = tree2$tip.label)) %>%
    arrange(.taxon) %>%
    mutate(.taxon = as.character(.taxon))
  
  list(tree = tree2, dat = dat2)
}

# =============================================================
# 3) Stage-1 collapse: prune monomorphic clades for BOTH x and y
#    (reduces pseudoreplication; retains only “informative” clades)
# =============================================================
build_children_list <- function(tree) {
  n_tot <- Ntip(tree) + Nnode(tree)
  kids <- vector("list", n_tot)
  for (i in seq_len(nrow(tree$edge))) {
    p <- tree$edge[i, 1]
    c <- tree$edge[i, 2]
    kids[[p]] <- c(kids[[p]], c)
  }
  kids
}

get_root_node <- function(tree) {
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  r <- setdiff(parents, children)
  if (length(r) == 0) stop("Could not determine root node.")
  r[1]
}

collapse_monomorphic_once_pair <- function(tree, dat, seed = 1L) {
  set.seed(seed)
  
  tree <- ape::reorder.phylo(tree, order = "postorder")
  ntip <- Ntip(tree)
  n_tot <- ntip + Nnode(tree)
  
  # bitmask encoding of tip states:
  #   0 -> 1 (bit 1), 1 -> 2 (bit 2)
  # so a node OR’s child bits; monomorphic if mask in {1,2}
  encode_bit <- function(v) ifelse(v == 0, 1L, ifelse(v == 1, 2L, 0L))
  
  x_mask <- integer(n_tot)
  y_mask <- integer(n_tot)
  tip_count <- integer(n_tot)
  rep_tip_index <- integer(n_tot)
  
  x_mask[1:ntip] <- encode_bit(dat$x)
  y_mask[1:ntip] <- encode_bit(dat$y)
  tip_count[1:ntip] <- 1L
  rep_tip_index[1:ntip] <- seq_len(ntip)
  
  kids <- build_children_list(tree)
  nodes_postorder <- unique(tree$edge[, 1])
  
  # Postorder accumulation: compute each internal node’s (x,y) variation
  for (node in nodes_postorder) {
    if (node <= ntip) next
    ch <- kids[[node]]
    if (length(ch) == 0) next
    
    xm <- 0L; ym <- 0L; tc <- 0L
    for (c in ch) {
      xm <- bitwOr(xm, x_mask[c])
      ym <- bitwOr(ym, y_mask[c])
      tc <- tc + tip_count[c]
    }
    x_mask[node] <- xm
    y_mask[node] <- ym
    tip_count[node] <- tc
    
    # If we collapse a monomorphic clade, pick ONE representative tip
    # weighted by descendant tip counts (so big clades aren't underweighted).
    w <- tip_count[ch]
    s <- sum(w)
    chosen_child <- if (s > 0) sample(ch, 1, prob = w / s) else ch[1]
    rep_tip_index[node] <- rep_tip_index[chosen_child]
  }
  
  is_mono <- function(mask) mask %in% c(1L, 2L)
  mono_pair <- is_mono(x_mask) & is_mono(y_mask)
  
  root <- get_root_node(tree)
  
  kept_tip_indices <- integer(0)
  stack <- c(root)
  
  # Traverse from root: if node is monomorphic for BOTH traits => keep 1 rep tip
  while (length(stack) > 0) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (node <= ntip) {
      kept_tip_indices <- c(kept_tip_indices, node)
    } else {
      if (mono_pair[node]) {
        kept_tip_indices <- c(kept_tip_indices, rep_tip_index[node])
      } else {
        stack <- c(stack, kids[[node]])
      }
    }
  }
  
  kept_tip_indices <- sort(unique(kept_tip_indices))
  kept_tips <- tree$tip.label[kept_tip_indices]
  
  tree_out <- drop.tip(tree, setdiff(tree$tip.label, kept_tips))
  dat_out <- dat %>%
    filter(.taxon %in% tree_out$tip.label) %>%
    mutate(.taxon = factor(.taxon, levels = tree_out$tip.label)) %>%
    arrange(.taxon) %>%
    mutate(.taxon = as.character(.taxon))
  
  list(tree = tree_out, dat = dat_out,
       n_before = ntip, n_after = Ntip(tree_out))
}

# =============================================================
# 4) Define Aculeata tips (your approach; uses the FULL tree)
# =============================================================
get_aculeata_tips <- function(insect.tree, pheno.df) {
  aculeate.families <- unique(subset(pheno.df, order=="Hymenoptera" & family!="")$family)
  
  tips_for_mrca <- row.names(subset(pheno.df, family %in% aculeate.families))
  if (is.null(tips_for_mrca) || length(tips_for_mrca) < 2) stop("Too few row.names tips for Aculeata MRCA.")
  
  tips_for_mrca <- intersect(tips_for_mrca, insect.tree$tip.label)
  if (length(tips_for_mrca) < 2) stop("Too few overlapping tips between row.names(pheno.df) and insect.tree tip labels.")
  
  node <- as.numeric(findMRCA(insect.tree, tips = tips_for_mrca))
  caper::clade.members(node, insect.tree, tip.labels=TRUE, include.nodes=TRUE)$tips
}


# =============================================================
# 5) CTMC machinery (THIS is where CTMC is defined/used)
# =============================================================

# -------------------------
# 5A) CTMC definition: build Q and convert to P(t)
# -------------------------
# This is the 2-state CTMC generator Q for eusociality:
#   state 0 = solitary, state 1 = eusocial
#
# Q = [[-q01, q01],
#      [ q10,-q10]]
#
# Then we compute the transition probability matrix:
#   P(t) = expm(Q*t)
#
# This is the fundamental CTMC object used in the likelihood.
P_2state <- function(q01, q10, t) {
  Q <- matrix(c(-q01, q01,
                q10, -q10), 2, 2, byrow = TRUE)
  ape::matexpo(Q * t)  # matrix exponential = CTMC transition probabilities
}

# -------------------------
# 5B) Use CTMC on each branch segment, with regime-specific q01
# -------------------------
# A simmap edge is split into segments with mapped HD states,
# e.g. names(edge_map) might be "DD","HD" or "0","1".
#
# For each segment we choose the regime-specific q01:
#   if segment is DD: q01 = q01_DD
#   if segment is HD and inside Acu: q01 = q01_Acu
#   if segment is HD and outside Acu: q01 = q01_HDother
#
# Then we multiply segment transition matrices along the edge:
#   P_edge = Π_k P_segment_k
edge_P_from_segments <- function(edge_map, is_in_acu_child,
                                 q01_DD, q01_HDother, q01_Acu, q10) {
  P <- diag(2)  # start with identity
  
  for (k in seq_along(edge_map)) {
    st <- names(edge_map)[k]
    t  <- as.numeric(edge_map[k])
    
    # interpret mapped state as HD or DD
    hd1 <- !(st %in% c("0", "DD", "dd", "Diploid", "diploid"))
    
    if (!hd1) {
      q01 <- q01_DD
    } else {
      q01 <- if (isTRUE(is_in_acu_child)) q01_Acu else q01_HDother
    }
    
    # CTMC step: compute P(t) for this segment and multiply
    Pseg <- P_2state(q01 = q01, q10 = q10, t = t)
    P <- P %*% Pseg
  }
  P
}


# -------------------------
# 5C) CTMC likelihood via pruning algorithm
# -------------------------
# Standard Felsenstein pruning for discrete traits:
# - Tip likelihood vectors are 1-hot for observed states
# - Propagate likelihoods up the tree using P_edge matrices
# - Combine at root with a root prior pi_root
#
# CTMC is used here because P_edge is derived from expm(Q*t).
loglik_eusocial_simmap <- function(tree_map, y_tip01, in_acu_tip,
                                   q01_DD, q01_HDother, q01_Acu, q10,
                                   pi_root = c(0.9, 0.1)) {
  tree <- ape::reorder.phylo(tree_map, order = "postorder")
  
  ntip <- Ntip(tree)
  nnode <- Nnode(tree)
  n_tot <- ntip + nnode
  
  # Likelihood vectors L[node, state] with states = {0,1}
  L <- matrix(1, n_tot, 2)
  
  # tips: 1-hot encoding
  for (i in seq_len(ntip)) {
    tip <- tree$tip.label[i]
    yi <- y_tip01[tip]
    if (is.na(yi)) stop("Missing y at tip in likelihood.")
    L[i, ] <- if (yi == 0) c(1, 0) else c(0, 1)
  }
  
  kids <- build_children_list(tree)
  
  # Determine which internal nodes are fully inside Aculeata
  # (used to classify edges by child membership)
  inA_mask <- logical(n_tot)
  tipcount <- integer(n_tot)
  
  inA_mask[1:ntip] <- as.logical(in_acu_tip[tree$tip.label])
  tipcount[1:ntip] <- 1L
  
  nodes_postorder <- unique(tree$edge[, 1])
  for (node in nodes_postorder) {
    if (node <= ntip) next
    ch <- kids[[node]]
    if (length(ch) == 0) next
    tipcount[node] <- sum(tipcount[ch])
    inA_mask[node] <- all(inA_mask[ch])  # strict “inside Acu” definition
  }
  
  # Pruning: compute likelihoods at internal nodes
  for (node in nodes_postorder) {
    if (node <= ntip) next
    ch <- kids[[node]]
    if (length(ch) == 0) next
    
    for (cnode in ch) {
      # identify edge parent=node -> child=cnode
      e_idx <- which(tree$edge[,1] == node & tree$edge[,2] == cnode)
      if (length(e_idx) != 1) stop("Edge index lookup failed.")
      
      edge_map <- tree$maps[[e_idx]]  # simmap segment lengths
      
      # edge regime depends on whether CHILD is inside Aculeata
      is_in_acu_child <- inA_mask[cnode]
      
      # CTMC on this edge: build P_edge from segment-wise expm(Q*t)
      P <- edge_P_from_segments(edge_map, is_in_acu_child,
                                q01_DD=q01_DD, q01_HDother=q01_HDother, q01_Acu=q01_Acu, q10=q10)
      
      # propagate: parent likelihood multiplied by sum over child states
      L_child <- L[cnode, ]
      contrib <- as.numeric(P %*% L_child)
      L[node, ] <- L[node, ] * contrib
    }
  }
  
  # root prior * root likelihood
  root <- get_root_node(tree)
  lik_root <- sum(pi_root * L[root, ])
  if (!is.finite(lik_root) || lik_root <= 0) return(-Inf)
  log(lik_root)
}


# =============================================================
# 6) CTMC parameter fitting (MLE)
# =============================================================
# Here we estimate q01 and q10 (and root state frequency) by maximizing
# the CTMC likelihood above.
#
# Multi-start optim is used because likelihoods can have local maxima.
make_negloglik <- function(model, tree_map, y_tip01, in_acu_tip) {
  force(model)
  function(par_free) {
    if (model == "M0") {
      # free params: q01, q10, pi_root_eus
      q01 <- exp(par_free[1]); q10 <- exp(par_free[2])
      pi_root <- c(1 - plogis(par_free[3]), plogis(par_free[3]))
      ll <- loglik_eusocial_simmap(tree_map, y_tip01, in_acu_tip,
                                   q01_DD=q01, q01_HDother=q01, q01_Acu=q01, q10=q10,
                                   pi_root=pi_root)
      return(-ll)
    }
    if (model == "Mploidy") {
      # q01 differs for DD vs HD (HDother + Acu share q01_HD)
      q01_DD <- exp(par_free[1]); q01_HD <- exp(par_free[2]); q10 <- exp(par_free[3])
      pi_root <- c(1 - plogis(par_free[4]), plogis(par_free[4]))
      ll <- loglik_eusocial_simmap(tree_map, y_tip01, in_acu_tip,
                                   q01_DD=q01_DD, q01_HDother=q01_HD, q01_Acu=q01_HD, q10=q10,
                                   pi_root=pi_root)
      return(-ll)
    }
    if (model == "Mfull") {
      # q01 differs for DD vs HDother vs Acu
      q01_DD <- exp(par_free[1]); q01_HDother <- exp(par_free[2]); q01_Acu <- exp(par_free[3]); q10 <- exp(par_free[4])
      pi_root <- c(1 - plogis(par_free[5]), plogis(par_free[5]))
      ll <- loglik_eusocial_simmap(tree_map, y_tip01, in_acu_tip,
                                   q01_DD=q01_DD, q01_HDother=q01_HDother, q01_Acu=q01_Acu, q10=q10,
                                   pi_root=pi_root)
      return(-ll)
    }
    stop("Unknown model")
  }
}


fit_one_model <- function(model, tree_map, y_tip01, in_acu_tip, nstarts = 50) {
  negll <- make_negloglik(model, tree_map, y_tip01, in_acu_tip)
  
  best <- list(value = Inf, par = NULL, conv = NA, logLik = NA)
  
  for (s in seq_len(nstarts)) {
    # random starts in log-rate space; pi_root_eus in (0.01, 0.5)
    if (model == "M0") {
      p0 <- c(runif(1, -6, 2), runif(1, -6, 2), qlogis(runif(1, 0.01, 0.5)))
      k <- 3
    } else if (model == "Mploidy") {
      p0 <- c(runif(1, -6, 2), runif(1, -6, 2), runif(1, -6, 2), qlogis(runif(1, 0.01, 0.5)))
      k <- 4
    } else {
      p0 <- c(runif(1, -6, 2), runif(1, -6, 2), runif(1, -6, 2), runif(1, -6, 2), qlogis(runif(1, 0.01, 0.5)))
      k <- 5
    }
    
    # MLE for CTMC parameters:
    fit <- tryCatch(
      optim(p0, negll, method = "Nelder-Mead",
            control = list(maxit = 5000, reltol = 1e-10)),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    if (is.finite(fit$value) && fit$value < best$value) {
      best$value <- fit$value
      best$par <- fit$par
      best$conv <- fit$convergence
      best$logLik <- -fit$value
      best$k <- k
    }
  }
  
  if (!is.finite(best$logLik)) stop(paste("All starts failed for model", model))
  
  # decode estimates into comparable columns
  if (model == "M0") {
    q01 <- exp(best$par[1]); q10 <- exp(best$par[2]); pi_e <- plogis(best$par[3])
    est <- tibble(q01_DD=q01, q01_HDother=q01, q01_Acu=q01, q10=q10, pi_root_eus=pi_e)
  } else if (model == "Mploidy") {
    q01_DD <- exp(best$par[1]); q01_HD <- exp(best$par[2]); q10 <- exp(best$par[3]); pi_e <- plogis(best$par[4])
    est <- tibble(q01_DD=q01_DD, q01_HDother=q01_HD, q01_Acu=q01_HD, q10=q10, pi_root_eus=pi_e)
  } else {
    q01_DD <- exp(best$par[1]); q01_HDother <- exp(best$par[2]); q01_Acu <- exp(best$par[3]); q10 <- exp(best$par[4]); pi_e <- plogis(best$par[5])
    est <- tibble(q01_DD=q01_DD, q01_HDother=q01_HDother, q01_Acu=q01_Acu, q10=q10, pi_root_eus=pi_e)
  }
  
  list(model=model, logLik=best$logLik, k=best$k,
       AIC = 2*best$k - 2*best$logLik,
       est=est, conv=best$conv)
}

# =============================================================
# 7) Pipeline: align -> collapse -> map HD -> fit CTMC repeatedly
# =============================================================

cat("\n============================================================\n")
cat("Aculeata vs non-Aculeata HD vs DD (simmap + CTMC)\n")
cat("============================================================\n")

# A) Align + collapse
base <- align_pair(insect.tree, pheno.df,
                   x_col="HD.arrhenotoky", y_col="Eusocial")
st1  <- collapse_monomorphic_once_pair(base$tree, base$dat, seed=1L)

treeC <- st1$tree
datC  <- st1$dat

cat("Stage-1 collapsed tips:", Ntip(treeC), " (", st1$n_before, " -> ", st1$n_after, ")\n", sep="")
cat("Counts x (DD/HD): ", sum(datC$x==0), "/", sum(datC$x==1), "\n", sep="")
cat("Counts y (sol/eus): ", sum(datC$y==0), "/", sum(datC$y==1), "\n", sep="")

saveRDS(treeC, "CTMC_Aculeata_collapsed_tree.rds")
saveRDS(datC,  "CTMC_Aculeata_collapsed_dat.rds")

# B) Define Aculeata on full tree, intersect with collapsed tree tips
aculeates_all <- get_aculeata_tips(insect.tree, pheno.df)
aculeates_C <- intersect(aculeates_all, treeC$tip.label)

cat("Aculeata tips on collapsed tree:", length(aculeates_C), "\n")

# indicator for each collapsed tip: inside Aculeata clade or not
in_acu_tip <- setNames(treeC$tip.label %in% aculeates_C, treeC$tip.label)

# C) Prepare tip states for mapping + CTMC likelihood
x_tip <- setNames(datC$x, datC$.taxon)  # 0/1 DD/HD
y_tip <- setNames(datC$y, datC$.taxon)  # 0/1 S/E

# D) Stochastic map x (DD<->HD) on the collapsed tree
# This generates trees where each edge is split into segments with an HD state.
x_fac <- factor(ifelse(x_tip[treeC$tip.label] == 1, "HD", "DD"), levels=c("DD","HD"))
names(x_fac) <- treeC$tip.label

cat("\nSimulating HD maps...\n")
nsim_maps <- 50
hd_maps <- phytools::make.simmap(treeC, x_fac,
                                 model="ER",      # simple mapping model for HD
                                 nsim=nsim_maps,
                                 pi="estimated",
                                 message=FALSE)


# E) For each mapped HD history, fit the eusociality CTMC models
# This is where the CTMC is repeatedly fit on different inferred regimes.
nstarts <- 50
fit_rows <- list()

pb <- txtProgressBar(min=0, max=nsim_maps, style=3)
for (i in seq_len(nsim_maps)) {
  trm <- hd_maps[[i]]  # a simmap tree with branch segments labeled DD/HD
  
  # Fit CTMC variants
  f0 <- fit_one_model("M0",     trm, y_tip01 = y_tip, in_acu_tip = in_acu_tip, nstarts = nstarts)
  fp <- fit_one_model("Mploidy",trm, y_tip01 = y_tip, in_acu_tip = in_acu_tip, nstarts = nstarts)
  ff <- fit_one_model("Mfull",  trm, y_tip01 = y_tip, in_acu_tip = in_acu_tip, nstarts = nstarts)
  
  # store results (one row per model per map)
  fit_rows[[length(fit_rows)+1]] <- tibble(
    map = i, model = f0$model, logLik = f0$logLik, k=f0$k, AIC=f0$AIC,
    q01_DD = f0$est$q01_DD, q01_HDother = f0$est$q01_HDother, q01_Acu = f0$est$q01_Acu,
    q10 = f0$est$q10, pi_root_eus = f0$est$pi_root_eus
  )
  fit_rows[[length(fit_rows)+1]] <- tibble(
    map = i, model = fp$model, logLik = fp$logLik, k=fp$k, AIC=fp$AIC,
    q01_DD = fp$est$q01_DD, q01_HDother = fp$est$q01_HDother, q01_Acu = fp$est$q01_Acu,
    q10 = fp$est$q10, pi_root_eus = fp$est$pi_root_eus
  )
  fit_rows[[length(fit_rows)+1]] <- tibble(
    map = i, model = ff$model, logLik = ff$logLik, k=ff$k, AIC=ff$AIC,
    q01_DD = ff$est$q01_DD, q01_HDother = ff$est$q01_HDother, q01_Acu = ff$est$q01_Acu,
    q10 = ff$est$q10, pi_root_eus = ff$est$pi_root_eus
  )
  
  setTxtProgressBar(pb, i)
}
close(pb)

fits_df <- bind_rows(fit_rows)


# F) Summarize across maps (medians)
summ_tbl <- fits_df %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    logLik_med = median(logLik),
    AIC_med    = median(AIC),
    q01_DD_med = median(q01_DD),
    q01_HDother_med = median(q01_HDother),
    q01_Acu_med = median(q01_Acu),
    q10_med = median(q10),
    .groups="drop"
  ) %>%
  dplyr::arrange(AIC_med) %>%
  dplyr::mutate(deltaAIC_med = AIC_med - min(AIC_med),
                wAIC_med = exp(-0.5*deltaAIC_med)/sum(exp(-0.5*deltaAIC_med)))

cat("\n===== Across-simmap summary (medians) =====\n")
print(summ_tbl)

setwd(res)
write.csv(
  summ_tbl,
  file = "FigC_model_comparison_summary_2020.csv",
  row.names = FALSE
)

# G) LRT per map for nested models (CTMC likelihood comparison)
lrt_df <- fits_df %>%
  dplyr::select(map, model, logLik, k) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(logLik, k)) %>%
  dplyr::mutate(
    LR_M0_vs_Mploidy = 2 * (logLik_Mploidy - logLik_M0),
    df_M0_vs_Mploidy = (k_Mploidy - k_M0),
    p_M0_vs_Mploidy  = stats::pchisq(LR_M0_vs_Mploidy, df = df_M0_vs_Mploidy, lower.tail = FALSE),
    
    LR_Mploidy_vs_Mfull = 2 * (logLik_Mfull - logLik_Mploidy),
    df_Mploidy_vs_Mfull = (k_Mfull - k_Mploidy),
    p_Mploidy_vs_Mfull  = stats::pchisq(LR_Mploidy_vs_Mfull, df = df_Mploidy_vs_Mfull, lower.tail = FALSE)
  )

lrt_summ <- tibble::tibble(
  comparison = c("M0 vs Mploidy", "Mploidy vs Mfull"),
  LR_median = c(stats::median(lrt_df$LR_M0_vs_Mploidy, na.rm = TRUE),
                stats::median(lrt_df$LR_Mploidy_vs_Mfull, na.rm = TRUE)),
  p_median  = c(stats::median(lrt_df$p_M0_vs_Mploidy, na.rm = TRUE),
                stats::median(lrt_df$p_Mploidy_vs_Mfull, na.rm = TRUE)),
  p_prop_lt_0.05 = c(mean(lrt_df$p_M0_vs_Mploidy < 0.05, na.rm = TRUE),
                     mean(lrt_df$p_Mploidy_vs_Mfull < 0.05, na.rm = TRUE))
)

cat("\n===== LRT summary across simmaps =====\n")
print(lrt_summ)

# H) Rate ratios (from Mfull) across maps
ratios_df <- fits_df %>%
  dplyr::filter(model == "Mfull") %>%
  dplyr::transmute(
    map,
    q01_Acu, q01_HDother, q01_DD,
    ratio_Acu_vs_HDother = q01_Acu / q01_HDother,
    ratio_Acu_vs_DD      = q01_Acu / q01_DD,
    ratio_HDother_vs_DD  = q01_HDother / q01_DD
  )

rat_summ <- ratios_df %>%
  dplyr::summarise(
    ratio_Acu_vs_HDother_med  = stats::median(ratio_Acu_vs_HDother, na.rm = TRUE),
    ratio_Acu_vs_DD_med       = stats::median(ratio_Acu_vs_DD, na.rm = TRUE),
    ratio_HDother_vs_DD_med   = stats::median(ratio_HDother_vs_DD, na.rm = TRUE),
    ratio_Acu_vs_HDother_q2.5  = stats::quantile(ratio_Acu_vs_HDother, 0.025, na.rm = TRUE),
    ratio_Acu_vs_HDother_q97.5 = stats::quantile(ratio_Acu_vs_HDother, 0.975, na.rm = TRUE),
    ratio_Acu_vs_DD_q2.5       = stats::quantile(ratio_Acu_vs_DD, 0.025, na.rm = TRUE),
    ratio_Acu_vs_DD_q97.5      = stats::quantile(ratio_Acu_vs_DD, 0.975, na.rm = TRUE),
    ratio_HDother_vs_DD_q2.5   = stats::quantile(ratio_HDother_vs_DD, 0.025, na.rm = TRUE),
    ratio_HDother_vs_DD_q97.5  = stats::quantile(ratio_HDother_vs_DD, 0.975, na.rm = TRUE)
  )

cat("\n===== Mfull rate-ratio summary across simmaps =====\n")
print(rat_summ)

cat("\nDONE.\n")



library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(phytools)
library(ape)

# Aculeata tip indicator
in_acu_tip <- setNames(treeC$tip.label %in% aculeates_C, treeC$tip.label)

# Tip annotation for ring
ann <- datC %>%
  filter(.taxon %in% treeC$tip.label) %>%
  transmute(
    label    = as.character(.taxon),
    ploidy01 = as.integer(x),   # 0 = DD, 1 = HD
    social01 = as.integer(y)    # 0 = Solitary, 1 = Eusocial
  ) %>%
  distinct(label, .keep_all = TRUE) %>%
  mutate(
    Sociality = factor(social01, levels = c(0, 1), labels = c("Solitary", "Eusocial"))
  )

# -----------------------------
# Helpers
# -----------------------------
build_children_list <- function(tree) {
  n_tot <- Ntip(tree) + Nnode(tree)
  kids <- vector("list", n_tot)
  for (i in seq_len(nrow(tree$edge))) {
    p <- tree$edge[i, 1]
    c <- tree$edge[i, 2]
    kids[[p]] <- c(kids[[p]], c)
  }
  kids
}

get_root_node <- function(tree) {
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  r <- setdiff(parents, children)
  if (length(r) == 0) stop("Could not determine root node.")
  r[1]
}

compute_node_summaries <- function(tree, tip_ploidy01, in_acu_tip) {
  tree <- ape::reorder.phylo(tree, "postorder")
  ntip <- Ntip(tree)
  n_tot <- ntip + Nnode(tree)
  
  kids <- build_children_list(tree)
  nodes_postorder <- unique(tree$edge[, 1])
  
  encode_bit <- function(v) ifelse(v == 0, 1L, ifelse(v == 1, 2L, 0L))
  
  ploidy_mask <- integer(n_tot)
  hd_ct <- integer(n_tot)
  dd_ct <- integer(n_tot)
  
  inA_mask <- logical(n_tot)
  
  # tips
  ploidy_mask[1:ntip] <- encode_bit(tip_ploidy01)
  hd_ct[1:ntip] <- as.integer(tip_ploidy01 == 1)
  dd_ct[1:ntip] <- as.integer(tip_ploidy01 == 0)
  
  inA_mask[1:ntip] <- as.logical(in_acu_tip[tree$tip.label])
  
  # postorder
  for (node in nodes_postorder) {
    if (node <= ntip) next
    ch <- kids[[node]]
    if (length(ch) == 0) next
    
    pm <- 0L; h <- 0L; d <- 0L
    inA_all <- TRUE
    
    for (c in ch) {
      pm <- bitwOr(pm, ploidy_mask[c])
      h  <- h + hd_ct[c]
      d  <- d + dd_ct[c]
      inA_all <- inA_all & inA_mask[c]
    }
    
    ploidy_mask[node] <- pm
    hd_ct[node] <- h
    dd_ct[node] <- d
    inA_mask[node] <- inA_all
  }
  
  list(tree = tree, ploidy_mask = ploidy_mask, hd_ct = hd_ct, dd_ct = dd_ct, inA_mask = inA_mask)
}

assign_edge_regimes_Acu_HD_DD <- function(tree, ann, in_acu_tip) {
  tip_x <- ann$ploidy01
  names(tip_x) <- ann$label
  tip_x <- tip_x[tree$tip.label]
  if (any(is.na(tip_x))) stop("Missing ploidy01 for some tips in ann.")
  
  ns <- compute_node_summaries(tree, tip_x, in_acu_tip)
  tree2 <- ns$tree
  
  child <- tree2$edge[, 2]
  edge_regime <- character(length(child))
  
  for (e in seq_along(child)) {
    nd <- child[e]
    
    # ploidy at node: monomorphic => exact; polymorphic => majority rule
    m <- ns$ploidy_mask[nd]
    is_HD <- if (m %in% c(1L, 2L)) {
      (m == 2L)
    } else {
      ns$hd_ct[nd] >= ns$dd_ct[nd]
    }
    
    if (!is_HD) {
      edge_regime[e] <- "DD"
    } else {
      edge_regime[e] <- if (isTRUE(ns$inA_mask[nd])) "Acu" else "HD_other"
    }
  }
  
  list(tree = tree2, edge_regime = edge_regime)
}

# -----------------------------
# 1) Assign regimes to edges
# -----------------------------
reg_out <- assign_edge_regimes_Acu_HD_DD(treeC, ann, in_acu_tip)
treeC <- reg_out$tree
edge_regime <- reg_out$edge_regime

# -----------------------------
# 2) OPTION A FIX: Force basal edges to DD (visual-only)
# -----------------------------
root <- get_root_node(treeC)
basal_edges <- which(treeC$edge[, 1] == root)   # edges out of root
edge_regime[basal_edges] <- "DD"

# -----------------------------
# 3) Heatmap (ring) data
# -----------------------------
heatmap_data <- ann %>%
  filter(label %in% treeC$tip.label) %>%
  transmute(label = as.character(label),
            Sociality = as.factor(Sociality)) %>%
  distinct(label, .keep_all = TRUE)

heatmap_mat <- heatmap_data["Sociality", drop = FALSE]
rownames(heatmap_mat) <- heatmap_data$label

# -----------------------------
# 4) Plot
# -----------------------------
COL_ACU <- "#B22222"  # muted brick red
COL_HD  <- "#2C7FB8"  # muted steel blue
COL_DD  <- "grey20"

reg_df <- tibble(
  node = treeC$edge[, 2],
  regime = factor(edge_regime, levels = c("Acu", "HD_other", "DD"))
)

RING_OFFSET <- 10
RING_WIDTH  <- 0.04

p <- ggtree(treeC, layout = "fan", open.angle = 5)

p$data <- p$data %>% left_join(reg_df, by = "node")

p <- p +
  geom_tree(aes(color = regime),
            linewidth = 0.8,
            lineend = "round",
            linejoin = "round",
            na.rm = TRUE) +
  scale_color_manual(
    name   = "Branch regime",
    values = c(Acu = COL_ACU, HD_other = COL_HD, DD = COL_DD),
    labels = c("Aculeate haplodiploids", "Non-aculeate haplodiploids", "Diploids")
  )

p <- gheatmap(
  p, heatmap_mat,
  offset = RING_OFFSET,
  width  = RING_WIDTH,
  colnames = FALSE,
  legend_title = "Tip sociality"
)

p <- p +
  scale_fill_manual(
    name = "Tip sociality",
    values = c(Solitary = "grey85", Eusocial = "grey45"),
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.6)),
    fill  = guide_legend(order = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(6, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(10, "mm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

print(p)

setwd(fig)
list.files()
# Save (edit path/filename as you like)
svg("CTMC_phylogeny_2020.svg", width = 10, height = 6)
print(p)
dev.off()



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

fits_df <- bind_rows(fit_rows)

# -----------------------------
# 1) Keep only Mfull (this is the model with separate q01_DD, q01_HDother, q01_Acu)
# -----------------------------
mfull_df <- fits_df %>%
  dplyr::filter(model == "Mfull") %>%
  dplyr::select(map, q01_DD, q01_HDother, q01_Acu)

# -----------------------------
# 2) Summarize across simmaps (median + 95% interval)
#    (Use quantiles because distributions can be skewed)
# -----------------------------
figB_summ <- mfull_df %>%
  pivot_longer(cols = c(q01_DD, q01_HDother, q01_Acu),
               names_to = "regime", values_to = "q01") %>%
  mutate(
    regime = recode(
      regime,
      q01_DD      = "Diploid",
      q01_HDother = "Non-aculeate haplodiploids",
      q01_Acu     = "Aculeate haplodiploids"
    ),
    regime = factor(regime, levels = c("Diploid",
                                       "Non-aculeate haplodiploids",
                                       "Aculeate haplodiploids"))
  ) %>%
  group_by(regime) %>%
  summarise(
    q01_med = median(q01, na.rm = TRUE),
    q01_lo  = quantile(q01, 0.025, na.rm = TRUE),
    q01_hi  = quantile(q01, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

figB_summ

# (Optional) save for reproducibility
write.csv(figB_summ, "FigB_rates_summary_Mfull_acu_nonacu_DD.csv", row.names = FALSE)

# -----------------------------
# 3) Plot (no title, with border, consistent styling)
# -----------------------------
COL_DD  <- "grey20"
COL_NON <- "#2C7FB8"  # your muted blue
COL_ACU <- "#B22222"   # your muted brick red

pB <- ggplot(figB_summ, aes(x = regime, y = q01_med, color = regime)) +
  geom_errorbar(aes(ymin = q01_lo, ymax = q01_hi),
                width = 0.12, linewidth = 0.6) +
  geom_point(
    size = 1.5) +
  scale_color_manual(
    values = c(
      "Diploid" = COL_DD,
      "Non-aculeate haplodiploids" = COL_NON,
      "Aculeate haplodiploids" = COL_ACU
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1),
    limits = c(0, 4),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "q01 (solitary → eusocial per branch-length unit)"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 11, color = "black"),
    axis.ticks   = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin  = margin(2, 2, 2, 2, unit = "mm")
  )

print(pB)
# Save
svg("CTMC_eus_rate_2020_v2.svg", width = 6.5, height = 4.5)
print(pB)
dev.off()


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Fig C: ΔAIC across models
# (uses Aculeata-based fits_df from the simmap CTMC script)
# -----------------------------

# 1) Summarize AIC across simmaps (median AIC per model)
figC_summ <- fits_df %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    AIC_med = median(AIC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(deltaAIC = AIC_med - min(AIC_med, na.rm = TRUE))

# 2) Set explicit order (bottom -> top after coord_flip)
model_order <- c("Mfull", "Mploidy", "M0")

figC_summ <- figC_summ %>%
  dplyr::filter(model %in% model_order) %>%
  dplyr::mutate(model = factor(model, levels = model_order))

# 3) Plot
pC <- ggplot(figC_summ, aes(x = model, y = deltaAIC)) +
  geom_col(width = 0.35, fill = "#8A7A7A") +   # thinner bars
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)) 
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  ) +
  labs(
    x = NULL,
    y = expression(Delta*AIC)
  )

print(pC)

# Optional save
svg("CTMC_model_comparison_2020.svg", width = 6.5, height = 4.5)
print(pC)
dev.off()




