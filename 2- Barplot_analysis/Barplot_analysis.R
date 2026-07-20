# ======================
# load libraries
# ======================

library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)

# ================
# function
# ================

prepare_resolution_counts <- function(win_salmonidae, win_salmoninae, win_NAME11, win_NAME13) {
  
  ############################
  # 1 AORE Salmonidae
  ############################
  
  AORE_salmonidae <- sum(
    win_salmonidae$total[win_salmonidae$state == "AORE"]
  )
  
  
  ############################
  # 2 LORE salmonidae - LORE salmoninae
  ############################
  
  dae_nae <- merge(win_salmonidae, win_salmoninae,
                   by = c("chr","win_start","win_end"),
                   suffixes = c("_dae","_nae"))
  
  # keep LORE in dae but not in nae
  salmonidae_specific_windows <- dae_nae[
    dae_nae$state_dae == "LORE" &
      dae_nae$state_nae != "LORE",
  ]
  
  LORE_dae_resolution <- sum(
    salmonidae_specific_windows$total_dae
  )
  
  
  ############################
  # 3 LORE salmoninae - LORE NAME11
  # LORE NAME11 = LORE NAME11 + NAME12
  ############################
  
  nae_nm11 <- merge(win_salmoninae, win_NAME11,
                    by = c("chr","win_start","win_end"),
                    suffixes = c("_nae","_nm11"))
  
  salmoninae_specific_windows <- nae_nm11[
    nae_nm11$state_nae == "LORE" &
      nae_nm11$state_nm11 != "LORE",
  ]
  
  LORE_nae_resolution <- sum(
    salmoninae_specific_windows$total_nae
  )
  
  
  ############################
  # 4 LORE NAME11 - LORE NAME13
  # LORE NAME13 = LORE NAME13 + NAME14
  ############################
  
  nm11_nm13 <- merge(win_NAME11, win_NAME13,
                     by = c("chr","win_start","win_end"),
                     suffixes = c("_nm11","_nm13"))
  
  NAME11_specific_windows <- nm11_nm13[
    nm11_nm13$state_nm11 == "LORE" &
      nm11_nm13$state_nm13 != "LORE",
  ]
  
  LORE_nm11_resolution <- sum(
    NAME11_specific_windows$total_nm11
  )
  
  ##########################
  # 5 LORE NAME13 
  ##########################
  
  LORE_nm13 <- sum(
    win_NAME13$total[win_NAME13$state == "LORE"]
  )
  
  
  
  ############################
  # final results
  ############################
  
  return(data.frame(
    category = c("AORE_Salmonidae",
                 "LORE_resolution_Salmonidae",
                 "LORE_resolution_Salmoninae",
                 "LORE_resolution_NAME11",
                 "LORE_NAME13"),
    
    total = c(AORE_salmonidae,
              LORE_dae_resolution,
              LORE_nae_resolution,
              LORE_nm11_resolution,
              LORE_nm13)
  ))
}

get_resolved_esoxIDs <- function(win_salmonidae, win_salmoninae, win_NAME11, win_NAME13) {
  
  # ================================================
  # 1. esoxID AORE in salmonidae
  # ================================================
  
  aore_windows <- win_salmonidae[win_salmonidae$state == "AORE", ]
  
  esoxID_AORE_salmonidae <- unique(unlist(aore_windows$esoxIDs))
  
  
  # ================================================
  # 2. esoxID LORE resolved salmonidae → salmoninae
  #    = LORE in dae but not LORE in nae
  # ================================================
  
  dae_nae <- merge(win_salmonidae, win_salmoninae,
                   by = c("chr", "win_start", "win_end"),
                   suffixes = c("_dae", "_nae"))
  
  salmonidae_specific <- dae_nae[
    dae_nae$state_dae == "LORE" & dae_nae$state_nae != "LORE", 
  ]
  
  # extract and flatten the esoxID of windows
  esoxID_LORE_salmonidae <- unique(unlist(salmonidae_specific$esoxIDs_dae))
  
  
  # ================================================
  # 3. esoxID LORE resolved salmoninae → NAME11
  #    = LORE in nae but not LORE in NAME11
  # ================================================
  
  nae_nm11 <- merge(win_salmoninae, win_NAME11,
                    by = c("chr", "win_start", "win_end"),
                    suffixes = c("_nae", "_nm11"))
  
  salmoninae_specific <- nae_nm11[
    nae_nm11$state_nae == "LORE" & nae_nm11$state_nm11 != "LORE", 
  ]
  
  esoxID_LORE_nae <- unique(unlist(salmoninae_specific$esoxIDs_nae))
  
  
  # ===============================================
  # 4 esoxID LORE resolved NAME11 -> NAME13
  #    = LORE in NAME11 but not in NAME13
  # ===============================================
  
  nm11_nm13 <- merge(win_NAME11, win_NAME13,
                     by = c("chr", "win_start", "win_end"),
                     suffixes = c("_nm11", "_nm13"))
  
  NAME11_specific <- nm11_nm13[
    nm11_nm13$state_nm11 == "LORE" & nm11_nm13$state_nm13 != "LORE", 
  ]
  
  esoxID_LORE_nm11 <- unique(unlist(NAME11_specific$esoxIDs_nm11))
  
  
  # ===============================================
  # 5 esoxID LORE NAME13 
  #    = LORE in NAME13 
  # ===============================================
  
  NAME13_windows <- win_NAME13[win_NAME13$state == "LORE", ]
  
  esoxID_LORE_nm13 <- unique(unlist(NAME13_windows$esoxIDs))
  
  
  # ================================================
  # Résultat 
  # ================================================
  
  return(list(
    AORE_Salmonidae            = esoxID_AORE_salmonidae,
    LORE_resolution_Salmonidae = esoxID_LORE_salmonidae,
    LORE_resolution_Salmoninae = esoxID_LORE_nae,
    LORE_resolution_NAME11     = esoxID_LORE_nm11,
    LORE_NAME13                = esoxID_LORE_nm13
  ))
}

get_resolved_esoxIDs_circos <- function(genes_salmonidae, genes_salmoninae, genes_NAME11, genes_NAME13) {
  # ================================================
  # 1. esoxID AORE in salmonidae
  # ================================================
  
  esoxID_AORE_salmonidae <- unique(genes_salmonidae$esoxID[genes_salmonidae$class == "AORE"])
  esoxID_LORE_salmonidae <- unique(genes_salmonidae$esoxID[genes_salmonidae$class == "LORE"])
  
  # ================================================
  # 2. esoxID LORE resolved salmonidae → salmoninae
  #    = LORE in dae but not LORE in nae
  # ================================================
  
  esoxID_LORE_salmoninae <- unique(genes_salmoninae$esoxID[genes_salmoninae$class == "LORE"])
  
  resolved_dae <- setdiff(esoxID_LORE_salmonidae, esoxID_LORE_salmoninae)
  
  # ================================================
  # 3. esoxID LORE résolu salmoninae → NAME11
  #    = LORE in nae but not LORE in NAME11
  # ================================================
  
  esoxID_LORE_NAME11 <- unique(genes_NAME11$esoxID[genes_NAME11$class == "LORE"])
  
  resolved_nae <- setdiff(esoxID_LORE_salmoninae, esoxID_LORE_NAME11)
  
  # ================================================
  # 4. esoxID LORE résolu NAME11 → NAME13
  #    = LORE in NAME11 but not LORE in NAME13
  # ================================================
  
  esoxID_LORE_NAME13 <- unique(genes_NAME13$esoxID[genes_NAME13$class == "LORE"])
  
  resolved_nm11 <- setdiff(esoxID_LORE_NAME11, esoxID_LORE_NAME13)
  
  
  # ================================================
  # 6. esoxID LORE NAME13
  #    = LORE in NAME13 
  # ================================================
  
  #esoxID_LORE_NAME13
  
  
  # ================================================
  # Résultat : liste nommée des esoxID par catégorie
  # ================================================
  
  return(list(
    AORE_Salmonidae            = esoxID_AORE_salmonidae,
    LORE_resolution_Salmonidae = esoxID_LORE_salmonidae,
    LORE_resolution_Salmoninae = esoxID_LORE_salmoninae,
    LORE_resolution_NAME11     = esoxID_LORE_NAME11,
    LORE_NAME13                = esoxID_LORE_NAME13
  ))
}

# ======================
# Chargement des fichier
# ======================

win_salmonidae <- readRDS("../Ancestor_Circos/win_salmonidae.rds")
win_salmoninae <- readRDS("../Ancestor_Circos/win_salmoninae.rds")
win_NAME11 <- readRDS("../Ancestor_Circos/win_NAME11.rds")
win_NAME13 <- readRDS("../Ancestor_Circos/win_NAME13.rds")

gene_salmonidae <- readRDS("../Ancestor_Circos/gene_salmonidae.rds")
gene_salmoninae <- readRDS("../Ancestor_Circos/gene_salmoninae.rds")
gene_NAME11 <- readRDS("../Ancestor_Circos/gene_NAME11.rds")
gene_NAME13 <- readRDS("../Ancestor_Circos/gene_NAME13.rds")

# ======================
# PIPELINE
# ======================

barplot <- prepare_resolution_counts(win_salmonidae, win_salmoninae, win_NAME11, win_NAME13)

# give the names of values of x on barplot 

barplot$category <- factor(
  barplot$category,
  levels = c(
    "AORE_Salmonidae",
    "LORE_resolution_Salmonidae",
    "LORE_resolution_Salmoninae",
    "LORE_resolution_NAME11",
    "LORE_NAME13"
  )
)

# barplot

p1 <- ggplot(barplot, aes(x = category, y = total, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey70", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3")) +
  ggtitle("Bar plot of rediploidization resolution") +
  ylab("Number of resolution") +
  theme(legend.position = "none")

print(p1)
ggsave("barplot.pdf", p1, width = 12, height = 8)
ggsave("barplot.svg", p1, width = 12, height = 8)

# recover the esoxID
resolved_IDs <- get_resolved_esoxIDs(win_salmonidae, win_salmoninae, win_NAME11, win_NAME13)

# see the esoxID LORE resolved in salmonidae
head(resolved_IDs$LORE_resolution_Salmonidae)

# see the esoxID LORE resolved in salmoninae  
head(resolved_IDs$LORE_resolution_Salmoninae)

# count by categories
sapply(resolved_IDs, length)

resolved_IDs_circos <- get_resolved_esoxIDs_circos(gene_salmonidae, gene_salmoninae, gene_NAME11, gene_NAME13)


# save for other script
saveRDS(resolved_IDs, "resolved_IDs.rds")
saveRDS(resolved_IDs_circos, "resolved_IDs_circos.rds")
