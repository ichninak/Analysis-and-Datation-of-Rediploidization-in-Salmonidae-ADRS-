# ======================
# load libraries
# ======================
library(dplyr)
library(circlize)

# ================
# function
# ================

prepare_data_circos <- function(coord, pair_file, name_output, filter_ids) {
  pairs = read.table(pair_file, col.names = c("esoxID", "geneA", "geneB", "flag"), stringsAsFactors = FALSE)
  
  # Filter esoxID
  pairs <- pairs %>%
    filter(esoxID %in% filter_ids)
  
  pairs <- pairs %>%
    dplyr::filter(esoxID %in% filter_ids)
  
  mA <- pairs %>%
    dplyr::inner_join(coord, by = c("geneA" = "gene")) %>%
    dplyr::rename(chrA = chr, startA = start, endA = end)
  
  mAB <- mA %>%
    dplyr::inner_join(coord, by = c("geneB" = "gene")) %>%
    dplyr::rename(chrB = chr, startB = start, endB = end)
  
  links <- mAB %>%
    dplyr::select(chrA, startA, endA, chrB, startB, endB, flag)
  
  write.table(links, name_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("OK ->", name_output, ":", nrow(links), "liens\n")
}

circos_precise_resolved <- function(karyo, links_file, name_output,
                                    windows_ratio, target = c("LORE", "AORE")) {
  library(circlize)
  library(GenomicRanges)
  target <- match.arg(target)
  
  links <- read.table(links_file, header = TRUE, sep = "\t",
                      fill = TRUE, comment.char = "")
  links$chrA <- as.character(links$chrA)
  links$chrB <- as.character(links$chrB)
  windows_ratio$chr <- as.character(windows_ratio$chr)
  karyo$chr <- as.character(karyo$chr)
  
  links <- links[links$chrA %in% karyo$chr & links$chrB %in% karyo$chr, ]
  links$class <- ifelse(links$flag == 1, "LORE", "AORE")
  col_map <- c(AORE = "#2C7FB8", LORE = "#D7301F")
  
  # state of windows
  windows_ratio$state <- "neutral"
  windows_ratio$state[windows_ratio$ratio_AORE > 0.6] <- "AORE"
  windows_ratio$state[windows_ratio$ratio_AORE < 0.4] <- "LORE"
  
  state_col <- c(AORE = "#2C7FB8", LORE = "#D7301F", neutral = "grey85")
  
  # ---- filterage : links `target` with 2 ends in windows `target` ----
  win_target <- windows_ratio[windows_ratio$state == target, ]
  
  if (nrow(win_target) > 0 && any(links$class == target)) {
    win_gr <- GRanges(win_target$chr,
                      IRanges(win_target$win_start, win_target$win_end))
    gA <- GRanges(links$chrA,
                  IRanges(pmin(links$startA, links$endA),
                          pmax(links$startA, links$endA)))
    gB <- GRanges(links$chrB,
                  IRanges(pmin(links$startB, links$endB),
                          pmax(links$startB, links$endB)))
    keep <- links$class == target &
      overlapsAny(gA, win_gr) &
      overlapsAny(gB, win_gr)
    links <- links[keep, ]
  } else {
    links <- links[0, ]
  }
  # -------------------------------------------------------------------------
  # print circos
  
  svg(name_output, width = 14, height = 14)
  circos.clear()
  circos.par(start.degree = 90,
             gap.after = rep(2, nrow(karyo)),
             cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = karyo$chr, xlim = cbind(karyo$start, karyo$end))
  
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr, facing = "clockwise",
                  niceFacing = TRUE, cex = 1.5)
    }
  )
  
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.08, bg.border = "grey20",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- windows_ratio[windows_ratio$chr == chr, ]
      if (nrow(sub) == 0) return()
      for (i in seq_len(nrow(sub))) {
        circos.rect(sub$win_start[i], 0, sub$win_end[i], 1,
                    col = state_col[sub$state[i]], border = NA)
      }
    }
  )
  
  for (i in seq_len(nrow(links))) {
    p1 <- (links$startA[i] + links$endA[i]) / 2
    p2 <- (links$startB[i] + links$endB[i]) / 2
    circos.link(links$chrA[i], p1, links$chrB[i], p2,
                col = adjustcolor(col_map[links$class[i]]),
                border = NA, lwd = 1)
  }
  
  dev.off()
  circos.clear()
}

circos_precise_resolved_pdf <- function(karyo, links_file, name_output,
                                    windows_ratio, target = c("LORE", "AORE")) {
  library(circlize)
  library(GenomicRanges)
  target <- match.arg(target)
  
  links <- read.table(links_file, header = TRUE, sep = "\t",
                      fill = TRUE, comment.char = "")
  links$chrA <- as.character(links$chrA)
  links$chrB <- as.character(links$chrB)
  windows_ratio$chr <- as.character(windows_ratio$chr)
  karyo$chr <- as.character(karyo$chr)
  
  links <- links[links$chrA %in% karyo$chr & links$chrB %in% karyo$chr, ]
  links$class <- ifelse(links$flag == 1, "LORE", "AORE")
  col_map <- c(AORE = "#2C7FB8", LORE = "#D7301F")
  
  # state of windows
  windows_ratio$state <- "neutral"
  windows_ratio$state[windows_ratio$ratio_AORE > 0.6] <- "AORE"
  windows_ratio$state[windows_ratio$ratio_AORE < 0.4] <- "LORE"
  
  state_col <- c(AORE = "#2C7FB8", LORE = "#D7301F", neutral = "grey85")
  
  # ---- filterage : links `target` with 2 ends in windows `target` ----
  win_target <- windows_ratio[windows_ratio$state == target, ]
  
  if (nrow(win_target) > 0 && any(links$class == target)) {
    win_gr <- GRanges(win_target$chr,
                      IRanges(win_target$win_start, win_target$win_end))
    gA <- GRanges(links$chrA,
                  IRanges(pmin(links$startA, links$endA),
                          pmax(links$startA, links$endA)))
    gB <- GRanges(links$chrB,
                  IRanges(pmin(links$startB, links$endB),
                          pmax(links$startB, links$endB)))
    keep <- links$class == target &
      overlapsAny(gA, win_gr) &
      overlapsAny(gB, win_gr)
    links <- links[keep, ]
  } else {
    links <- links[0, ]
  }
  # -------------------------------------------------------------------------
  
  pdf(name_output, width = 14, height = 14)
  circos.clear()
  circos.par(start.degree = 90,
             gap.after = rep(2, nrow(karyo)),
             cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = karyo$chr, xlim = cbind(karyo$start, karyo$end))
  
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr, facing = "clockwise",
                  niceFacing = TRUE, cex = 1.5)
    }
  )
  
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.08, bg.border = "grey20",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- windows_ratio[windows_ratio$chr == chr, ]
      if (nrow(sub) == 0) return()
      for (i in seq_len(nrow(sub))) {
        circos.rect(sub$win_start[i], 0, sub$win_end[i], 1,
                    col = state_col[sub$state[i]], border = NA)
      }
    }
  )
  
  for (i in seq_len(nrow(links))) {
    p1 <- (links$startA[i] + links$endA[i]) / 2
    p2 <- (links$startB[i] + links$endB[i]) / 2
    circos.link(links$chrA[i], p1, links$chrB[i], p2,
                col = adjustcolor(col_map[links$class[i]]),
                border = NA, lwd = 1)
  }
  
  dev.off()
  circos.clear()
}

# ======================
# load files 
# ======================

karyo <- read.table("../Ancestor_Circos/karyo_mykiss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
karyo$chr <- as.character(karyo$chr)
coord <- read.delim(
  "../Ancestor_Circos/coord_mykiss.tsv",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("chr", "start", "end", "gene")
)

win_salmonidae <- readRDS("../Ancestor_Circos/win_salmonidae.rds")
win_salmoninae <- readRDS("../Ancestor_Circos/win_salmoninae.rds")
win_NAME11 <- readRDS("../Ancestor_Circos/win_NAME11.rds")
win_NAME13 <- readRDS("../Ancestor_Circos/win_NAME13.rds")

resolved_IDs <- readRDS("../Barplot_analysis/resolved_IDs.rds")
resolved_IDs_circos <- readRDS("../Barplot_analysis/resolved_IDs_circos.rds")

# ======================
# PIPELINE
# ======================


# preparation of dataframe (links) (link for build circos)

prepare_data_circos(coord, "../Ancestor_Circos/lore_gar_mykiss1_Salmonidae.tsv", "link_cas_1.tsv", resolved_IDs_circos$AORE_Salmonidae)
prepare_data_circos(coord, "../Ancestor_Circos/lore_gar_mykiss1_Salmonidae.tsv", "link_cas_2.tsv", resolved_IDs_circos$LORE_resolution_Salmonidae)
prepare_data_circos(coord, "../Ancestor_Circos/lore_gar_mykiss1_Salmoninae.tsv", "link_cas_3.tsv", resolved_IDs_circos$LORE_resolution_Salmoninae)
prepare_data_circos(coord, "../Ancestor_Circos/lore_gar_mykiss1_NAME11_combine.tsv", "link_cas_4.tsv", resolved_IDs_circos$LORE_resolution_NAME11)
prepare_data_circos(coord, "../Ancestor_Circos/lore_gar_mykiss1_NAME13_combine.tsv", "link_cas_5.tsv", resolved_IDs_circos$LORE_NAME13)


# preparation of dataframe

# aore salmonidae (AORe)
win_cas_1 <- win_salmonidae %>%
  filter(state == "AORE",
         sapply(esoxIDs, function(ids) any(ids %in% resolved_IDs$AORE_Salmonidae)))


# lore salmonidae (LORe1)
dae_nae <- merge(win_salmonidae, win_salmoninae,
                 by = c("chr", "win_start", "win_end"),
                 suffixes = c("_dae", "_nae"))

win_cas_2 <- dae_nae %>%
  filter(
    state_dae == "LORE",
    state_nae != "LORE",
    sapply(esoxIDs_dae, function(ids) any(ids %in% resolved_IDs$LORE_resolution_Salmonidae))
  ) %>%
  # Rename for circos_precise work
  dplyr::rename(
    state = state_dae,
    esoxIDs = esoxIDs_dae,
    ratio_AORE = ratio_AORE_dae,
    total = total_dae
  ) %>%
  select(chr, win_start, win_end, state, ratio_AORE, esoxIDs, total)


# lore salmoninae (LORe2)

nae_nm11 <- merge(win_salmoninae, win_NAME11,
                  by = c("chr","win_start","win_end"),
                  suffixes = c("_nae","_nm11"))

win_cas_3 <- nae_nm11 %>%
  filter(
    state_nae == "LORE",
    state_nm11 != "LORE",
    sapply(esoxIDs_nae, function(ids) any(ids %in% resolved_IDs$LORE_resolution_Salmoninae))
  ) %>%
  # Rename for circos_precise work
  dplyr::rename(state = state_nae,
                esoxIDs = esoxIDs_nae,
                ratio_AORE = ratio_AORE_nae,
                total = total_nae
  ) %>%
  select(chr, win_start, win_end, state, ratio_AORE, esoxIDs, total)

# lore NAME11(LORe3)

nm11_nm13 <- merge(win_NAME11, win_NAME13,
                   by = c("chr","win_start","win_end"),
                   suffixes = c("_nm11","_nm13"))

win_cas_4 <- nm11_nm13 %>%
  filter(
    state_nm11 == "LORE",
    state_nm13 != "LORE",
    sapply(esoxIDs_nm11, function(ids) any(ids %in% resolved_IDs$LORE_resolution_NAME11))
  ) %>%
  # Rename for circos_precise work
  dplyr::rename(state = state_nm11,
                esoxIDs = esoxIDs_nm11,
                ratio_AORE = ratio_AORE_nm11,
                total = total_nm11
  ) %>%
  select(chr, win_start, win_end, state, ratio_AORE, esoxIDs, total)


# lore NAME13 (LORe4)

win_cas_5 <- win_NAME13 %>%
  filter(state == "LORE",
         sapply(esoxIDs, function(ids) any(ids %in% resolved_IDs$LORE_NAME13)))



# plot circos

## AORE salmonidae 
circos_precise_resolved(karyo, "link_cas_1.tsv", "circos_AORe.svg", win_cas_1, target = "AORE")
circos_precise_resolved_pdf(karyo, "link_cas_1.tsv", "circos_AORe.pdf", win_cas_1, target = "AORE")

## LORE resolved salmonidae
circos_precise_resolved(karyo, "link_cas_2.tsv", "circos_LORe1.svg", win_cas_2, target = "LORE")
circos_precise_resolved_pdf(karyo, "link_cas_2.tsv", "circos_LORe1.pdf", win_cas_2, target = "LORE")

## LORE resolved salmoninae
circos_precise_resolved(karyo, "link_cas_3.tsv", "circos_LORe2.svg", win_cas_3)
circos_precise_resolved_pdf(karyo, "link_cas_3.tsv", "circos_LORe2.pdf", win_cas_3)

## LORE resolved NAME11
circos_precise_resolved(karyo, "link_cas_4.tsv", "circos_LORe3.svg", win_cas_4)
circos_precise_resolved_pdf(karyo, "link_cas_4.tsv", "circos_LORe3.pdf", win_cas_4)

## LORE NAME13
circos_precise_resolved(karyo, "link_cas_5.tsv", "circos_LORe4.svg", win_cas_5)
circos_precise_resolved_pdf(karyo, "link_cas_5.tsv", "circos_LORe4.pdf", win_cas_5)



## export file 

saveRDS(win_cas_1, "win_cas_1.rds") # AORE
saveRDS(win_cas_2, "win_cas_2.rds") # LORe1
saveRDS(win_cas_3, "win_cas_3.rds") # LORe2
saveRDS(win_cas_4, "win_cas_4.rds") # LORe3
saveRDS(win_cas_5, "win_cas_5.rds") # LORe4




