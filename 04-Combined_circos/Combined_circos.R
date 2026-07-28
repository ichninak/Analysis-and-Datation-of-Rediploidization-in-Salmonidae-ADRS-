# ================
# load libraries
# ===================

library(circlize)
library(dplyr)
library(pdftools)

# =================
# load files 
# =================

karyo <- read.table("karyo_mykiss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
karyo$chr <- as.character(karyo$chr)
coord <- read.delim(
  "coord_mykiss.tsv",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("chr", "start", "end", "gene")
)
write.table(coord, "coord.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

links_all <- rbind(
  read_links("../03-Cas_circos/link_cas_1.tsv",    "AORe"),
  read_links("../03-Cas_circos/link_cas_2.tsv",    "LORe1"),
  read_links("../03-Cas_circos/link_cas_3.tsv",    "LORe2"),
  read_links("../03-Cas_circos/link_cas_4.tsv",    "LORe3"),
  read_links("../03-Cas_circos/link_cas_5.tsv",    "LORe4")
)

win_cas_1 <- readRDS("../03-Cas_circos/win_cas_1.rds")
win_cas_2 <- readRDS("../03-Cas_circos/win_cas_2.rds")
win_cas_3 <- readRDS("../03-Cas_circos/win_cas_3.rds")
win_cas_4 <- readRDS("../03-Cas_circos/win_cas_4.rds")
win_cas_5 <- readRDS("../03-Cas_circos/win_cas_5.rds")

common_cols <- c("chr", "win_start", "win_end", "state", "ratio_AORE", "esoxIDs")

win_all <- rbind(
  mutate(win_cas_1[, common_cols],  cas = "AORe"),
  mutate(win_cas_2[, common_cols],  cas = "LORe1"),
  mutate(win_cas_3[, common_cols],  cas = "LORe2"),
  mutate(win_cas_4[, common_cols],  cas = "LORe3"),
  mutate(win_cas_5[, common_cols],  cas = "LORe4")
)

resolved_IDs_circos <- readRDS("../02-Barplot_analysis/resolved_IDs_circos.rds")
resolved_IDs <- readRDS("../02-Barplot_analysis/resolved_IDs.rds")

pair_files <- c(
  "../01-Ancestor_circos/lore_gar_mykiss1_Salmonidae.tsv",
  "../01-Ancestor_circos/lore_gar_mykiss1_Salmoninae.tsv",
  "../01-Ancestor_circos/lore_gar_mykiss1_NAME11_combine.tsv",
  "../01-Ancestor_circoss/lore_gar_mykiss1_NAME13_combine.tsv"
)

ks_per_win    <- readRDS("../07-Ks_value/ks_per_win.rds")



# ===============
# function
# ==============

read_links <- function(file, cas_label) {
  df <- read.table(file, header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
  df$cas <- cas_label
  return(df)
}

circos_combined <- function(karyo, links_all, win_all, col_cas, output_file) {
  
  links_all$chrA <- as.character(links_all$chrA)
  links_all$chrB <- as.character(links_all$chrB)
  win_all$chr    <- as.character(win_all$chr)
  karyo$chr      <- as.character(karyo$chr)
  
  valid_chr <- karyo$chr
  links_all <- links_all[
    links_all$chrA %in% valid_chr &
      links_all$chrB %in% valid_chr, ]
  
  cas_levels = colnames(col_cas)
  links_all <- links_all[order(match(links_all$cas, cas_levels)), ]
  
  pdf(output_file, width = 14, height = 14)
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.after    = rep(2, nrow(karyo)),
    cell.padding = c(0, 0, 0, 0)
  )
  
  circos.initialize(
    factors = karyo$chr,
    xlim    = cbind(karyo$start, karyo$end)
  )
  
  # Track labels chromosomes
  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.06,
    bg.border    = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl  <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr,
                  facing = "clockwise", niceFacing = TRUE, cex = 1.5)
    }
  )
  
  # Track windows colored by categories
  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.08,
    bg.border    = "grey30",
    bg.col       = "grey95",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- win_all[win_all$chr == chr, ]
      
      if (nrow(sub) == 0) return()
      
      for (i in seq_len(nrow(sub))) {
        cas_i <- as.character(sub$cas[i])
        col_i <- col_cas[cas_i]
        if (is.na(col_i)) col_i <- "grey50"
        
        circos.rect(
          xleft   = sub$win_start[i],
          ybottom = 0,
          xright  = sub$win_end[i],
          ytop    = 1,
          col     = col_i,
          border  = NA
        )
      }
    }
  )
  
  # links
  for (i in seq_len(nrow(links_all))) {
    r  <- links_all[i, ]
    p1 <- (r$startA + r$endA) / 2
    p2 <- (r$startB + r$endB) / 2
    
    circos.link(
      r$chrA, p1,
      r$chrB, p2,
      col    = adjustcolor(col_cas[as.character(r$cas)]),
      border = NA,
      lwd    = 1
    )
  }
  
  
  legend("bottomright",
         legend = names(col_cas),
         fill   = col_cas,
         border = NA,
         bty    = "n",
         cex    = 0.8
  )
  
  circos.clear()
  dev.off()
}

circos_combined_ks <- function(karyo, links_all, win_all, ks_per_win,
                               col_cas, output_file) {
  
  # preparation of types
  links_all$chrA  <- as.character(links_all$chrA)
  links_all$chrB  <- as.character(links_all$chrB)
  win_all$chr     <- as.character(win_all$chr)
  karyo$chr       <- as.character(karyo$chr)
  ks_per_win$chr  <- as.character(ks_per_win$chr)
  
  valid_chr  <- karyo$chr
  links_all  <- links_all[links_all$chrA %in% valid_chr &
                            links_all$chrB %in% valid_chr, ]
  ks_per_win <- ks_per_win[ks_per_win$chr %in% valid_chr, ]
  cas_levels <- names(col_cas)
  links_all  <- links_all[order(match(links_all$cas, cas_levels)), ]
  
  # ylim  
  ks_ylim_by_chr <- do.call(rbind, lapply(valid_chr, function(ch) {
    vals <- ks_per_win$Ks_mean[ks_per_win$chr == ch]
    vals <- vals[!is.na(vals)]
    ymax <- if (length(vals) == 0) 1 else max(vals) * 1.05
    data.frame(chr = ch, ymax = ymax, stringsAsFactors = FALSE)
  }))
  
  # LOESS 
  loess_by_chr <- lapply(valid_chr, function(ch) {
    sub <- ks_per_win[ks_per_win$chr == ch, ]
    sub <- sub[order(sub$mid), ]
    sub <- sub[!is.na(sub$Ks_mean), ]
    if (nrow(sub) < 5) return(NULL)
    tryCatch({
      lo    <- loess(Ks_mean ~ mid, data = sub, span = 0.4)
      x_seq <- seq(min(sub$mid), max(sub$mid), length.out = 200)
      y_seq <- predict(lo, newdata = data.frame(mid = x_seq))
      data.frame(x = x_seq, y = y_seq)
    }, error = function(e) NULL)
  })
  names(loess_by_chr) <- valid_chr
  
  # opening PDF 
  pdf(output_file, width = 14, height = 14)
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.after    = rep(2, nrow(karyo)),
    cell.padding = c(0, 0, 0, 0)
  )
  
  circos.initialize(
    factors = karyo$chr,
    xlim    = cbind(karyo$start, karyo$end)
  )
  
  # Track 1 : labels chromosomes 
  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.06,
    bg.border    = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl  <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr,
                  facing = "clockwise", niceFacing = TRUE, cex = 1.2)
    }
  )
  
  # Track 2 : curve Ks par chromosome
  circos.trackPlotRegion(
    ylim         = c(0, 1),          
    track.height = 0.18,
    bg.border    = "grey60",
    bg.col       = "grey98",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      
      # ymax local for this chromosome
      ymax_local <- ks_ylim_by_chr$ymax[ks_ylim_by_chr$chr == chr]
      if (length(ymax_local) == 0 || is.na(ymax_local)) ymax_local <- 1
      
      # function of rescaling Ks -> [0, 1]
      rsc <- function(v) pmax(0, pmin(v, ymax_local)) / ymax_local
      
      sub <- ks_per_win[ks_per_win$chr == chr, ]
      if (nrow(sub) == 0) return()
      sub <- sub[order(sub$mid), ]
      sub <- sub[!is.na(sub$Ks_mean), ]
      
      # line grey 
      if (nrow(sub) >= 2) {
        circos.lines(sub$mid, rsc(sub$Ks_mean),
                     col = "grey50", lwd = 0.8)
      }
      
      # points colored by categories
      cols <- col_cas[as.character(sub$cas)]
      cols[is.na(cols)] <- "grey50"
      circos.points(sub$mid, rsc(sub$Ks_mean),
                    col = cols, pch = 16, cex = 0.4)
      
      # curve LOESS black
      lo_data <- loess_by_chr[[chr]]
      if (!is.null(lo_data)) {
        lo_data <- lo_data[!is.na(lo_data$y), ]
        if (nrow(lo_data) >= 2) {
          circos.lines(lo_data$x, rsc(lo_data$y),
                       col = adjustcolor("black", alpha.f = 0.5), lwd = 1.2)
        }
      }
      
      # axis Y with value of ks
      at_real     <- round(c(0,
                             ymax_local * 0.25,
                             ymax_local * 0.50,
                             ymax_local * 0.75,
                             ymax_local), 2)
      at_rescaled <- rsc(at_real)
      circos.yaxis(
        side        = "left",
        at          = at_rescaled,
        labels      = at_real,
        tick.length = convert_x(0.5, "mm", chr),
        labels.cex  = 0.35
      )
    }
  )
  
  # Track 3 : windows colored by categories 
  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.06,
    bg.border    = "grey30",
    bg.col       = "grey95",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- win_all[win_all$chr == chr, ]
      if (nrow(sub) == 0) return()
      
      for (i in seq_len(nrow(sub))) {
        cas_i <- as.character(sub$cas[i])
        col_i <- col_cas[cas_i]
        if (is.na(col_i)) col_i <- "grey50"
        circos.rect(
          xleft   = sub$win_start[i],
          ybottom = 0,
          xright  = sub$win_end[i],
          ytop    = 1,
          col     = col_i,
          border  = NA
        )
      }
    }
  )
  
  # links ohnologs
  for (i in seq_len(nrow(links_all))) {
    r  <- links_all[i, ]
    p1 <- (r$startA + r$endA) / 2
    p2 <- (r$startB + r$endB) / 2
    circos.link(
      r$chrA, p1,
      r$chrB, p2,
      col    = adjustcolor(col_cas[as.character(r$cas)], alpha.f = 0.5),
      border = NA,
      lwd    = 1
    )
  }
  
  # legends
  legend("bottomright",
         legend = names(col_cas),
         fill   = col_cas,
         border = NA, bty = "n", cex = 0.9)
  
  circos.clear()
  dev.off()
  message("Circos saved → ", output_file)
}

# same for make svg
circos_combined_ks_svg <- function(karyo, links_all, win_all, ks_per_win,
                               col_cas, output_file) {
  
  # preparation of types 
  links_all$chrA  <- as.character(links_all$chrA)
  links_all$chrB  <- as.character(links_all$chrB)
  win_all$chr     <- as.character(win_all$chr)
  karyo$chr       <- as.character(karyo$chr)
  ks_per_win$chr  <- as.character(ks_per_win$chr)
  
  valid_chr  <- karyo$chr
  links_all  <- links_all[links_all$chrA %in% valid_chr &
                            links_all$chrB %in% valid_chr, ]
  ks_per_win <- ks_per_win[ks_per_win$chr %in% valid_chr, ]
  cas_levels <- names(col_cas)
  links_all  <- links_all[order(match(links_all$cas, cas_levels)), ]
  
  # ylim 
  ks_ylim_by_chr <- do.call(rbind, lapply(valid_chr, function(ch) {
    vals <- ks_per_win$Ks_mean[ks_per_win$chr == ch]
    vals <- vals[!is.na(vals)]
    ymax <- if (length(vals) == 0) 1 else max(vals) * 1.05
    data.frame(chr = ch, ymax = ymax, stringsAsFactors = FALSE)
  }))
  
  # LOESS 
  loess_by_chr <- lapply(valid_chr, function(ch) {
    sub <- ks_per_win[ks_per_win$chr == ch, ]
    sub <- sub[order(sub$mid), ]
    sub <- sub[!is.na(sub$Ks_mean), ]
    if (nrow(sub) < 5) return(NULL)
    tryCatch({
      lo    <- loess(Ks_mean ~ mid, data = sub, span = 0.4)
      x_seq <- seq(min(sub$mid), max(sub$mid), length.out = 200)
      y_seq <- predict(lo, newdata = data.frame(mid = x_seq))
      data.frame(x = x_seq, y = y_seq)
    }, error = function(e) NULL)
  })
  names(loess_by_chr) <- valid_chr
  
  # opening PDF
  svg(output_file, width = 14, height = 14)
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.after    = rep(2, nrow(karyo)),
    cell.padding = c(0, 0, 0, 0)
  )
  
  circos.initialize(
    factors = karyo$chr,
    xlim    = cbind(karyo$start, karyo$end)
  )
  
 
  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.06,
    bg.border    = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl  <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr,
                  facing = "clockwise", niceFacing = TRUE, cex = 1.2)
    }
  )
  
  
  circos.trackPlotRegion(
    ylim         = c(0, 1),          
    track.height = 0.18,
    bg.border    = "grey60",
    bg.col       = "grey98",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      
     
      ymax_local <- ks_ylim_by_chr$ymax[ks_ylim_by_chr$chr == chr]
      if (length(ymax_local) == 0 || is.na(ymax_local)) ymax_local <- 1
      
      
      rsc <- function(v) pmax(0, pmin(v, ymax_local)) / ymax_local
      
      sub <- ks_per_win[ks_per_win$chr == chr, ]
      if (nrow(sub) == 0) return()
      sub <- sub[order(sub$mid), ]
      sub <- sub[!is.na(sub$Ks_mean), ]
      
     
      if (nrow(sub) >= 2) {
        circos.lines(sub$mid, rsc(sub$Ks_mean),
                     col = "grey50", lwd = 0.8)
      }
      
      
      cols <- col_cas[as.character(sub$cas)]
      cols[is.na(cols)] <- "grey50"
      circos.points(sub$mid, rsc(sub$Ks_mean),
                    col = cols, pch = 16, cex = 0.4)
      
   
      lo_data <- loess_by_chr[[chr]]
      if (!is.null(lo_data)) {
        lo_data <- lo_data[!is.na(lo_data$y), ]
        if (nrow(lo_data) >= 2) {
          circos.lines(lo_data$x, rsc(lo_data$y),
                       col = adjustcolor("black", alpha.f = 0.5), lwd = 1.2)
        }
      }
      
  
      at_real     <- round(c(0,
                             ymax_local * 0.25,
                             ymax_local * 0.50,
                             ymax_local * 0.75,
                             ymax_local), 2)
      at_rescaled <- rsc(at_real)
      circos.yaxis(
        side        = "left",
        at          = at_rescaled,
        labels      = at_real,
        tick.length = convert_x(0.5, "mm", chr),
        labels.cex  = 0.35
      )
    }
  )
  

  circos.trackPlotRegion(
    ylim         = c(0, 1),
    track.height = 0.06,
    bg.border    = "grey30",
    bg.col       = "grey95",
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- win_all[win_all$chr == chr, ]
      if (nrow(sub) == 0) return()
      
      for (i in seq_len(nrow(sub))) {
        cas_i <- as.character(sub$cas[i])
        col_i <- col_cas[cas_i]
        if (is.na(col_i)) col_i <- "grey50"
        circos.rect(
          xleft   = sub$win_start[i],
          ybottom = 0,
          xright  = sub$win_end[i],
          ytop    = 1,
          col     = col_i,
          border  = NA
        )
      }
    }
  )
  

  for (i in seq_len(nrow(links_all))) {
    r  <- links_all[i, ]
    p1 <- (r$startA + r$endA) / 2
    p2 <- (r$startB + r$endB) / 2
    circos.link(
      r$chrA, p1,
      r$chrB, p2,
      col    = adjustcolor(col_cas[as.character(r$cas)], alpha.f = 0.5),
      border = NA,
      lwd    = 1
    )
  }
  

  legend("bottomright",
         legend = names(col_cas),
         fill   = col_cas,
         border = NA, bty = "n", cex = 0.9)
  
  circos.clear()
  dev.off()
  message("Circos saved → ", output_file)
}

prepare_ohno_links <- function(pair_files, coord_file, resolved_ids_circos) {
  
  coord <- read.table(coord_file, header = TRUE, stringsAsFactors = FALSE)
  coord$gene <- as.character(coord$gene)
  
  res <- resolved_ids_circos
  
  assign_cas <- function(esoxid) {
    if (esoxid %in% res$AORE_Salmonidae)            return("AORe")
    if (esoxid %in% res$LORE_resolution_Salmonidae) return("LORe1")
    if (esoxid %in% res$LORE_resolution_Salmoninae) return("LORe2")
    if (esoxid %in% res$LORE_resolution_NAME11)     return("LORe3")
    if (esoxid %in% res$LORE_NAME13)                return("LORe4")
    return(NA)
  }
  
  # read and combine all files of pairs 
  all_pairs <- do.call(rbind, lapply(pair_files, function(f) {
    read.table(f, col.names = c("esoxID", "gene1", "gene2", "flag"),
               stringsAsFactors = FALSE)
  }))
  
  # duplicaiton
  all_pairs <- unique(all_pairs)
  all_pairs$gene1 <- as.character(all_pairs$gene1)
  all_pairs$gene2 <- as.character(all_pairs$gene2)
  
  # Assign category
  all_pairs$cas <- sapply(all_pairs$esoxID, assign_cas)
  all_pairs <- all_pairs[!is.na(all_pairs$cas), ]
  
  # Join gene1
  all_pairs <- merge(all_pairs, coord, by.x = "gene1", by.y = "gene", all.x = FALSE)
  names(all_pairs)[names(all_pairs) == "chr"]   <- "chrA"
  names(all_pairs)[names(all_pairs) == "start"] <- "startA"
  names(all_pairs)[names(all_pairs) == "end"]   <- "endA"
  
  # Join gene2
  all_pairs <- merge(all_pairs, coord, by.x = "gene2", by.y = "gene", all.x = FALSE)
  names(all_pairs)[names(all_pairs) == "chr"]   <- "chrB"
  names(all_pairs)[names(all_pairs) == "start"] <- "startB"
  names(all_pairs)[names(all_pairs) == "end"]   <- "endB"
  
  links <- all_pairs[, c("esoxID", "cas", "chrA", "startA", "endA", "chrB", "startB", "endB")]
  
  return(links)
}

# ===============
# PIPELINE
# ==============


links_ohno <- prepare_ohno_links(pair_files, "coord.tsv", resolved_IDs_circos)
links_ohno_bin <- prepare_ohno_links(pair_files, "coord.tsv", resolved_IDs)

links_ohno <- links_ohno %>%
  mutate(cas = recode(cas,
                      "cas1" = "AORe",
                      "cas2" = "LORe1",
                      "cas3" = "LORe2",
                      "cas4" = "LORe3",
                      "cas5" = "LORe4"))

links_ohno_bin <- links_ohno_bin %>%
  mutate(cas = recode(cas,
                      "cas1" = "AORe",
                      "cas2" = "LORe1",
                      "cas3" = "LORe2",
                      "cas4" = "LORe3",
                      "cas5" = "LORe4"))


# Colors by category 


col_cas <- c(
  AORe  = "grey70",   
  LORe1 = "#e41a1c",   
  LORe2 = "#377eb8",
  LORe3 = "#4daf4a",
  LORe4 = "#984ea3"
)

win_all_plot <- win_all[, c("chr", "win_start", "win_end", "ratio_AORE", "cas")]



# 
# circos_combined(
#   karyo       = karyo,
#   links_all   = links_all,
#   win_all     = win_all_plot,
#   col_cas     = col_cas,
#   output_file = "circos_combined4.pdf"
# )
circos_combined(
  karyo       = karyo,
  links_all   = links_ohno_bin,
  win_all     = win_all,
  col_cas     = col_cas,
  output_file = "circos_combined.pdf"
)

circos_combined_ks(
  karyo       = karyo,
  links_all   = links_ohno_bin,
  win_all     = win_all,
  ks_per_win  = ks_per_win,     
  col_cas     = col_cas,
  output_file = "circos_with_ks.pdf"
)

circos_combined_ks_svg(
  karyo       = karyo,
  links_all   = links_ohno_bin,
  win_all     = win_all,
  ks_per_win  = ks_per_win,     
  col_cas     = col_cas,
  output_file = "circos_with_ks.svg"
)



