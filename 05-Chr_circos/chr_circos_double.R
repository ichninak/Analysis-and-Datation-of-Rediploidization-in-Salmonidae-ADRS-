# ====================================================
# Circos for only double chromosomes or triple chromosomes
# ====================================================

# ===============
# load libraries
# ===============

library(dplyr)
library(tidyr)
library(circlize)

# ===============
# function
# ===============

plot_circos_links <- function(
    links_data,
    karyo_data,
    col_cas,
    win_all,
    title        = NULL,
    alpha        = 1,
    lwd          = 1.6,
    bdr          = NA,
    track.height = 0.02,
    canvas.lim   = 1.3,
    legend.x     = 1.1,
    legend.y     = 0.5,
    show.legend  = TRUE,
    output.file  = NULL,
    width        = 14,
    height       = 14
) {
  
  if (!is.null(output.file)) {
    ext <- tools::file_ext(output.file)
    if (ext == "svg") svg(output.file, width = width, height = height)
    else if (ext == "pdf") pdf(output.file, width = width, height = height)
    else if (ext == "png") png(output.file, width = width, height = height,
                               units = "in", res = 300)
    else stop("Format non supporté : svg, pdf ou png")
  }
  
  
  cas_levels <- names(col_cas)
  links_data$cas <- as.character(links_data$cas)
  links_data     <- links_data[order(match(links_data$cas, cas_levels)), ]
  
  par(mar = c(4, 4, 4, 8))
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after    = 2,
    cell.padding = c(0, 0, 0, 0),
    #track.margin = c(0.005, 0.005),
    canvas.xlim  = c(-canvas.lim, canvas.lim),
    canvas.ylim  = c(-canvas.lim, canvas.lim)
  )
  
  circos.initializeWithIdeogram(karyo_data, plotType = NULL)
  
  
  # labels
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr, facing = "clockwise", niceFacing = TRUE, cex = 1.5)
    }
  )
  
  # # --- Piste chromosome rectangle fin ---
  # circos.genomicTrack(
  #   karyo_data,
  #   ylim         = c(0, 1),
  #   track.height = track.height,
  #   bg.border    = "grey40",
  #   bg.col       = "grey80",
  #   panel.fun = function(region, value, ...) {
  #     circos.genomicRect(region, value,
  #                        ytop    = 1,
  #                        ybottom = 0,
  #                        col     = "grey80",
  #                        border  = NA, ...)
  #   }
  # )
  
  # Track fenêtres colorées par cas
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
  
  
  
  # --- Liens ---
  for (i in seq_len(nrow(links_data))) {
    cl <- links_data$cas[i]
    p1 <- (links_data$startA[i] + links_data$endA[i]) / 2
    p2 <- (links_data$startB[i] + links_data$endB[i]) / 2
    
    
    
    circos.link(
      links_data$chrA[i], p1,
      links_data$chrB[i], p2,
      col    = adjustcolor(col_cas[as.character(cl)], alpha.f = alpha),
      border = bdr,
      lwd    = lwd
    )
  }
  
  if (!is.null(title)) title(title, cex.main = 1.2)
  
  if (show.legend) {
    cas_presents <- intersect(names(col_cas), unique(links_data$cas))
    legend(
      x      = legend.x,
      y      = legend.y,
      legend = cas_presents,
      fill   = col_cas[cas_presents],
      border = NA,
      bty    = "n",
      cex    = 1,
      title  = "type of rediplodization",
      xpd    = TRUE
    )
  }
  
  circos.clear()
  
  if (!is.null(output.file)) dev.off()
}

plot_circos_links_ks <- function(
    links_data,
    karyo_data,
    col_cas,
    win_all,
    ks_per_win,                       
    title        = NULL,
    alpha        = 1,
    lwd          = 1.6,
    bdr          = NA,
    track.height = 0.02,
    canvas.lim   = 1.3,
    legend.x     = 1.1,
    legend.y     = 0.5,
    show.legend  = TRUE,
    output.file  = NULL,
    width        = 14,
    height       = 14
) {
  
  if (!is.null(output.file)) {
    ext <- tools::file_ext(output.file)
    if (ext == "svg") svg(output.file, width = width, height = height)
    else if (ext == "pdf") pdf(output.file, width = width, height = height)
    else if (ext == "png") png(output.file, width = width, height = height,
                               units = "in", res = 300)
    else stop("Format non supporté : svg, pdf ou png")
  }
  
  cas_levels <- names(col_cas)
  links_data$cas <- as.character(links_data$cas)
  links_data     <- links_data[order(match(links_data$cas, cas_levels)), ]
  
  # ── préparation Ks ──────────────────────────────────────────────────────
  ks_per_win$chr <- as.character(ks_per_win$chr)
  valid_chr      <- as.character(karyo_data[[1]])
  ks_per_win     <- ks_per_win[ks_per_win$chr %in% valid_chr, ]
  
  ks_ylim_by_chr <- do.call(rbind, lapply(valid_chr, function(ch) {
    vals <- ks_per_win$Ks_mean[ks_per_win$chr == ch]
    vals <- vals[!is.na(vals)]
    ymax <- if (length(vals) == 0) 1 else max(vals) * 1.05
    data.frame(chr = ch, ymax = ymax, stringsAsFactors = FALSE)
  }))
  
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
  
  par(mar = c(4, 4, 4, 8))
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after    = 2,
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim  = c(-canvas.lim, canvas.lim),
    canvas.ylim  = c(-canvas.lim, canvas.lim)
  )
  
  circos.initializeWithIdeogram(karyo_data, plotType = NULL)
  
  # ── Track 1 : labels chromosomes ────────────────────────────────────────
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      xl  <- get.cell.meta.data("xlim")
      circos.text(mean(xl), 0.5, chr,
                  facing = "clockwise", niceFacing = TRUE, cex = 1.5)
    }
  )
  
  # ── Track 2 : courbe Ks par chromosome (Y libre) ────────────────────────
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
  
  # ── Track 3 : fenêtres colorées par cas ─────────────────────────────────
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
  
  # ── Liens ──────────────────────────────────────────────────────────────
  for (i in seq_len(nrow(links_data))) {
    cl <- links_data$cas[i]
    p1 <- (links_data$startA[i] + links_data$endA[i]) / 2
    p2 <- (links_data$startB[i] + links_data$endB[i]) / 2
    
    circos.link(
      links_data$chrA[i], p1,
      links_data$chrB[i], p2,
      col    = adjustcolor(col_cas[as.character(cl)], alpha.f = alpha),
      border = bdr,
      lwd    = lwd
    )
  }
  
  if (!is.null(title)) title(title, cex.main = 1.2)
  
  if (show.legend) {
    cas_presents <- intersect(names(col_cas), unique(links_data$cas))
    legend(
      x      = legend.x,
      y      = legend.y,
      legend = cas_presents,
      fill   = col_cas[cas_presents],
      border = NA,
      bty    = "n",
      cex    = 1,
      title  = "type of rediplodization",
      xpd    = TRUE
    )
  }
  
  circos.clear()
  
  if (!is.null(output.file)) dev.off()
}

plot_circos <- function(links_data, karyo_data, col_liens, title = NULL) {
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after    = 2,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.005, 0.005)
  )
  
  circos.initializeWithIdeogram(karyo_data, plotType = c("axis", "labels"))
  
  for (i in seq_len(nrow(links_data))) {
    cl <- links_data$cas[i]
    p1 <- (links_data$startA[i] + links_data$endA[i]) / 2
    p2 <- (links_data$startB[i] + links_data$endB[i]) / 2
    circos.link(
      links_data$chrA[i], p1,
      links_data$chrB[i], p2,
      col    = adjustcolor(col_liens[cl], alpha.f = 0.5),
      border = 2,
      lwd    = 0.8
    )
  }
  
  if (!is.null(title)) title(title, cex.main = 1.2)
  circos.clear()
}

prepare_ohno_links <- function(pair_files, coord_file, resolved_ids_circos) {
  
  coord <- read.table(coord_file, header = TRUE, stringsAsFactors = FALSE)
  coord$gene <- as.character(coord$gene)
  
  res <- resolved_ids_circos
  
  assign_cas <- function(esoxid) {
    if (esoxid %in% res$AORE_Salmonidae)            return("cas1")
    if (esoxid %in% res$LORE_resolution_Salmonidae) return("cas2")
    if (esoxid %in% res$LORE_resolution_Salmoninae) return("cas3")
    if (esoxid %in% res$LORE_resolution_NAME11)     return("cas4")
    if (esoxid %in% res$LORE_NAME13)                return("cas5")
    return(NA)
  }
  
  # Lire et combiner tous les fichiers de paires
  all_pairs <- do.call(rbind, lapply(pair_files, function(f) {
    read.table(f, col.names = c("esoxID", "gene1", "gene2", "flag"),
               stringsAsFactors = FALSE)
  }))
  
  # Dédupliquer
  all_pairs <- unique(all_pairs)
  all_pairs$gene1 <- as.character(all_pairs$gene1)
  all_pairs$gene2 <- as.character(all_pairs$gene2)
  
  # Assigner le cas
  all_pairs$cas <- sapply(all_pairs$esoxID, assign_cas)
  all_pairs <- all_pairs[!is.na(all_pairs$cas), ]
  
  # Jointure gene1
  all_pairs <- merge(all_pairs, coord, by.x = "gene1", by.y = "gene", all.x = FALSE)
  names(all_pairs)[names(all_pairs) == "chr"]   <- "chrA"
  names(all_pairs)[names(all_pairs) == "start"] <- "startA"
  names(all_pairs)[names(all_pairs) == "end"]   <- "endA"
  
  # Jointure gene2
  all_pairs <- merge(all_pairs, coord, by.x = "gene2", by.y = "gene", all.x = FALSE)
  names(all_pairs)[names(all_pairs) == "chr"]   <- "chrB"
  names(all_pairs)[names(all_pairs) == "start"] <- "startB"
  names(all_pairs)[names(all_pairs) == "end"]   <- "endB"
  
  links <- all_pairs[, c("esoxID", "cas", "chrA", "startA", "endA", "chrB", "startB", "endB")]
  
  return(links)
}


make_track <- function(cas_name) {
  rbind(
    links_ohno %>% filter(cas == cas_name) %>%
      select(chr = chrA, start = startA, end = endA),
    links_ohno %>% filter(cas == cas_name) %>%
      select(chr = chrB, start = startB, end = endB)
  ) %>%
    distinct() %>%
    mutate(value = 1)
}

# ===============
# load file
# ===============

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
write.table(coord, "coord.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pair_files <- c(
  "../Ancestor_Circos/lore_gar_mykiss1_Salmonidae.tsv",
  "../Ancestor_Circos/lore_gar_mykiss1_Salmoninae.tsv",
  "../Ancestor_Circos/lore_gar_mykiss1_NAME11_combine.tsv",
  "../Ancestor_Circos/lore_gar_mykiss1_NAME13_combine.tsv"
)

resolved_IDs_circos <- readRDS("../Barplot_analysis/resolved_IDs_circos.rds")
resolved_IDs <- readRDS("../Barplot_analysis/resolved_IDs.rds")

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

track_cas1 <- make_track("AORe")
track_cas2 <- make_track("LORe1")
track_cas3 <- make_track("LORe2")
track_cas4 <- make_track("LORe3")
track_cas5 <- make_track("LORe4")


col_test <- c(
  AORe = "grey70",
  LORe1 = "#e41a1c",
  LORe2 = "#377eb8",
  LORe3 = "#4daf4a",
  LORe4 = "#984ea3"
)


win_cas_1 <- readRDS("../cas_Circos/win_cas_1.rds")
win_cas_2 <- readRDS("../cas_Circos/win_cas_2.rds")
win_cas_3 <- readRDS("../cas_Circos/win_cas_3.rds")
win_cas_4 <- readRDS("../cas_Circos/win_cas_4.rds")
win_cas_5 <- readRDS("../cas_Circos/win_cas_5.rds")

common_cols <- c("chr", "win_start", "win_end", "state", "ratio_AORE", "esoxIDs")

win_all <- rbind(
  mutate(win_cas_1[, common_cols],  cas = "AORe"),
  mutate(win_cas_2[, common_cols],  cas = "LORe1"),
  mutate(win_cas_3[, common_cols],  cas = "LORe2"),
  mutate(win_cas_4[, common_cols],  cas = "LORe3"),
  mutate(win_cas_5[, common_cols],  cas = "LORe4")
)

ks_per_win    <- readRDS("../Ks_value/ks_per_win.rds")

# =============
# PIPELINE
# =============

# Define interest chromosome
# you need to change if you want to have other double
chr_sel <- c("18", "7")

# Filter the links
links_filtered <- links_ohno_bin %>%
  filter(
    (chrA == "18" & chrB == "7") |
      (chrA == "7" & chrB == "18")   
  )


# Filter karyotype
karyo_filtered <- karyo %>% filter(chr %in% chr_sel)  

# Filter of tracks
filter_track <- function(track) track %>% filter(chr %in% chr_sel)

track_cas1_f <- filter_track(track_cas1)
track_cas2_f <- filter_track(track_cas2)
track_cas3_f <- filter_track(track_cas3)
track_cas4_f <- filter_track(track_cas4)
track_cas5_f <- filter_track(track_cas5)



# all categories
plot_circos_links(
  links_data  = links_filtered,
  karyo_data  = karyo_filtered,
  col_cas     = col_test,
  win_all = win_all,
  title       = "",
  output.file = "circos_chr18_17_all.svg"
)

# all categories with curve of Ks
plot_circos_links_ks(
  links_data  = links_filtered,
  karyo_data  = karyo_filtered,
  col_cas     = col_test,
  win_all = win_all,
  ks_per_win = ks_per_win,
  title       = "",
  output.file = "circos_chr19_10_tout_ks.pdf"
)

# # all categories separately in one file 
# # 4 column x 2 lines, format large
# svg("circos_chr24_27_bycategories.svg", width = 20, height = 10)
# par(mfrow = c(2, 4))  # 2 lines x 4 column 
# 
# for (cas_name in paste0("cas", 1:7)) {
#   links_cas <- links_filtered %>% filter(cas == cas_name)
#   if (nrow(links_cas) == 0) next
#   plot_circos(links_cas, karyo_filtered, col_cas, title = cas_name)
# }
# dev.off()





# Filter karyo and links on triple chromosomes
# you need to change the chr if you want to test other triple
chrs_interesse <- c("12", "13", "17")

karyo_filtered <- karyo %>% filter(chr %in% chrs_interesse)

links_filtered <- links_ohno_bin %>% 
  filter(chrA %in% chrs_interesse & chrB %in% chrs_interesse)

# Plot
plot_circos_links(
  links_data   = links_filtered,
  karyo_data   = karyo_filtered,
  col_cas      = col_test,
  win_all = win_all,
  title        = "",
  show.legend  = TRUE,
  output.file  = "circos_chr12_13_17.svg"
)

plot_circos_links_ks(
  links_data   = links_filtered,
  karyo_data   = karyo_filtered,
  col_cas      = col_test,
  win_all = win_all,
  ks_per_win = ks_per_win,
  title        = "",
  show.legend  = TRUE,
  output.file  = "circos_chr21_9_15_ks.pdf"
)

# svg("circos_chr15_21_9_bycategories.svg", width = 20, height = 10)
# par(mfrow = c(2, 4))
# 
# for (cas_name in paste0("cas", 1:7)) {
#   links_cas <- links_filtered %>% filter(cas == cas_name)
#   if (nrow(links_cas) == 0) next
#   plot_circos(
#     links_data = links_cas,
#     karyo_data = karyo_filtered,
#     col_liens  = col_test,
#     title      = cas_name
#   )
# }
# 
# dev.off()






#library(pdftools)

# change pdf to png
#pdf_convert(pdf = "results/circos_chr12_13_17.pdf", format = "png", dpi = 300)
