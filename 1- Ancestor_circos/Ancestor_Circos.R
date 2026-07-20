# ======================
# load library
# ======================

library(dplyr)
library(circlize)
library(tidyr)


# ================
# function
# ================

prepare_data <- function(pairs_file, coord_file) {
  
  pairs <- read.table(pairs_file, col.names = c("esoxID", "gene1", "gene2", "flag"), stringsAsFactors = FALSE)
  
  coord <- read.delim(coord_file, header = TRUE, sep = "\t", col.names = c("chr", "start", "end", "gene"),
                      stringsAsFactors = FALSE)
  # remove scaffolds
  coord <- coord[coord$chr != "Unplaced Scaffold", ]
  # put classes/categories
  pairs$class <- ifelse(pairs$flag == 1, "LORE", "AORE")
  
  genes1 <- data.frame(
    gene = pairs$gene1,
    esoxID = pairs$esoxID,
    class = pairs$class
  )
  genes2 <- data.frame(
    gene = pairs$gene2,
    esoxID = pairs$esoxID,
    class = pairs$class
  )
  
  # have genes and their classes/categories separately 
  genes_flag <- rbind(genes1, genes2)
  
  # avoid double 
  genes_flag <- unique(genes_flag)
  
  # merge both table of flags with coord by genes 
  genes_full <- merge(coord, genes_flag, by ="gene")
  
  # avoid numeric 
  genes_full$chr <- as.character(genes_full$chr)
  
  # keep only genes with classes 
  genes_full <- genes_full[genes_full$class %in% c("AORE", "LORE"), ]
  
  return(genes_full)
}

# to have windows of 2000kb on each chromosomes to study more precisely rediploidization
compute <- function(karyo, genes_full, win_size = 2000000, step = 2000000) {
  res_list <- list()
  
  for (j in seq_len(nrow(karyo))) {
    
    chr_j <- karyo$chr[j]
    chr_start <- karyo$start[j]
    chr_end <- karyo$end[j]
    
    g_chr <- genes_full[genes_full$chr == chr_j, ]
    
    if (nrow(g_chr) == 0) next
    if ((chr_end - chr_start + 1) < win_size) next
    
    win_starts <- seq(chr_start, chr_end - win_size + 1, by = step)
    win_ends <- win_starts + win_size - 1
    n_wins <- length(win_starts)
    
    chr_res <- data.frame(
      chr = chr_j,
      win_start = win_starts,
      win_end = win_ends,
      n_AORE = integer(n_wins),
      n_LORE = integer(n_wins),
      total = integer(n_wins),
      ratio_AORE = numeric(n_wins),
      esoxIDs = I(vector("list", n_wins)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(win_starts)) {
      s <- win_starts[i]
      e <- win_ends[i]
      
      sub <- g_chr[g_chr$start <= e & g_chr$end >= s, ]
      
      n_AORE <- sum(sub$class == "AORE")
      n_LORE <- sum(sub$class == "LORE")
      total <- n_AORE + n_LORE
      
      chr_res$n_AORE[i] <- n_AORE
      chr_res$n_LORE[i] <- n_LORE
      chr_res$total[i] <- total
      chr_res$ratio_AORE[i] <- (n_AORE + 0.1) / (n_AORE + n_LORE + 0.2)
      chr_res$esoxIDs[[i]] <- unique(sub$esoxID)
      
    }
    
    chr_res$state <- "neutral"
    chr_res$state[chr_res$ratio_AORE > 0.55] <- "AORE"
    chr_res$state[chr_res$ratio_AORE < 0.45] <- "LORE"
    
    res_list[[chr_j]] <- chr_res
  }
  
  windows_ratio <- do.call(rbind, res_list)
  
  # pour les boxplot
  windows_ratio$mid <- (windows_ratio$win_start + windows_ratio$win_end) /2
  
  return(windows_ratio)
}

# to generate circos
circos_precise_svg <- function(karyo, links_file, name_output, windows_ratio) {
  library(circlize)
  links <- read.table(links_file, header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
  links$chrA <- as.character(links$chrA)
  links$chrB <- as.character(links$chrB)
  windows_ratio$chr <- as.character(windows_ratio$chr)
  karyo$chr <- as.character(karyo$chr)
  valid_chr <- karyo$chr
  
  links <- links[
    links$chrA %in% valid_chr &
      links$chrB %in% valid_chr,
  ]
  
  links$class <- ifelse(links$flag == 1, "LORE", "AORE")
  col_map <- c(AORE="#2C7FB8", LORE="#D7301F")
  
  svg(name_output, width = 14, height = 14)
  
  circos.clear()
  
  # circle
  circos.par(
    start.degree = 90,
    gap.after = rep(2, nrow(karyo)),
    #track.margin = c(0.001, 0.001),
    cell.padding = c(0,0,0,0)
  )
  
  circos.initialize(
    factors = karyo$chr,
    xlim = cbind(karyo$start, karyo$end)
  )
  
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
  
  # colors in chromosomes
  windows_ratio$chr <- as.character(windows_ratio$chr)
  
  windows_ratio$state <- "neutral"
  windows_ratio$state[windows_ratio$ratio_AORE > 0.6] <- "AORE"
  windows_ratio$state[windows_ratio$ratio_AORE < 0.4] <- "LORE"
  
  state_col <- c(
    AORE = "#2C7FB8",   # bleu
    LORE = "#D7301F",   # rouge
    neutral = "grey85"
  )
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.08,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- windows_ratio[windows_ratio$chr == chr, ]
      
      if(nrow(sub) == 0) return()
      
      for(i in seq_len(nrow(sub))) {
        circos.rect(
          xleft   = sub$win_start[i],
          ybottom = 0,
          xright  = sub$win_end[i],
          ytop    = 1,
          col     = state_col[sub$state[i]],
          border  = NA
        )
      }
    }
  )
  
  # links
  for(i in seq_len(nrow(links))){
    cl <- links$class[i]
    p1 <- (links$startA[i] + links$endA[i]) / 2
    p2 <- (links$startB[i] + links$endB[i]) / 2
    
    circos.link(
      links$chrA[i], p1,
      links$chrB[i], p2,
      col = adjustcolor(col_map[cl]),
      border = NA,
      lwd = 1
    )
  }
  
  
  dev.off()
  circos.clear()
}

circos_precise_pdf <- function(karyo, links_file, name_output, windows_ratio) {
  library(circlize)
  links <- read.table(links_file, header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
  links$chrA <- as.character(links$chrA)
  links$chrB <- as.character(links$chrB)
  windows_ratio$chr <- as.character(windows_ratio$chr)
  karyo$chr <- as.character(karyo$chr)
  valid_chr <- karyo$chr
  
  links <- links[
    links$chrA %in% valid_chr &
      links$chrB %in% valid_chr,
  ]
  
  links$class <- ifelse(links$flag == 1, "LORE", "AORE")
  col_map <- c(AORE="#2C7FB8", LORE="#D7301F")
  
  pdf(name_output, width = 14, height = 14)
  
  circos.clear()
  
  
  circos.par(
    start.degree = 90,
    gap.after = rep(2, nrow(karyo)),
    #track.margin = c(0.001, 0.001),
    cell.padding = c(0,0,0,0)
  )
  
  circos.initialize(
    factors = karyo$chr,
    xlim = cbind(karyo$start, karyo$end)
  )
  
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
  
  windows_ratio$chr <- as.character(windows_ratio$chr)
  
  windows_ratio$state <- "neutral"
  windows_ratio$state[windows_ratio$ratio_AORE > 0.6] <- "AORE"
  windows_ratio$state[windows_ratio$ratio_AORE < 0.4] <- "LORE"
  
  state_col <- c(
    AORE = "#2C7FB8",   # bleu
    LORE = "#D7301F",   # rouge
    neutral = "grey85"
  )
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.08,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- get.cell.meta.data("sector.index")
      sub <- windows_ratio[windows_ratio$chr == chr, ]
      
      if(nrow(sub) == 0) return()
      
      for(i in seq_len(nrow(sub))) {
        circos.rect(
          xleft   = sub$win_start[i],
          ybottom = 0,
          xright  = sub$win_end[i],
          ytop    = 1,
          col     = state_col[sub$state[i]],
          border  = NA
        )
      }
    }
  )
  
  for(i in seq_len(nrow(links))){
    cl <- links$class[i]
    p1 <- (links$startA[i] + links$endA[i]) / 2
    p2 <- (links$startB[i] + links$endB[i]) / 2
    
    circos.link(
      links$chrA[i], p1,
      links$chrB[i], p2,
      col = adjustcolor(col_map[cl]),
      border = NA,
      lwd = 1
    )
  }
  
  
  dev.off()
  circos.clear()
}

prepare_data_circos <- function(coord, pair_file, name_output) {
  pairs = read.table(pair_file, sep = "\t", col.names = c("esoxID", "geneA", "geneB", "flag"), stringsAsFactors = FALSE)
  
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

# ======================
# load files
# ======================

# karyotype
karyo <- read.table("karyo_mykiss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
karyo$chr <- as.character(karyo$chr)
# coordinate of karyotype
coord <- read.delim(
  "coord_mykiss.tsv",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("chr", "start", "end", "gene")
)
write.table(coord, "coord.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# preparation of links for the pipeline

prepare_data_circos(coord, "lore_gar_mykiss1_Salmonidae.tsv", "link_salmonidae.tsv")

prepare_data_circos(coord, "lore_gar_mykiss1_Salmoninae.tsv", "link_salmoninae.tsv")

prepare_data_circos(coord, "lore_gar_mykiss1_NAME11_combine.tsv", "link_NAME11.tsv")

prepare_data_circos(coord, "lore_gar_mykiss1_NAME13_combine.tsv", "link_NAME13.tsv")

# ======================
# PIPELINE
# ======================

# circos de chaque ancetre

## salmonidae ancestor

genes_salmonidae <- prepare_data("lore_gar_mykiss1_salmonidae.tsv", "coord.tsv")
win_salmonidae <- compute(karyo, genes_salmonidae) 

# circos
circos_precise_svg(karyo, "link_salmonidae.tsv", "circos_salmonidae.svg", win_salmonidae)
circos_precise_pdf(karyo, "link_salmonidae.tsv", "circos_salmonidae.pdf", win_salmonidae)

## salmoninae ancestor 

genes_salmoninae <- prepare_data("lore_gar_mykiss1_Salmoninae.tsv", "coord.tsv")
win_salmoninae <- compute(karyo, genes_salmoninae)

# circos
circos_precise_svg(karyo, "link_salmoninae.tsv", "circos_salmoninae.svg", win_salmoninae)
circos_precise_pdf(karyo, "link_salmoninae.tsv", "circos_salmoninae.pdf", win_salmoninae)

## NAME11 ancestor
### NAME11 = NAME11 + NMAE12


genes_NAME11 <- prepare_data("lore_gar_mykiss1_NAME11_combine.tsv", "coord.tsv")
win_NAME11 <- compute(karyo, genes_NAME11)

# circos

circos_precise_svg(karyo, "link_NAME11.tsv", "circos_NAME11.svg", win_NAME11)
circos_precise_pdf(karyo, "link_NAME11.tsv", "circos_NAME11.pdf", win_NAME11)

## NAME13 ancestor
### NAME13 = NAME13 + NAME14

genes_NAME13 <- prepare_data("lore_gar_mykiss1_NAME13_combine.tsv", "coord.tsv")
win_NAME13 <- compute(karyo, genes_NAME13)

# circos

circos_precise_svg(karyo, "link_NAME13.tsv", "circos_NAME13.svg", win_NAME13)
circos_precise_pdf(karyo, "link_NAME13.tsv", "circos_NAME13.pdf", win_NAME13)

message("end of pipeline")

# export file


saveRDS(win_salmonidae,   "win_salmonidae.rds")
saveRDS(genes_salmonidae, "gene_salmonidae.rds")

saveRDS(win_salmoninae,   "win_salmoninae.rds")
saveRDS(genes_salmoninae, "gene_salmoninae.rds")

saveRDS(win_NAME11,   "win_NAME11.rds")
saveRDS(genes_NAME11, "gene_NAME11.rds")

saveRDS(win_NAME13,   "win_NAME13.rds")
saveRDS(genes_NAME13, "gene_NAME13.rds")


#write.table(win_salmonidae, "win_salmonidae.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#write.table(genes_salmonidae, "gene_salmonidae.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#write.table(win_salmoninae, "win_salmoninae.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#write.table(genes_salmoninae, "gene_salmoninae.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#write.table(win_NAME11, "win_NAME11.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#write.table(genes_NAME11, "gene_NAME11.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#write.table(win_NAME13, "win_NAME13.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#write.table(genes_NAME13, "gene_NAME13.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

