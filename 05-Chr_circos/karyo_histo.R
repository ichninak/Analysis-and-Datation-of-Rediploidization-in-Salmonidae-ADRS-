# to have karyotype with color of categories in histogram 
# and also plot of density

# ==============
# load library
# ==============

library(circlize)
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)


# ==============
# function
# ==============

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
  
  # read and combine all file of pairs 
  all_pairs <- do.call(rbind, lapply(pair_files, function(f) {
    read.table(f, col.names = c("esoxID", "gene1", "gene2", "flag"),
               stringsAsFactors = FALSE)
  }))
  
  # duplication 
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


plot_density_chr <- function(data, 
                             chr_name,
                             karyo,
                             colours = c("navy", "blue", "cyan", "green", "yellow", "orange", "red"),
                             scale = 2,
                             alpha = 0.8) {
  
  # keep start and end since karyotype
  chr_info <- karyo[karyo$chr == chr_name, ]
  
  if (nrow(chr_info) == 0) {
    stop(paste("Chromosome", chr_name, "non trouvé dans karyo"))
  }
  
  chr_start <- chr_info$start
  chr_end   <- chr_info$end
  chr_length <- chr_end - chr_start
  
  ggplot(data, aes(x = position, y = cas, fill = after_stat(x))) +
    geom_density_ridges_gradient(
      scale = scale,
      alpha = alpha,
      rel_min_height = 0.01
    ) +
    scale_x_continuous(
      limits = c(chr_start, chr_end),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0, 0)
    ) +
    scale_fill_gradientn(
      colours = colours,
      limits = c(chr_start, chr_end),   # 🔑 gradient calé sur start→end
      name = "Position (Mb)",
      labels = scales::label_number(scale = 1e-6)
    ) +
    theme_ridges() +
    labs(
      title = paste("Densité de gènes - Chromosome", chr_name),
      subtitle = paste0("Longueur : ", round(chr_length/1e6, 1), " Mb"),
      x = "Position chromosomique",
      y = "Cas"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(face = "italic")
    ) +
    guides(fill = guide_colorbar(
      barwidth = 15, barheight = 0.8,
      title.position = "top", title.hjust = 0.5
    ))
}

plot_histo_chr <- function(data, 
                           chr_name,
                           karyo,
                           colours = c("navy", "blue", "cyan", "green", "yellow", "orange", "red"),
                           binwidth = 800000,
                           scale = 2,
                           spacing = 2,   
                           alpha = 0.9) {
  
  chr_info <- karyo[karyo$chr == chr_name, ]
  if (nrow(chr_info) == 0) stop(paste("Chromosome", chr_name, "non trouvé"))
  chr_start  <- chr_info$start
  chr_end    <- chr_info$end
  chr_length <- chr_end - chr_start
  
  bins <- seq(chr_start, chr_end, by = binwidth)
  
  data_binned <- data %>%
    mutate(bin = cut(position, breaks = bins, include.lowest = TRUE,
                     labels = bins[-length(bins)])) %>%
    mutate(bin = as.numeric(as.character(bin))) %>%
    count(cas, bin) %>%
    complete(cas, bin = bins[-length(bins)], fill = list(n = 0)) %>%
    group_by(cas) %>%
    mutate(height = n / max(n, na.rm = TRUE) * scale) %>%
    ungroup() %>%
    mutate(cas = factor(cas),
           y_base = as.numeric(cas) * spacing,   
           y_top  = y_base + height)
  
  ggplot(data_binned) +
    geom_rect(aes(xmin = bin, xmax = bin + binwidth,
                  ymin = y_base, ymax = y_top,
                  fill = bin),
              alpha = alpha, colour = "black", linewidth = 0.1) +
    scale_x_continuous(
      limits = c(chr_start, chr_end),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(data_binned$cas)) * spacing,
      labels = levels(data_binned$cas),
      expand = expansion(add = c(0.1, scale + 0.2))
    ) +
    scale_fill_gradientn(
      colours = colours,
      limits = c(chr_start, chr_end),
      name = "Position (Mb)",
      labels = scales::label_number(scale = 1e-6)
    ) +
    labs(
      title = paste("Densité de gènes - Chromosome", chr_name),
      subtitle = paste0("Longueur : ", round(chr_length/1e6, 1), " Mb"),
      x = "Position chromosomique",
      y = "Cas"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(face = "italic"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_colorbar(
      barwidth = 15, barheight = 0.8,
      title.position = "top", title.hjust = 0.5
    ))
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

extract_link_ohno <- function(chr_nbr, links_ohno){
  chr_data <- rbind(
    links_ohno %>% filter(chrA == chr_nbr) %>% select(cas, position = startA),
    links_ohno %>% filter(chrB == chr_nbr) %>% select(cas, position = startB)
  ) %>% distinct()
  
  chr_data$cas <- factor(chr_data$cas, levels = c("LORe4", "LORe3", "LORe2", "LORe1", "AORe"))
  
  return(chr_data)
}

# ==============
# load file
# ==============

karyo <- read.table("../01-Ancestor_circos/karyo_mykiss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
karyo$chr <- as.character(karyo$chr)
coord <- read.delim(
  "../01-Ancestor_circos/coord_mykiss.tsv",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  col.names = c("chr", "start", "end", "gene")
)
write.table(coord, "coord.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pair_files <- c(
  "../01-Ancestor_circos/lore_gar_mykiss1_Salmonidae.tsv",
  "../01-Ancestor_circos/lore_gar_mykiss1_Salmoninae.tsv",
  "../01-Ancestor_circos/lore_gar_mykiss1_NAME11_combine.tsv",
  "../01-Ancestor_circos/lore_gar_mykiss1_NAME13_combine.tsv"
)

resolved_IDs_circos <- readRDS("../02-Barplot_analysis/resolved_IDs_circos.rds")
resolved_IDs <- readRDS("../02-Barplot_analysis/resolved_IDs.rds")

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

ks_per_win    <- readRDS("../07-Ks_value/ks_per_win.rds")

# ==============
# PIPELINE
# ==============

############# density plot ############

# extract for each chromosomes

chr1_data <- extract_link_ohno("1", links_ohno)
chr2_data <- extract_link_ohno("2", links_ohno)
chr3_data <- extract_link_ohno("3", links_ohno)
chr6_data <- extract_link_ohno("6", links_ohno)
chr7_data <- extract_link_ohno("7", links_ohno)
chr10_data <- extract_link_ohno("10", links_ohno)
chr12_data <- extract_link_ohno("12", links_ohno)
chr13_data <- extract_link_ohno("13", links_ohno)
chr14_data <- extract_link_ohno("14", links_ohno)
chr15_data <- extract_link_ohno("15", links_ohno)
chr17_data <- extract_link_ohno("17", links_ohno)
chr18_data <- extract_link_ohno("18", links_ohno)
chr19_data <- extract_link_ohno("19", links_ohno)
chr21_data <- extract_link_ohno("21", links_ohno)
chr23_data <- extract_link_ohno("23", links_ohno)
chr24_data <- extract_link_ohno("24", links_ohno)
chr26_data <- extract_link_ohno("26", links_ohno)
chr27_data <- extract_link_ohno("27", links_ohno)
chr31_data <- extract_link_ohno("31", links_ohno)


# Plot of density for chromosome
# try the spacing you want

plot_histo_chr(chr21_data, "21", karyo, spacing = 2)    
plot_histo_chr(chr21_data, "21", karyo, spacing = 2.5)
plot_histo_chr(chr15_data, "15", karyo, spacing = 2.5)
plot_histo_chr(chr21_data, "21", karyo, spacing = 4)    

# another type of density
plot_density_chr(chr2_data, "2", karyo)

chr_list <- list(
  "1" = chr1_data,
  "2" = chr2_data,
  "3" = chr3_data,
  "6" = chr6_data,
  "7" = chr7_data,
  "10" = chr10_data,
  "12" = chr12_data,
  "13" = chr13_data,
  "14" = chr14_data,
  "15" = chr15_data,
  "17" = chr17_data,
  "18" = chr18_data,
  "19" = chr19_data,
  "21" = chr21_data,
  "23" = chr23_data,
  "24" = chr24_data,
  "26" = chr26_data,
  "27" = chr27_data,
  "31" = chr31_data
)

# see all density plot in same time

plots <- imap(chr_list, ~ plot_density_chr(.x, chr_name = .y, karyo = karyo))
plots2 <- imap(chr_list, ~ plot_histo_chr(.x, chr_name = .y, karyo = karyo, spacing = 3))

plots[["21"]]
plots2[["21"]]

# save all
iwalk(plots, ~ ggsave(paste0("density_chr", .y, ".pdf"),
                      .x, width = 8, height = 5, dpi = 300))

iwalk(plots2, ~ ggsave(paste0("histo_chr", .y, ".pdf"),
                       .x, width = 8, height = 5, dpi = 300))


wrap_plots(plots, ncol = 2)



######## version karyo histogram #########


cas_list <- list(AORe = win_cas_1, LORe1 = win_cas_2, LORe2 = win_cas_3,
                 LORe3 = win_cas_4, LORe4 = win_cas_5)

all_cas <- bind_rows(cas_list, .id = "cas")

# 1. Define the levels uniquely with chromosomes present in data
chrs_presents <- unique(all_cas$chr)
chrs_ordonnes <- c(as.character(1:32), "X", "Y")
chrs_final <- chrs_ordonnes[chrs_ordonnes %in% chrs_presents]


all_cas <- all_cas %>%
  mutate(chr = factor(chr, levels = chrs_final)) %>%
  arrange(chr, cas) %>%
  mutate(track = factor(paste(chr, cas, sep = "_"),
                        levels = unique(paste(chr, cas, sep = "_"))))

# 2. Chromosome outlines (also filtered)
h <- 0.35

chr_contour <- karyo %>%
  filter(chr %in% chrs_final) %>%
  group_by(chr) %>%
  summarise(xmin = 0, xmax = max(end)) %>%
  mutate(chr = factor(chr, levels = chrs_final),
         ymin = as.numeric(chr) - h,
         ymax = as.numeric(chr) + h)

# 3. Labels Y
chr_breaks <- seq_along(chrs_final)
chr_labels <- chrs_final

# 4. Plot
g <- ggplot(all_cas, aes(xmin = win_start, xmax = win_end,
                         ymin = as.numeric(chr) - h,
                         ymax = as.numeric(chr) + h,
                         fill = cas)) +
  geom_rect() +
  geom_vline(xintercept = seq(0, 100e6, by = 25e6),
             linetype = "dashed", color = "grey40",
             alpha = 0.6, linewidth = 0.3) +
  geom_rect(data = chr_contour,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linewidth = 0.3,
            inherit.aes = FALSE) +
  scale_y_continuous(breaks = chr_breaks,
                     labels = chr_labels,
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(0, 100e6, by = 25e6),
                     labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  scale_fill_manual(values = col_test, name = "Categories") +
  labs(title = "Karyotype with rediploidization categories", x = "Position", y = "Chromosome") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

print(g)
ggsave("karyotype_5cas.pdf", plot = g, width = 14, height = 10)





######## version karyo histogram + Ks curve #########


cas_list <- list(AORe = win_cas_1, LORe1 = win_cas_2, LORe2 = win_cas_3,
                 LORe3 = win_cas_4, LORe4 = win_cas_5)
all_cas <- bind_rows(cas_list, .id = "cas") %>%
  filter(chr %in% chrs_final) %>%
  mutate(chr = factor(chr, levels = chrs_final))

ks_per_win <- readRDS("../07-Ks_value/ks_per_win.rds")
ks_df <- ks_per_win %>%
  filter(chr %in% chrs_final, !is.na(Ks_mean)) %>%
  mutate(chr = factor(chr, levels = chrs_final))

karyo <- karyo %>% mutate(chr = as.character(chr))

ks_max_df <- ks_df %>%
  group_by(chr) %>%
  summarise(ks_max = max(Ks_mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(chr = as.character(chr))

all_cas <- all_cas %>%
  mutate(chr = as.character(chr)) %>%
  left_join(ks_max_df, by = "chr") %>%
  mutate(chr = factor(chr, levels = chrs_final),
         ymin = -ks_max * 0.12, ymax = 0)

chr_contour <- karyo %>%
  filter(chr %in% chrs_final) %>%
  left_join(ks_max_df, by = "chr") %>%
  mutate(chr = factor(chr, levels = chrs_final),
         xmin = 0, xmax = end,
         ymin = -ks_max * 0.12, ymax = 0)

chr_limits <- karyo %>%
  filter(chr %in% chrs_final) %>%
  mutate(chr = factor(chr, levels = chrs_final))

g <- ggplot() +
  geom_blank(data = chr_limits, aes(x = 0, y = 0)) +
  geom_blank(data = chr_limits, aes(x = end / 1e6, y = 0)) +
  
  geom_line(data = ks_df,
            aes(x = mid / 1e6, y = Ks_mean),
            color = "grey50", linewidth = 0.4) +
  geom_point(data = ks_df,
             aes(x = mid / 1e6, y = Ks_mean, color = cas),
             size = 1.2, alpha = 0.9) +
  geom_smooth(data = ks_df,
              aes(x = mid / 1e6, y = Ks_mean),
              method = "loess", span = 0.4, se = FALSE,
              color = "black", linewidth = 0.5) +
  
  geom_rect(data = chr_contour,
            aes(xmin = xmin / 1e6, xmax = xmax / 1e6,
                ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 0.3) +
  
  geom_rect(data = all_cas,
            aes(xmin = win_start / 1e6, xmax = win_end / 1e6,
                ymin = ymin, ymax = ymax,
                fill = cas)) +
  
  facet_wrap(~ chr, scales = "free", ncol = 6) +
  scale_color_manual(values = col_test, name = "Catégorie") +
  scale_fill_manual(values = col_test, name = "Catégorie") +
  labs(title = "Profil Ks et catégories de rediploïdisation par chromosome",
       x = "Position (Mb)", y = "Ks moyen") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

print(g)
ggsave("karyotype_ks_facets.pdf", plot = g, width = 16, height = 12)
ggsave("karyotype_ks_facets.png", plot = g, width = 16, height = 12, dpi = 300)
