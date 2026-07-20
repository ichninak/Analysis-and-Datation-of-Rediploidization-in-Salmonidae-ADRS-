# ======================
# load libraries
# ======================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(ggridges)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

# ================
# function
# ================

process_density <- function(windows_ratio, gff_file, output_density) {
  
  # Classification AORE / LORE / equal

  df <- windows_ratio %>%
    mutate(class = case_when(
      ratio_AORE > 0.55 ~ "AORE",
      ratio_AORE < 0.45 ~ "LORE",
      TRUE              ~ "equal"
    )) %>%
    filter(total > 0)
  
  cat("--- Classification ---\n")
  print(df %>% group_by(class) %>%
          summarise(n_windows = n(),
                    total_relations = sum(total, na.rm = TRUE),
                    .groups = "drop"))
  
  # cleaning of column of type list
  list_cols <- sapply(df, is.list)
  if (any(list_cols)) {
    cat("Colonnes 'list' supprimées :",
        paste(names(df)[list_cols], collapse = ", "), "\n")
    df <- df[, !list_cols]
  }
  
  # Construction of GRanges (windows)
  windows_gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$win_start, end = df$win_end)
  )
  
  # Import GFF + rename chromosomes
  gff <- import(gff_file)
  
  regions <- gff[gff$type == "region" & !is.na(gff$chromosome)]
  map <- setNames(as.character(regions$chromosome),
                  as.character(seqnames(regions)))
  map <- map[map != "Unknown"]
  
  gff2     <- renameSeqlevels(gff, map)
  genes_gr <- gff2[gff2$type == "gene"]
  
  # Harmonization of chromosomes
  common_chr <- intersect(
    unique(as.character(seqnames(windows_gr))),
    unique(as.character(seqnames(genes_gr)))
  )
  
  windows_gr <- keepSeqlevels(windows_gr, common_chr, pruning.mode = "coarse")
  genes_gr   <- keepSeqlevels(genes_gr,   common_chr, pruning.mode = "coarse")
  

  # count and density 
  df$gene_count <- 0L
  hits <- countOverlaps(
    GRanges(seqnames = df$chr,
            ranges   = IRanges(start = df$win_start, end = df$win_end)),
    genes_gr
  )
  df$gene_count   <- hits
  df$window_size  <- df$win_end - df$win_start + 1
  df$gene_density <- df$gene_count / (df$window_size / 1e6)  # genes/Mb
  
  # save file
  write.table(df, output_density,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("✓ writing file :", output_density, "\n")
  invisible(df)
}

# ======================
# load files
# ======================
gff_file <- "GCF_013265735.2_USDA_OmykA_1.1_genomic.gff"

win_salmonidae = readRDS("../Ancestor_Circos/win_salmoninae.rds")
win_NAME13 = readRDS("../Ancestor_Circos/win_NAME13.rds")

win_cas_1 <- readRDS("../cas_Circos/win_cas_1.rds")
win_cas_2 <- readRDS("../cas_Circos/win_cas_2.rds")
win_cas_3 <- readRDS("../cas_Circos/win_cas_3.rds")
win_cas_4 <- readRDS("../cas_Circos/win_cas_4.rds")
win_cas_5 <- readRDS("../cas_Circos/win_cas_5.rds")

# ======================
# PIPELINE
# ======================

density_salmonidae <- process_density(win_salmonidae, gff_file, "density_salmonidae.tsv")
density_NAME13     <- process_density(win_NAME13,     gff_file, "density_NAME13.tsv")

# boxplot density 
max_val <- quantile(density_AORe$gene_count, 0.98, na.rm = TRUE)

ggplot(density_salmonidae, aes(x = class, y = gene_count, fill = class)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = c(AORE = "blue", LORE = "red")) +
  coord_cartesian(ylim = c(0, max_val)) +
  theme_bw() +
  labs(
    x = "Classe",
    y = "Number total of genes"
  )

# without limit and with p-value 
ggplot(density_salmonidae, aes(x = class, y = gene_count, fill = class)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = c(AORE = "blue", LORE = "red")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("AORE", "LORE")),
                     label = "p.signif") +
  theme_bw() +
  labs(
    x = "Classe",
    y = "Number total of genes"
  )

# with only AORE and LORE categories 
density_salmonidae %>%
  filter(class %in% c("AORE", "LORE")) %>%
  ggplot(aes(x = class, y = gene_count, fill = class)) +
  geom_jitter(width = 0.4, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("AORE", "LORE")),
                     label = "p.signif") +
  coord_cartesian(ylim = c(0, max_val)) +
  scale_fill_manual(values = c("AORE" = "grey70", "LORE" = "#e41a1c")) +
  theme_bw() +
  labs(
    title = 'Density of genes to compare AORE vs LORE',
    x = "Categories",
    y = "Number total of genes"
  ) 
ggsave("density_of_gene_AOREvsLORE.png", width = 8, height = 8, dpi = 300)

wilcox.test(gene_count ~ class,
            data = density_salmonidae %>% filter(class %in% c("AORE", "LORE")),
            alternative = "two.sided")
wilcox.test(gene_count ~ class,
            data = density_salmonidae %>% filter(class %in% c("AORE", "LORE")),
            alternative = "greater")
wilcox.test(gene_count ~ class,
            data = density_salmonidae %>% filter(class %in% c("AORE", "LORE")),
            alternative = "less")

wilcox.test(gene_count ~ class,
            data = density_NAME13 %>% filter(class %in% c("AORE", "LORE")),
            alternative = "two.sided")
wilcox.test(gene_count ~ class,
            data = density_NAME13 %>% filter(class %in% c("AORE", "LORE")),
            alternative = "greater")
wilcox.test(gene_count ~ class,
            data = density_NAME13 %>% filter(class %in% c("AORE", "LORE")),
            alternative = "less")

