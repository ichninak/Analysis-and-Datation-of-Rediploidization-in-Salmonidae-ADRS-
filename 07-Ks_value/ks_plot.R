
# ============================================================
# Analyse Ka/Ks des paires d'ohnologues - O. mykiss
# ============================================================

################
# load libraries
################

library(ggplot2)
library(dplyr)
library(patchwork)  
library(scales)
library(tidyr)
library(zoo)

############
# load files 
############

data <- read.table("all_mykiss.axt.kaks", 
                   header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE,
                   na.strings = c("NA", "", "-0"))

cat("Nombre de paires brutes :", nrow(data), "\n")

# load merged.tsv (all ohnologs)
merged <- read.table("merged.tsv",
                     header = FALSE,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     col.names = c("esoxID", "gene1", "gene2", "flag"))

cat("merged.tsv :", nrow(merged), "lignes\n")
head(merged)

# Construction of "gene1&gene2" 
merged$pair1 <- paste(merged$gene1, merged$gene2, sep = "&")
merged$pair2 <- paste(merged$gene2, merged$gene1, sep = "&")

# Dictionary 
pair_dict <- c(
  setNames(merged$flag, merged$pair1),
  setNames(merged$flag, merged$pair2)
)

# load the file .kaks (Kaks_calculator)
kaks <- read.table("all_mykiss.axt.kaks",
                   header = TRUE,        
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   fill = TRUE,
                   quote = "")

cat("kaks :", nrow(kaks), "lignes\n")
head(kaks)

# first column "gene1&gene2"
colnames(kaks)[1] <- "pair"

# add AORE/LORE
kaks$flag <- pair_dict[kaks$pair]

# check pairs 
n_unmatched <- sum(is.na(kaks$flag))
cat("Paires non trouvées dans merged.tsv :", n_unmatched, "\n")

# Separation 
aore <- kaks %>% filter(flag == 0)
lore <- kaks %>% filter(flag == 1)
unmatched <- kaks %>% filter(is.na(flag))

cat("AORE :", nrow(aore), "lignes\n")
cat("LORE :", nrow(lore), "lignes\n")
cat("Non matché :", nrow(unmatched), "lignes\n")

# Save 
write.table(aore, "AORE.kaks", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lore, "LORE.kaks", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(unmatched, "unmatched.kaks", sep = "\t", row.names = FALSE, quote = FALSE)




############
# PIPELINE
############



####### lore aore


cat("Nombre de paires brutes :", nrow(lore), "\n")
cat("Nombre de paires brutes :", nrow(aore), "\n")

# clean data
lore <- lore %>%
  rename(KaKs = `Ka.Ks`,
         Pvalue = `P.Value.Fisher.`) %>%
  mutate(across(c(Ka, Ks, KaKs, Pvalue, Length), 
                ~ as.numeric(as.character(.))))

aore <- aore %>%
  rename(KaKs = `Ka.Ks`,
         Pvalue = `P.Value.Fisher.`) %>%
  mutate(across(c(Ka, Ks, KaKs, Pvalue, Length), 
                ~ as.numeric(as.character(.))))

# keep quality
lore_clean <- lore %>%
  filter(!is.na(Ka), !is.na(Ks), !is.na(KaKs)) %>%   
  filter(is.finite(KaKs)) %>%                         
  filter(Ks >= 0.01) %>%                              
  filter(Ks <= 2) %>%                                 
  filter(Pvalue < 0.05)                               

aore_clean <- aore %>%
  filter(!is.na(Ka), !is.na(Ks), !is.na(KaKs)) %>%   
  filter(is.finite(KaKs)) %>%                         
  filter(Ks >= 0.01) %>%                              
  filter(Ks <= 2) %>%                                 
  filter(Pvalue < 0.05)                               

cat("Number of pairs after filterage:", nrow(lore_clean), "\n")
cat("Pairs remove:", nrow(lore) - nrow(lore_clean), "\n\n")

cat("Number of pairs after filterage:", nrow(aore_clean), "\n")
cat("Pairs remove:", nrow(aore) - nrow(aore_clean), "\n\n")

# Stats
cat("Median Ka/Ks  :", round(median(aore_clean$KaKs), 3), "\n")
cat("Mean Ka/Ks  :", round(mean(aore_clean$KaKs), 3), "\n")
cat("Median Ks     :", round(median(aore_clean$Ks), 3), "\n")
cat("% pairs Ka/Ks < 1:", round(mean(aore_clean$KaKs < 1) * 100, 1), "%\n") # Purifying selection
cat("% paires Ka/Ks > 1:", round(mean(aore_clean$KaKs > 1) * 100, 1), "%\n") # Positive selection

# for colors 
lore_clean <- lore_clean %>%
  mutate(selection = case_when(
    KaKs < 0.1 ~ "Strong Purifying selection",
    KaKs < 1   ~ "Purifying selection",
    KaKs == 1  ~ "Neutral",
    KaKs > 1   ~ "Positive selection"
  ))

aore_clean <- aore_clean %>%
  mutate(selection = case_when(
    KaKs < 0.1 ~ "Strong Purifying selection",
    KaKs < 1   ~ "Purifying selection",
    KaKs == 1  ~ "Neutral",
    KaKs > 1   ~ "Positive selection"
  ))

# ============================================================
# PLOT A — Histogram Ka/Ks
# ============================================================
p1_lore <- ggplot(lore_clean, aes(x = KaKs)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median(data_clean$KaKs), 
             color = "darkgreen", linetype = "dotted", linewidth = 1) +
  annotate("text", x = 1.05, y = Inf, label = "Ka/Ks = 1\n(neutral)", 
           hjust = 0, vjust = 1.5, color = "red", size = 3) +
  labs(title = "A1. Distribution Ka/Ks of ohnologs LORE",
       x = "Ka/Ks", y = "Number of pairs") +
  theme_bw(base_size = 12)

p1_aore <- ggplot(aore_clean, aes(x = KaKs)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median(data_clean$KaKs), 
             color = "darkgreen", linetype = "dotted", linewidth = 1) +
  annotate("text", x = 1.05, y = Inf, label = "Ka/Ks = 1\n(neutral)", 
           hjust = 0, vjust = 1.5, color = "red", size = 3) +
  labs(title = "A2. Distribution Ka/Ks of ohnologs AORE",
       x = "Ka/Ks", y = "Number of pairs") +
  theme_bw(base_size = 12)

print(p1_lore)
print(p1_aore)

# ============================================================
# PLOT B — Scatter Ka vs Ks
# ============================================================
p2_lore <- ggplot(lore_clean, aes(x = Ks, y = Ka, color = selection)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_abline(slope = 0.5, intercept = 0, color = "orange", linetype = "dotted") +
  scale_color_manual(values = c("Strong Purifying selection" = "darkblue",
                                "Purifying selection" = "steelblue",
                                "Positive selection" = "red")) +
  labs(title = "B1. Ka vs Ks by pairs of ohnologs LORE",
       x = "Ks (synonymous)", y = "Ka (non-synonymous)",
       color = "Régime") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p2_aore <- ggplot(aore_clean, aes(x = Ks, y = Ka, color = selection)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_abline(slope = 0.5, intercept = 0, color = "orange", linetype = "dotted") +
  scale_color_manual(values = c("Strong Purifying selection" = "darkblue",
                                "Purifying selection" = "steelblue",
                                "Positive selection" = "red")) +
  labs(title = "B2. Ka vs Ks by pairs of ohnologs AORE",
       x = "Ks (synonymous)", y = "Ka (non-synonymous)",
       color = "Régime") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p2_lore)
print(p2_aore)
# ============================================================
# PLOT C — Density of Ks 
# ============================================================

p3_lore <- ggplot(lore_clean, aes(x = Ks)) +
  geom_density(fill = "lightgreen", alpha = 0.6, color = "darkgreen") +
  geom_vline(xintercept = median(data_clean$Ks), 
             color = "darkred", linetype = "dashed") +
  labs(title = "C1. Distribution of Ks - signature of the duplication of LORE",
       x = "Ks", y = "Density") +
  xlim(0, 1) +
  theme_bw(base_size = 12)

p3_aore <- ggplot(aore_clean, aes(x = Ks)) +
  geom_density(fill = "lightgreen", alpha = 0.6, color = "darkgreen") +
  geom_vline(xintercept = median(data_clean$Ks), 
             color = "darkred", linetype = "dashed") +
  labs(title = "C2. Distribution of Ks - signature of the duplication of AORE",
       x = "Ks", y = "Density") +
  xlim(0, 1) +
  theme_bw(base_size = 12)

print(p3_lore)
print(p3_aore)
# ============================================================
# combined plot
# ============================================================
final_plot <- (p1_lore | p1_aore) / (p2_lore | p2_aore) / (p3_lore | p3_aore) +
  plot_annotation(title = "Evolution of ohnologs in O. mykiss",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

print(final_plot)

# --- Sauvegarde ---
ggsave("ohnologues_kaks_analysis.pdf", final_plot, 
       width = 12, height = 10)
ggsave("ohnologues_kaks_analysis_lore_aore.png", final_plot, 
       width = 12, height = 10, dpi = 300)



#############
# split by rediploidization categories 
##############


# files of categories
cas_file <- read.table("GeneID_par_cas.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

# combine NAME11 + NAME12 and NAME13 + NAME14
cas_file <- cas_file %>%
  mutate(cas = case_when(
    cas == "cas4" | cas == "cas5" ~ "cas4",
    cas == "cas6" | cas == "cas7" ~ "cas5",
    TRUE ~ cas # keep other same 
  ))

cat("GeneID_par_cas.tsv :", nrow(cas_file), "lignes\n")
cat("Catégories trouvées :", paste(sort(unique(cas_file$cas)), collapse = ", "), "\n\n")

# ============================================================
# 2. Dictioniary geneID -> cas (categories)
# ============================================================
gene_to_cas <- setNames(cas_file$cas, as.character(cas_file$gene))

# ============================================================
# 3. annotation of each pairs of merged.tsv with his 'cas' == categories
# ============================================================
merged$cas_gene1 <- gene_to_cas[as.character(merged$gene1)]
merged$cas_gene2 <- gene_to_cas[as.character(merged$gene2)]

# check if the both gene have a categories 
# else keep the one not NA
merged$cas <- ifelse(!is.na(merged$cas_gene1), merged$cas_gene1, merged$cas_gene2)

# Diagnostic
cat("Paires annotées :", sum(!is.na(merged$cas)), "/", nrow(merged), "\n")

# Check the conflict 
conflicts <- merged %>% 
  filter(!is.na(cas_gene1) & !is.na(cas_gene2) & cas_gene1 != cas_gene2)
cat("Conflits (gene1 et gene2 dans des cas différents) :", nrow(conflicts), "\n\n")

# ============================================================
# Build the key "gene1&gene2" to match with kaks
# ============================================================
merged$pair1 <- paste(merged$gene1, merged$gene2, sep = "&")
merged$pair2 <- paste(merged$gene2, merged$gene1, sep = "&")

# Dictioniary of pair -> cas
pair_to_cas <- c(
  setNames(merged$cas, merged$pair1),
  setNames(merged$cas, merged$pair2)
)

# Annotation of file kaks
kaks$cas <- pair_to_cas[kaks$pair]

cat("Lignes kaks annotées :", sum(!is.na(kaks$cas)), "/", nrow(kaks), "\n\n")

# ============================================================
# loop create a file by categories 
# ============================================================
for (c in paste0("cas", 1:5)) {
  
  pairs_cas <- merged %>% filter(cas == c)
  
  kaks_cas <- kaks %>% filter(cas == c)

  write.table(pairs_cas,
              file = paste0("pairs_", c, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(kaks_cas,
              file = paste0("kaks_", c, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(sprintf("%s : %d paires | %d lignes kaks\n", 
              c, nrow(pairs_cas), nrow(kaks_cas)))
}


# clean and filter kaks
kaks_clean <- kaks %>%
  filter(!is.na(cas)) %>%
  mutate(Ka.Ks = as.numeric(Ka.Ks),
         Ks = as.numeric(Ks),
         Ka = as.numeric(Ka)) %>%
  filter(!is.na(Ka.Ks), Ks > 0.01, Ks < 2) %>%
  filter(!is.na(Ka), Ka >= 0)

# color for categories
col_cas <- c(AORe = "grey60", 
             LORe1 = "#e41a1c", 
             LORe2 = "#377eb8",
             LORe3 = "#4daf4a", 
             LORe4 = "#984ea3"
)


kaks_clean <- kaks_clean %>%
  mutate(cas = recode(cas,
                      "cas1" = "AORe",
                      "cas2" = "LORe1",
                      "cas3" = "LORe2",
                      "cas4" = "LORe3",
                      "cas5" = "LORe4"))

# boxplot ka by categories
ggplot(kaks_clean, aes(x = cas, y = Ka, fill = cas)) +
  geom_jitter(alpha = 0.1, col = "grey") +
  scale_fill_manual(values = col_cas) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(title = "Ka by categories",
       x = "", y = "Ka") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Distribution_ka_by_categories.pdf", width = 12, height = 8, dpi = 300)

#stat
pairwise.wilcox.test(kaks_clean$Ka.Ks, kaks_clean$cas, p.adjust.method = "bonferroni", alternative = "less")

# Boxplot Ka/Ks by categories
ggplot(kaks_clean, aes(x = cas, y = Ka.Ks, fill = cas)) +
  geom_jitter(alpha = 0.1, col = "grey") +
  scale_fill_manual(values = col_cas) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(title = "Ka/Ks by categories",
       x = "", y = "Ka/Ks") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Distribution_ka_ks_by_categories.pdf", width = 12, height = 8, dpi = 300)

# stat
pairwise.wilcox.test(kaks_clean$Ka.Ks, kaks_clean$cas, p.adjust.method = "bonferroni", alternative = "greater")





# Distribution Ks
ggplot(kaks_clean, aes(x = cas, y = Ks, fill = cas)) +
  geom_jitter(alpha = 0.1, col = "grey") +
  scale_fill_manual(values = col_cas) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Distribution of Ks by rediploidization catégories",
       x = NULL,
       fill = "categories") +
  theme_minimal()

# stats
pairwise.wilcox.test(kaks_clean$Ks, kaks_clean$cas, p.adjust.method = "bonferroni", alternative = "less")
#save
ggsave("Distribution_ks_by_cas.pdf", width = 12, height = 8, dpi = 300)


ggplot(kaks_clean, aes(x = Ks, color = cas, fill = cas)) +
  geom_density(alpha = 0.15, linewidth = 1) +
  labs(title = "Distribution Ks par catégorie cas",
       x = "Ks", y = "Densité",
       color = "Catégorie", fill = "Catégorie") +
  theme_minimal() +
  theme(legend.position = "right")

###############
# ks par chromosome
#################

# ============
# load files 
# ============

# files of pairs
pairs <- read.table("merged.tsv",
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("esoxID", "gene1", "gene2", "flag"))

# categories file 
gene_chr <- read.table("GeneID_par_cas.tsv", header = TRUE, 
                       sep = "\t", stringsAsFactors = FALSE)

# file KaKs
kaks <- read.table("all_mykiss.axt.kaks", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)


# extract gene1 and gene2 

kaks <- kaks %>%
  separate(Sequence, into = c("gene1", "gene2"), sep = "&", remove = FALSE) %>%
  select(gene1, gene2, Ka, Ks, Ka.Ks = `Ka.Ks`) %>%
  mutate(Ks = as.numeric(Ks),
         Ka = as.numeric(Ka),
         Ka.Ks = as.numeric(Ka.Ks))

# annotate each pairs 

gene_chr$gene <- as.character(gene_chr$gene)

kaks_annotated <- kaks %>%
  left_join(gene_chr, by = c("gene1" = "gene")) %>%
  rename(chr1 = chr) %>%
  left_join(gene_chr, by = c("gene2" = "gene")) %>%
  rename(chr2 = chr) %>%
  filter(!is.na(chr1), !is.na(chr2), !is.na(Ks))

# filter

kaks_filtered <- kaks_annotated %>%
  filter(Ks > 0.01,    
         Ks < 2)        

write.table(kaks_filtered, sep = "\t", quote = FALSE, row.names = FALSE,)


# Duplicate each pairs 

chr_ks <- bind_rows(
  kaks_filtered %>% select(chr = chr1, partner_chr = chr2, gene = gene1, Ks),
  kaks_filtered %>% select(chr = chr2, partner_chr = chr1, gene = gene2, Ks)
)

# order by chr (choose what you want)
chr_ks$chr <- factor(chr_ks$chr, 
                     levels = unique(chr_ks$chr)[order(as.numeric(
                       gsub("[^0-9]", "", unique(chr_ks$chr))))])

# order by ks (choose what you want)
chr_order <- chr_ks %>%
  group_by(chr) %>%
  summarise(ks_median = median(Ks, na.rm = TRUE)) %>%
  arrange(ks_median) %>%
  pull(chr)

chr_ks$chr <- factor(chr_ks$chr, levels = chr_order)


# PLOT 1 : Violin + boxplot — distribution Ks by chromosome

p1 <- ggplot(chr_ks, aes(x = chr, y = Ks, fill = chr)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Distribution of Ks by chromosome",
       subtitle = "",
       x = "Chromosome", y = "Ks")

ggsave("ks_violin_by_chromosome.png", p1, width = 14, height = 6, dpi = 300)

# PLOT 2 : Density for comparison ohnologs chr

# adapt with each chr you want to test 
p2 <- chr_ks %>%
  filter(chr %in% c("6", "26")) %>%   # here
  ggplot(aes(x = Ks, color = chr, fill = chr)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  theme_bw() +
  labs(title = "Distribution Ks : chr6 vs chr26",
       subtitle = "",
       x = "Ks", y = "Density")

ggsave("ks_density_chr6_vs_chr26.png", p2, width = 8, height = 5, dpi = 300)

# 8. PLOT 3 : Toutes les densités par chromosome (facettes)

p3 <- ggplot(chr_ks, aes(x = Ks, fill = chr)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ chr, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Distribution Ks par chromosome",
       x = "Ks", y = "Densité")

ggsave("ks_density_facets.png", p3, width = 14, height = 10, dpi = 300)


# Stats by chromosome

stats_chr <- chr_ks %>%
  group_by(chr) %>%
  summarise(n_paires = n(),
            ks_median = median(Ks),
            ks_mean = mean(Ks),
            ks_sd = sd(Ks)) %>%
  arrange(chr)

print(stats_chr)
write.table(stats_chr, "stats_ks_par_chromosome.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#############
# ks gradient like histogram
#############


# Ks mean by gene
gene_ks <- bind_rows(
  kaks_annotated %>% select(gene = gene1, chr = chr1, Ks),
  kaks_annotated %>% select(gene = gene2, chr = chr2, Ks)
) %>%
  group_by(gene, chr) %>%
  summarise(Ks_mean = mean(Ks, na.rm = TRUE), .groups = "drop")

# add position and calculate the rank
gene_ks <- gene_ks %>%
  left_join(gene_chr %>% select(gene, start), by = "gene") %>%  
  filter(!is.na(start)) %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# mean slicing 
window_size <- 10
gene_ks_smooth <- gene_ks %>%
  arrange(chr, rank) %>%
  group_by(chr) %>%
  mutate(Ks_smooth = rollmean(Ks_mean, k = window_size, 
                              fill = NA, align = "center")) %>%
  ungroup() %>%
  filter(!is.na(Ks_smooth))

# order the chr (descending)
chr_order <- gene_ks_smooth %>%
  group_by(chr) %>% summarise(n = n()) %>%
  arrange(desc(n)) %>% pull(chr)

gene_ks_smooth$chr <- factor(gene_ks_smooth$chr, levels = rev(chr_order))

# Order chr by ks ascending 
chr_order <- gene_ks_smooth %>%
  group_by(chr) %>%
  summarise(Ks_median = median(Ks_smooth, na.rm = TRUE)) %>%
  arrange(Ks_median) %>%
  pull(chr)

# Ks lower is up in plot
gene_ks_smooth$chr <- factor(gene_ks_smooth$chr, levels = rev(chr_order))

# Plot
p_gradient <- ggplot(gene_ks_smooth, aes(x = rank, y = chr, fill = Ks_smooth)) +
  geom_tile(height = 0.9) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "log10",
    name = "Mean Ks",
    direction = -1
  ) +
  labs(
    title = "Gradients of Ks along chromosomes",
    subtitle = paste0("windows = ", window_size, " genes"),
    x = "Rank",
    y = "Chromosome"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p_gradient)
ggsave("gradients_Ks_chromosomes.png", p_gradient, width = 14, height = 8, dpi = 300)


#####################################################################
# this below its for other script to plot ks on circos or other plot
#####################################################################

# load file 
win_cas_1 <- readRDS("../cas_Circos/win_cas_1.rds")
win_cas_2 <- readRDS("../cas_Circos/win_cas_2.rds")
win_cas_3 <- readRDS("../cas_Circos/win_cas_3.rds")
win_cas_4 <- readRDS("../cas_Circos/win_cas_4.rds")
win_cas_5 <- readRDS("../cas_Circos/win_cas_5.rds")

all_windows <- bind_rows(
  win_cas_1 %>% select(chr, win_start, win_end) %>% mutate(cas = "AORe"),
  win_cas_2 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe1"),
  win_cas_3 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe2"),
  win_cas_4 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe3"),
  win_cas_5 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe4")
) %>% mutate(chr = as.character(chr),
             mid = (win_start + win_end) / 2)

# ks by gene with position 
gene_ks_pos <- bind_rows(
  kaks_annotated %>% select(gene = gene1, chr = chr1, Ks),
  kaks_annotated %>% select(gene = gene2, chr = chr2, Ks)
) %>%
  filter(!is.na(Ks)) %>%
  left_join(gene_chr %>% select(gene, start), by = "gene") %>%
  filter(!is.na(start)) %>%
  mutate(chr = as.character(chr))

# assign each gene to windows
gene_in_win <- gene_ks_pos %>%
  inner_join(all_windows, by = "chr", relationship = "many-to-many") %>%
  filter(start >= win_start, start <= win_end)

# ks mean by windows
ks_per_win <- gene_in_win %>%
  group_by(chr, win_start, win_end, mid, cas) %>%
  summarise(Ks_mean   = mean(Ks, na.rm = TRUE),
            Ks_median = median(Ks, na.rm = TRUE),
            n_genes   = n(),
            se        = sd(Ks, na.rm = TRUE) / sqrt(n()),
            .groups   = "drop") %>%
  filter(n_genes >= 3)

# PLOT each chromosomes in one plot 

ks_per_win$chr <- factor(ks_per_win$chr,
                         levels = unique(ks_per_win$chr)[order(as.numeric(
                           gsub("[^0-9]", "", unique(ks_per_win$chr))))])

p_all <- ggplot(ks_per_win, aes(x = mid / 1e6, y = Ks_mean)) +
  geom_line(color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = cas), size = 1.2, alpha = 0.8) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE,
              color = "black", linewidth = 0.5) +
  scale_color_manual(values = col_cas, name = "Catégorie") +
  facet_wrap(~ chr, scales = "free") +   # <-- X et Y libres par facet
  labs(title = "Profil Ks par chromosome",
       subtitle = "Bins de 2 Mb",
       x = "Position (Mb)", y = "Ks moyen") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")


print(p_all)
ggsave("ks_curves_all_chr_windows.pdf", p_all,
       width = 16, height = 10, dpi = 300)

# save table 
###### very important ###########
write.table(ks_per_win, "ks_per_window.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# important also 
saveRDS(ks_per_win, "ks_per_win.rds")

# Optional
saveRDS(gene_ks_smooth, "gene_ks_smooth.rds")

