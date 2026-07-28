library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# ============================================
# load file bed
# ============================================
bed <- read.table("genes.Oncorhynchus.mykiss1.bed", sep = "\t", header = FALSE,
                  col.names = c("chr", "start", "end", "gene_id"),
                  colClasses = c("character", "integer", "integer", "character"))
bed <- bed[!grepl("^NW_", bed$chr), ]

gene_to_chr   <- setNames(bed$chr,   bed$gene_id)
gene_to_start <- setNames(bed$start, bed$gene_id)
mykiss_genes  <- bed$gene_id

# ============================================
# read ancgenes: familly of double genes (≥2 mykiss)
# ============================================
# we need to keep each line of O.mykiss genes that is ohnologs so minimum 2 genes 

lines <- readLines("../08-ohnologue_analysis/ancGenes.Protacanthopterygii.list")

double_rows <- list()
for (line in lines) {
  fields <- strsplit(line, "[ \t]+")[[1]]
  if (length(fields) < 2) next
  family <- fields[1]
  myk    <- fields[-1][fields[-1] %in% mykiss_genes]
  if (length(myk) >= 2) {
    double_rows[[length(double_rows) + 1]] <- data.frame(
      family       = family,
      mykiss_genes = paste(myk, collapse = ","),
      n_mykiss     = length(myk),
      stringsAsFactors = FALSE
    )
  }
}
df_double <- do.call(rbind, double_rows)
cat("Familles doubles (≥2 mykiss) :", nrow(df_double), "\n")

# ============================================
# load Tau
# ============================================
tau_data <- read.csv("ohno_ts.csv", stringsAsFactors = FALSE)
tau_data$gene <- as.character(tau_data$gene)

# ============================================
# modify the table 
#    and after that calculate a tau-range by family 
# ============================================
fam_long <- df_double %>%
  separate_rows(mykiss_genes, sep = ",") %>%
  rename(gene_id = mykiss_genes)

# keep tau-range of each copy (for calculate the tau-range "family")
fam_with_tau <- fam_long %>%
  left_join(tau_data %>% select(gene, tau_range, tau_mean, tau_med, Tau),
            by = c("gene_id" = "gene"))

# keep one value for each family 
fam_tau <- fam_with_tau %>%
  group_by(family) %>%
  summarise(
    tau_range_fam = mean(tau_range, na.rm = TRUE),
    tau_mean_fam  = mean(tau_mean,  na.rm = TRUE),
    tau_med_fam   = mean(tau_med,   na.rm = TRUE),
    Tau_fam       = mean(Tau,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.nan(tau_range_fam))

cat("Familles avec tau_range défini :", nrow(fam_tau), "\n")

# re-assign for each copy the value tau-range of his family 
fam_long <- fam_long %>%
  inner_join(fam_tau, by = "family") %>%
  mutate(
    chr   = unname(gene_to_chr[gene_id]),
    start = unname(gene_to_start[gene_id])
  ) %>%
  filter(!is.na(chr)) %>%
  # rename 
  rename(tau_range = tau_range_fam,
         tau_mean  = tau_mean_fam,
         tau_med   = tau_med_fam,
         Tau       = Tau_fam)


# ============================================
# load windows
# ============================================
win_cas_1 <- readRDS("../03-Cas_circos/win_cas_1.rds")
win_cas_2 <- readRDS("../03-Cas_circos/win_cas_2.rds")
win_cas_3 <- readRDS("../03-Cas_circos/win_cas_3.rds")
win_cas_4 <- readRDS("../03-Cas_circos/win_cas_4.rds")
win_cas_5 <- readRDS("../03-Cas_circos/win_cas_5.rds")

all_windows <- bind_rows(
  win_cas_1 %>% select(chr, win_start, win_end) %>% mutate(cas = "AORe"),
  win_cas_2 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe1"),
  win_cas_3 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe2"),
  win_cas_4 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe3"),
  win_cas_5 %>% select(chr, win_start, win_end) %>% mutate(cas = "LORe4")
) %>% mutate(chr = as.character(chr))

# ============================================
# assign each copy to one windows 
# ============================================
genes_by_win <- fam_long %>%
  inner_join(all_windows, by = "chr", relationship = "many-to-many") %>%
  filter(start >= win_start, start <= win_end)

print(genes_by_win %>% count(cas))

# ============================================
# tau-range by windows
# ============================================
window_taurange <- genes_by_win %>%
  group_by(cas, chr, win_start, win_end) %>%
  summarise(
    n_genes        = n(),
    tau_range_mean = mean(tau_range, na.rm = TRUE),
    tau_range_med  = median(tau_range, na.rm = TRUE),
    tau_mean_mean  = mean(tau_mean,  na.rm = TRUE),
    Tau_mean       = mean(Tau,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 2)

print(window_taurange %>% count(cas))

# ============================================
# Boxplots : one point = one windows
# ============================================
col_test <- c(AORe = "grey70", LORe1 = "#e41a1c", LORe2 = "#377eb8",
              LORe3 = "#4daf4a", LORe4 = "#984ea3")

p_mean <- ggplot(window_taurange, aes(x = cas, y = tau_range_mean, fill = cas)) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4, color = "grey20") +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Mean tau_range of ancestral ohnologs by category",
       subtitle = "One point = one window",
       x = NULL, y = "Mean tau_range per window") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_mean)
ggsave("tau_range_mean_anc_per_cas.png", p_mean, width = 10, height = 8, dpi = 300)

p_med <- ggplot(window_taurange, aes(x = cas, y = tau_range_med, fill = cas)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "grey20") +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Median tau_range of ancestral ohnologs by category",
       subtitle = "One point = one window",
       x = NULL, y = "Median tau_range per window") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_med)
ggsave("tau_range_median_anc_per_cas.svg", p_med, width = 10, height = 8, dpi = 300)

(p_mean | p_med)
ggsave("tau_range_anc_combined.pdf", width = 16, height = 7, dpi = 300)

# ============================================
# statistics
# ============================================
# tau-range-mean
print(kruskal.test(tau_range_mean ~ cas, data = window_taurange))
print(pairwise.wilcox.test(window_taurange$tau_range_mean, window_taurange$cas,
                           p.adjust.method = "bonferroni", alternative = "less"))

#tau-range-med
print(kruskal.test(tau_range_med ~ cas, data = window_taurange))
print(pairwise.wilcox.test(window_taurange$tau_range_med, window_taurange$cas,
                           p.adjust.method = "bonferroni"))
