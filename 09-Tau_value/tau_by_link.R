library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================
# Load GeneID_par_cas.tsv
# ============================================
gene_cas <- read.table("GeneID_par_cas.tsv", sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "esoxID_bed", "cas",
                                     "chr", "start", "end"),
                       colClasses = c("character", "character", "character",
                                      "character", "integer", "integer"))

# Fusion NAME11+12 -> cas4, NAME13+14 -> cas5
gene_cas <- gene_cas %>%
  mutate(cas = case_when(
    cas %in% c("cas4", "cas5") ~ "cas4",
    cas %in% c("cas6", "cas7") ~ "cas5",
    TRUE                       ~ cas
  ))

cas_map <- c(cas1 = "AORe", cas2 = "LORe1", cas3 = "LORe2",
             cas4 = "LORe3", cas5 = "LORe4")
gene_cas <- gene_cas %>% mutate(cas = unname(cas_map[cas]))

print(table(gene_cas$cas, useNA = "ifany"))

# ============================================
# load merged.tsv pairs of ohnologs
# ============================================
merged <- read.table("../08-ohnologue_analysis/merged.tsv", header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE,
                     colClasses = c("character", "character", "character", "character"))
colnames(merged) <- c("esoxID", "gene1", "gene2", "flag")

# ============================================
# load Tau 
# ============================================
tau_data <- read.csv("ohno_ts.csv", stringsAsFactors = FALSE)

# ============================================
# join tau_range to each pairs
# ============================================
merged_tau <- merged %>%
  inner_join(tau_data %>% select(esoxID, tau_range, tau_mean, tau_med),
             by = "esoxID") %>%
  mutate(pair_id = paste0("pair_", row_number()))


pairs_cas <- merged_tau %>%
  left_join(gene_cas %>% select(gene1 = gene_id, cas1 = cas), by = "gene1") %>%
  left_join(gene_cas %>% select(gene2 = gene_id, cas2 = cas), by = "gene2") %>%
  mutate(cas = case_when(
    !is.na(cas1) & !is.na(cas2) & cas1 == cas2 ~ cas1,   
    !is.na(cas1) & is.na(cas2)                  ~ cas1,   
    is.na(cas1) & !is.na(cas2)                  ~ cas2,
    !is.na(cas1) & !is.na(cas2) & cas1 != cas2  ~ NA_character_, 
    TRUE                                         ~ NA_character_
  ))

pairs_cas <- pairs_cas %>% filter(!is.na(cas))
print(pairs_cas %>% count(cas))

# ============================================
# Boxplot : 1 point = 1 pairs of ohnologs
# ============================================
ord <- c("AORe", "LORe1", "LORe2", "LORe3", "LORe4")
pairs_cas$cas <- factor(pairs_cas$cas, levels = ord)

col_test <- c(AORe = "grey70", LORe1 = "#e41a1c", LORe2 = "#377eb8",
              LORe3 = "#4daf4a", LORe4 = "#984ea3")

p_pair <- ggplot(pairs_cas, aes(x = cas, y = tau_range, fill = cas)) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.15, color = "grey20") +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "tau_range of ohnolog pairs by rediploidization category",
       subtitle = "One point = one ohnolog pair",
       x = NULL, y = "tau_range") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_pair)
ggsave("tau_range_per_pair.svg", p_pair, width = 10, height = 8, dpi = 300)

# ============================================
# statistics
# ============================================
# tau-range by pairs
print(kruskal.test(tau_range ~ cas, data = pairs_cas))
print(pairwise.wilcox.test(pairs_cas$tau_range, pairs_cas$cas,
                           p.adjust.method = "bonferroni", alternative = "less"))
