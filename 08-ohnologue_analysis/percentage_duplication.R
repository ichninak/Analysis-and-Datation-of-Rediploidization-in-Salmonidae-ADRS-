library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# ============================================
# load bed file 
# ============================================
bed <- read.table("genes.Oncorhynchus.mykiss1.bed", sep = "\t", header = FALSE,
                  col.names = c("chr", "start", "end", "gene_id"),
                  colClasses = c("character", "integer", "integer", "character"))

bed <- bed[!grepl("^NW_", bed$chr), ]
write.table(bed, "genes_filtered.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print(unique(bed$chr))

gene_to_chr  <- setNames(bed$chr, bed$gene_id)
mykiss_genes <- bed$gene_id

# ============================================
# read ancestor file 
# ============================================
lines <- readLines("ancGenes.Protacanthopterygii.list")

double_rows  <- list()
all_anc_myk  <- c()   # all genes mykiss see on ancestor file 

for (line in lines) {
  fields <- strsplit(line, "[ \t]+")[[1]]
  if (length(fields) < 2) next
  
  family    <- fields[1]
  all_genes <- fields[-1]
  
  myk <- all_genes[all_genes %in% mykiss_genes]
  if (length(myk) == 0) next
  
  all_anc_myk <- c(all_anc_myk, myk)
  
  # Ohnologs : family with >= 2 copies mykiss
  if (length(myk) >= 2) {
    double_rows[[length(double_rows) + 1]] <- data.frame(
      family       = family,
      mykiss_genes = paste(myk, collapse = ","),
      chrs         = paste(unname(gene_to_chr[myk]), collapse = ","),
      n_mykiss     = length(myk),
      type         = "double",
      stringsAsFactors = FALSE
    )
  }
}

df_double <- do.call(rbind, double_rows)
all_anc_myk <- unique(all_anc_myk)

# ohnologs gene 
ohno_genes <- unique(unlist(strsplit(df_double$mykiss_genes, ",")))

ohno__anc_genes <- unique(unlist(strsplit(df_double$mykiss_genes, ",")))
length(ohno__anc_genes)

# Genes "simple" = all genes mykiss not ohnologs
simple_genes <- setdiff(all_anc_myk, ohno_genes)

df_simple <- data.frame(
  family       = paste0("single_", seq_along(simple_genes)),
  mykiss_genes = simple_genes,
  chrs         = unname(gene_to_chr[simple_genes]),
  n_mykiss     = 1,
  type         = "simple",
  stringsAsFactors = FALSE
)

df <- bind_rows(df_double, df_simple)

# ============================================
# table per chr 
# ============================================
per_gene <- df %>%
  separate_rows(mykiss_genes, chrs, sep = ",") %>%
  rename(gene_id = mykiss_genes, chr = chrs) %>%
  select(gene_id, chr, type)

table_chr <- per_gene %>%
  group_by(chr, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  mutate(
    double = if ("double" %in% names(.)) double else 0,
    simple = if ("simple" %in% names(.)) simple else 0,
    total  = double + simple,
    ratio_ohnology = double / total
  ) %>%
  arrange(desc(ratio_ohnology))

print(table_chr)

# ============================================
# load windows
# ============================================
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
) %>% mutate(chr = as.character(chr))

print(table(all_windows$cas))

# ============================================
# status of ohnologs per genes 
# ============================================
gene_status <- per_gene %>%
  mutate(is_ohno = ifelse(type == "double", 1, 0)) %>%
  select(gene_id, is_ohno) %>%
  distinct()

bed_status <- bed %>%
  inner_join(gene_status, by = "gene_id")

# ============================================
# Ratio per windows
# ============================================
window_ratios <- bed_status %>%
  inner_join(all_windows, by = "chr", relationship = "many-to-many") %>%
  filter(start >= win_start, start <= win_end) %>%
  group_by(cas, chr, win_start, win_end) %>%
  summarise(
    n_genes    = n(),
    n_ohno     = sum(is_ohno),
    n_simple   = sum(is_ohno == 0),
    ohno_ratio = mean(is_ohno),
    simple_ratio = mean(is_ohno == 0),
    .groups    = "drop"
  ) %>%
  filter(n_genes >= 2)

print(window_ratios %>% count(cas))

# ============================================
# Boxplots
# ============================================
col_test <- c(AORe = "grey70", LORe1 = "#e41a1c", LORe2 = "#377eb8",
              LORe3 = "#4daf4a", LORe4 = "#984ea3")

p <- ggplot(window_ratios, aes(x = cas, y = ohno_ratio, fill = cas)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey30") +
  scale_fill_manual(values = col_test, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(title = "Proportion of ohnologs by rediploidization categories",
       subtitle = "Windows of 2 Mb", x = NULL, y = "Proportion of ohnologs") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p)
ggsave("ohnology_ratio_per_cas.pdf", p, width = 12, height = 9, dpi = 300)

p_ohno <- ggplot(window_ratios, aes(x = cas, y = n_ohno, fill = cas)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey30") +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Number of ohnologs by rediploidization category",
       subtitle = "Window of 2 Mb", x = NULL, y = "Number of ohnologs") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_ohno)
ggsave("ohnology_count_per_cas.pdf", p_ohno, width = 12, height = 9, dpi = 300)

p_total <- ggplot(window_ratios, aes(x = cas, y = n_genes, fill = cas)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey30") +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Total number of genes per window by rediploidization category",
       subtitle = "Window of 2 Mb", x = NULL, y = "Number of genes") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_total)
ggsave("gene_count_per_cas.pdf", p_total, width = 12, height = 9, dpi = 300)

p_simple <- ggplot(window_ratios, aes(x = cas, y = n_simple, fill = cas)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey30") +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Total number of simple genes per window by rediploidization category",
       subtitle = "Window of 2 Mb", x = NULL, y = "Number of simple genes") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_simple)
ggsave("gene_count_per_cas.svg", p_total, width = 12, height = 9, dpi = 300)

p_s <- ggplot(window_ratios, aes(x = cas, y = simple_ratio, fill = cas)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey30") +
  scale_fill_manual(values = col_test, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(title = "Pourcentage of singletons by rediploidization categories",
       subtitle = "Windows of 2 Mb", x = NULL, y = "Proportion of singletons") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_s)

(p | p_ohno | p_total | p_simple)
ggsave("combined_boxplots.pdf", width = 48, height = 9, dpi = 300)

# ============================================
# total per categories 
# ============================================
totals_per_cas <- window_ratios %>%
  group_by(cas) %>%
  summarise(total_genes = sum(n_genes),
            total_ohno  = sum(n_ohno), .groups = "drop")
print(totals_per_cas)

p_tot_genes <- ggplot(totals_per_cas, aes(x = cas, y = total_genes, fill = cas)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = total_genes), vjust = -0.4, size = 4) +
  scale_fill_manual(values = col_test, guide = "none") +
  labs(title = "Number total of genes per categories",
       x = NULL, y = "Number total of genes") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_tot_genes)
ggsave("total_genes_per_cas.png", p_tot_genes, width = 8, height = 6, dpi = 300)

# ============================================
# statistics
# ============================================
kw <- kruskal.test(ohno_ratio ~ cas, data = window_ratios)
print(kw)

pairwise <- pairwise.wilcox.test(window_ratios$ohno_ratio,
                                 window_ratios$cas,
                                 p.adjust.method = "bonferroni", alternative = "greater")
print(pairwise)

pairwise <- pairwise.wilcox.test(window_ratios$n_ohno,
                                 window_ratios$cas,
                                 p.adjust.method = "bonferroni", alternative = "greater")
print(pairwise)
