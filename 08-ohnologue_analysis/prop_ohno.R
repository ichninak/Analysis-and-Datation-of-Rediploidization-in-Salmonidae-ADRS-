# script to have a mean of ratio of ohnologs per chr 

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ============================================
# load bed file 
# ============================================
bed <- read.table("genes.Oncorhynchus.mykiss1.bed", sep = "\t", header = FALSE,
                  col.names = c("chr", "start", "end", "gene_id"),
                  colClasses = c("character", "integer", "integer", "character"))

# remove nw scaffold
bed <- bed[!grepl("^NW_", bed$chr), ]

# (optionnal) 
write.table(bed, "genes_filtered.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print(unique(bed$chr))

# dico gene_id -> chr
gene_to_chr <- setNames(bed$chr, bed$gene_id)
mykiss_genes <- bed$gene_id

# ============================================
# read ancestor file 
# ============================================
lines <- readLines("ancGenome.Protacanthopterygii.list")

cars_to_keep <- paste0("CAR_", 1:62)

results <- list()
excluded_expansions <- 0

for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) < 5) next
  
  car_chr <- parts[1]
  if (!(car_chr %in% cars_to_keep)) next
  
  fields <- strsplit(parts[5], " ")[[1]]
  family <- fields[1]
  all_genes <- fields[-1]
  
  # keep only mykiss genes 
  myk <- all_genes[all_genes %in% mykiss_genes]
  
  # Filter
  if (length(myk) == 0) next               
  if (length(myk) > 2) {                  
    excluded_expansions <- excluded_expansions + 1
    next
  }
  
  gene_type <- if (length(myk) == 2) "double" else "simple"
  chrs <- unname(gene_to_chr[myk])
  
  results[[length(results) + 1]] <- data.frame(
    CAR = car_chr,
    family = family,
    mykiss_genes = paste(myk, collapse = ","),
    chrs = paste(chrs, collapse = ","),
    n_mykiss = length(myk),
    type = gene_type,
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, results)

# ============================================
# table per chr one gene = one line
# ============================================
per_gene <- df %>%
  separate_rows(chrs, sep = ",") %>%
  rename(chr = chrs) %>%
  select(chr, type)

table_chr <- per_gene %>%
  group_by(chr, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  mutate(
    double = if ("double" %in% names(.)) double else 0,
    simple = if ("simple" %in% names(.)) simple else 0,
    total = double + simple,
    ratio_ohnology = double / total
  ) %>%
  arrange(desc(ratio_ohnology))

print(table_chr)
write.table(table_chr, "ohnology_per_chr.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

# ============================================
# plot 
# ============================================
# order the chr descending with mean ohno
table_chr$chr <- factor(table_chr$chr,
                        levels = rev(table_chr$chr))  


p <- ggplot(table_chr, aes(x = ratio_ohnology, y = chr)) +
  geom_col(fill = "#7BA7C7", width = 0.85) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Mean ohnolog gene ratio per chromosome",
    x = "Ohnolog gene percentage",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 9)
  )

print(p)

ggsave("ohnology_ratio_per_chr.pdf", p, width = 10, height = 7, dpi = 200)

