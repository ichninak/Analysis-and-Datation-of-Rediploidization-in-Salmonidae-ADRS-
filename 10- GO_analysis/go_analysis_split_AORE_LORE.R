################
# Load libraries
################

library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

# ===========
# Load files 
# ===========

AORE <- read.delim(
  "DAVIDFunctAnnotClusterReport_zebra_A_2026-06-24.csv",
  header = TRUE, sep = ","
)

LORE <- read.delim(
  "DAVIDFunctAnnotClusterReport_zebra_L_2026-06-24.csv",
  header = TRUE, sep = ","
)

AORE$Group <- "AORE"
LORE$Group <- "LORE"

david <- rbind(AORE, LORE)

# Table GO

go_map <- AnnotationDbi::select(
  GO.db,
  keys = keys(GO.db, keytype = "GOID"),
  columns = c("GOID", "TERM", "ONTOLOGY"),
  keytype = "GOID"
)

# ADD GOID / ONTOLOGY

david_go <- david %>%
  left_join(go_map, by = c("Term" = "TERM"))

# Prepare the data

prepare_go <- function(data, ontology_type) {
  data %>%
    filter(ONTOLOGY == ontology_type) %>%
    filter(P.Value < 0.05) %>%
    mutate(
      logP = -log10(P.Value),
      Term_label = paste0(GOID, " - ", Term)
    )
}

# Function for plot 
# order by p-value

plot_go_group <- function(go_data, group_name, ontology_name, y_label, top_n = NULL) {
  
  df <- go_data %>%
    filter(Group == group_name) %>%
    # keep one line by GO 
    group_by(Term_label) %>%
    slice_max(order_by = logP, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(logP))
  
  if (!is.null(top_n)) {
    df <- df %>% slice_head(n = top_n)
  }
  
  # order axis y 
  df <- df %>%
    mutate(Term_label = factor(Term_label, levels = rev(unique(Term_label))))
  
  ggplot(df, aes(x = logP, y = Term_label)) +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_gradient(low = "skyblue", high = "red") +
    theme_bw() +
    labs(
      title = paste0("GO enrichment ", ontology_name, " - ", group_name),
      subtitle = "Terms ordered by p-value",
      x = "-log10(PValue)",
      y = y_label,
      size = "Gene count",
      color = "-log10(PValue)"
    ) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold")
    )
}


# Prepare BP, CC, MF

go_bp <- prepare_go(david_go, "BP")
go_cc <- prepare_go(david_go, "CC")
go_mf <- prepare_go(david_go, "MF")

# plots
# split AORE LORE 

## BP 
plot_bp_AORE <- plot_go_group(go_bp, "AORE", "BP", "GO biological process")
plot_bp_LORE <- plot_go_group(go_bp, "LORE", "BP", "GO biological process")

## CC 
plot_cc_AORE <- plot_go_group(go_cc, "AORE", "CC", "GO cellular component")
plot_cc_LORE <- plot_go_group(go_cc, "LORE", "CC", "GO cellular component")

## MF 
plot_mf_AORE <- plot_go_group(go_mf, "AORE", "MF", "GO molecular function")
plot_mf_LORE <- plot_go_group(go_mf, "LORE", "MF", "GO molecular function")


# Display plot and save plot 

# BP
plot_bp_AORE
ggsave("Go_enr_BP_AORE.pdf", plot_bp_AORE, width = 10, height = 8, dpi = 300)
plot_bp_LORE
ggsave("Go_enr_BP_LORE.pdf", plot_bp_LORE, width = 10, height = 8, dpi = 300)

# CC
plot_cc_AORE
ggsave("Go_enr_CC_AORE.pdf", plot_cc_AORE, width = 10, height = 8, dpi = 300)
plot_cc_LORE
ggsave("Go_enr_CC_LORE.pdf", plot_cc_LORE, width = 10, height = 8, dpi = 300)

# MF
plot_mf_AORE
ggsave("Go_enr_MF_AORE.pdf", plot_mf_AORE, width = 10, height = 8, dpi = 300)
plot_mf_LORE
ggsave("Go_enr_MF_LORE.pdf", plot_mf_LORE, width = 10, height = 8, dpi = 300)
