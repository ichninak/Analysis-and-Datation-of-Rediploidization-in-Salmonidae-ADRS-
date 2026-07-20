###################
# load libraries
###################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("GO.db")
BiocManager::install("AnnotationDbi")

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
  header = TRUE,
  sep = ","
)

LORE <- read.delim(
  "DAVIDFunctAnnotClusterReport_zebra_L_2026-06-24.csv",
  header = TRUE,
  sep = ","
)

AORE$Group <- "AORE"
LORE$Group <- "LORE"

david <- rbind(AORE, LORE)

#############
# PIPELINE
#############

# Table GO

go_map <- AnnotationDbi::select(
  GO.db,
  keys = keys(GO.db, keytype = "GOID"),
  columns = c("GOID", "TERM", "ONTOLOGY"),
  keytype = "GOID"
)

# Add GOID / ONTOLOGY

david_go <- david %>%
  left_join(go_map, by = c("Term" = "TERM"))

# Function to prepare the data 

prepare_go <- function(data, ontology_type) {
  
  go_data <- data %>%
    filter(ONTOLOGY == ontology_type) %>%
    filter(P.Value < 0.05) %>%
    mutate(
      logP = -log10(P.Value),
      Term_label = paste0(GOID, " - ", Term)
    )
  
  return(go_data)
}

# Function for order GO 

order_go <- function(go_data) {
  
  go_order <- go_data %>%
    select(Term_label, Group, logP) %>%
    group_by(Term_label, Group) %>%
    summarise(
      logP = max(logP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = logP,
      values_fill = list(logP = 0)
    ) %>%
    mutate(
      AORE = ifelse(is.na(AORE), 0, AORE),
      LORE = ifelse(is.na(LORE), 0, LORE),
      diff_logP = AORE - LORE,
      dominant_group = ifelse(diff_logP >= 0, "AORE", "LORE")
    ) %>%
    arrange(
      factor(dominant_group, levels = c("AORE", "LORE")),
      desc(abs(diff_logP))
    )
  
  return(go_order)
}

# Function for plot 

plot_go <- function(go_data, go_order, ontology_name, y_label) {
  
  go_plot <- go_data %>%
    left_join(
      go_order %>% select(Term_label, diff_logP, dominant_group),
      by = "Term_label"
    ) %>%
    mutate(
      Term_label = factor(Term_label, levels = rev(go_order$Term_label))
    )
  
  p <- ggplot(go_plot,
              aes(x = Group,
                  y = Term_label)) +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_gradient(low = "skyblue", high = "red") +
    theme_bw() +
    labs(
      title = paste0("GO enrichment comparison ", ontology_name),
      subtitle = "AORE enriched terms first, then LORE enriched terms",
      x = "",
      y = y_label,
      size = "Gene count",
      color = "-log10(PValue)"
    ) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# Prepare BP, CC, MF

go_bp <- prepare_go(david_go, "BP")
go_cc <- prepare_go(david_go, "CC")
go_mf <- prepare_go(david_go, "MF")


# Calculate the order

go_order_bp <- order_go(go_bp)
go_order_cc <- order_go(go_cc)
go_order_mf <- order_go(go_mf)


# plots

plot_bp <- plot_go(
  go_data = go_bp,
  go_order = go_order_bp,
  ontology_name = "BP",
  y_label = "GO biological process"
)

plot_cc <- plot_go(
  go_data = go_cc,
  go_order = go_order_cc,
  ontology_name = "CC",
  y_label = "GO cellular component"
)

plot_mf <- plot_go(
  go_data = go_mf,
  go_order = go_order_mf,
  ontology_name = "MF",
  y_label = "GO molecular function"
)


# display plots

plot_bp
ggsave("Go_enr_BP.pdf", plot_bp, width = 10, height = 8, dpi = 300)
plot_cc
ggsave("Go_enr_CC.png", plot_cc, width = 10, height = 8, dpi = 300)
plot_mf
ggsave("Go_enr_MF.png", plot_mf, width = 10, height = 8, dpi = 300)
