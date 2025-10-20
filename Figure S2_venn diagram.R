# =========================
# GWAS SNP Overlap Plots
# - Venn per trait (≤5 models)
# - UpSet per trait (>5 models)
# - Summary tables of totals & overlaps
# =========================

# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)   # Publication-quality Venns
  library(ComplexUpset)    # UpSet plots when many sets
  library(ggplot2)
})

# ---- Inputs ----
infile  <- "methods comparision_2.csv"   # path to your file
out_dir <- "gwas_overlap_outputs"        # folder for outputs
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read + normalize column names ----
raw <- readr::read_csv(infile, show_col_types = FALSE)

# Try to standardize to Trait, Model, SNP
df <- raw %>%
  janitor::clean_names() %>% 
  rename(
    Trait = tidyselect::any_of(c("trait", "phenotype", "trait_name")),
    Model = tidyselect::any_of(c("model", "method", "gwas_model", "algorithm")),
    SNP   = tidyselect::any_of(c("snp", "marker", "snp_id", "marker_id"))
  ) %>%
  mutate(
    Trait = as.character(Trait),
    Model = as.character(Model),
    SNP   = as.character(SNP)
  ) %>%
  filter(!is.na(Trait), !is.na(Model), !is.na(SNP)) %>%
  mutate(
    SNP = str_trim(SNP),
    Trait = str_trim(Trait),
    Model = str_trim(Model)
  ) %>%
  distinct(Trait, Model, SNP) # ensure unique triplets

# ---- Summary: totals per trait × model ----
totals_tbl <- df %>%
  count(Trait, Model, name = "Total_SNPs") %>%
  arrange(Trait, desc(Total_SNPs))

readr::write_csv(totals_tbl, file.path(out_dir, "totals_per_trait_model.csv"))

# ---- Helper: all intersections per trait (pairwise, k-wise) ----
# For reporting overlaps "between and across methods"
compute_intersections <- function(dfi) {
  # dfi is filtered to one trait
  sets <- dfi %>% split(.$Model) %>% lapply(\(x) unique(x$SNP))
  models <- names(sets)
  # Pairwise overlaps
  pairs <- combn(models, 2, simplify = FALSE)
  pw <- purrr::map_dfr(pairs, function(mm) {
    tibble(
      Trait = dfi$Trait[1],
      Models = paste(mm, collapse = " ∩ "),
      k = 2L,
      Overlap = length(intersect(sets[[mm[1]]], sets[[mm[2]]]))
    )
  })
  # k-wise overlaps for k >= 3 up to all models
  higher <- list()
  if (length(models) >= 3) {
    for (k in 3:length(models)) {
      combs <- combn(models, k, simplify = FALSE)
      tmp <- purrr::map_dfr(combs, function(mm) {
        intr <- Reduce(intersect, sets[mm])
        tibble(
          Trait = dfi$Trait[1],
          Models = paste(mm, collapse = " ∩ "),
          k = k,
          Overlap = length(intr)
        )
      })
      higher[[as.character(k)]] <- tmp
    }
  }
  bind_rows(pw, bind_rows(higher))
}

# ---- Plotters ----
plot_venn <- function(dfi, title) {
  # dfi is one trait; up to 5 models recommended for Venn
  sets <- dfi %>%
    split(.$Model) %>%
    lapply(\(x) unique(x$SNP))
  
  p <- ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient(low = "#e6e6e6", high = "#9fbad6") +
    labs(title = title, subtitle = "Numbers indicate overlapping SNP counts") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  p
}

plot_upset <- function(dfi, title) {
  # Build presence/absence matrix for models per SNP
  wide <- dfi %>%
    distinct(SNP, Model) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = Model, values_from = value, values_fill = 0L)
  
  set_cols <- setdiff(colnames(wide), "SNP")
  ggplot(wide, aes(x = after_stat(intersection))) +
    upset(set_cols, min_size = 1, width_ratio = 0.65) +
    labs(title = title, subtitle = "Intersection sizes of SNP sets across models") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# ---- Iterate by trait: compute overlaps and plot ----
traits <- sort(unique(df$Trait))

all_overlap <- list()

for (tr in traits) {
  dft <- df %>% filter(Trait == tr)
  n_models <- n_distinct(dft$Model)
  
  # Intersections table
  inter_tbl <- compute_intersections(dft)
  all_overlap[[tr]] <- inter_tbl
  
  # Save intersection table per trait
  readr::write_csv(inter_tbl, file.path(out_dir, paste0("overlaps_", make.names(tr), ".csv")))
  
  # Plot choice
  title_txt <- paste0(tr)
  if (n_models <= 5) {
    p <- plot_venn(dft, title_txt)
    ggsave(
      filename = file.path(out_dir, paste0("venn_", make.names(tr), ".pdf")),
      plot = p, width = 7, height = 6, device = cairo_pdf, dpi = 300
    )
    ggsave(
      filename = file.path(out_dir, paste0("venn_", make.names(tr), ".png")),
      plot = p, width = 7, height = 6, dpi = 300
    )
  } else {
    p <- plot_upset(dft, title_txt)
    ggsave(
      filename = file.path(out_dir, paste0("upset_", make.names(tr), ".pdf")),
      plot = p, width = 9, height = 6, device = cairo_pdf, dpi = 300
    )
    ggsave(
      filename = file.path(out_dir, paste0("upset_", make.names(tr), ".png")),
      plot = p, width = 9, height = 6, dpi = 300
    )
  }
}

# ---- Combined overlap summary across all traits ----
overlap_tbl <- bind_rows(all_overlap) %>%
  arrange(Trait, k, desc(Overlap))
readr::write_csv(overlap_tbl, file.path(out_dir, "overlaps_all_traits.csv"))

# ---- Also export a wide table: total per trait×model + all-models common count ----
all_models_common <- df %>%
  group_by(Trait) %>%
  summarize(
    Models = list(unique(Model)),
    Common_across_all = {
      lst <- split(SNP, Model)
      if (length(lst) >= 2) length(Reduce(intersect, lst)) else NA_integer_
    },
    .groups = "drop"
  )

totals_wide <- totals_tbl %>%
  pivot_wider(names_from = Model, values_from = Total_SNPs, values_fill = 0) %>%
  left_join(all_models_common %>% select(Trait, Common_across_all), by = "Trait")

readr::write_csv(totals_wide, file.path(out_dir, "totals_wide_with_common.csv"))

# ---- Console messages ----
message("Done. Outputs written to: ", normalizePath(out_dir))
message("- totals_per_trait_model.csv")
message("- overlaps_all_traits.csv and overlaps_<Trait>.csv")
message("- venn_<Trait>.pdf/png (≤5 models) or upset_<Trait>.pdf/png (>5 models)")
