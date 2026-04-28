############################################################
# Variant calling analysis summary
############################################################

library(dplyr)
library(tidyr)
library(pheatmap)

# Read filtered variant-calling RDS files
wtop50_vc <- bind_rows(
  readRDS("wtop50_rep1_filtered.rds"),
  readRDS("wtop50_rep2_filtered.rds")
)

wtpa14_vc <- bind_rows(
  readRDS("wtpa14_rep1_filtered.rds"),
  readRDS("wtpa14_rep2_filtered.rds")
)

# Check number of SAILOR sites detected in each condition
wtop50_vc %>% distinct(chr_pos) %>% nrow()
wtpa14_vc %>% distinct(chr_pos) %>% nrow()

# Check number of sites per replicate
wtop50_vc %>% group_by(rep) %>% tally()
wtpa14_vc %>% group_by(rep) %>% tally()

# Identify sites shared between OP50 and PA14
wtop50_wtpa14_sites <- inner_join(
  wtop50_vc %>% distinct(chr_pos),
  wtpa14_vc %>% distinct(chr_pos),
  by = "chr_pos"
) %>%
  pull(chr_pos)

length(wtop50_wtpa14_sites)

# Combine OP50 and PA14 data, keeping only overlapping sites
vc <- bind_rows(wtop50_vc, wtpa14_vc) %>%
  filter(chr_pos %in% wtop50_wtpa14_sites)

# Keep only sites with coverage >= 10,
# both replicates per condition,
# and both OP50 and PA14 represented
vc_2rep_list <- vc %>%
  filter(coverage >= 10) %>%
  group_by(condition, chr_pos) %>%
  filter(n_distinct(rep) == 2) %>%
  ungroup() %>%
  group_by(chr_pos) %>%
  filter(n_distinct(condition) == 2) %>%
  ungroup() %>%
  distinct(chr_pos) %>%
  pull(chr_pos)

# Final variant-calling dataset
vc_wtop50_wtpa14_final <- vc %>%
  filter(chr_pos %in% vc_2rep_list,
         coverage >= 10)

# Final number of sites used for downstream analysis
vc_wtop50_wtpa14_final %>%
  distinct(chr_pos) %>%
  nrow()

############################################################
# Heatmap of percent variant / editing
############################################################

plot_variant_heatmap <- function(vc) {
  
  vc_wide <- vc %>%
    mutate(sample = paste(geno, condition, rep, sep = "_")) %>%
    select(chr_pos, sample, per_variant) %>%
    group_by(chr_pos, sample) %>%
    summarize(per_variant = mean(per_variant, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = sample,
                values_from = per_variant)
  
  vc_wide <- vc_wide %>%
    rowwise() %>%
    mutate(OP50_mean = mean(c_across(contains("OP50")),
                            na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(OP50_mean))
  
  mat <- vc_wide %>%
    select(-chr_pos, -OP50_mean) %>%
    as.data.frame()
  
  rownames(mat) <- vc_wide$chr_pos
  
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE
  )
}

plot_variant_heatmap(vc_wtop50_wtpa14_final)

############################################################
# Identify sites with >=10% change between OP50 and PA14
############################################################

editing_10pct_changes <- vc_wtop50_wtpa14_final %>%
  group_by(chr_pos, condition) %>%
  summarize(mean_editing = mean(per_variant, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = condition,
              values_from = mean_editing) %>%
  mutate(
    diff_PA14_minus_OP50 = PA14 - OP50,
    change_direction = case_when(
      diff_PA14_minus_OP50 >= 10 ~ "Increased by ≥10%",
      diff_PA14_minus_OP50 <= -10 ~ "Decreased by ≥10%",
      TRUE ~ "Changed <10%"
    )
  )

# Count sites by change direction
editing_10pct_changes %>%
  count(change_direction)
