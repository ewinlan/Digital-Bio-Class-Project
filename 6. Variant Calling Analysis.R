cd /Class_Project/sailor/outputs
module load r
R

############################################################
# Variant calling
############################################################

# Because these files are large, use fread from data.table
library(data.table)
library(pheatmap)

# Set working directory to wherever your variant CSV files are
# setwd("/N/slate/ewinlan/Class_Project/sailor")

# Renaming some of the columns and add a column denoting replicate

wtop50_rep1_vc <- fread("WTOP50_rep1_rmdup_variant.csv") %>%
  dplyr::rename("per_variant" = `%variant`) %>%
  mutate(geno = "WTOP50") %>%
  mutate(rep = "rep1") %>%
  dplyr::rename("num_variant" = `#variant_reads`)

wtop50_rep2_vc <- fread("WTOP50_rep2_rmdup_variant.csv") %>%
  dplyr::rename("per_variant" = `%variant`) %>%
  mutate(geno = "WTOP50") %>%
  mutate(rep = "rep2") %>%
  dplyr::rename("num_variant" = `#variant_reads`)

wtpa14_rep1_vc <- fread("WTPA14_rep1_rmdup_variant.csv") %>%
  dplyr::rename("per_variant" = `%variant`) %>%
  mutate(geno = "WTPA14") %>%
  mutate(rep = "rep1") %>%
  dplyr::rename("num_variant" = `#variant_reads`)

wtpa14_rep2_vc <- fread("WTPA14_rep2_rmdup_variant.csv") %>%
  dplyr::rename("per_variant" = `%variant`) %>%
  mutate(geno = "WTPA14") %>%
  mutate(rep = "rep2") %>%
  dplyr::rename("num_variant" = `#variant_reads`)

# Example of what table looks like
head(wtop50_rep1_vc)

# Combine the conditions into one table
wtop50_vc = rbind(wtop50_rep1_vc, wtop50_rep2_vc) %>%
  mutate(chr_pos = paste(chr, pos, sep = "_")) %>%
  dplyr::select(-chr, -pos)

wtpa14_vc = rbind(wtpa14_rep1_vc, wtpa14_rep2_vc) %>%
  mutate(chr_pos = paste(chr, pos, sep = "_")) %>%
  dplyr::select(-chr, -pos)

# Pulling out the variant calling information for the high confidence sites
# Use your final SAILOR site list here

wtop50_vc %>%
  filter(chr_pos %in% n2_sites_final) %>%
  distinct(chr_pos) %>%
  nrow()

# But, each replicate may not have a particular site, so group by replicate

wtop50_vc %>%
  filter(chr_pos %in% n2_sites_list_final) %>%
  group_by(rep) %>%
  tally()

# Create the table with only high confidence sites

wtop50_vc_hc_sites = wtop50_vc %>%
  filter(chr_pos %in% n2_sites_final)

# Repeat for PA14; the site list stays the same

wtpa14_vc %>%
  filter(chr_pos %in% n2_sites_final) %>%
  distinct(chr_pos) %>%
  nrow()

wtpa14_vc %>%
  filter(chr_pos %in% n2_sites_final) %>%
  group_by(rep) %>%
  tally()

wtpa14_vc_hc_sites = wtpa14_vc %>%
  filter(chr_pos %in% n2_sites_final)

# Do the variant sites found in OP50 and PA14 overlap?
# We can only compare sites that appear in both conditions.

c(overlap =
    nrow(inner_join(wtop50_vc_hc_sites %>% distinct(chr_pos),
                    wtpa14_vc_hc_sites %>% distinct(chr_pos),
                    by = "chr_pos")),
  unique_wtop50 =
    nrow(anti_join(wtop50_vc_hc_sites %>% distinct(chr_pos),
                   wtpa14_vc_hc_sites %>% distinct(chr_pos),
                   by = "chr_pos")),
  unique_wtpa14 =
    nrow(anti_join(wtpa14_vc_hc_sites %>% distinct(chr_pos),
                   wtop50_vc_hc_sites %>% distinct(chr_pos),
                   by = "chr_pos")))

# Create a list of the overlapping sites

wtop50_wtpa14_sites = inner_join(wtop50_vc_hc_sites %>% distinct(chr_pos),
                                 wtpa14_vc_hc_sites %>% distinct(chr_pos),
                                 by = "chr_pos") %>%
  pull(chr_pos)

# Confirm the number of overlapping sites matches above
length(wtop50_wtpa14_sites)

# Combine OP50 and PA14 variant dataframe and pull out only overlapping sites

vc = rbind(wtop50_vc, wtpa14_vc) %>%
  filter(chr_pos %in% wtop50_wtpa14_sites)

# Confirm overlapping sites are present
vc %>%
  distinct(chr_pos) %>%
  nrow()

# We only want sites that have both replicates for both OP50 and PA14

vc_2rep_list = vc %>%
  group_by(geno, chr_pos) %>%
  filter(n_distinct(rep) == 2) %>%
  ungroup() %>%
  group_by(chr_pos) %>%
  filter(n_distinct(geno) == 2) %>%
  ungroup() %>%
  distinct(chr_pos) %>%
  pull(chr_pos)

length(vc_2rep_list)

# Pull these sites out of the variant data

vc_wtop50_wtpa14_final = vc %>%
  filter(chr_pos %in% vc_2rep_list)

vc_wtop50_wtpa14_final %>%
  distinct(chr_pos) %>%
  nrow()

# Example of data frame to use in plotting

head(vc_wtop50_wtpa14_final)

plot_variant_heatmap <- function(vc) {
  
  vc_wide <- vc %>%
    mutate(sample = paste(geno, condition, rep, sep = "_")) %>%
    dplyr::select(chr_pos, sample, per_variant) %>%
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
    dplyr::select(-chr_pos, -OP50_mean) %>%
    as.data.frame()
  
  rownames(mat) <- vc_wide$chr_pos
  
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE
  )
}

options(repr.plot.width = 5, repr.plot.height = 8)

plot_variant_heatmap(vc_op50_pa14_final)

####### to identify sites that either increased or decreased by 10% #######

editing_10pct_changes <- vc_op50_pa14_final %>%
+   group_by(chr_pos, condition) %>%
+   summarize(mean_editing = mean(per_variant, na.rm = TRUE),
+             .groups = "drop") %>%
+   pivot_wider(names_from = condition,
+               values_from = mean_editing) %>%
+   mutate(
+     diff_PA14_minus_OP50 = PA14 - OP50,
+     change_direction = case_when(
+       diff_PA14_minus_OP50 >= 10 ~ "Increased by ≥10%",
+       diff_PA14_minus_OP50 <= -10 ~ "Decreased by ≥10%",
+       TRUE ~ "Changed <10%"
+     )
+   )
> editing_10pct_changes %>%
+   count(change_direction)
# A tibble: 3 × 2
  change_direction      n
  <chr>             <int>
1 Changed <10%        265
2 Decreased by ≥10%    47
3 Increased by ≥10%   123
> editing_10pct_changes %>%
+   filter(abs(diff_PA14_minus_OP50) >= 10) %>%
+   nrow()
[1] 170
> changed_10pct_sites <- editing_10pct_changes %>%
+   filter(abs(diff_PA14_minus_OP50) >= 10)
> 
> write.csv(changed_10pct_sites,
+           "sites_changed_by_10percent.csv",
+           row.names = FALSE)
> changed_10pct_sites <- editing_10pct_changes %>%
+   filter(abs(diff_PA14_minus_OP50) >= 10)
> 
> write.csv(changed_10pct_sites,
+           "sites_changed_by_10percent.csv",
+           row.names = FALSE)
> 
