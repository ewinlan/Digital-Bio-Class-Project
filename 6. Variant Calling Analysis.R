############################################################
# SAILOR RNA Editing Site Filtering Workflow
############################################################

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)

############################################################
# Publication theme for plots
############################################################

pub_theme <- function() {
  theme_bw() +
    theme(
      axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 18, color = "black", face = "bold"),
      legend.text = element_text(size = 20, color = "black"),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white", linewidth = 1.8),
      strip.text = element_text(size = 25, color = "black", face = "bold"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
      axis.ticks = element_line(colour = "black", linewidth = 2)
    )
}

############################################################
# Read in annotated SAILOR output files
############################################################

n2OP50_rep1_sail <- read.csv("WT_OP50_SRR23261284.annotated.sites.csv", sep = ",") %>%
  mutate(genotype = "N2OP50", rep = "rep1")

n2OP50_rep2_sail <- read.csv("WT_OP50_SRR23261288.annotated.sites.csv", sep = ",") %>%
  mutate(genotype = "N2OP50", rep = "rep2")

n2PA14_rep1_sail <- read.csv("WT_PA14_SRR23261295.annotated.sites.csv", sep = ",") %>%
  mutate(genotype = "N2PA14", rep = "rep1")

n2PA14_rep2_sail <- read.csv("WT_PA14_SRR23261298.annotated.sites.csv", sep = ",") %>%
  mutate(genotype = "N2PA14", rep = "rep2")

adr2_sail <- read.csv("adr2.annotated.sites.csv", sep = ",") %>%
  mutate(genotype = "adr2", rep = "comb")

############################################################
# Combine tables
############################################################

edit_sail <- bind_rows(
  n2OP50_rep1_sail,
  n2OP50_rep2_sail,
  n2PA14_rep1_sail,
  n2PA14_rep2_sail,
  adr2_sail
)

############################################################
# Separate coverage column and convert columns to numeric
############################################################

edit_sail <- edit_sail %>%
  separate(
    coverage,
    into = c("coverage", "ref_alt", "per_editing"),
    sep = "\\|",
    remove = TRUE
  ) %>%
  mutate(
    coverage = as.numeric(coverage),
    per_editing = as.numeric(per_editing)
  )

############################################################
# Extract gene and feature information
############################################################

feature_gene <- edit_sail %>%
  unite(col = "chr_pos", c("chr", "pos"), remove = FALSE) %>%
  dplyr::select(chr_pos, feature, wbID)

############################################################
# Number of variants by genotype and replicate
############################################################

edit_sail %>%
  group_by(genotype, rep) %>%
  tally()

# Expected output:
# genotype  rep        n
# N2OP50    rep1   76339
# N2OP50    rep2  118730
# N2PA14    rep1   49634
# N2PA14    rep2   54625
# adr2      comb  153159

############################################################
# Number of genes represented by genotype and replicate
############################################################

edit_sail %>%
  distinct(genotype, wbID, rep) %>%
  group_by(genotype, rep) %>%
  tally()

# Expected output:
# genotype  rep       n
# N2OP50    rep1  10356
# N2OP50    rep2  11901
# N2PA14    rep1   8513
# N2PA14    rep2   8910
# adr2      comb  12815

############################################################
# Function to count sites found in each replicate
############################################################

venn_counts_2reps <- function(df, genotype_value, rep_col = "rep", chr_col = "chr", pos_col = "pos") {
  
  df2 <- df %>%
    filter(genotype == genotype_value) %>%
    mutate(
      chr = as.character(.data[[chr_col]]),
      pos = as.integer(.data[[pos_col]]),
      rep = as.character(.data[[rep_col]]),
      site = paste0(chr, "_", pos)
    )
  
  site_reps <- df2 %>%
    group_by(site) %>%
    reframe(
      rep_key = paste(sort(unique(rep)), collapse = "_")
    )
  
  venn_counts <- site_reps %>%
    count(rep_key, name = "n") %>%
    arrange(rep_key) %>%
    mutate(
      category = case_when(
        rep_key == "rep1" ~ "rep1 only",
        rep_key == "rep2" ~ "rep2 only",
        rep_key == "rep1_rep2" ~ "rep1 & rep2",
        TRUE ~ "other"
      )
    )
  
  return(venn_counts)
}

# Optional checks:
venn_counts_2reps(edit_sail, genotype_value = "N2OP50")
venn_counts_2reps(edit_sail, genotype_value = "N2PA14")

############################################################
# OP50: identify sites found in both replicates
############################################################

op50_sites <- edit_sail %>%
  filter(genotype == "N2OP50") %>%
  mutate(chr_pos = paste(chr, pos, sep = "_"))

op50_sites_reps <- op50_sites %>%
  group_by(chr_pos) %>%
  summarise(
    reps_present = list(unique(rep)),
    n_reps = n_distinct(rep),
    .groups = "drop"
  ) %>%
  filter(n_reps >= 2) %>%
  distinct(chr_pos) %>%
  pull(chr_pos)

length(op50_sites_reps)

# Expected output:
# 9043

############################################################
# PA14: identify sites found in both replicates
############################################################

pa14_sites <- edit_sail %>%
  filter(genotype == "N2PA14") %>%
  mutate(chr_pos = paste(chr, pos, sep = "_"))

pa14_sites_reps <- pa14_sites %>%
  group_by(chr_pos) %>%
  summarise(
    reps_present = list(unique(rep)),
    n_reps = n_distinct(rep),
    .groups = "drop"
  ) %>%
  filter(n_reps >= 2) %>%
  distinct(chr_pos) %>%
  pull(chr_pos)

length(pa14_sites_reps)

# Expected output:
# 5386

############################################################
# OP50: keep sites in both replicates with confidence >= 0.75
############################################################

op50_sites %>%
  filter(chr_pos %in% op50_sites_reps) %>%
  filter(conf >= 0.75) %>%
  distinct(chr_pos) %>%
  nrow()

# Expected output before removing sites with only one remaining replicate:
# 4070

op50_sites_list <- op50_sites %>%
  filter(chr_pos %in% op50_sites_reps) %>%
  filter(conf >= 0.75) %>%
  group_by(chr_pos) %>%
  tally() %>%
  filter(n != 1) %>%
  pull(chr_pos)

length(op50_sites_list)

# Expected output:
# 3532

############################################################
# PA14: keep sites in both replicates with confidence >= 0.75
############################################################

pa14_sites %>%
  filter(chr_pos %in% pa14_sites_reps) %>%
  filter(conf >= 0.75) %>%
  distinct(chr_pos) %>%
  nrow()

# Expected output before removing sites with only one remaining replicate:
# 2023

pa14_sites_list <- pa14_sites %>%
  filter(chr_pos %in% pa14_sites_reps) %>%
  filter(conf >= 0.75) %>%
  group_by(chr_pos) %>%
  tally() %>%
  filter(n != 1) %>%
  pull(chr_pos)

length(pa14_sites_list)

# Expected output:
# 1816

############################################################
# Create list of sites detected in adr-2(-)
############################################################

adr2_sites_list <- edit_sail %>%
  filter(genotype == "adr2") %>%
  mutate(chr_pos = paste(chr, pos, sep = "_")) %>%
  distinct(chr_pos) %>%
  pull(chr_pos)

############################################################
# OP50: remove adr-2(-) sites as likely false positives
############################################################

length(op50_sites_list)

# Expected output before adr-2 filtering:
# 3532

setdiff(op50_sites_list, adr2_sites_list) %>%
  length()

# Expected output after adr-2 filtering:
# 3446

op50_sites_list_final <- setdiff(op50_sites_list, adr2_sites_list)

op50_sites_final <- op50_sites %>%
  filter(chr_pos %in% op50_sites_list_final)

op50_sites_final %>%
  distinct(chr_pos) %>%
  nrow()

# Expected output:
# 3446

op50_sites_final_position <- op50_sites_final %>%
  distinct(chr_pos) %>%
  dplyr::select(chr_pos)

############################################################
# PA14: remove adr-2(-) sites as likely false positives
############################################################

length(pa14_sites_list)

# Expected output before adr-2 filtering:
# 1816

setdiff(pa14_sites_list, adr2_sites_list) %>%
  length()

# Expected output after adr-2 filtering:
# 1764

pa14_sites_list_final <- setdiff(pa14_sites_list, adr2_sites_list)

pa14_sites_final <- pa14_sites %>%
  filter(chr_pos %in% pa14_sites_list_final)

pa14_sites_final %>%
  distinct(chr_pos) %>%
  nrow()

# Expected output:
# 1764

pa14_sites_final_position <- pa14_sites_final %>%
  distinct(chr_pos) %>%
  dplyr::select(chr_pos)

############################################################
# Optional: save final tables
############################################################

write.csv(op50_sites_final, "op50_sites_final.csv", row.names = FALSE)
write.csv(pa14_sites_final, "pa14_sites_final.csv", row.names = FALSE)

write.csv(op50_sites_final_position, "op50_sites_final_positions.csv", row.names = FALSE)
write.csv(pa14_sites_final_position, "pa14_sites_final_positions.csv", row.names = FALSE)

write.csv(feature_gene, "feature_gene_table.csv", row.names = FALSE)
