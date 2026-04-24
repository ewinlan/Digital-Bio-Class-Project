
module load r
R
> library(dplyr)

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

> library(tidyr)
> library(ggplot2)
> 
> pub_theme <- function() {
+   theme_bw() +
+   theme(axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
+           axis.text.y = element_text(size = 15, color = "black"),
+           axis.title = element_text(size = 18, color = "black", face = "bold"), #panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
+           legend.text = element_text(size = 20, color = "black"), 
+           legend.title = element_blank(), strip.background = element_rect(fill = "white", size=1.8),
+           strip.text = element_text(size = 25, color = "black", face = "bold"), legend.position = "none", 
+           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
+         axis.ticks = element_line(colour = "black", linewidth = 2)
+   )
+ }


> n2OP50_rep1_sail = read.csv("WT_OP50_SRR23261284.annotated.sites.csv", sep = ",") %>% mutate(genotype = "N2OP50") %>% mutate(rep = "rep1")
> n2OP50_rep2_sail = read.csv("WT_OP50_SRR23261288.annotated.sites.csv", sep = ",") %>% mutate(genotype = "N2OP50") %>% mutate(rep = "rep2")
> n2PA14_rep1_sail = read.csv("WT_PA14_SRR23261295.annotated.sites.csv", sep = ",") %>% mutate(genotype = "N2PA14") %>% mutate(rep = "rep1")
> n2PA14_rep2_sail = read.csv("WT_PA14_SRR23261298.annotated.sites.csv", sep = ",") %>% mutate(genotype = "N2PA14") %>% mutate(rep = "rep2")
> adr2_sail = read.csv("adr2.annotated.sites.csv", sep = ",") %>% mutate(genotype = "adr2") %>% mutate(rep = "comb")


> # Combine the 4 tables into one
> edit_sail = rbind(n2OP50_rep1_sail, n2OP50_rep2_sail, n2PA14_rep1_sail, n2PA14_rep2_sail, adr2_sail)


> # Separate the poorly-formated "coverage" column and calculate
> # Convert the "per_editing" and "coverage" columns to numeric
> edit_sail = edit_sail %>% separate(coverage, c("coverage", "ref_alt", "per_editing"), sep = "\\|") %>% 
+     mutate(per_editing = as.numeric(per_editing)) %>% mutate(coverage = as.numeric(coverage))


> # Extract the genes and their feature type for later use
> feature_gene = edit_sail %>% unite(col = "chr_pos", c("chr", "pos"), remove = FALSE) %>% 


+     dplyr::select(chr_pos, feature, wbID)
> # Number of variants by genotype and replicate (if applicable)
> edit_sail %>% group_by(genotype, rep) %>% tally()
Warning: ‘timedatectl’ indicates the non-existent timezone name ‘n/a’
Warning: Your system is mis-configured: ‘/etc/localtime’ is not a symlink
Warning: It is strongly recommended to set envionment variable TZ to ‘America/New_York’ (or equivalent)
# A tibble: 5 × 3
# Groups:   genotype [3]
  genotype rep        n
  <chr>    <chr>  <int>
1 N2OP50   rep1   76339
2 N2OP50   rep2  118730
3 N2PA14   rep1   49634
4 N2PA14   rep2   54625
5 adr2     comb  153159


> # Number of genes represented
> edit_sail %>% 
+     distinct(genotype, wbID, rep) %>% group_by(genotype, rep) %>% tally()
# A tibble: 5 × 3
# Groups:   genotype [3]
  genotype rep       n
  <chr>    <chr> <int>
1 N2OP50   rep1  10356
2 N2OP50   rep2  11901
3 N2PA14   rep1   8513
4 N2PA14   rep2   8910
5 adr2     comb  12815


> # Create a fucntion to assign each site as being in 1, 2, or 3 replicates
> venn_counts_2reps <- function(df, genotype_value = "N2", rep_col = "rep", chr_col = "chr", pos_col = "pos") {
+   
+   df2 <- df %>%
+     filter(genotype == genotype_value) %>%
+     mutate(
+       chr = as.character(.data[[chr_col]]),
+       pos = as.integer(.data[[pos_col]]),
+       rep = as.character(.data[[rep_col]]),
+       site = paste0(chr, "_", pos)
+     )
+   
+   # Use reframe to safely return one row per site
+   site_reps <- df2 %>%
+     group_by(site) %>%
+     reframe(
+       rep_key = paste(sort(unique(rep)), collapse = "_")
+     )
+   
+   venn_counts <- site_reps %>%
+     count(rep_key, name = "n") %>%
+     arrange(rep_key) %>%
+     mutate(
+       category = case_when(
+         rep_key == "1"     ~ "rep1 only",
+         rep_key == "2"     ~ "rep2 only",
+         rep_key == "1_2"   ~ "rep1 & rep2",
+         TRUE               ~ "other"
+       )
+     )
+   
+   return(venn_counts)
+ }

> 
> # OP50 only
> op50_sites <- edit_sail %>%
+   filter(genotype == "N2OP50") %>%
+   mutate(chr_pos = paste(chr, pos, sep = "_"))
> 
> op50_sites_reps <- op50_sites %>%
+   group_by(chr_pos) %>%
+   summarise(
+     reps_present = list(unique(rep)),
+     n_reps = n_distinct(rep),
+     .groups = "drop"
+   ) %>%
+   filter(n_reps >= 2) %>%
+   distinct(chr_pos) %>%
+   pull(chr_pos)

> length(op50_sites_reps)
[1] 9043


> # PA14 only
> pa14_sites <- edit_sail %>%
+   filter(genotype == "N2PA14") %>%
+   mutate(chr_pos = paste(chr, pos, sep = "_"))
> 
> pa14_sites_reps <- pa14_sites %>%
+   group_by(chr_pos) %>%
+   summarise(
+     reps_present = list(unique(rep)),
+     n_reps = n_distinct(rep),
+     .groups = "drop"
+   ) %>%
+   filter(n_reps >= 2) %>%
+   distinct(chr_pos) %>%
+   pull(chr_pos)
> 
> length(pa14_sites_reps)
[1] 5386


> # OP50: pull out sites in 2 reps, apply conf >= 0.75,
> # then remove sites that only have 1 remaining replicate after filtering
> op50_sites %>%
+   filter(chr_pos %in% op50_sites_reps) %>%
+   filter(conf >= 0.75) %>%
+   distinct(chr_pos) %>%
+   nrow()
[1] 4070
> 
> op50_sites_list = op50_sites %>%
+   filter(chr_pos %in% op50_sites_reps) %>%
+   filter(conf >= 0.75) %>%
+   group_by(chr_pos) %>%
+   tally() %>%
+   filter(n != 1) %>%
+   pull(chr_pos)

> # PA14: same thing
> 
> pa14_sites %>%
+   filter(chr_pos %in% pa14_sites_reps) %>%
+   filter(conf >= 0.75) %>%
+   distinct(chr_pos) %>%
+   nrow()
[1] 2023
> 
> pa14_sites_list = pa14_sites %>%
+   filter(chr_pos %in% pa14_sites_reps) %>%
+   filter(conf >= 0.75) %>%
+   group_by(chr_pos) %>%
+   tally() %>%
+   filter(n != 1) %>%
+   pull(chr_pos)
> length(op50_sites_list)
[1] 3532
> length(pa14_sites_list)
[1] 1816


> # Create a list of any sites in adr2(-)
> adr2_sites_list = edit_sail %>%
+   filter(genotype == "adr2") %>%
+   mutate(chr_pos = paste(chr, pos, sep = "_")) %>%
+   distinct(chr_pos) %>%
+   pull(chr_pos)
al PA14 list with adr2 sites removed
pa14_sites_list_final = setdiff(pa14_sites_list, adr2_sites_list)

# Create final PA14 table
pa14_sites_final = pa14_sites %>%
  filter(chr_pos %in% pa14_sites_list_final)

pa14_sites_final %>%
  distinct(chr_pos) %>%
  nrow()

# Create PA14 table with just positions
pa14_sites_final_position = pa14_sites_final %>%
  distinct(chr_pos) %>%
  dplyr::select(chr_pos)> 
> # Number of OP50 sites before removing adr2 false positives
> length(op50_sites_list)
[1] 3532



> # Number of OP50 sites after removing adr2 false positives
> setdiff(op50_sites_list, adr2_sites_list) %>% length()
[1] 3446


> # Save final OP50 list with adr2 sites removed
> op50_sites_list_final = setdiff(op50_sites_list, adr2_sites_list)

> # Create final OP50 table
> op50_sites_final = op50_sites %>%
+   filter(chr_pos %in% op50_sites_list_final)

> op50_sites_final %>%
+   distinct(chr_pos) %>%
+   nrow()
[1] 3446

> # Create OP50 table with just positions
> op50_sites_final_position = op50_sites_final %>%
+   distinct(chr_pos) %>%
+   dplyr::select(chr_pos)


> # Number of PA14 sites before removing adr2 false positives
> length(pa14_sites_list)
[1] 1816

> # Number of PA14 sites after removing adr2 false positives
> setdiff(pa14_sites_list, adr2_sites_list) %>% length()
[1] 1764

> # Save final PA14 list with adr2 sites removed
> pa14_sites_list_final = setdiff(pa14_sites_list, adr2_sites_list)

> # Create final PA14 table
> pa14_sites_final = pa14_sites %>%
+   filter(chr_pos %in% pa14_sites_list_final)
> pa14_sites_final %>%
+   distinct(chr_pos) %>%
+   nrow()
[1] 1764

> # Create PA14 table with just positions
> pa14_sites_final_position = pa14_sites_final %>%
+   distinct(chr_pos) %>%
+   dplyr::select(chr_pos)
