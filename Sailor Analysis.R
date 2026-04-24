R

Example of SAILOR and variant calling analysis

April 2026
Load packages

library(dplyr)
library(tidyr)
library(ggplot2)
Set plot theme

pub_theme <- function() {
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 18, color = "black", face = "bold"), #panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
          legend.text = element_text(size = 20, color = "black"), 
          legend.title = element_blank(), strip.background = element_rect(fill = "white", size=1.8),
          strip.text = element_text(size = 25, color = "black", face = "bold"), legend.position = "none", 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2)
  )
}
SAILOR analysis

Running SAILOR on the wildtype and adr-2(-) replicates to identify edit sites and remove false positives. These edit sites will be used later in the comparison of editing frequency between wildtype and glh-1(ok439).

Import data

setwd("/geode2/home/u050/lhkelley/Quartz/Desktop/GSF4254_FLARE_results/")

# Import N2 replicates and merged adr-2(-) 
# Adding "genotype" and "replicate" columns
n2_rep1_sail = read.csv("GSF4254_N2_rep1_anno_sites.csv", sep = ",") %>% mutate(genotype = "N2") %>% mutate(rep = "rep1")
n2_rep2_sail = read.csv("GSF4254_N2_rep2_anno_sites.csv", sep = ",") %>% mutate(genotype = "N2") %>% mutate(rep = "rep2")
n2_rep3_sail = read.csv("GSF4254_N2_rep3_anno_sites.csv", sep = ",") %>% mutate(genotype = "N2") %>% mutate(rep = "rep3")
adr2_sail = read.csv("adr2.merged.FLAREannotated.sites.csv", sep = ",") %>% mutate(genotype = "adr2") %>% mutate(rep = "comb")
# How the data frame should look
head(n2_rep1_sail)
A data.frame: 6 × 12
X	chr	pos.1	pos	coverage	conf	strand	feature	wbID	geneID	genotype	rep
<int>	<chr>	<int>	<int>	<chr>	<dbl>	<chr>	<chr>	<chr>	<lgl>	<chr>	<chr>
1	0	X	17625436	17625437	108|T>C|0.009259259	0.341166062	-	CDS	WBGene00019211	NA	N2	rep1
2	1	I	   15004	   15005	775|A>G|0.001290323	0.000418447	+	CDS	WBGene00022276	NA	N2	rep1
3	2	X	17609974	17609975	18|T>C|0.055555556  	0.842943193	-	CDS	WBGene00015868	NA	N2	rep1
4	3	X	17470644	17470645	26|T>C|0.038461538  	0.777821359	-	CDS	WBGene00006602	NA	N2	rep1
5	4	I	   15083	   15084	1036|A>G|0.000965251	0.000030369	+	CDS	WBGene00022276	NA	N2	rep1
6	5	I	   15049	   15050	1117|A>G|0.000895255	0.000013455	+	CDS	WBGene00022276	NA	N2	rep1
# Combine the 4 tables into one
edit_sail = rbind(n2_rep1_sail, n2_rep2_sail, n2_rep3_sail, adr2_sail)

# Separate the poorly-formated "coverage" column and calculate
# Convert the "per_editing" and "coverage" columns to numeric
edit_sail = edit_sail %>% separate(coverage, c("coverage", "ref_alt", "per_editing"), sep = "\\|") %>% 
    mutate(per_editing = as.numeric(per_editing)) %>% mutate(coverage = as.numeric(coverage))
# Extract the genes and their feature type for later use
feature_gene = edit_sail %>% unite(col = "chr_pos", c("chr", "pos"), remove = FALSE) %>% 
    dplyr::select(chr_pos, feature, wbID)
Sanity checks

# Number of variants by genotype and replicate (if applicable)
edit_sail %>% group_by(genotype, rep) %>% tally()
A grouped_df: 4 × 3
genotype	rep	n
<chr>	<chr>	<int>
N2	rep1	80614
N2	rep2	71746
N2	rep3	159211
adr2	comb	96080
# Number of genes represented
edit_sail %>% 
    distinct(genotype, wbID, rep) %>% group_by(genotype, rep) %>% tally()
A grouped_df: 4 × 3
genotype	rep	n
<chr>	<chr>	<int>
N2	rep1	8338
N2	rep2	8277
N2	rep3	10065
adr2	comb	6855
Implementing confidence and replicate coverage cutoffs
For now, this code is taking sites with a confidence score of 0.75 or greater and sites must be in 2+ replicates.

Determining replicate coverage
# Create a fucntion to assign each site as being in 1, 2, or 3 replicates
venn_counts_N2 <- function(df, genotype_value = "N2", rep_col = "rep", chr_col = "chr", pos_col = "pos") {
  
  df2 <- df %>%
    filter(genotype == genotype_value) %>%
    mutate(
      chr = as.character(.data[[chr_col]]),
      pos = as.integer(.data[[pos_col]]),
      rep = as.character(.data[[rep_col]]),
      site = paste0(chr, "_", pos)
    )
  
  # Use reframe to safely return one row per site
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
        rep_key == "1"     ~ "rep1 only",
        rep_key == "2"     ~ "rep2 only",
        rep_key == "3"     ~ "rep3 only",
        rep_key == "1_2"   ~ "rep1 & rep2",
        rep_key == "1_3"   ~ "rep1 & rep3",
        rep_key == "2_3"   ~ "rep2 & rep3",
        rep_key == "1_2_3" ~ "rep1 & rep2 & rep3",
        TRUE               ~ "other"
      )
    )
  
  return(venn_counts)
}
Run the function with the edit_sail dataframe

# Compute Venn counts for N2
venn_n2 <- venn_counts_N2(edit_sail)
# View results
venn_n2
A tibble: 7 × 3
rep_key	n	category
<chr>	<int>	<chr>
rep1	54732	other
rep1_rep2	3529	other
rep1_rep2_rep3	10534	other
rep1_rep3	11819	other
rep2	47890	other
rep2_rep3	9793	other
rep3	127065	other
From the edit_sail dataframe, we're pulling out the wildtype sites because we don't care about the adr-2(-) data anymore. Now we are working with the dataframe n2_sites.

# Currently, the edit dataframe has "chr" and "pos" as separate columns - in the long run, it's easier to combine those columsn into one column/variable that we'll call "chr_pos"
n2_sites <- edit_sail %>%
  filter(genotype == "N2") %>%
  mutate(chr_pos = paste(chr, pos, sep = "_"))
# Find sites present in 2 or 3 replicates of wildtype and pull those out into a vector
n2_sites_reps <- n2_sites %>%
  group_by(chr_pos) %>%
  summarise(reps_present = list(unique(rep)), n_reps = n_distinct(rep), .groups = "drop") %>%
  filter(n_reps == 3 | n_reps == 2) %>% distinct(chr_pos) %>%
  pull(chr_pos)
# Number of sites here should equal the sum of the 2 and 3 rep groups in the table above
length(n2_sites_reps)
35675
Now pull out these sites that are in 2 or 3 replicates and implement the 0.75 confidence cutoff
# How many sites are left with these cutoffs?
n2_sites %>% filter(chr_pos %in% n2_sites_reps) %>% filter(conf >= 0.75) %>% distinct(chr_pos) %>% nrow()
5317
# But, because of the confidence cutoff, sites that were originally in 2+ reps could now be removed
# So, we want to remove sites now only in 1 rep
# Create a vector of the list of sites that are actually conf >= 0.75 and in 2+ reps
n2_sites_list = n2_sites %>% filter(chr_pos %in% n2_sites_reps) %>% filter(conf >= 0.75) %>% 
    group_by(chr_pos) %>% tally() %>% filter(n != 1) %>% pull(chr_pos)
Lastly, we need to remove any false positives from the adr-2(-) dataset
# Creating a list of any "sites" in adr-2(-)
adr2_sites_list = adr2_sail %>% filter(genotype == "adr2") %>%
  mutate(chr_pos = paste(chr, pos, sep = "_")) %>% distinct(chr_pos) %>% pull(chr_pos)
# Number of sites in N2 list
length(n2_sites_list)

# Number of sites after removing the false positives
setdiff(n2_sites_list, adr2_sites_list) %>% length()

# 89 false positives were removed

# Save list with removed false positives
n2_sites_list_final = setdiff(n2_sites_list, adr2_sites_list)
4742
4653
# Create the a table of final sites
# When we implement this criteria, how many sites are there?
n2_sites_final = n2_sites %>% filter(chr_pos %in% n2_sites_list_final) 
n2_sites_final %>% distinct(chr_pos) %>% nrow()
4653
# Create a version of the table with just the position
n2_sites_final_position = n2_sites_final %>% distinct(chr_pos) %>% dplyr::select(chr_pos)
Variant calling
Because these files are so large, it's better to use the fread command to read in the files. This is from the data.table package.

library(data.table)
library(pheatmap)
Set working directory and import data

setwd("/geode2/home/u050/lhkelley/Quartz/Desktop/GSF4254_VC_27Aug2025/")
# Renaming some of the columns and add a column denoting replicate
n2_rep1_vc <- fread("GSF4254-N2-rep1_S10_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "N2") %>% mutate(rep = "rep1") %>% dplyr::rename("num_variant" = `#variant_reads`)
n2_rep2_vc <- fread("GSF4254-N2-rep2_S11_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "N2") %>% mutate(rep = "rep2") %>% dplyr::rename("num_variant" = `#variant_reads`)
n2_rep3_vc <- fread("GSF4254-N2-rep3_S12_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "N2") %>% mutate(rep = "rep3") %>% dplyr::rename("num_variant" = `#variant_reads`)
glh1_rep1_vc <- fread("GSF4254-712-rep1_S1_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "glh1") %>% mutate(rep = "rep1") %>% dplyr::rename("num_variant" = `#variant_reads`)
glh1_rep2_vc <- fread("GSF4254-712-rep2_S2_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "glh1") %>% mutate(rep = "rep2") %>% dplyr::rename("num_variant" = `#variant_reads`)
glh1_rep3_vc <- fread("GSF4254-712-rep3_S3_R1_001_Aligned.sortedByCoord.out_27Aug_rmdup_variant.csv") %>% dplyr::rename("per_variant" = `%variant`) %>% mutate(geno = "glh1") %>% mutate(rep = "rep3") %>% dplyr::rename("num_variant" = `#variant_reads`)
# Example of what table looks like
head(n2_rep1_vc)
A data.table: 6 × 9
chr	pos	ref	alt	coverage	per_variant	num_variant	geno	rep
<chr>	<int>	<chr>	<chr>	<int>	<dbl>	<int>	<chr>	<chr>
I	31065	T	C	2	100.00	2	N2	rep1
I	71840	T	G	2	100.00	2	N2	rep1
I	71843	C	T	2	100.00	2	N2	rep1
I	71844	C	G	2	100.00	2	N2	rep1
I	111034	C	T	12	66.67	8	N2	rep1
I	111035	C	G	12	83.33	10	N2	rep1
# Combine the genotypes into one table
n2_vc = rbind(n2_rep1_vc, n2_rep2_vc, n2_rep3_vc) %>% mutate(chr_pos = paste(chr, pos, sep = "_")) %>% dplyr::select(-chr, -pos)
glh1_vc = rbind(glh1_rep1_vc, glh1_rep2_vc, glh1_rep3_vc) %>% mutate(chr_pos = paste(chr, pos, sep = "_")) %>% dplyr::select(-chr, -pos)
With the known edits sites found in the SAILOR section, pull those out of the variant calling data frames. Use the list version of the n2_sites_final.

# Pulling out the variant calling information for the 4,653 high confidence sites
# See how many distinct sites there are total (all 4,653 sites might not have variant data)
n2_vc %>% filter(chr_pos %in% n2_sites_list_final) %>% distinct(chr_pos) %>% nrow()

# Only 38% of our high confidence sites were in our variant data
1771
# But, each replicate may not have a particular site, so group by replicate and see how many sites fall into each group
n2_vc %>% filter(chr_pos %in% n2_sites_list_final) %>% group_by(rep) %>% tally
A tibble: 3 × 2
rep	n
<chr>	<int>
rep1	1298
rep2	1183
rep3	1334
# Create the table with the 1,803 sites
n2_vc_hc_sites = n2_vc %>% filter(chr_pos %in% n2_sites_list_final)
# Repeat for mutant (the site list stays the same)
glh1_vc %>% filter(chr_pos %in% n2_sites_list_final) %>% distinct(chr_pos) %>% nrow()
glh1_vc %>% filter(chr_pos %in% n2_sites_list_final) %>% group_by(rep) %>% tally
glh1_vc_hc_sites = glh1_vc %>% filter(chr_pos %in% n2_sites_list_final)
1759
A tibble: 3 × 2
rep	n
<chr>	<int>
rep1	1224
rep2	1317
rep3	1218
Do the variant sites found in N2 and the mutant overlap? We can only compare sites that appear in both genotypes.

c(overlap = nrow(inner_join(n2_vc_hc_sites %>% distinct(chr_pos), glh1_vc_hc_sites %>% distinct(chr_pos), by="chr_pos")),
    unique_n2 = nrow(anti_join(n2_vc_hc_sites %>% distinct(chr_pos), glh1_vc_hc_sites %>% distinct(chr_pos), by="chr_pos")),
    unique_adr2 = nrow(anti_join(glh1_vc_hc_sites %>% distinct(chr_pos), n2_vc_hc_sites %>% distinct(chr_pos), by="chr_pos")))

# 1,440 sites overlap between the two genotypes; extract those sites to use for analysis
overlap1440unique_n2331unique_adr2319
# Create a list of the overlapping sites
n2_glh1_sites = inner_join(n2_vc_hc_sites %>% distinct(chr_pos), glh1_vc_hc_sites %>% distinct(chr_pos), by="chr_pos") %>% pull(chr_pos)

# Confirm the number of overlapping sites matches above
length(n2_glh1_sites)
1440
# Combine N2 and mutant variant dataframe and pull out only the 1,440 sites
vc = rbind(n2_vc, glh1_vc) %>% filter(chr_pos %in% n2_glh1_sites)

# Confirm 1,440 sites are present
vc %>% distinct(chr_pos) %>% nrow()
1440
# But, we only want sites that have 3 replicates for both N2 and the mutant; so remove sites with only 1-2 reps
vc_3rep_list = vc %>% group_by(geno, chr_pos) %>% filter(n_distinct(rep) == 3) %>% ungroup() %>%
    group_by(chr_pos) %>% filter(n_distinct(geno) == 2) %>% ungroup() %>% distinct(chr_pos) %>% pull(chr_pos)
length(vc_3rep_list)
# 640 sites meet this requirement
640
# Pull these 640 sites out of the variant data
vc_n2_glh1_final = vc %>% filter(chr_pos %in% vc_3rep_list) 
vc_n2_glh1_final %>% distinct(chr_pos) %>% nrow()
# There are 640 sites that have 3 replicates in both N2 and the mutant
640
# Example of data frame to use in plotting
head(vc_n2_glh1_final)
A data.table: 6 × 8
ref	alt	coverage	per_variant	num_variant	geno	rep	chr_pos
<chr>	<chr>	<int>	<dbl>	<int>	<chr>	<chr>	<chr>
T	C	15	46.67	7	N2	rep1	I_318466
T	C	20	40.00	8	N2	rep1	I_629796
T	C	38	52.63	20	N2	rep1	I_629832
T	C	38	71.05	27	N2	rep1	I_629833
T	C	43	51.16	22	N2	rep1	I_629841
T	C	43	37.21	16	N2	rep1	I_629842
What is the editing frequency between wildtype and mutant by site?
plot_variant_heatmap <- function(vc) {
  
  # Make a combined column name like "N2_rep1"
  vc_wide <- vc_n2_glh1_final %>%
    mutate(sample = paste(geno, rep, sep = "_")) %>%
    dplyr::select(chr_pos, sample, per_variant) %>%
    pivot_wider(names_from = sample, values_from = per_variant)
  
  # Ensure desired column order
  desired_cols <- c(
    "N2_rep1", "glh1_rep1", "N2_rep2", "glh1_rep2", "N2_rep3", "glh1_rep3"
  )
  
  vc_wide <- vc_wide %>%
    dplyr::select(chr_pos, all_of(desired_cols))
  
  # Compute N2 mean for ordering
  vc_wide <- vc_wide %>%
    rowwise() %>%
    mutate(N2_mean = mean(c_across(starts_with("N2_")), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(N2_mean))
  
  # Convert to matrix
  mat <- vc_wide %>%
    dplyr::select(-chr_pos, -N2_mean) %>%
    as.data.frame()
  
  rownames(mat) <- vc_wide$chr_pos
  
  # Plot heatmap
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE
  )
}
options(repr.plot.width = 5, repr.plot.height = 8)
plot_variant_heatmap(vc_n2_glh1_final)
