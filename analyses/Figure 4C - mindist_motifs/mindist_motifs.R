rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(glue)
library(scales)

library(cowplot)
library(magick)
library(grImport)

# library(broom)
# library(ggrepel)
# library(magrittr)
# library(naturalsort)
# library(progress)
# library(reshape2)

setwd("~/Documents/Projects/RBPome/mindist_motifs")

# rbp1 <- "FMR1"; rbp2 <- "FXR2" # FXR2 near FMR1
# rbp1 <- "FXR2"; rbp2 <- "FMR1" # FMR1 near FXR2

# q <- read.table(paste0("output-motif_ratios-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1-closest_is_distant.txt"), header=T, sep="\t", quote="")
# q <- read.table(paste0("output-motif_ratios-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1-solo.txt"), header=T, sep="\t", quote="")
q <- read.table(paste0("output-motif_ratios-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1.txt"), header=T, sep="\t", quote="")
q <- as_tibble(q)
q$pair <- paste0(q$symbol1,'|',q$symbol2)

# Add motif & motif source columns
q %<>% 
  separate(motifpair, c("motif1", "motif2"), sep=" & ", remove=F) %>%
  separate(motif1, c("source1", NA), sep="\\|", remove=F) %>%
  separate(motif2, c("source2", NA), sep="\\|", remove=F) %>%
  mutate(motif1=str_replace_all(motif1, "\\|", "_"), motif2=str_replace_all(motif2, "\\|", "_"))
q



# Get counts too (just for diagnostic reasons)
counts <- read.table(paste0("output-motif_counts-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1.txt"), header=T, sep="\t", quote="")
# counts <- read.table(paste0("tmp/output-motif_counts-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1.txt"), header=T, sep="\t", quote="")
counts <- as_tibble(counts)
counts$pair <- paste0(counts$symbol1,'|',counts$symbol2)
counts
counts %>% group_by(bindtype) %>% summarise(mean=mean(n), sd=sd(n), median=median(n), mad=mad(n))
# counts %>% group_by(bindtype, motifpair) %>% summarise(mean=mean(n), sd=sd(n), median=median(n), mad=mad(n), min=min(n), max=max(n)) %>% summarise(mean_of_means=mean(mean), median_of_medians=median(median), mean_of_mins=mean(min), mean_of_maxs=mean(max))
# >> Both types of control (closest_is_distant, i.e. "ratio", and solo, i.e. "ratio_solo") have very similar sample sizes behind their resamples, so it doesn't matter much which one I use.
# "closest_is_distant" means motif occurrences in genes where the other RBP binds, but is farther away than 54 nt.
# "solo" means motif occurrences in genes where the other RBP does not bind at all.
# >> Conceptually, "solo" seems like the better control.
# >> Go with "ratio_solo".


# Get RBP pairs of interest (the 30 from the "strict" Figure 3A network, i.e. the ones that have significant evidence of co-binding in proximity)
# Strictly filtered for cobinding (resampling p-value and resampling Wilcoxon p-value significant in one orientation at least)
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (21: the "strict" Cytoscape network, Figure 5a)
rbp_pairs <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
str(rbp_pairs)
# #DEBUG LIMIT 3
# q <- Query("SELECT protein_a_unique_pairs AS symbol1, protein_b_unique_pairs AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b!='N/A' AND resampling_p_value_b_vs_a!='N/A' AND resampling_wilcoxon_p_value_a_vs_b!='N/A' AND resampling_wilcoxon_p_value_b_vs_a!='N/A') AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a_unique_pairs, protein_b_unique_pairs ORDER BY protein_a_unique_pairs, protein_b_unique_pairs LIMIT 3")
# #END DEBUG
pairs_of_interest <-
  rbp_pairs %>% 
  mutate(pair=glue("{symbol1}|{symbol2}")) %>% 
  mutate(pair_inv=glue("{symbol2}|{symbol1}")) %>% 
  pivot_longer(cols=c(pair, pair_inv)) %>% 
  pull(value) %>%
  unique %>%
  sort
str(pairs_of_interest)

# Without adding inverse pairs
pairs_of_interest_noinv <-
  rbp_pairs %>% 
  mutate(pair=glue("{symbol1}|{symbol2}")) %>% 
  # mutate(pair_inv=glue("{symbol2}|{symbol1}")) %>% 
  # pivot_longer(cols=c(pair, pair_inv)) %>% 
  pull(pair) %>%
  unique %>%
  sort
str(pairs_of_interest_noinv)





# Ratio (with pseudocounts)
mypair <- "EWSR1|HNRNPK"
mypair <- "HNRNPK|PCBP1"
# for (mypair in c("EWSR1|HNRNPK", "FXR2|FMR1", "HNRNPK|PCBP1", "NONO|EWSR1", "NONO|SFPQ")) {
for (mypair in c("EWSR1|NONO", "NONO|EWSR1")) {
  q %>%
    filter(pair==mypair) %>%
    ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed") + ggtitle(mypair)
  ggsave(paste0("output-motifpairs-ratio-boxplots-",str_replace(mypair, '\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5*max((length(unique(subset(q, pair==mypair)$motifpair)) / 54), 1), units="mm", limitsize=F)
}
for (mypair in c("EWSR1|NONO", "NONO|EWSR1")) {
  q %>%
    filter(pair==mypair) %>%
    ggplot(aes(x=fct_reorder(motifpair, ratio_solo), y=ratio_solo)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed") + ggtitle(mypair)
  ggsave(paste0("output-motifpairs-ratio-boxplots-solo-",str_replace(mypair, '\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5*max((length(unique(subset(q, pair==mypair)$motifpair)) / 54), 1), units="mm", limitsize=F)
}
# Conclusion: These are asymmetric (pair orientation matters)
# q %>%
#   filter(pair=="FXR2|FMR1") %>%
#   ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed")
# q %>%
#   filter(pair=="HNRNPK|PCBP1") %>%
#   ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed")
# q %>%
#   filter(pair=="NONO|EWSR1") %>%
#   ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed")
# q %>%
#   filter(pair=="NONO|SFPQ") %>%
#   ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed")
q




# # Mean ratio
# q %>%
#   group_by(symbol1, symbol2, pair, motifpair) %>%
#   summarise(mean_ratio=mean(ratio), median_ratio=median(ratio), sd_ratio=sd(ratio), median_raw_ratio=median(raw_ratio), sd_raw_ratio=sd(raw_ratio), .groups="drop_last") %>%
#   arrange(-mean_ratio) %>%
#   filter(mean_ratio - sd_ratio >= 4) %>%
#   # head(n=100) %>%
#   ggplot(aes(x=fct_reorder(motifpair, mean_ratio), y=mean_ratio)) + coord_flip() + geom_boxplot(notch=T) + geom_errorbar(aes(ymin=mean_ratio-sd_ratio, ymax=mean_ratio+sd_ratio)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=4, linetype="dashed")
# # ggsave(paste0("output-motifpairs-ratio-boxplots-top100.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")
# 
# 
# # Raw ratio (no pseudocounts)
# q %>% pull(pair) %>% unique %>% length
# q %>%
#   # head(n=100000) %>%
#   ggplot(aes(x=ratio, y=fct_reorder(pair, ratio, .fun = mean))) + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(paste0("output-motifpairs-raw_ratio-boxplots.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")
# 
# # Ratio (with pseudocounts)
# q %>%
#   head(n=1000) %>%
#   ggplot(aes(x=factor(motifpair), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(paste0("output-motifpairs-raw_ratio-boxplots.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")






# Set ratio threshold (motif pairs of interest)
ratio_threshold <- 2


# Try: motif pairs of interest based on the resampled median of ratio_solo
q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(pair) %>% unique %>% length
# >> 8-fold leads to 9 RBP pairs of interest


# Try: motif pairs of interest based on the resampled median of ratio_solo, for individual sources
mysource <- 'attract'
mysource <- 'cisbp'
mysource <- 'dominguez'
mysource <- 'rbpdb'
mysource <- 'rbpmap'
mysource <- 'rnacompete'
q %>%
  filter(source1 == mysource & source2==mysource) %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(pair) %>% unique %>% length
# >> ATTRACT: 8 RBP pairs of interest with 8-fold, 14 with ≥4-fold, 16 with ≥2-fold
# >> CISBP: Just 2 RBP pairs of interest with 8-fold, 4 with ≥2 fold (FMR1|FXR2 & vice versa) (no gain)
# >> Dominguez: 0 RBP pairs of interest, even with ≥2 fold (0) (no gain)
# >> RBPDB: 0 RBP pairs of interest, even with ≥2 fold (0) (no gain)
# >> RBPMAP: 2 RBP pairs of interest, 2 with ≥2 fold (FMR1|FXR2 & vice versa) (no gain)
# >> RNAcompete: Just 2 RBP pairs of interest, even with ≥2 fold (FMR1|FXR2 & vice versa) (no gain)


q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(motifpair) %>% unique %>% length
# >> 8-fold leads to 253 motif pairs of interest

# Assign motif pairs of interest for ratio
motifpairs_of_interest_ratio <- q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio=median(ratio), sd_ratio=sd(ratio)) %>%
  arrange(-median_ratio) %>%
  filter(median_ratio >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(motifpair) %>% unique %>% sort
str(motifpairs_of_interest_ratio)
# >> Leads to 147 motif pairs of interest










# PLOT A: Show each pair only in one orientation, and only one top motif

# Assign motif pairs of interest for ratio_solo (top by median ratio_solo)
q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest_noinv) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest_ratio_solo
str(motifpairs_of_interest_ratio_solo)
# >> Leads to 675 motif pairs of interest


# # Overlap between motif pairs of interest selected by ratio or ratio_solo:
# intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo) %>% length
# # >> 78 pairs are robust
# motifpairs_of_interest_both <- intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo)


# Final:
# Only get the top ones for ratio_solo (randomly choosing a motif combination where there's a tie, which happens for a few RBPs)
# DEBUG: Should investigate if these motifs are identical (with_ties=T)
# q %>%
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
#   group_by(pair, motifpair) %>%
#   summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
#   slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
#   pull(motifpair) %>% 
#   unique %>% 
#   sort -> 
#   motifpairs_of_interest
# END DEBUG
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
  pull(motifpair) %>% unique %>% sort %>% tibble %>% print(n=1000)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
  select(pair, motifpair) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest
motifpairs_of_interest




# Filter for RBP pairs and motif pairs of interest
# qf <- q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% arrange(desc(ratio_solo))
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo)
# %>% arrange(desc(ratio_solo))
# q
# max(q$ratio_solo)
# q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)))
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair))) %>% summarise(max(n))
# max(q$ratio_solo) was 100, which I found suspicious, but it's correct (it results from the pseudocount ratio for a count of 99 vs. 0: 100/1)

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio: 4 RBP pairs, 147 motif pairs

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio_solo: 9 RBP pairs, 253 motif pairs

# After filtering for both: 4 RBP pairs, 78 motif pairs

# # Remove "CCCCCCC" pairs
# q %>% 
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
#   group_by(motif1, motif2, motifpair) %>% 
#   summarise(median_ratio_solo=median(ratio_solo)) %>%
#   arrange(desc(median_ratio_solo)) %>% 
#   select(motif1, motif2) %>% 
#   mutate(motifpair=glue("{motif1} & {motif2}")) %>%
#   ungroup %>%
#   select(motifpair) %>%
#   unique %>%
#   head(4) %>% pull(motifpair) -> motifpairs_to_remove

# Do NOT remove "CCCCCCC" pairs (where both RBPs have "CCCCCCC" as their top motif)
motifpairs_to_remove <- NULL

motifpairs_to_remove
motifpairs_of_interest

str(motifpairs_to_remove)
str(motifpairs_of_interest)

motifpairs_of_interest %>% tibble %>% mutate(tmp=str_replace_all(., "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(.) %>% length
motifpairs_of_interest %<>% tibble %>% rename(orig=".") %>% mutate(tmp=str_replace_all(orig, "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(orig)




# Ratio (with pseudocounts)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>%
  ggplot(aes(x=ratio_solo, y=fct_reorder(pair, ratio_solo, .fun = median))) + 
  geom_boxplot(notch=T) + 
  geom_vline(xintercept=1, linetype="dashed") +
  scale_x_continuous(limits=c(0.5, NA), trans=log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x)) +
  # labels = trans_format("log2", math_format(2^.x))) +
  xlab("Motif pair pseudocount enrichment at cobound sites") +
  ylab(NULL) +
  # theme_bw() + 
  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin=unit(c(5.5, 75.5, 5.5, 55.5), "points")) ->
  # theme(plot.margin=unit(c(5.5, 125.5, 5.5, 5.5), "points")) ->
  p
p

# Add logos
# motifpairs_of_interest
q %>% 
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
  group_by(motif1, motif2, motifpair) %>% 
  summarise(median_ratio_solo=median(ratio_solo)) %>%
  arrange(desc(median_ratio_solo)) %>% 
  select(motif1, motif2) %>% 
  unique ->
  logo_pairs_to_plot
logo_pairs_to_plot

# Left
lr <- 1
i <- 1
# Left/right hand side
for (lr in 1:2) {
  print(glue(" >> {lr}"))
  
  # Cycle through the motifs
  for (i in 1:nrow(logo_pairs_to_plot[lr])) {
    print(glue("   >> {i}"))
    
    mymotif <- logo_pairs_to_plot[i, lr]
    mymotif
    
    
    if (nrow(logo_pairs_to_plot[lr]) == 9) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.9, hjust=0, vjust=0, width=0.13, height=0.06)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.11, hjust=0, vjust=0, width=0.13, height=0.06)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.9, hjust=0, vjust=0, width=0.13, height=0.06)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.11, hjust=0, vjust=0, width=0.13, height=0.06)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.895 - (0.74/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.06) -> p
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.7 + 0.14*(lr-1), y=0.92 - (0.84/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.84, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.04, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.84, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.04, hjust=0, vjust=0, width=0.13, height=0.2)
      
      # Dynamic
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2/(nrow(logo_pairs_to_plot[lr])/9)) -> p
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2) -> p
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2) -> p
    }
  }
}
# p # Saving to PDF is much quicker than plotting in RStudio here:
ggsave(glue("output-ratio_solo-boxplots-motifs-ratio_threshold{ratio_threshold}-{nrow(logo_pairs_to_plot[lr])}_pairs_noinv.pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*1, units="mm")

















# PLOT B: Show only one top motif for each RBP

# Assign motif pairs of interest for ratio_solo (top by median ratio_solo)
q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest_ratio_solo
str(motifpairs_of_interest_ratio_solo)
# >> Leads to 253 motif pairs of interest


# # Overlap between motif pairs of interest selected by ratio or ratio_solo:
# intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo) %>% length
# # >> 78 pairs are robust
# motifpairs_of_interest_both <- intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo)


# Final:
# Only get the top ones for ratio_solo (randomly choosing a motif combination where there's a tie, which happens for a few RBPs)
# DEBUG: Should investigate if these motifs are identical (with_ties=T)
# q %>%
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
#   group_by(pair, motifpair) %>%
#   summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
#   slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
#   pull(motifpair) %>% 
#   unique %>% 
#   sort -> 
#   motifpairs_of_interest
# END DEBUG
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
  pull(motifpair) %>% unique %>% sort %>% tibble %>% print(n=1000)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
  select(pair, motifpair) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest
motifpairs_of_interest




# Filter for RBP pairs and motif pairs of interest
# qf <- q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% arrange(desc(ratio_solo))
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo)
# %>% arrange(desc(ratio_solo))
# q
# max(q$ratio_solo)
# q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)))
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair))) %>% summarise(max(n))
# max(q$ratio_solo) was 100, which I found suspicious, but it's correct (it results from the pseudocount ratio for a count of 99 vs. 0: 100/1)

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio: 4 RBP pairs, 147 motif pairs

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio_solo: 9 RBP pairs, 253 motif pairs

# After filtering for both: 4 RBP pairs, 78 motif pairs

# # Remove "CCCCCCC" pairs
# q %>% 
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
#   group_by(motif1, motif2, motifpair) %>% 
#   summarise(median_ratio_solo=median(ratio_solo)) %>%
#   arrange(desc(median_ratio_solo)) %>% 
#   select(motif1, motif2) %>% 
#   mutate(motifpair=glue("{motif1} & {motif2}")) %>%
#   ungroup %>%
#   select(motifpair) %>%
#   unique %>%
#   head(4) %>% pull(motifpair) -> motifpairs_to_remove

# Do NOT remove "CCCCCCC" pairs (where both RBPs have "CCCCCCC" as their top motif)
motifpairs_to_remove <- NULL

motifpairs_to_remove
motifpairs_of_interest

str(motifpairs_to_remove)
str(motifpairs_of_interest)

motifpairs_of_interest %>% tibble %>% mutate(tmp=str_replace_all(., "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(.) %>% length
motifpairs_of_interest %<>% tibble %>% rename(orig=".") %>% mutate(tmp=str_replace_all(orig, "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(orig)




# Ratio (with pseudocounts)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>%
  ggplot(aes(x=ratio_solo, y=fct_reorder(pair, ratio_solo, .fun = median))) + 
  geom_boxplot(notch=T) + 
  geom_vline(xintercept=1, linetype="dashed") +
  scale_x_continuous(limits=c(0.5, NA), trans=log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x)) +
  # labels = trans_format("log2", math_format(2^.x))) +
  xlab("Motif pair pseudocount enrichment at cobound sites") +
  ylab(NULL) +
  # theme_bw() + 
  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin=unit(c(5.5, 75.5, 5.5, 55.5), "points")) ->
  # theme(plot.margin=unit(c(5.5, 125.5, 5.5, 5.5), "points")) ->
  p
p

# Add logos
# motifpairs_of_interest
q %>% 
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
  group_by(motif1, motif2, motifpair) %>% 
  summarise(median_ratio_solo=median(ratio_solo)) %>%
  arrange(desc(median_ratio_solo)) %>% 
  select(motif1, motif2) %>% 
  unique ->
  logo_pairs_to_plot
logo_pairs_to_plot

# Left
lr <- 1
i <- 1
# Left/right hand side
for (lr in 1:2) {
  print(glue(" >> {lr}"))
  
  # Cycle through the motifs
  for (i in 1:nrow(logo_pairs_to_plot[lr])) {
    print(glue("   >> {i}"))
    
    mymotif <- logo_pairs_to_plot[i, lr]
    mymotif
    
    
    if (nrow(logo_pairs_to_plot[lr]) == 18) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.94, hjust=0, vjust=0, width=0.13, height=0.04)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.09, hjust=0, vjust=0, width=0.13, height=0.04)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.94, hjust=0, vjust=0, width=0.13, height=0.04)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.09, hjust=0, vjust=0, width=0.13, height=0.04)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.94 - (0.85/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.04) -> p
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.7 + 0.14*(lr-1), y=0.92 - (0.84/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else if (nrow(logo_pairs_to_plot[lr]) == 15) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.92, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.08, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.92, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.08, hjust=0, vjust=0, width=0.13, height=0.07)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.92 - (0.84/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else if (nrow(logo_pairs_to_plot[lr]) == 14) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.92, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.08, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.92, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.08, hjust=0, vjust=0, width=0.13, height=0.07)
      
      # Dynamic
      # Not yet fine-tuned
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.92 - (0.84/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else if (nrow(logo_pairs_to_plot[lr]) == 11) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.91, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.09, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.91, hjust=0, vjust=0, width=0.13, height=0.07)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.09, hjust=0, vjust=0, width=0.13, height=0.07)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.91 - (0.82/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else if (nrow(logo_pairs_to_plot[lr]) == 5) {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.79, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.09, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.79, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.09, hjust=0, vjust=0, width=0.13, height=0.2)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.79 - (0.7/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.07) -> p
    } else {
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.84, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.04, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.84, hjust=0, vjust=0, width=0.13, height=0.2)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.04, hjust=0, vjust=0, width=0.13, height=0.2)
      
      # Dynamic
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2/(nrow(logo_pairs_to_plot[lr])/9)) -> p
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2) -> p
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84*(lr-1), y=0.84 - (0.8/(nrow(logo_pairs_to_plot[lr])-1)) * (i-1), hjust=0, vjust=0, width=0.13, height=0.2) -> p
    }
  }
}
# p # Saving to PDF is much quicker than plotting in RStudio here:
ggsave(glue("output-ratio_solo-boxplots-motifs-ratio_threshold{ratio_threshold}-{nrow(logo_pairs_to_plot[lr])}_pairs.pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*1.5, units="mm")

















# PLOT C: Show all motifs for each RBP

# Assign motif pairs of interest for ratio_solo (top by median ratio_solo)
q %>%
  group_by(symbol1, symbol2, pair, motifpair) %>%
  summarise(median_ratio_solo=median(ratio_solo), sd_ratio_solo=sd(ratio_solo)) %>%
  arrange(-median_ratio_solo) %>%
  filter(median_ratio_solo >= ratio_threshold) %>%
  filter(pair %in% pairs_of_interest) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest_ratio_solo
str(motifpairs_of_interest_ratio_solo)
# >> Leads to 1319 motif pairs of interest


# # Overlap between motif pairs of interest selected by ratio or ratio_solo:
# intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo) %>% length
# # >> 78 pairs are robust
# motifpairs_of_interest_both <- intersect(motifpairs_of_interest_ratio, motifpairs_of_interest_ratio_solo)


# Final:
# Only get the top ones for ratio_solo (retaining all motif combinations where there's a tie, which happens for a few RBPs)
# DEBUG: Should investigate if these motifs are identical (with_ties=T)
# q %>%
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
#   group_by(pair, motifpair) %>%
#   summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
#   slice_max(motifpair_median_ratio_solo, with_ties=F) %>%
#   pull(motifpair) %>% 
#   unique %>% 
#   sort -> 
#   motifpairs_of_interest
# END DEBUG
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=T) %>%
  pull(motifpair) %>% unique %>% sort %>% tibble %>% print(n=1000)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>%
  group_by(pair, motifpair) %>%
  summarise(motifpair_median_ratio_solo=median(ratio_solo)) %>%
  slice_max(motifpair_median_ratio_solo, with_ties=T) %>%
  select(pair, motifpair) %>%
  pull(motifpair) %>%
  unique %>%
  sort ->
  motifpairs_of_interest
motifpairs_of_interest




# Filter for RBP pairs and motif pairs of interest
# qf <- q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% arrange(desc(ratio_solo))
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo)
# %>% arrange(desc(ratio_solo))
# q
# max(q$ratio_solo)
# q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair)))
# counts %>% filter(motifpair %in% (q %>% filter(ratio_solo==max(ratio_solo)) %>% pull(motifpair))) %>% summarise(max(n))
# max(q$ratio_solo) was 100, which I found suspicious, but it's correct (it results from the pseudocount ratio for a count of 99 vs. 0: 100/1)

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio: 4 RBP pairs, 147 motif pairs

q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(pair) %>% unique %>% length
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest_ratio_solo) %>% pull(motifpair) %>% unique %>% length
# After filtering for ratio_solo: 9 RBP pairs, 253 motif pairs

# After filtering for both: 4 RBP pairs, 78 motif pairs

# # Remove "CCCCCCC" pairs
# q %>% 
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
#   group_by(motif1, motif2, motifpair) %>% 
#   summarise(median_ratio_solo=median(ratio_solo)) %>%
#   arrange(desc(median_ratio_solo)) %>% 
#   select(motif1, motif2) %>% 
#   mutate(motifpair=glue("{motif1} & {motif2}")) %>%
#   ungroup %>%
#   select(motifpair) %>%
#   unique %>%
#   head(4) %>% pull(motifpair) -> motifpairs_to_remove

# Do NOT remove "CCCCCCC" pairs (where both RBPs have "CCCCCCC" as their top motif)
motifpairs_to_remove <- NULL

motifpairs_to_remove
motifpairs_of_interest

str(motifpairs_to_remove)
str(motifpairs_of_interest)

motifpairs_of_interest %>% tibble %>% mutate(tmp=str_replace_all(., "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(.) %>% length
motifpairs_of_interest %<>% tibble %>% rename(orig=".") %>% mutate(tmp=str_replace_all(orig, "\\|", "_")) %>% filter(!(tmp %in% motifpairs_to_remove)) %>% pull(orig)




# Ratio (with pseudocounts)
q %>%
  filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>%
  ggplot(aes(x=ratio_solo, y=fct_reorder(pair, ratio_solo, .fun = median))) + 
  geom_boxplot(notch=T) + 
  geom_vline(xintercept=1, linetype="dashed") +
  scale_x_continuous(limits=c(0.5, NA), trans=log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x)) +
  # labels = trans_format("log2", math_format(2^.x))) +
  xlab("Motif pair pseudocount enrichment at cobound sites") +
  ylab(NULL) +
  # theme_bw() + 
  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin=unit(c(5.5, 165.5, 5.5, 165.5), "points")) ->
  p
p

# Get RBP pairs of interest (the ones being plotted here)
# q %>% 
#   filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% 
#   group_by(pair, motifpair) %>% 
#   summarise(median_ratio_solo=median(ratio_solo)) %>%
#   arrange(desc(median_ratio_solo)) %>%
#   pull(pair) %>% 
#   unique ->
#     my_rbp_pairs_of_interest
# Get them directly from the plot:
my_rbp_pairs_of_interest <- rev(ggplot_build(p)$layout$panel_params[[1]]$y$breaks)


# Get maximum number of motif pairs per RBP pair
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% select(pair, motifpair) %>% unique %>%
  group_by(pair) %>% tally %>% arrange(desc(n))
q %>% filter(pair %in% pairs_of_interest & motifpair %in% motifpairs_of_interest) %>% select(pair, motifpair) %>% unique %>%
  group_by(pair) %>% tally %>% arrange(desc(n)) %>% pull(n) %>% max ->
  max_motifpairs


# Left
lr <- 1
rbp <- 1
i <- 1
# Left/right hand side
for (lr in 1:2) {
  print(glue(" >> lr {lr}"))
  
  # Cycle through RBP pairs of interest (the ones being plotted here)
  for (rbp in 1:length(my_rbp_pairs_of_interest)) {
    print(glue("   >> rbp pair {rbp} ({my_rbp_pairs_of_interest[rbp]})"))

    # Get their motifs
    q %>% 
      filter(pair==my_rbp_pairs_of_interest[[rbp]] & motifpair %in% motifpairs_of_interest) %>%
      group_by(motif1, motif2, motifpair) %>%
      summarise(median_ratio_solo=median(ratio_solo)) %>%
      arrange(desc(median_ratio_solo)) %>% 
      select(motif1, motif2) %>% 
      unique ->
      logo_pairs_to_plot
    logo_pairs_to_plot
    
    # Cycle through their motifs
    for (i in 1:nrow(logo_pairs_to_plot[lr])) {
      print(glue("     >> {i}"))
      
      mymotif <- logo_pairs_to_plot[i, lr]
      mymotif
      
      
      # # Top left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=(((i-1)/(max_motifpairs-1)) * 0.2), y=0.94, hjust=0, vjust=0, width=0.1, height=0.03)
      # # Bottom left
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0, y=0.09, hjust=0, vjust=0, width=0.1, height=0.03)
      # # Top right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=(((i-1)/(max_motifpairs-1)) * 0.2) + 0.67, y=0.94, hjust=0, vjust=0, width=0.1, height=0.03)
      # # Bottom right
      # ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=0.84, y=0.09, hjust=0, vjust=0, width=0.1, height=0.03)
      
      # Dynamic
      ggdraw(p) + draw_image(glue("output_logos_png_trim/{mymotif}.png"), x=((0.25 - (((i-1)/(max_motifpairs-1)) * 0.25)) * abs((lr-2))) + (((((i-1)/(max_motifpairs-1)) * 0.25) + 0.67) * (lr-1)), y=0.94 - (0.85/(length(my_rbp_pairs_of_interest)-1)) * (rbp-1), hjust=0, vjust=0, width=0.07, height=0.025) -> p
    }
  }
}
# p # Saving to PDF is much quicker than plotting in RStudio here:
# Normal plot
ggsave(glue("output-ratio_solo-boxplots-motifs-ratio_threshold{ratio_threshold}-{length(my_rbp_pairs_of_interest)}_pairs_multi.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")
# # Huge long plot
# ggsave(glue("output-ratio_solo-boxplots-motifs-ratio_threshold{ratio_threshold}-{nrow(logo_pairs_to_plot[lr])}_pairs.pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*3, units="mm")



