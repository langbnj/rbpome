rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(glue)
library(scales)

# library(broom)
# library(ggrepel)
library(magrittr)
# library(naturalsort)
# library(progress)
# library(reshape2)

# install.packages("foreach")
# install.packages("doParallel")
# install.packages("GenomicRanges")
library(foreach)
library(doParallel)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
library(GenomicRanges)

setwd("~/Documents/Projects/RBPome/Motif distances")

mindist_motifs <- 54
# motif_source <- "all"
# motif_source <- "attract"
# motif_source <- "dominguez"
# motif_source <- "rnacompete"
motif_sources <- c("attract", "dominguez", "rnacompete")




# "The authors should also analyze in detail the RBPs which are known to bind to motifs located proximally in target RNAs,
#  but are not detected by Y2H. What fraction of all pairs with proximal binding motifs are actually detected as binders?"

q <- Query("SELECT * FROM rbp_motifs_extend_fimo_eclip_encode_12")
q

set.seed(1)
# qs <- slice_head(q, n=10)
# qs <- slice_sample(q, n=1000)
# qs <- slice_sample(q, n=10000)
# qs <- slice_sample(q, prop=0.1) # 10% of data
qs <- q # All data
qs

qs %<>% mutate(motifcentre = motifstart + ((motifstop - motifstart) / 2))
qs %<>% select(-acc, -species, -method, -extend, -type, -pvalue, -qvalue, -psig, -qsig, -score)
qs

# qs %>%
#   group_by(symbol, chr, strand, id) %>%
#   # rowwise() %>%
#   # mutate(mindist = min(abs(motifcentre - qs[qs$symbol==symbol & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre)))
#   # Avoid warnings by setting the motifcentre vector to NA if it has length 0 (if there are no motif hits on the same chromosome)
#   mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre))))
# %>%
#   mutate(mindist = na_if(mindist, Inf))
# 
# ifelse(length(qs[qs$symbol=="TIA1" & qs$chr=="chr20" & qs$strand=="-" & qs$id!=856018,]$motifcentre) == 0


# Screen hit pairs
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course)
qtmp <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND sumis>=7.1 GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# qtmp %>% select(symbol1, symbol2) %>% unique %>% transmute(pair = glue("{symbol1} {symbol2}")) -> qtmppairs
qtmp %>% select(symbol1, symbol2) %>% unique -> screenpairs
print(screenpairs, n=50)
# Add inverse pairs
screenpairs %<>% bind_rows(screenpairs %>% transmute(newsymbol1=symbol2, newsymbol2=symbol1) %>% transmute(symbol1=newsymbol1, symbol2=newsymbol2)) %>% unique
screenpairs
screenpairs %<>% mutate(pair=as.character(glue("{symbol1}|{symbol2}"))) %>% mutate(screenpair=as.integer(1))
screenpairs

pairs <- tibble(symbol1=character(), symbol2=character())
# pairs
for (symbol1 in unique(qs$symbol)) {
  for (symbol2 in unique(qs$symbol)) {
    if (symbol1 != symbol2) {
      pairs %<>% bind_rows(tibble(symbol1=symbol1, symbol2=symbol2))
    }
  }
}
pairs
pairs %<>% left_join(screenpairs)
pairs
pairs %<>% 
  mutate(pair=as.character(glue("{symbol1}|{symbol2}"))) %>% 
  mutate(screenpair = as.integer(replace_na(screenpair, 0)))
pairs
pairs %>% summary

# # For performance: Sample a smaller number of non-screenpairs (random pairs)
# pairs %>% filter(screenpair == 1) %>% nrow -> n_screenpairs
# bind_rows(pairs %>% filter(screenpair == 0) %>% slice_sample(n=n_screenpairs),
#           pairs %>% filter(screenpair == 1)) -> pairs_sample
# pairs_sample
# pairs_sample %>% mutate(screenpair=as.factor(screenpair)) %>% summary

# Use all pairs
pairs_sample <- pairs

symbol1 = "FMR1"
symbol2 = "FXR2"
symbol1 = "FXR1"
symbol2 = "FXR2"
symbol1 = "FXR2"
symbol2 = "FXR1"

# for (symbol1 in unique(qs$symbol)) {
#   for (symbol2 in unique(qs$symbol)) {
#     if (symbol1 != symbol2) {
pairs_sample %<>% arrange(symbol1, symbol2)
pairs_sample
# suppressWarnings(try(rm(qfracs)))

# # Single-thread
# for (i in 1:nrow(pairs_sample)) {
#   symbol1 <- pairs_sample[i,]$symbol1
#   symbol2 <- pairs_sample[i,]$symbol2
#   screenpair <- pairs_sample[i,]$screenpair
#   cat(glue(" >> {symbol1}|{symbol2} >> screenpair {screenpair}\n\n"))
#   system.time(
#     qs %>%
#     filter(symbol == symbol2) %>%
#     group_by(symbol, chr, strand, id) %>%
#     # This is the fastest method (not checking the centre positions)
#     mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre)))) ->
#     # mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id & abs(qs$motifcentre - motifcentre) <= 1000,]$motifcentre)))) ->
#     # mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id & qs$motifcentre >= motifcentre-1000 & qs$motifcentre <= motifcentre+1000,]$motifcentre)))) ->
#       # %>%
#     # filter(!is.na(mindist)) %>%
#     # filter(mindist <= mindist_motifs) ->
#     qsp)
#   # qsp
#   qsp %>% group_by(symbol) %>% summarise(n=n(), n1000=sum(mindist<=1000), n54=sum(mindist<=54), frac54=n54/n1000, frac54_total=n54/n) ->
#     qfrac
#   qfrac %<>% mutate(pair=glue("{symbol1}|{symbol2}")) %>% select(-symbol) %>% relocate(pair)
#   qfrac
#   if (!exists("qfracs")) {
#     cat("   >> New\n\n")
#     print(qfrac)
#     qfracs <- qfrac
#   } else {
#     cat("   >> Adding\n\n")
#     print(qfrac)
#     qfracs %<>% bind_rows(qfrac)
#   }
# }




# # Parallelised:
# cores = detectCores()
# cores
# cl <- makeForkCluster(cores[1] - 2) # Leave 2 cores free
# registerDoParallel(cl)
# 
# # #DEBUG
# # pairs_sample %>% filter(pair %in% c("FMR1|FXR2", "FXR1|FXR2", "FXR2|FXR1")) -> pairs_sample_tmp
# # pairs_sample
# # pairs_sample_tmp
# #END DEBUG
# pairs_sample
# qfracs <- foreach(i = 1:nrow(pairs_sample), .combine = 'bind_rows') %dopar% {
#   symbol1 <- pairs_sample[i,]$symbol1
#   symbol2 <- pairs_sample[i,]$symbol2
#   screenpair <- pairs_sample[i,]$screenpair
#   # cat(glue(" >> {symbol1}|{symbol2} >> screenpair {screenpair}\n\n"))
# 
#   qs %>%
#     filter(symbol == symbol2) %>%
#     group_by(symbol, chr, strand, id) %>%
#     # This is the fastest method (not checking the centre positions)
#     mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre)))) ->
#     # mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id & abs(qs$motifcentre - motifcentre) <= 1000,]$motifcentre)))) ->
#     # mutate(mindist = min(abs(motifcentre - ifelse(length(qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id,]$motifcentre) == 0, NA, qs[qs$symbol==symbol1 & qs$chr==chr & qs$strand==strand & qs$id!=id & qs$motifcentre >= motifcentre-1000 & qs$motifcentre <= motifcentre+1000,]$motifcentre)))) ->
#     # %>%
#     # filter(!is.na(mindist)) %>%
#     # filter(mindist <= mindist_motifs) ->
#     qsp
#   # qsp
#   qsp %>% group_by(symbol) %>% summarise(n=n(), n1000=sum(mindist<=1000), n54=sum(mindist<=54), frac54=n54/n1000, frac54_total=n54/n) ->
#     qfrac
#   qfrac %<>% mutate(pair=glue("{symbol1}|{symbol2}")) %>% select(-symbol) %>% relocate(pair)
#   qfrac
# }
# stopCluster(cl)
# write_rds(qfracs, "qfracs.rds")







# Parallelised with GenomicRanges (much much faster):
suppressWarnings(try(rm(qaccompanied)))



# #DEBUG
# pairs_sample %>% filter(pair %in% c("FMR1|FXR2", "FXR1|FXR2", "FXR2|FXR1")) -> pairs_sample_tmp
# pairs_sample
# pairs_sample_tmp
# symbol1 = "FMR1"
# symbol2 = "FXR2"
# # symbol1 = "FXR1"
# # symbol2 = "FXR2"
# # symbol1 = "FXR2"
# # symbol2 = "FXR1"
# # END DEBUG
pairs_sample
# pairs_sample_tmp

get_qfracs <- function(motif_source) {
  suppressWarnings(try(rm(my_qaccompanied)))
  
  cores = detectCores()
  cores
  cl <- makeForkCluster(cores[1] - 2) # Leave 2 cores free
  registerDoParallel(cl)
  
  my_qaccompanied <- foreach(i = 1:nrow(pairs_sample), .combine = 'bind_rows') %dopar% {
    symbol1 <- pairs_sample[i,]$symbol1
    symbol2 <- pairs_sample[i,]$symbol2
    screenpair <- pairs_sample[i,]$screenpair
    # cat(glue(" >> {symbol1}|{symbol2} >> screenpair {screenpair}\n\n"))
    
    # Get motif coordinates for symbol1
    qs %>%
      filter(symbol == symbol1) %>%
      filter(source == motif_source) %>%
      transmute(coords=glue("{chr}:{motifstart}-{motifstop}:{strand}")) %>%
      pull(coords) %>%
      as("GRanges") -> coords1
    
    # Get motif coordinates for symbol2 (which we are relative to in symbol1|symbol2)
    qs %>%
      filter(symbol == symbol2) %>%
      filter(source == motif_source) %>%
      transmute(coords=glue("{chr}:{motifstart}-{motifstop}:{strand}")) %>%
      pull(coords) %>%
      as("GRanges") -> coords2
    coords2
    # coords2[1:2]
    # resize(coords2[1:2], width = width(coords2[1:2]) + (54 * 2), fix = "center")
    
    
    # Expand symbol2 coordinates (which we are relative to in symbol1|symbol2) by 54 nt in each direction.
    # No need to expand symbol1 coordinates, we'll simply check if they're in this range.
    # coords1
    # coords2
    # width(coords2)
    coords2 <- resize(coords2, width = width(coords2) + (mindist_motifs * 2), fix = "center")
    # coords2
    # width(coords2)
    
    # hit_regions_agreed_grange <- reduce(suppressWarnings(GenomicRanges::intersect(hit_regions_before_grange, hit_regions_after_grange)))
    # GenomicRanges::intersect(coords1, coords2)
    # reduce(GenomicRanges::intersect(coords1, coords2))
    # # reduce merges overlapping regions
    # 
    # coords2[1:2]
    # GenomicRanges::intersect(coords1, coords2[1:2])
    # GenomicRanges::intersect(coords2[1:2], coords1)
    
    # # ?countOverlaps(query, subject)
    # length(countOverlaps(coords1, coords2))
    # length(countOverlaps(coords2, coords1))
    # coords2
    # countOverlaps(coords1, coords2)
    # countOverlaps(coords2, coords1)
    # width(coords1)
    # width(coords2)
    
    # This is what we want (number of symbol1 motifs near symbol2 motifs, with symbol2 motifs being the reference)
    # countOverlaps(coords2, coords1)
    symbol2_motifs_total <- length(coords2)
    # symbol2_motifs_solo <- sum(countOverlaps(coords2, coords1) == 0)
    symbol2_motifs_accompanied <- sum(countOverlaps(coords2, coords1) > 0)
    symbol2_motifs_accompanied_fraction <- symbol2_motifs_accompanied / symbol2_motifs_total
    # symbol2_motifs_accompanied_fraction
    # symbol2_motifs_solo
    # symbol2_motifs_company <- sum(countOverlaps(coords2, coords1) > 0)
    # symbol2_motifs_company
    # symbol2_motifs_solo + symbol2_motifs_company
    
    qfrac <- tibble_row(symbol1=symbol1, symbol2=symbol2, pair=glue("{symbol1}|{symbol2}"), motif_source=motif_source, symbol2_motifs_accompanied=symbol2_motifs_accompanied, symbol2_motifs_total=symbol2_motifs_total, accompanied_fraction=symbol2_motifs_accompanied_fraction)
    qfrac
  }
  stopCluster(cl)
  write_rds(my_qaccompanied, glue("qaccompanied-{motif_source}.rds"))
  
  return(my_qaccompanied)
}

qaccompanied_attract <- get_qfracs("attract")
qaccompanied_dominguez <- get_qfracs("dominguez")
qaccompanied_rnacompete <- get_qfracs("rnacompete")

qaccompanied_attract
qaccompanied_dominguez
qaccompanied_rnacompete

qaccompanied <- bind_rows(qaccompanied_attract, qaccompanied_dominguez, qaccompanied_rnacompete)





# # "What fraction of all pairs with proximal binding motifs are actually detected as binders?"
# qfracs2 <- read_rds("qfracs.rds")
# qfracs2
# 
# # Add screenpair information
# qfracs2 %>%
#   left_join(screenpairs) %>%
#   mutate(screenpair = replace_na(screenpair, 0)) %>%
#   mutate(screenpair = as.factor(screenpair)) %>%
#   ggplot(aes(x=sum(screenpair == 1) / length(screenpair))) +
#   geom_bar() +
#   facet_grid(cols = vars(frac54 >= 0.1)) +
#   ggsave("output-motif-distances.pdf") +
#   ggsave("output-motif-distances.png")



# "What fraction of all pairs with proximal binding motifs are actually detected as binders?"
# qaccompanied <- read_rds("qaccompanied.rds")
qaccompanied

# # Add screenpair information
# qaccompanied %<>%
#   left_join(screenpairs) %>%
#   mutate(screenpair = replace_na(screenpair, 0)) %>%
#   mutate(screenpair = as.factor(screenpair))

# # Plot
# qaccompanied %>%
#   left_join(screenpairs) %>%
#   mutate(screenpair = replace_na(screenpair, 0)) %>%
#   mutate(screenpair = as.factor(screenpair)) %>%
#   ggplot(aes(x=sum(screenpair == 1) / length(screenpair))) +
#   geom_bar() +
#   facet_grid(cols = vars(accompanied_fraction >= 0.1)) +
#   ggsave("output-motif-distances.pdf") +
#   ggsave("output-motif-distances.png")

# Plot
qaccompanied %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  mutate(screenpair = as.factor(screenpair)) %>%
  # filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
  ggplot(aes(x=screenpair, y=accompanied_fraction)) +
  geom_boxplot(notch=T) +
  # facet_grid(cols = vars(accompanied_fraction>=0.1)) +
  theme_minimal()

# Plot ≥1000 motif hits
qaccompanied %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  mutate(screenpair = as.factor(screenpair)) %>%
  mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
  mutate(accnum = as.character(symbol2_motifs_accompanied>=1000)) %>%
  mutate(accnum = recode(accnum, "FALSE"="Adjacent motifs\n<1000", "TRUE"="Adjacent motifs\n≥1000")) %>%
  mutate(motif_source=recode(motif_source, "attract"="ATtRACT\ndatabase", "dominguez"="Dominguez et al.\nRNA Bind-N-Seq", "rnacompete"="RNAcompete")) %>%
  # filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
  ggplot(aes(x=screenpair)) +
  geom_bar() +
  # facet_wrap(~(symbol2_motifs_accompanied>=1000), scales="free_y") +
  facet_wrap(motif_source ~ accnum, ncol=2, scales="free") +
  # facet_grid(rows = vars(motif_source), cols = vars(accnum), scales="free") +
  theme_minimal() +
  xlab("") +
  ylab("Number of RBP pairs") +
  ggsave("Figure R10A output-motif_distances-count1000.pdf", width=5, height=7, device=cairo_pdf) +
  ggsave("output-motif_distances-count1000.pdf", width=5, height=7, device=cairo_pdf) +
  ggsave("output-motif-distances-count1000.png", width=5, height=7)

# # Plot 20%
# qaccompanied %>%
#   left_join(screenpairs) %>%
#   mutate(screenpair = replace_na(screenpair, 0)) %>%
#   mutate(screenpair = as.factor(screenpair)) %>%
#   mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
#   mutate(accfrac = as.character(accompanied_fraction>=0.2)) %>%
#   mutate(accfrac = recode(accfrac, "FALSE"="Adjacent motifs\n<20%", "TRUE"="Adjacent motifs\n≥20%")) %>%
#   mutate(motif_source=recode(motif_source, "attract"="ATtRACT\ndatabase", "dominguez"="Dominguez et al.\nRNA Bind-N-Seq", "rnacompete"="RNAcompete")) %>%
#   filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
#   ggplot(aes(x=screenpair)) +
#   geom_bar() +
#   # facet_wrap(~accfrac, scales="free_y") +
#   facet_wrap(motif_source ~ accfrac, ncol=2, scales="free") +
#   # facet_grid(rows = vars(motif_source), cols = vars(accfrac), scales="free") +
#   theme_minimal() +
#   xlab("") +
#   ylab("Number of RBP pairs") +
#   # ggsave("Figure R10B output-motif_distances-20percent.pdf", width=5, height=7, device=cairo_pdf) +
#   ggsave("output-motif_distances-20percent.pdf", width=5, height=7, device=cairo_pdf) +
#   ggsave("output-motif_distances-20percent.png", width=5, height=7)
# 
# # Plot 15%
# qaccompanied %>%
#   left_join(screenpairs) %>%
#   mutate(screenpair = replace_na(screenpair, 0)) %>%
#   mutate(screenpair = as.factor(screenpair)) %>%
#   mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
#   mutate(accfrac = as.character(accompanied_fraction>=0.15)) %>%
#   mutate(accfrac = recode(accfrac, "FALSE"="Adjacent motifs\n<15%", "TRUE"="Adjacent motifs\n≥15%")) %>%
#   mutate(motif_source=recode(motif_source, "attract"="ATtRACT\ndatabase", "dominguez"="Dominguez et al.\nRNA Bind-N-Seq", "rnacompete"="RNAcompete")) %>%
#   filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
#   ggplot(aes(x=screenpair)) +
#   geom_bar() +
#   # facet_wrap(~accfrac, scales="free_y") +
#   facet_wrap(motif_source ~ accfrac, ncol=2, scales="free") +
#   # facet_grid(rows = vars(motif_source), cols = vars(accfrac), scales="free") +
#   theme_minimal() +
#   xlab("") +
#   ylab("Number of RBP pairs") +
#   # ggsave("Figure R10B output-motif_distances-15percent.pdf", width=5, height=7, device=cairo_pdf) +
#   ggsave("output-motif_distances-15percent.pdf", width=5, height=7, device=cairo_pdf) +
#   ggsave("output-motif_distances-15percent.png", width=5, height=7)

# Plot ≥100 motif hits and 10%
qaccompanied %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  mutate(screenpair = as.factor(screenpair)) %>%
  mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
  mutate(accfrac = as.character(accompanied_fraction>=0.1)) %>%
  mutate(accfrac = recode(accfrac, "FALSE"="Adjacent motifs\n<10%", "TRUE"="Adjacent motifs\n≥10%")) %>%
  mutate(motif_source=recode(motif_source, "attract"="ATtRACT\ndatabase", "dominguez"="Dominguez et al.\nRNA Bind-N-Seq", "rnacompete"="RNAcompete")) %>%
  filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
  ggplot(aes(x=screenpair)) +
  geom_bar() +
  # facet_wrap(~accfrac, scales="free_y") +
  facet_wrap(motif_source ~ accfrac, ncol=2, scales="free") +
  # facet_grid(rows = vars(motif_source), cols = vars(accfrac), scales="free") +
  theme_minimal() +
  xlab("") +
  ylab("Number of RBP pairs") +
  ggsave("Figure R10B output-motif_distances-10percent.pdf", width=5, height=7, device=cairo_pdf) +
  ggsave("output-motif_distances-10percent.pdf", width=5, height=7, device=cairo_pdf) +
  ggsave("output-motif_distances-10percent.png", width=5, height=7)

# Best pairs (RNAcompete):
qaccompanied %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  mutate(screenpair = as.factor(screenpair)) %>%
  mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
  # mutate(accfrac = as.character(accompanied_fraction>=0.1)) %>%
  # mutate(accfrac = recode(accfrac, "FALSE"="Adjacent motifs\n<10%", "TRUE"="Adjacent motifs\n≥10%")) %>%
  # filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
  filter(motif_source == "rnacompete") %>%
  # filter(accfrac == "Adjacent motifs\n≥10%" & screenpair == "Screen hits\n(interacting RBPs)") %>%
  # filter(accfrac == "Adjacent motifs\n≥10%") %>%
  filter(accompanied_fraction >= 0.2) %>%
  arrange(-accompanied_fraction) %>%
  print(n=100)

# Best pairs (RNAcompete):
qaccompanied %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  mutate(screenpair = as.factor(screenpair)) %>%
  # mutate(screenpair = recode(screenpair, "0"="Random\nRBP pairs", "1"="Screen hits\n(interacting RBPs)")) %>%
  # mutate(accfrac = as.character(accompanied_fraction>=0.1)) %>%
  # mutate(accfrac = recode(accfrac, "FALSE"="Adjacent motifs\n<10%", "TRUE"="Adjacent motifs\n≥10%")) %>%
  # filter(symbol2_motifs_total >= 100) %>% # At least 100 motif hits for symbol2
  # filter(motif_source == "rnacompete") %>%
  # filter(accfrac == "Adjacent motifs\n≥10%" & screenpair == "Screen hits\n(interacting RBPs)") %>%
  # filter(accfrac == "Adjacent motifs\n≥10%") %>%
  # filter(accompanied_fraction >= 0.3) %>%
  arrange(-accompanied_fraction) %>%
  print(n=100) %>%
  write_tsv("output-table-S9.tsv")
