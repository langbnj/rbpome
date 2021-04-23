rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(reshape2)
library(glue)
library(broom)
library(scales)
library(magrittr)


eclip_rbps <- Query("SELECT DISTINCT symbol FROM clip_raw_gene WHERE type='eclip_encode_12'")
eclip_rbps

# All interactions that include at least one eCLIP RBP (77 pairs)
q1 <- Query("SELECT DISTINCT symbol1 FROM rbpome WHERE scoretype='SUM_IS' AND avg_is>=7.1 AND homodimer=0")
q2 <- Query("SELECT DISTINCT symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND avg_is>=7.1 AND homodimer=0")
q1 %<>% transmute(symbol = symbol1)
q2 %<>% transmute(symbol = symbol2)
q1
q2
q <- bind_rows(q1, q2)
q %<>% unique
q
q %<>% filter(symbol %in% eclip_rbps$symbol) %>% unique
q
q$symbol %>% unique %>% cat(sep="\n")
# Get unique RBPs in these pairs (50 RBPs)
rbps_both <- q$symbol %>% unique
rbps_both
rbps_both %>% nrow

# eCLIP-eCLIP interactions only (71 pairs, as in pipeline/rbpome_analysis/main.pl)
q <- Query("SELECT DISTINCT symbol1, symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND avg_is>=7.1 AND homodimer=0")
q %<>% filter(symbol1 %in% eclip_rbps$symbol & symbol2 %in% eclip_rbps$symbol) %>% unique
q
# c(q$symbol1, q$symbol2) %>% unique %>% length
# Get unique RBPs in these pairs (50 RBPs)
rbps_both <- bind_rows(q %>% select(symbol = symbol1), q%>% select(symbol = symbol2)) %>% unique
rbps_both
rbps_both %>% nrow

# Get target gene lengths for these 50 RBPs
q

# write_tsv(q, "~/Dropbox/RBPome_paper/#FINAL/NAR/Revision/Table R2/Table R2.tsv")



# fit_vs_random.pl uses this (50 RBPs):
q <- Query("SELECT symbol1, symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2")
# Get unique RBPs in these pairs (50 RBPs)
rbps_both <- bind_rows(q %>% select(symbol = symbol1), q%>% select(symbol = symbol2)) %>% unique
rbps_both
rbps_both %>% nrow



# Get target gene lengths for these 50 RBPs
print(paste0("SELECT c.symbol, (g.stop-g.start)+1 AS len FROM clip_raw_gene c JOIN gencode_gff3_gene g ON g.ensgv=c.ensgv WHERE c.type='eclip_encode_12' AND c.symbol IN ('",paste0(rbps_both$symbol, collapse="', '"),"') GROUP BY c.symbol, g.ensgv"))
# # Group by RBP, ENSGV
# q <- Query(paste0("SELECT c.symbol, (g.stop-g.start)+1 AS len FROM clip_raw_gene c JOIN gencode_gff3_gene g ON g.ensgv=c.ensgv WHERE c.type='eclip_encode_12' AND c.symbol IN ('",paste0(rbps_both$symbol, collapse="', '"),"') GROUP BY c.symbol, g.ensgv"))
# # Group by ENSGV
# q <- Query(paste0("SELECT (g.stop-g.start)+1 AS len FROM clip_raw_gene c JOIN gencode_gff3_gene g ON g.ensgv=c.ensgv WHERE c.type='eclip_encode_12' AND c.symbol IN ('",paste0(rbps_both$symbol, collapse="', '"),"') GROUP BY g.ensgv"))
# Group by RBP, ENSGV and require more than one peak per RBP-ENSGV (so that there is actually a distance to be calculated)
q <- Query(paste0("SELECT * FROM (SELECT c.symbol, g.ensgv, (g.stop-g.start)+1 AS len, COUNT(DISTINCT c.id) AS peaks FROM clip_raw_gene c JOIN gencode_gff3_gene g ON g.ensgv=c.ensgv AND g.species='human' WHERE c.type='eclip_encode_12' AND c.symbol IN ('",paste0(rbps_both$symbol, collapse="', '"),"') GROUP BY c.symbol, g.ensgv) t WHERE peaks>1"))
q
# Plot a histogram
q %>% 
  ggplot(aes(x=len)) + 
  geom_histogram() + 
  # scale_x_log10() +
  scale_x_log10(labels = scales::comma) +
  # scale_x_log10(labels = scales::comma, breaks= breaks_log(8)) +
  annotation_logticks(sides = "b") +
  # geom_vline(xintercept = median(q$len), linetype="dashed") +
  # geom_vline(xintercept = 1000, linetype="solid") +
  geom_vline(xintercept = 1000, linetype="dashed") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  xlab("Gene length (of target RNAs)") +
  ylab("RBP-Gene interactions") +
  ggsave("Figure R1.pdf", width=4, height=4)

median(q$len)
