rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
setwd("~/Documents/Projects/RBPome/RBP-RNA interactor counts/")

library(tidyverse)
library(scales)


# # q <- Query("SELECT a.celltype, a.ensgv, COUNT(DISTINCT a.rep) AS reps, COUNT(DISTINCT c.symbol) AS rbps FROM encode_abundance a, clip_raw_gene c WHERE a.expressed=1 AND a.ensgv=c.ensgv AND c.type='eclip_encode_12' AND a.celltype=c.celltype GROUP BY a.celltype, a.ensgv")
# q <- Query("SELECT a.ensgv, COUNT(DISTINCT a.celltype) AS celltypes, COUNT(DISTINCT a.rep) AS reps, COUNT(DISTINCT c.symbol) AS rbps FROM encode_abundance a, clip_raw_gene c WHERE a.expressed=1 AND a.ensgv=c.ensgv AND c.type='eclip_encode_12' GROUP BY a.ensgv")
# q %<>% filter(reps==2 & celltypes==2)

# qa <- Query("SELECT a.celltype, a.ensgv, COUNT(DISTINCT a.rep) AS reps FROM encode_abundance a WHERE a.expressed=1 GROUP BY a.celltype, a.ensgv")
# qc <- Query("SELECT c.celltype, c.ensgv, COUNT(DISTINCT c.symbol) AS rbps FROM clip_raw_gene c WHERE c.type='eclip_encode_12' GROUP BY c.celltype, c.ensgv")

# qa %<>% filter(reps==2)
# qa %<>% filter(reps==2 & celltypes==2)

qa <- Query("SELECT a.ensgv, COUNT(DISTINCT a.celltype, a.rep) AS celltypereps, COUNT(DISTINCT a.celltype) AS celltypes, COUNT(DISTINCT a.rep) AS reps FROM encode_abundance a WHERE a.expressed=1 GROUP BY a.ensgv")
qc <- Query("SELECT c.ensgv, COUNT(DISTINCT c.symbol) AS rbps FROM clip_raw_gene c WHERE c.type='eclip_encode_12' GROUP BY c.ensgv")

qa %<>% filter(celltypereps==4)
qa
qc
qj <- inner_join(qa, qc)
qj

# all_equal(q, qj)
q <- qj

q %>% 
  ggplot(aes(x=rbps)) + 
  geom_histogram(binwidth=10, boundary=0) +
  # facet_grid(cols = vars(celltype)) +
  geom_vline(xintercept = median(q$rbps), linetype="dashed") +
  scale_x_continuous(breaks=breaks_pretty(8), expand = c(0,0)) +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  xlab("RBPs bound per RNA gene") +
  ylab("RNA genes") +
  ggsave("output-rbp_rna_interactor_counts-histogram.pdf", width=5, height=4)

# q %>% 
#   ggplot(aes(x="hi", y=rbps)) + 
#   geom_violin(draw_quantiles = c(0.5)) +
#   # facet_grid(cols = vars(celltype)) +
#   # geom_vline(xintercept = median(q$rbps), linetype="dashed") +
#   # scale_x_continuous(breaks=breaks_pretty(8), expand = c(0,0)) +
#   theme_minimal()

q %>% 
  # ggplot(aes(x=rbps, colour=celltype, fill=celltype)) + 
  ggplot(aes(x=rbps)) + 
  geom_density() +
  # facet_grid(cols = vars(celltype)) +
  geom_vline(xintercept = median(q$rbps), linetype="dashed") +
  theme_minimal()

q %>% 
  ggplot(aes(y=celltype, x=rbps, colour=celltype, fill=celltype)) + 
  geom_violin(alpha=0.3, draw_quantiles = c(0.5)) +
  # facet_grid(cols = vars(celltype)) +
  # geom_vline(xintercept = median(q$rbps), linetype="dashed") +
  guides(colour=F, fill=F) +
  theme_minimal() +
  xlab("RBPs bound per RNA gene") +
  ylab("Cell line")

median(q$rbps)
