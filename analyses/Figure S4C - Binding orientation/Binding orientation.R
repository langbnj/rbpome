rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(reshape2)
library(glue)
library(broom)
library(scales)
library(forcats)
library(magrittr)
library(ggbeeswarm)
library(ggrepel)

setwd("~/Documents/Projects/RBPome/Binding orientation/")


# Make overview box plot (close binding events in the real data vs. in the random distribution)


table <- tibble(read.table("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed.txt", header=T, sep="\t", quote=""))
table

# q <- rbind(tibble(class="Observed", value=table$fraction_close_real), tibble(class="Randomised", value=table$fraction_close_random))
# str(q)
# head(q)
# 
# wilcox.test(table$fraction_close_real, table$fraction_close_random)
# 
# # Plot difference between median_dist_real and median_dist_random
# table$median_dist_delta <- table$median_dist_random - table$median_dist_real
# q <- tibble(value=table$median_dist_delta)
# str(q)
# head(q)



# Boxplot
ggplot(table, aes(x=median_dist_real)) + geom_boxplot(notch=T) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot, filtered to significant pairs
ggplot(subset(table, resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit "close"
ggplot(subset(table, abs(median_dist_real) <= 54), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-54, 54)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-limit54.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit "close", filtered to significant pairs
ggplot(subset(table, abs(median_dist_real) <= 54 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-54, 54)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-limit54-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 100, filtered to significant pairs
ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-limit100-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 1000
ggplot(subset(table, abs(median_dist_real) <= 1000), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-1000, 1000)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-limit1000.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 1000, filtered to significant pairs
ggplot(subset(table, abs(median_dist_real) <= 1000 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-1000, 1000)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-limit1000-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)



# Violin limit 100, filtered to significant pairs
# ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=NA)) + geom_violin() + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_density() + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=0)) + geom_violin() + geom_beeswarm(groupOnX=T) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=factor(""))) + geom_violin() + geom_point(position="jitter", alpha=0.3) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")

# Set labels
# table$label <- NULL
# table$label <- NA


# Beeswarm limit 54, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt54) < 10 | abs(table$median_dist_real_lt54) > 54,]$label <- NA
# table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_random) < 15 | abs(table$median_dist_random) > 100,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt54, median_dist_random_lt54)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit54-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt54) %>% arrange(median_dist_real_lt54)
# >> RBFOX2|SUPV3L1 -10



# Beeswarm limit 100, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt100) < 10 | abs(table$median_dist_real_lt100) > 100,]$label <- NA
# table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_random) < 15 | abs(table$median_dist_random) > 100,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt100, median_dist_random_lt100)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit100-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt100) %>% arrange(median_dist_real_lt100)
# >> RBFOX2|SUPV3L1 -11.5, RBFOX2|DDX21 +11


# Beeswarm limit 117, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt117) < 10 | abs(table$median_dist_real_lt117) > 117,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt117, median_dist_random_lt117)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  # coord_cartesian(xlim=c(-117, 117)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit117-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt117) %>% arrange(median_dist_real_lt117)
# >> RBFOX2|SUPV3L1 -11.5, RBFOX2|DDX21 +15


# Beeswarm limit 250, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt250) < 10 | abs(table$median_dist_real_lt250) > 250,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt250, median_dist_random_lt250)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  # coord_cartesian(xlim=c(-250, 250)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit250-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt250) %>% arrange(median_dist_real_lt250)
# >> RBM15|SRSF9 -13.5, RBFOX2|SUPV3L1 -12, RBFOX2|DDX21 +13


# Beeswarm limit 400, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt400) < 10 | abs(table$median_dist_real_lt400) > 400,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt400, median_dist_random_lt400)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  # coord_cartesian(xlim=c(-400, 400)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit400-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt400) %>% arrange(median_dist_real_lt400)
# >> RBM15|SRSF9 -17, RBFOX2|SUPV3L1 -11, RBFOX2|DDX21 +13


# Beeswarm limit 1000, filtered to significant pairs
table %<>% unite(label, c(symbol1, symbol2), sep="|", remove=F)
table
table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_real_lt1000) < 10 | abs(table$median_dist_real_lt1000) > 1000,]$label <- NA
# table[table$resampling_p >= 0.05 | table$resampling_wilcox >= 0.05 | abs(table$median_dist_random) < 15 | abs(table$median_dist_random) > 100,]$label <- NA

table %>% 
  pivot_longer(cols=c(median_dist_real_lt1000, median_dist_random_lt1000)) %>%
  filter(resampling_p < 0.05 & resampling_wilcox < 0.05) %>%
  ggplot(aes(x=value, y=name)) + 
  # geom_violin() +
  geom_beeswarm(groupOnX=F) +
  geom_text_repel(aes(label=label), nudge_y=0.5, force = 20) +
  # coord_cartesian(xlim=c(-54, 54)) +
  # coord_cartesian(xlim=c(-100, 100)) +
  # coord_cartesian(xlim=c(-1000, 1000)) +
  theme_minimal() + 
  xlab("Median distance [nt]") +
  ylab("") +
  ggsave("output-fit_vs_random-rbpome-eclip_encode_12-SUM_IS-54-signed-delta-beeswarm-limit1000-filtered.pdf", width=183, height=91.5, units="mm", device = cairo_pdf)
 table %>% filter(!is.na(label)) %>% select(label, median_dist_real_lt1000) %>% arrange(median_dist_real_lt1000)
# >> RBM15|SRSF9 -25.5, RBFOX2|SUPV3L1 -14, PCBP1|PTBP1 -10, NONO|RBM5 +15, TRA2A|SRSF9 +97


