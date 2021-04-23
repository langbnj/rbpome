rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(ggplot2)
library(dplyr)
# library(reshape2)
# library(tidyr)
# library(viridis)
# # library(plyr)
# library(scales)
# library(stringr)


setwd("~/Documents/Projects/RBPome/mindist_cdf")
# q <- read.table("output-fit_vs_random-txt-eclip_encode_12-SUM_IS.txt", quote="", sep="\t", header=T)
q <- read.table("output-fit_vs_random-txt-rbpome-eclip_encode_12-SUM_IS.txt", quote="", sep="\t", header=T)

q$pair <- paste0(q$symbol1,"|",q$symbol2)
q$randompair <- as.factor(q$randompair)
q$random <- as.factor(q$random)
q$pair <- as.factor(q$pair)
head(q)
str(q)

# Flag positive pairs
# pos <- Query("SELECT DISTINCT symbol1, symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND symbol1!=symbol2 AND bg_alldirect=1 AND homodimer=0 AND eclip1=1 AND eclip2=1 AND hc=1")
# Require at least two "direct" studies in BioGRID confirming the interaction (this leads to 5 positive pairs)
pos <- Query("SELECT DISTINCT r.symbol1, r.symbol2 FROM rbpome r, (SELECT gene1, gene2, COUNT(DISTINCT pmid) AS pmids FROM biogrid WHERE direct=1 GROUP BY gene1, gene2) g WHERE r.scoretype='SUM_IS' AND r.symbol1!=r.symbol2 AND r.bg_alldirect=1 AND r.homodimer=0 AND r.eclip1=1 AND r.eclip2=1 AND r.hc=1 AND g.gene1=r.symbol1 AND g.gene2=r.symbol2 AND g.pmids>=2")
posflip <- pos
pos$pair <- paste0(pos$symbol1,"|",pos$symbol2)
posflip$pair <- paste0(posflip$symbol2,"|",posflip$symbol1)
posflip$tmp <- posflip$symbol1
posflip$symbol1 <- posflip$symbol2
posflip$symbol2 <- posflip$tmp
posflip$tmp <- NULL
pos
posflip
pos <- rbind(pos, posflip)
pos
q$set <- NA
q[q$randompair==1,]$set <- "background_resamples"
q[q$randompair==0,]$set <- "screen_hits"
q[q$pair %in% pos$pair,]$set <- "positives"
# >> The positives turn up in the random set a lot!
# unique(subset(q, set=="positives")$pair)


# Flag negative pairs
q[q$pair == "FASTKD2|SLBP",]$set <- "negatives"
q[q$pair == "SLBP|FASTKD2",]$set <- "negatives"


# Add on a "superpositive" set (FMR1|FXR2, FXR1|FXR2, U2AF1|U2AF2, and no inverse pairs since they don't work as well)
# (Funnily, FXR1|FXR2 is NOT in the positives set - it doesn't have more than one direct study confirming it)
qs <- subset(q, set=="positives" & pair %in% c("FMR1|FXR2", "FXR1|FXR2", "U2AF1|U2AF2"))
# qs$symbol1 <- as.factor(as.character(qs$symbol1))
# qs$symbol2 <- as.factor(as.character(qs$symbol2))
# qs$pair <- as.factor(as.character(qs$pair))
qs$set <- "superpositives"
summary(qs)
q <- rbind(qs, q)

q$set <- as.factor(q$set)
summary(q)
summary(subset(q, set=="superpositives"))
summary(subset(q, set=="positives"))
summary(subset(q, set=="screen_hits"))
summary(subset(q, set=="background_resamples"))
summary(subset(q, set=="negatives"))

# Back up the original data frame
qo <- q


# Get list of RBPs used in these plots
q1 <- subset(q, random==0 & set=="superpositives") %>% select(symbol1, symbol2) %>% as_tibble %>% unique
q1
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 4 RBPs
q1 <- subset(q, random==0 & set=="positives") %>% select(symbol1, symbol2) %>% as_tibble %>% unique
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 9 RBPs
q1 <- subset(q, random==0 & set=="screen_hits") %>% select(symbol1, symbol2) %>% as_tibble %>% unique
q1
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 45 RBPs, 114 pairs
q1 <- subset(q, random==0 & (set=="screen_hits" | set=="positives")) %>% select(symbol1, symbol2) %>% as_tibble %>% unique
q1
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 47 RBPs, 124 pairs
q1 <- subset(q, random==0 & (set=="screen_hits" | set=="positives" | set=="superpositives")) %>% select(symbol1, symbol2) %>% as_tibble %>% unique
q1
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 47 RBPs, 124 pairs
q1 <- subset(q, random==0 & set=="background_resamples") %>% select(symbol1, symbol2) %>% as_tibble %>% unique
q1
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 52 RBPs, 278 pairs
q1 <- subset(q, random==0) %>% select(symbol1, symbol2) %>% as_tibble %>% unique
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique
# 52 RBPs
# Print list
bind_rows(q1 %>% transmute(symbol = symbol1), q1 %>% transmute(symbol = symbol2)) %>% unique %>% pull %>% cat(sep = "\n")


# myset <- "negatives"
# myset <- "screen_hits"
for (myset in c("negatives", "superpositives", "positives", "screen_hits", "background_resamples")) {
# for (myset in c("negatives")) {
# for (myset in c("superpositives", "positives")) {
  q <- qo
  
  real <- subset(q, random==0 & set==myset)
  random <- subset(q, random==1 & set==myset)
  
  # For the large sets:
  if (myset %in% c("screen_hits", "background_resamples")) {
    # Subsample the random set
    set.seed(2020)
    random_sub <- sample_n(random, nrow(real), replace=T)
  } else {
    # Change nothing
    random_sub <- random
  }
  
  # Recombine
  q <- rbind(real, random_sub)
  # str(q)

  
  
  
  # # All
  ggplot(real, aes(x=mindist)) + geom_histogram(binwidth=1) + xlim(0, 250)
  # ggplot(real, aes(x=mindist)) + stat_ecdf() + theme_minimal()
  # ggplot(real, aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + theme_minimal()
  # ggplot(real, aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + theme_minimal()
  # ggplot(real, aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 100)) + theme_minimal()
  # 
  # ggplot(subset(real, mindist<=1000), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggplot(subset(real, mindist<=1000), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggplot(subset(real, mindist<=1000), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 100)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  
  
  # Best plot - combined real & random
  # ≤1000 only
  ggplot(subset(q, mindist<=1000), aes(x=mindist, colour=random)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below1000.png"), width=5, height=5)
  # ≤400 only
  ggplot(subset(q, mindist<=400), aes(x=mindist, colour=random)) + stat_ecdf() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below400.png"), width=5, height=5)
  # ≤250 only
  ggplot(subset(q, mindist<=250), aes(x=mindist, colour=random)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below250.png"), width=5, height=5)
  # all data
  ggplot(q, aes(x=mindist, colour=random)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata.png"), width=5, height=5)
  
  
  
  
  # Difference: alldata
  realfun <- ecdf(subset(q, random==0)$mindist)
  randfun <- ecdf(subset(q, random==1)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  
  ggplot(dif, aes(x, y)) + geom_line() + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-difference.png"), width=5, height=5)
  ggplot(dif, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-difference-zoom1000.png"), width=5, height=5)
  
  
  
  # Difference: below1000
  realfun <- ecdf(subset(q, random==0 & mindist<=1000)$mindist)
  randfun <- ecdf(subset(q, random==1 & mindist<=1000)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  
  ggplot(dif, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below1000-difference.png"), width=5, height=5)
  
  
  
  # Difference: below400
  realfun <- ecdf(subset(q, random==0 & mindist<=400)$mindist)
  randfun <- ecdf(subset(q, random==1 & mindist<=400)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  
  ggplot(dif, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below400-difference.png"), width=5, height=5)
  
  
  
  # Difference: below250
  realfun <- ecdf(subset(q, random==0 & mindist<=250)$mindist)
  randfun <- ecdf(subset(q, random==1 & mindist<=250)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  
  ggplot(dif, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below250-difference.png"), width=5, height=5)
  
  
  
  # "Factor" (real divided by random): doesn't look great
  # fac <- matrix(ncol=2, nrow=max(q$mindist)+1)
  # 
  # fac[,1] <- 0:max(q$mindist)
  # fac[,2] <- apply(fac, 1, function(x){realfun(x[1]) / randfun(x[1])})
  # 
  # fac <- as.data.frame(fac)
  # colnames(fac) <- c("x", "y")
  # 
  # ggplot(fac, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-combined-below1000-factor.png"), width=5, height=5)
  # ggplot(fac, aes(x, y)) + geom_line() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-combined-below1000-factor.png"), width=5, height=5)
  
  
  
  
  
  # # Best plot - real
  # # ≤1000 only
  # ggplot(subset(real, mindist<=1000), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below1000-real.png"), width=5, height=5)
  # # ≤400 only
  # ggplot(subset(real, mindist<=400), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below400-real.png"), width=5, height=5)
  # # ≤250 only
  # ggplot(subset(real, mindist<=250), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below250-real.png"), width=5, height=5)
  # # all data
  # ggplot(real, aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-real.png"), width=5, height=5)
  # 
  # 
  # # Best plot - random
  # # ≤1000 only
  # ggplot(subset(random, mindist<=1000), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below1000-random.png"), width=5, height=5)
  # # ≤400 only
  # ggplot(subset(random, mindist<=400), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below400-random.png"), width=5, height=5)
  # # ≤250 only
  # ggplot(subset(random, mindist<=250), aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-below250-random.png"), width=5, height=5)
  # # all data
  # ggplot(random, aes(x=mindist)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  # ggsave(paste0("output-",myset,"-mindist_cdf-pooled-alldata-random.png"), width=5, height=5)
  
  
  # Split by pair - real
  # ≤1000 only
  ggplot(subset(real, mindist<=1000), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below1000-real.png"), width=5, height=5)
  # ≤400 only
  ggplot(subset(real, mindist<=400), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below400-real.png"), width=5, height=5)
  # ≤250 only
  ggplot(subset(real, mindist<=250), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-real.png"), width=5, height=5)
  # all data
  ggplot(real, aes(x=mindist, colour=pair)) + stat_ecdf() + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-real.png"), width=5, height=5)
  # all data - zoom1000
  ggplot(real, aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-real-zoom1000.png"), width=5, height=5)
  
  # all data - zoom1000 with legend, so I can find out what the two cool random pairs are:
  ggplot(real, aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-real-zoom1000-with_legend.png"), width=5+(ceiling(length(unique(real$pair))/18))*1.75, height=5)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-real-zoom1000-with_legend.pdf"), width=5+(ceiling(length(unique(real$pair))/18))*1.75, height=5)
  
  # ≤1000 only with legend
  ggplot(subset(real, mindist<=1000), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below1000-real-with_legend.png"), width=5+(ceiling(length(unique(subset(real, mindist<=1000)$pair))/18))*1.75, height=5)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below1000-real-with_legend.pdf"), width=5+(ceiling(length(unique(subset(real, mindist<=1000)$pair))/18))*1.75, height=5)

  
  
    
  # Split by pair - random
  # ≤1000 only
  ggplot(subset(random_sub, mindist<=1000), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below1000-random.png"), width=5, height=5)
  # ≤400 only
  ggplot(subset(random_sub, mindist<=400), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 400)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below400-random.png"), width=5, height=5)
  # ≤250 only
  ggplot(subset(random_sub, mindist<=250), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-random.png"), width=5, height=5)
  # all data
  ggplot(random_sub, aes(x=mindist, colour=pair)) + stat_ecdf() + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-random.png"), width=5, height=5)
  # all data - zoom1000
  ggplot(random_sub, aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 1000)) + theme_minimal() + guides(colour=F)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-alldata-random-zoom1000.png"), width=5, height=5)
  
  # ≤250 only with legend, so I can find out what the two high random pairs are:
  ggplot(subset(random_sub, mindist<=250), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-random-with_legend_random_sub.png"), width=5+(ceiling(length(unique(random_sub$pair))/18))*1.75, height=5)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-random-with_legend_random_sub.pdf"), width=5+(ceiling(length(unique(random_sub$pair))/18))*1.75, height=5)
  ggplot(subset(random, mindist<=250), aes(x=mindist, colour=pair)) + stat_ecdf() + coord_cartesian(xlim=c(0, 250)) + geom_vline(xintercept = 50, linetype="dashed") + theme_minimal()
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-random-with_legend_random.png"), width=5+(ceiling(length(unique(random$pair))/18))*1.75, height=5)
  ggsave(paste0("output-",myset,"-mindist_cdf-pairs-below250-random-with_legend_random.pdf"), width=5+(ceiling(length(unique(random$pair))/18))*1.75, height=5)
  

  
  
  
}






# Make two plots that contain superpositives, screen_hits, and background_resamples:
q <- qo

# Split into real & random, subsample random, then recombine
# Real
real <- subset(q, random==0)
# Random
if (exists("random")) { rm(random) }
for (myset in unique(q$set)) {
  cat(paste0(" >> ",myset,"\n"))
  
  # Subsample the random set
  set.seed(2020)
  random_sub <- sample_n(subset(q, random==1 & set==myset), nrow(subset(q, random==0 & set==myset)), replace=T)
  
  if (exists("random") == T) {
    random <- rbind(random, random_sub)
  } else {
    random <- random_sub
  }
}
# Recombine
summary(real)
summary(random)
q <- rbind(real, random)
# str(q)



cols <- c("negatives 0" = "#46ba1f",
          "negatives 1" = "#b5e3a5",
          "superpositives 0" = "#1e3d59",
          "superpositives 1" = "#a3afbb",
          "positives 0" = "#1e3d59",
          "positives 1" = "#a3afbb",
          "screen_hits 0" = "#ff7f00",
          "screen_hits 1" = "#ffcb96",
          "background_resamples 0" = "#bec1c0",
          "background_resamples 1" = "#e4e6e5")



# install.packages("plotROC")
# library(plotROC)
# + style_roc() + scale_x_continuous()
# >> Doesn't look right because of the x-axis, which obviously isn't 0-1 here.
#   >> Can fix this by adding + scale_x_continuous(), but weirdly it sets ylab() to "True positive fraction" and doesn't let me change it anymore.




# Difference: below1000
# Get graphs
if (exists("alldif")) { rm(alldif) }
myset <- "positives"
for (myset in c("superpositives", "positives", "screen_hits", "background_resamples")) {
  cat(paste0(" >> ",myset,"\n"))
  # try(print(str(alldif)))
  realfun <- ecdf(subset(q, set==myset & random==0 & mindist<=1000)$mindist)
  randfun <- ecdf(subset(q, set==myset & random==1 & mindist<=1000)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  dif$set <- myset
  
  if (exists("alldif") == T) {
    alldif <- rbind(alldif, dif)
  } else {
    alldif <- dif
  }
}
dif <- alldif
dif$set <- as.factor(dif$set)
str(dif)

# Get maxima (positives) (≤1000)
# subset(dif, set %in% c("positives", "screen_hits", "background_resamples"))
# str(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")))
# p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
#   geom_line() + 
#   geom_vline(xintercept = 50, linetype="dashed") + 
#   scale_colour_manual(values = cols) + 
#   coord_cartesian(xlim=c(0, 1000), expand=0) + 
#   theme_minimal() +
#   xlab("Δ binding site distance [nt]") +
#   ylab("Probability density") +
#   guides(colour=F)
# p

# Find maximum: superpositives (≤1000)
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max_superpositives
thresh_superpositives <- dif[dif$set=="superpositives" & dif$y==max_superpositives,]$x
thresh_superpositives

# Find maximum: positives (≤1000)
max_positives <- max(subset(dif, set=="positives")$y)
max_positives
thresh_positives <- dif[dif$set=="positives" & dif$y==max_positives,]$x
thresh_positives

# Find maximum: screen_hits (≤1000)
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max_screen_hits
thresh_screen_hits <- dif[dif$set=="screen_hits" & dif$y==max_screen_hits,]$x
thresh_screen_hits

# Find maximum: background_resamples
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max_background_resamples
thresh_background_resamples <- dif[dif$set=="background_resamples" & dif$y==max_background_resamples,]$x
thresh_background_resamples


# 90/95% of maximum, upper bound:
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max95_superpositives <- max_superpositives * 0.95
thresh95_superpositives <- max(dif[dif$set=="superpositives" & dif$y>=max95_superpositives,]$x)
thresh95_superpositives
max90_superpositives <- max_superpositives * 0.9
thresh90_superpositives <- max(dif[dif$set=="superpositives" & dif$y>=max90_superpositives,]$x)
thresh90_superpositives

# 90/95% of maximum, upper bound:
max_positives <- max(subset(dif, set=="positives")$y)
max95_positives <- max_positives * 0.95
thresh95_positives <- max(dif[dif$set=="positives" & dif$y>=max95_positives,]$x)
thresh95_positives
max90_positives <- max_positives * 0.9
thresh90_positives <- max(dif[dif$set=="positives" & dif$y>=max90_positives,]$x)
thresh90_positives

# 90/95% of maximum, upper bound:
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max95_screen_hits <- max_screen_hits * 0.95
thresh95_screen_hits <- max(dif[dif$set=="screen_hits" & dif$y>=max95_screen_hits,]$x)
thresh95_screen_hits
max90_screen_hits <- max_screen_hits * 0.9
thresh90_screen_hits <- max(dif[dif$set=="screen_hits" & dif$y>=max90_screen_hits,]$x)
thresh90_screen_hits

# 90/95% of maximum, upper bound:
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max95_background_resamples <- max_background_resamples * 0.95
thresh95_background_resamples <- max(dif[dif$set=="background_resamples" & dif$y>=max95_background_resamples,]$x)
thresh95_background_resamples
max90_background_resamples <- max_background_resamples * 0.9
thresh90_background_resamples <- max(dif[dif$set=="background_resamples" & dif$y>=max90_background_resamples,]$x)
thresh90_background_resamples


# # Use 95%
# thresh_superpositives <- thresh95_superpositives
# thresh_positives <- thresh95_positives
# thresh_screen_hits <- thresh95_screen_hits
# thresh_background_resamples <- thresh95_background_resamples




# Plot
ggplot(subset(dif, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  # geom_vline(xintercept = thresh95_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_positives)," nt"), colour=cols["positives 0"]) +
  # annotate("text", x=thresh95_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh95_positives)," nt"), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_screen_hits)," nt"), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_background_resamples)," nt"), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
ggsave(paste0("output-all-mindist_cdf-pooled-below1000-difference-superpositives.png"), width=5, height=5)

p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  # geom_vline(xintercept = thresh95_positives, linetype="dotted", colour=cols["positives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_positives)," nt"), colour=cols["positives 0"]) +
  # annotate("text", x=thresh95_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh95_positives)," nt"), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_screen_hits)," nt"), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_background_resamples)," nt"), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Δ binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below1000-difference-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4e.pdf"), width=91.5, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e.png"), width=91.5, height=91.5,units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4e-legend.pdf"), width=183, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-legend.png"), width=183, height=91.5,units="mm")






nrow(subset(q, mindist<=1000 & set=="positives"))
nrow(subset(q, mindist<=10000 & set=="positives"))
nrow(subset(q, mindist<=100000 & set=="positives"))
nrow(subset(q, set=="positives"))

# Best plot - combined real & random
# ≤1000 only
ggplot(subset(q, mindist<=1000 & set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  # geom_vline(xintercept = thresh95_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_positives)," nt"), colour=cols["positives 0"]) +
  # annotate("text", x=thresh95_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh95_positives)," nt"), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_screen_hits)," nt"), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_background_resamples)," nt"), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  # theme_minimal() + style_roc() + scale_x_continuous() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F) +
  ggsave(paste0("output-all-mindist_cdf-pooled-below1000-superpositives.png"), width=5, height=5)

p <- ggplot(subset(q, mindist<=1000 & set %in% c("positives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  # geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  # geom_vline(xintercept = thresh95_positives, linetype="dotted", colour=cols["positives 0"]) + 
  # geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  # geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  # annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_positives)," nt"), colour=cols["positives 0"]) +
  # annotate("text", x=thresh95_positives, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh95_positives)," nt"), colour=cols["positives 0"]) +
  # annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_screen_hits)," nt"), colour=cols["screen_hits 0"]) +
  # annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=paste0(round(thresh_background_resamples)," nt"), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below1000-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4d.pdf"), width=91.5, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d.png"), width=91.5, height=91.5, units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4d-legend.pdf"), width=183, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-legend.png"), width=183, height=91.5, units="mm")















# 4d ≤400

# Best plot - combined real & random
# ≤400 only
ggplot(subset(q, mindist<=400 & set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 400), expand=0) + 
  theme_minimal() +
  # theme_minimal() + style_roc() + scale_x_continuous() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F) +
  ggsave(paste0("output-all-mindist_cdf-pooled-below400-superpositives.png"), width=5, height=5)

p <- ggplot(subset(q, mindist<=400 & set %in% c("positives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 400), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below400-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4d-400.pdf"), width=91.5, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-400.png"), width=91.5, height=91.5, units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4d-400-legend.pdf"), width=183, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-400-legend.png"), width=183, height=91.5, units="mm")















# 4e ≤400
# Difference: below400
# Get graphs
if (exists("alldif")) { rm(alldif) }
for (myset in c("superpositives", "positives", "screen_hits", "background_resamples")) {
  cat(paste0(" >> ",myset,"\n"))
  # try(print(str(alldif)))
  realfun <- ecdf(subset(q, set==myset & random==0 & mindist<=400)$mindist)
  randfun <- ecdf(subset(q, set==myset & random==1 & mindist<=400)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  dif$set <- myset
  
  if (exists("alldif") == T) {
    alldif <- rbind(alldif, dif)
  } else {
    alldif <- dif
  }
}
dif <- alldif
dif$set <- as.factor(dif$set)
str(dif)

# Find maximum: superpositives (≤1000)
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max_superpositives
thresh_superpositives <- dif[dif$set=="superpositives" & dif$y==max_superpositives,]$x
thresh_superpositives

# Find maximum: positives (≤1000)
max_positives <- max(subset(dif, set=="positives")$y)
max_positives
thresh_positives <- dif[dif$set=="positives" & dif$y==max_positives,]$x
thresh_positives

# Find maximum: screen_hits (≤1000)
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max_screen_hits
thresh_screen_hits <- dif[dif$set=="screen_hits" & dif$y==max_screen_hits,]$x
thresh_screen_hits

# Find maximum: background_resamples
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max_background_resamples
thresh_background_resamples <- dif[dif$set=="background_resamples" & dif$y==max_background_resamples,]$x
thresh_background_resamples


# Plot
ggplot(subset(dif, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_superpositives, y=Inf, hjust=0, vjust=1, label=round(thresh_superpositives), colour=cols["superpositives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 400), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
ggsave(paste0("output-all-mindist_cdf-pooled-below400-difference-superpositives.png"), width=5, height=5)

p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=round(thresh_positives), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 400), expand=0) + 
  theme_minimal() +
  xlab("Δ binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below400-difference-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4e-400.pdf"), width=91.5, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-400.png"), width=91.5, height=91.5,units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4e-400-legend.pdf"), width=183, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-400-legend.png"), width=183, height=91.5,units="mm")




















# 4d ≤10,000

# Best plot - combined real & random
# ≤10,000 only
ggplot(subset(q, mindist<=10000 & set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 10000), expand=0) + 
  theme_minimal() +
  # theme_minimal() + style_roc() + scale_x_continuous() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F) +
  ggsave(paste0("output-all-mindist_cdf-pooled-below10000-superpositives.png"), width=5, height=5)

p <- ggplot(subset(q, mindist<=10000 & set %in% c("positives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 10000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below10000-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4d-10000.pdf"), width=91.5, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-10000.png"), width=91.5, height=91.5, units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4d-10000-legend.pdf"), width=183, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-10000-legend.png"), width=183, height=91.5, units="mm")















# 4e ≤10,000
# Difference: below10000
# Get graphs
if (exists("alldif")) { rm(alldif) }
for (myset in c("superpositives", "positives", "screen_hits", "background_resamples")) {
  cat(paste0(" >> ",myset,"\n"))
  # try(print(str(alldif)))
  realfun <- ecdf(subset(q, set==myset & random==0 & mindist<=10000)$mindist)
  randfun <- ecdf(subset(q, set==myset & random==1 & mindist<=10000)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  dif$set <- myset
  
  if (exists("alldif") == T) {
    alldif <- rbind(alldif, dif)
  } else {
    alldif <- dif
  }
}
dif <- alldif
dif$set <- as.factor(dif$set)
str(dif)

# Find maximum: superpositives (≤1000)
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max_superpositives
thresh_superpositives <- dif[dif$set=="superpositives" & dif$y==max_superpositives,]$x
thresh_superpositives

# Find maximum: positives (≤1000)
max_positives <- max(subset(dif, set=="positives")$y)
max_positives
thresh_positives <- dif[dif$set=="positives" & dif$y==max_positives,]$x
thresh_positives

# Find maximum: screen_hits (≤1000)
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max_screen_hits
thresh_screen_hits <- dif[dif$set=="screen_hits" & dif$y==max_screen_hits,]$x
thresh_screen_hits

# Find maximum: background_resamples
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max_background_resamples
thresh_background_resamples <- dif[dif$set=="background_resamples" & dif$y==max_background_resamples,]$x
thresh_background_resamples


# Plot
ggplot(subset(dif, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_superpositives, y=Inf, hjust=0, vjust=1, label=round(thresh_superpositives), colour=cols["superpositives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 10000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
ggsave(paste0("output-all-mindist_cdf-pooled-below10000-difference-superpositives.png"), width=5, height=5)

p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=round(thresh_positives), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 10000), expand=0) + 
  theme_minimal() +
  xlab("Δ binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below10000-difference-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4e-10000.pdf"), width=91.5, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-10000.png"), width=91.5, height=91.5,units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4e-10000-legend.pdf"), width=183, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-10000-legend.png"), width=183, height=91.5,units="mm")




















# 4d ≤100,000

# Best plot - combined real & random
# ≤100,000 only
ggplot(subset(q, mindist<=100000 & set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 100000), expand=0) + 
  theme_minimal() +
  # theme_minimal() + style_roc() + scale_x_continuous() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F) +
  ggsave(paste0("output-all-mindist_cdf-pooled-below100000-superpositives.png"), width=5, height=5)

p <- ggplot(subset(q, mindist<=100000 & set %in% c("positives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 100000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below100000-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4d-100000.pdf"), width=91.5, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-100000.png"), width=91.5, height=91.5, units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4d-100000-legend.pdf"), width=183, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-100000-legend.png"), width=183, height=91.5, units="mm")















# 4e ≤100,000
# Difference: below100000
# Get graphs
if (exists("alldif")) { rm(alldif) }
for (myset in c("superpositives", "positives", "screen_hits", "background_resamples")) {
  cat(paste0(" >> ",myset,"\n"))
  # try(print(str(alldif)))
  realfun <- ecdf(subset(q, set==myset & random==0 & mindist<=100000)$mindist)
  randfun <- ecdf(subset(q, set==myset & random==1 & mindist<=100000)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  dif$set <- myset
  
  if (exists("alldif") == T) {
    alldif <- rbind(alldif, dif)
  } else {
    alldif <- dif
  }
}
dif <- alldif
dif$set <- as.factor(dif$set)
str(dif)

# Find maximum: superpositives (≤1000)
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max_superpositives
thresh_superpositives <- dif[dif$set=="superpositives" & dif$y==max_superpositives,]$x
thresh_superpositives

# Find maximum: positives (≤1000)
max_positives <- max(subset(dif, set=="positives")$y)
max_positives
thresh_positives <- dif[dif$set=="positives" & dif$y==max_positives,]$x
thresh_positives

# Find maximum: screen_hits (≤1000)
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max_screen_hits
thresh_screen_hits <- dif[dif$set=="screen_hits" & dif$y==max_screen_hits,]$x
thresh_screen_hits

# Find maximum: background_resamples
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max_background_resamples
thresh_background_resamples <- dif[dif$set=="background_resamples" & dif$y==max_background_resamples,]$x
thresh_background_resamples


# Plot
ggplot(subset(dif, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_superpositives, y=Inf, hjust=0, vjust=1, label=round(thresh_superpositives), colour=cols["superpositives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 100000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
ggsave(paste0("output-all-mindist_cdf-pooled-below100000-difference-superpositives.png"), width=5, height=5)

p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=round(thresh_positives), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  coord_cartesian(xlim=c(0, 100000), expand=0) + 
  theme_minimal() +
  xlab("Δ binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-below100000-difference-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4e-100000.pdf"), width=91.5, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-100000.png"), width=91.5, height=91.5,units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4e-100000-legend.pdf"), width=183, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-100000-legend.png"), width=183, height=91.5,units="mm")














# 4d Unlimited
# Best plot - combined real & random
ggplot(subset(q, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  # coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  # theme_minimal() + style_roc() + scale_x_continuous() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F) +
  ggsave(paste0("output-all-mindist_cdf-pooled-all-superpositives.png"), width=5, height=5)

p <- ggplot(subset(q, set %in% c("positives", "screen_hits", "background_resamples")), aes(x=mindist, colour=paste0(set,' ',random), linetype=random)) +
  stat_ecdf() + 
  geom_vline(xintercept = 50, linetype="dashed") + 
  scale_colour_manual(values = cols) + 
  # coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F, linetype=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-all-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4d-unlimited.pdf"), width=91.5, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-unlimited.png"), width=91.5, height=91.5, units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4d-unlimited-legend.pdf"), width=183, height=91.5, units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4d-unlimited-legend.png"), width=183, height=91.5, units="mm")













# 4e Unlimited
# Get graphs
if (exists("alldif")) { rm(alldif) }
for (myset in c("superpositives", "positives", "screen_hits", "background_resamples")) {
  cat(paste0(" >> ",myset,"\n"))
  # try(print(str(alldif)))
  # realfun <- ecdf(subset(q, set==myset & random==0 & mindist<=1000)$mindist)
  # randfun <- ecdf(subset(q, set==myset & random==1 & mindist<=1000)$mindist)
  realfun <- ecdf(subset(q, set==myset & random==0)$mindist)
  randfun <- ecdf(subset(q, set==myset & random==1)$mindist)
  
  dif <- matrix(ncol=2, nrow=max(q$mindist)+1)
  dif[,1] <- 0:max(q$mindist)
  dif[,2] <- apply(dif, 1, function(x){realfun(x[1]) - randfun(x[1])})
  dif <- as.data.frame(dif)
  colnames(dif) <- c("x", "y")
  dif$set <- myset
  
  if (exists("alldif") == T) {
    alldif <- rbind(alldif, dif)
  } else {
    alldif <- dif
  }
}
dif <- alldif
dif$set <- as.factor(dif$set)
str(dif)

# Find maximum: superpositives (≤1000)
max_superpositives <- max(subset(dif, set=="superpositives")$y)
max_superpositives
thresh_superpositives <- dif[dif$set=="superpositives" & dif$y==max_superpositives,]$x
thresh_superpositives

# Find maximum: positives (≤1000)
max_positives <- max(subset(dif, set=="positives")$y)
max_positives
thresh_positives <- dif[dif$set=="positives" & dif$y==max_positives,]$x
thresh_positives

# Find maximum: screen_hits (≤1000)
max_screen_hits <- max(subset(dif, set=="screen_hits")$y)
max_screen_hits
thresh_screen_hits <- dif[dif$set=="screen_hits" & dif$y==max_screen_hits,]$x
thresh_screen_hits

# Find maximum: background_resamples
max_background_resamples <- max(subset(dif, set=="background_resamples")$y)
max_background_resamples
thresh_background_resamples <- dif[dif$set=="background_resamples" & dif$y==max_background_resamples,]$x
thresh_background_resamples

# Plot
ggplot(subset(dif, set %in% c("superpositives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_superpositives, linetype="dotted", colour=cols["superpositives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_superpositives, y=Inf, hjust=0, vjust=1, label=round(thresh_superpositives), colour=cols["superpositives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  # coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
ggsave(paste0("output-all-mindist_cdf-pooled-all-difference-superpositives.png"), width=5, height=5)

p <- ggplot(subset(dif, set %in% c("positives", "screen_hits", "background_resamples")), aes(x, y, colour=paste0(set,' 0'))) + 
  geom_line() + 
  geom_vline(xintercept = thresh_positives, linetype="dotted", colour=cols["positives 0"]) + 
  geom_vline(xintercept = thresh_screen_hits, linetype="dotted", colour=cols["screen_hits 0"]) + 
  geom_vline(xintercept = thresh_background_resamples, linetype="dotted", colour=cols["background_resamples 0"]) + 
  annotate("text", x=thresh_positives, y=Inf, hjust=0, vjust=1, label=round(thresh_positives), colour=cols["positives 0"]) +
  annotate("text", x=thresh_screen_hits, y=Inf, hjust=0, vjust=1, label=round(thresh_screen_hits), colour=cols["screen_hits 0"]) +
  annotate("text", x=thresh_background_resamples, y=Inf, hjust=0, vjust=1, label=round(thresh_background_resamples), colour=cols["background_resamples 0"]) +
  scale_colour_manual(values = cols) + 
  # coord_cartesian(xlim=c(0, 1000), expand=0) + 
  theme_minimal() +
  xlab("Δ binding site distance [nt]") +
  ylab("Probability density") +
  guides(colour=F)
p
ggsave(paste0("output-all-mindist_cdf-pooled-all-difference-positives.png"), width=5, height=5)

# Final for Figure 4
ggsave(paste0("output-figure4e-unlimited.pdf"), width=91.5, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-unlimited.png"), width=91.5, height=91.5,units="mm")

# With legend
p <- p + guides(colour="legend", linetype="legend")
p
ggsave(paste0("output-figure4e-unlimited-legend.pdf"), width=183, height=91.5,units="mm", device=cairo_pdf)
ggsave(paste0("output-figure4e-unlimited-legend.png"), width=183, height=91.5,units="mm")









# Plot random (originally dashed) lines only (to show how close we get to the diagonal with different distance cutoffs)
# 4d random only
q <- tibble(qo)
subset(q, set %in% c("positives", "screen_hits", "background_resamples"))
q %>% 
  # filter(set %in% c("positives", "screen_hits", "background_resamples")) %>% 
  filter(set == "background_resamples") %>% 
  filter(random == 1) -> qtmp
qtmp
bind_rows(qtmp %>% filter(mindist<=50) %>% mutate(mindist=mindist/50, maxdist=50),
          qtmp %>% filter(mindist<=100) %>% mutate(mindist=mindist/100, maxdist=100),
          qtmp %>% filter(mindist<=200) %>% mutate(mindist=mindist/200, maxdist=200),
          qtmp %>% filter(mindist<=400) %>% mutate(mindist=mindist/400, maxdist=400),
          qtmp %>% filter(mindist<=700) %>% mutate(mindist=mindist/700, maxdist=700),
          qtmp %>% filter(mindist<=1000) %>% mutate(mindist=mindist/1000, maxdist=1000),
          qtmp %>% filter(mindist<=2000) %>% mutate(mindist=mindist/2000, maxdist=2000),
          qtmp %>% filter(mindist<=4000) %>% mutate(mindist=mindist/4000, maxdist=4000),
          qtmp %>% filter(mindist<=7000) %>% mutate(mindist=mindist/7000, maxdist=7000),
          qtmp %>% filter(mindist<=10000) %>% mutate(mindist=mindist/10000, maxdist=10000),
          qtmp %>% filter(mindist<=20000) %>% mutate(mindist=mindist/20000, maxdist=20000),
          qtmp %>% filter(mindist<=100000) %>% mutate(mindist=mindist/100000, maxdist=100000)) %>%
  ggplot(aes(x=mindist, colour=as.factor(maxdist), linetype=(maxdist!=1000))) +
  stat_ecdf() +
  # stat_ecdf(geom="line") +
  # geom_vline(xintercept = 50, linetype="dashed") + 
  # geom_abline(linetype="dashed") +
  geom_abline() +
  # scale_colour_manual(values = cols) + 
  scale_colour_discrete(name = "Maximum distance") +
  # scale_linetype_manual(values = c("dotted", "solid")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  # scale_linetype_manual(values = c("solid", "dashed")) +
  coord_cartesian(xlim=c(0, 1), expand=0) +
  theme_minimal() +
  xlab("Binding site distance (relative to maximum distance)") +
  ylab("Cumulative probability density") +
  guides(linetype=F) +
  # ggsave(paste0("output-figure-R1B-maxdists.pdf"), width=137.25, height=91.5, units="mm") +
  # ggsave(paste0("output-figure-R1B-maxdists-cairo.pdf"), width=137.25, height=91.5, units="mm", device = cairo_pdf) +
  ggsave(paste0("output-figure-R1B-maxdists.png"), width=137.25, height=91.5, units="mm")











# Plot binding site distance distribution
q <- tibble(qo)

# Split into real & random, subsample random, then recombine
# Real
real <- subset(q, random==0)
# Random
if (exists("random")) { rm(random) }
for (myset in unique(q$set)) {
  cat(paste0(" >> ",myset,"\n"))
  
  # Subsample the random set
  set.seed(2020)
  random_sub <- sample_n(subset(q, random==1 & set==myset), nrow(subset(q, random==0 & set==myset)), replace=T)
  
  if (exists("random") == T) {
    random <- rbind(random, random_sub)
  } else {
    random <- random_sub
  }
}
# Recombine
summary(real)
summary(random)
q <- rbind(real, random)
# str(q)

# How many "distance=0" cases are there?
real %>% 
  filter(set == "screen_hits") %>% 
  filter(mindist<=1000) %>%
  group_by(pair) %>%
  summarise(mediandist = median(mindist)) %>%
  ggplot(aes(x=mediandist, y=pair)) + geom_bar(stat="identity")

real %>% 
  filter(set == "screen_hits") %>% 
  filter(mindist<=10000) %>%
  ggplot(aes(x=mindist)) + geom_histogram() + theme_minimal()

real %>% 
  filter(set == "screen_hits") %>% 
  filter(mindist<=54) %>%
  ggplot(aes(x=mindist)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=pretty_breaks(20), expand=c(0,0)) +
  theme_minimal() +
  xlab("Binding site distance") +
  ylab("Frequency") +
  # theme(panel.grid.minor = element_blank()) +
  # ggtitle("Screen hits") +
  ggsave(paste0("output-figure-R6-screen_hits.pdf"), width=6, height=4)

real %>% 
  filter(set == "positives") %>% 
  filter(mindist<=54) %>%
  ggplot(aes(x=mindist)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=pretty_breaks(20), expand=c(0,0)) +
  theme_minimal() +
  xlab("Binding site distance") +
  ylab("Frequency") +
  # theme(panel.grid.minor = element_blank()) +
  ggsave(paste0("output-figure-R6-positives.pdf"), width=6, height=4)

q %>% 
  filter(set %in% c("positives", "screen_hits", "background_resamples")) %>% 
  filter(mindist<=54) %>%
  ggplot(aes(x=mindist)) +
  geom_histogram(binwidth=1) +
  # geom_density() +
  scale_x_continuous(breaks=pretty_breaks(20), expand=c(0,0)) +
  facet_wrap(~fct_relevel(set, "positives", "screen_hits"), ncol=1, scales="free_y") +
  theme_minimal() +
  xlab("Binding site distance") +
  ylab("Frequency") +
  # theme(panel.grid.minor = element_blank()) +
  ggsave(paste0("output-figure-R6.pdf"), width=5, height=5)

q %>% 
  filter(set %in% c("positives", "screen_hits", "background_resamples")) %>% 
  filter(mindist<=54) %>%
  ggplot(aes(x=mindist, colour=fct_relevel(set, "positives", "screen_hits"))) +
  # geom_histogram(binwidth=1) +
  geom_density() +
  scale_x_continuous(breaks=pretty_breaks(20), expand=c(0,0)) +
  # facet_wrap(~fct_relevel(set, "positives", "screen_hits"), ncol=1, scales="free_y") +
  theme_minimal() +
  guides(fill=F, colour=F) +
  xlab("Binding site distance") +
  ylab("Frequency") +
  # theme(panel.grid.minor = element_blank()) +
  ggsave(paste0("output-figure-R6-density.pdf"), width=6, height=8)

# Which RBP pairs are they for?

# Strictly filtered for cobinding (resampling p-value and resampling Wilcoxon p-value significant in one orientation at least)
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (21: the "strict" Cytoscape network, Figure 5a)
qtmp <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b IS NOT NULL AND resampling_p_value_b_vs_a IS NOT NULL AND resampling_wilcoxon_p_value_a_vs_b IS NOT NULL AND resampling_wilcoxon_p_value_b_vs_a IS NOT NULL) AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# qtmp %>% select(symbol1, symbol2) %>% unique %>% transmute(pair = glue("{symbol1} {symbol2}")) -> qtmppairs
qtmp %>% select(symbol1, symbol2) %>% unique -> hcpairs
print(hcpairs, n=50)
# Add inverse pairs
hcpairs %<>% bind_rows(hcpairs %>% transmute(newsymbol1=symbol2, newsymbol2=symbol1) %>% transmute(symbol1=newsymbol1, symbol2=newsymbol2)) %>% unique
hcpairs
hcpairs %<>% transmute(pair=as.character(glue("{symbol1}|{symbol2}"))) %>% mutate(hc=1)
hcpairs

# Screen hit pairs
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course)
qtmp <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND sumis>=7.1 GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# qtmp %>% select(symbol1, symbol2) %>% unique %>% transmute(pair = glue("{symbol1} {symbol2}")) -> qtmppairs
qtmp %>% select(symbol1, symbol2) %>% unique -> screenpairs
print(screenpairs, n=50)
# # Add inverse pairs
# screenpairs %<>% bind_rows(screenpairs %>% transmute(newsymbol1=symbol2, newsymbol2=symbol1) %>% transmute(symbol1=newsymbol1, symbol2=newsymbol2)) %>% unique
# screenpairs
screenpairs %<>% transmute(pair=as.character(glue("{symbol1}|{symbol2}"))) %>% mutate(screenpair=1)
screenpairs

# Screen hit pairs (with inverse pairs)
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course)
qtmp <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND sumis>=7.1 GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# qtmp %>% select(symbol1, symbol2) %>% unique %>% transmute(pair = glue("{symbol1} {symbol2}")) -> qtmppairs
qtmp %>% select(symbol1, symbol2) %>% unique -> screenpairflip
print(screenpairflip, n=50)
# Add inverse pairs
screenpairflip %<>% bind_rows(screenpairflip %>% transmute(newsymbol1=symbol2, newsymbol2=symbol1) %>% transmute(symbol1=newsymbol1, symbol2=newsymbol2)) %>% unique
screenpairflip
screenpairflip %<>% transmute(pair=as.character(glue("{symbol1}|{symbol2}"))) %>% mutate(screenpairflip=1)
screenpairflip


real %>% 
  filter(set %in% c("positives", "screen_hits")) %>% 
  filter(mindist==0) %>%
  group_by(pair) %>%
  summarise(n0=n()) -> real0

real %>% 
  filter(set %in% c("positives", "screen_hits")) %>% 
  group_by(pair) %>%
  summarise(ntotal=n()) -> realtotal

real %>% 
  filter(set %in% c("positives", "screen_hits")) %>% 
  filter(mindist<=5) %>%
  group_by(pair) %>%
  summarise(n5=n()) -> real5

real %>% 
  filter(set %in% c("positives", "screen_hits")) %>% 
  filter(mindist<=54) %>%
  group_by(pair) %>%
  summarise(n54=n()) -> real54

left_join(real54, real0) %>%
  mutate(n0 = replace_na(n0, 0)) %>% 
  mutate(frac0 = n0/n54) %>%
  left_join(realtotal) %>%
  # mutate(ntotal = replace_na(ntotal, 0)) %>% 
  left_join(hcpairs) %>%
  mutate(hc = replace_na(hc, 0)) %>%
  filter(frac0>=0.1) %>%
  filter(n54>=10) %>%
  # arrange(-frac0)
  # filter(frac0>=0.25) %>%
  # filter(frac0>=0.5) %>%
  # filter(ntotal>=100) %>%
  # filter(ntotal>=100) %>%
  ggplot(aes(x=fct_reorder(pair, -frac0), y=frac0, fill=hc)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  xlab("RBP pair") +
  ylab("Fraction of distances = 0 nt") +
  ggsave(paste0("output-figure-R7.pdf"), width=6, height=4)

left_join(real54, real5) %>%
  mutate(n5 = replace_na(n5, 0)) %>% 
  mutate(frac5 = n5/n54) %>%
  left_join(real0) %>%
  mutate(n0 = replace_na(n0, 0)) %>% 
  mutate(frac0 = n0/n54) %>%
  left_join(realtotal) %>%
  # mutate(ntotal = replace_na(ntotal, 0)) %>% 
  left_join(hcpairs) %>%
  mutate(hc = replace_na(hc, 0)) %>%
  # filter(frac5>=0.2) %>%
  filter(frac5>=0.25) %>%
  filter(n54>=10) %>%
  arrange(-frac5) %>%
  # filter(frac5>=0.5) %>%
  # filter(ntotal>=100) %>%
  ggplot(aes(x=fct_reorder(pair, -frac5), y=frac5, fill=hc)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  xlab("RBP pair") +
  ylab("Fraction of distances <=5 nt") +
  ggsave(paste0("output-figure-R8.pdf"), width=10, height=4)


left_join(real54, real5) %>%
  mutate(n5 = replace_na(n5, 0)) %>% 
  mutate(frac5 = n5/n54) %>%
  left_join(hcpairs) %>%
  mutate(hc = replace_na(hc, 0)) %>%
  # filter(frac5>=0.2) %>%
  # filter(frac5>=0.25) %>%
  # filter(frac5>=0.5) %>%
  filter(frac5>=1) %>%
  left_join(real) %>%
  filter(mindist<1000) %>%
  # print(n=10000)
  group_by(pair) %>%
filter(mindist<=54) %>% print(n=100)
  # summarise(n=n())








# Random data: Do close distances mean we see more interactions, even if the proteins don't interact (RNA bridging)?
real %>% 
  filter(set %in% c("background_resamples")) %>% 
  group_by(set, pair) %>%
  summarise(ntotal=n()) -> randomtotal

real %>% 
  filter(set %in% c("background_resamples")) %>% 
  filter(mindist<=5) %>%
  group_by(set, pair) %>%
  summarise(n5=n()) -> random5

real %>% 
  filter(set %in% c("background_resamples")) %>% 
  filter(mindist<=54) %>%
  group_by(set, pair) %>%
  summarise(n54=n()) -> random54

real %>% 
  filter(set %in% c("background_resamples")) %>% 
  filter(mindist<=1000) %>%
  group_by(set, pair) %>%
  summarise(n1000=n()) -> random1000

real %>% 
  filter(set %in% c("screen_hits")) %>% 
  filter(mindist<=1000) %>%
  group_by(set, pair) %>%
  summarise(n1000=n()) -> real1000

real %>% 
  filter(set %in% c("screen_hits")) %>% 
  filter(mindist<=54) %>%
  group_by(set, pair) %>%
  summarise(n54=n()) -> real54

left_join(random1000, random54, by=c("set", "pair")) %>%
  bind_rows(left_join(real1000, real54, by=c("set", "pair"))) %>%
  mutate(n54 = replace_na(n54, 0)) %>% 
  mutate(frac54 = n54/n1000) %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  left_join(screenpairflip) %>%
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>%
  # filter(frac54>=0.2) %>%
  # filter(frac54>=0.25) %>%
  filter(n1000>=100) %>%
  # filter(set=="background_resamples" & screenpair==1) %>% pull(pair) %>% cat(sep="\n") 
  filter(!(set=="background_resamples" & screenpairflip==1)) %>%   # Filter out background_resample RBP pairs that are screen hits (including inverse pairs)
  # filter(!(set=="screen_hits")) %>%
  arrange(-frac54) %>%
  # filter(frac54>=0.5) %>%
  # filter(ntotal>=100) %>%
  ggplot(aes(x=fct_reorder(pair, -frac54), y=frac54, fill=screenpair)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  xlab("RBP pair") +
  ylab("Fraction of distances <=54 nt") +
  ggsave(paste0("output-figure-R7A.pdf"), width=12, height=4)

left_join(random1000, random54, by=c("set", "pair")) %>%
  bind_rows(left_join(real1000, real54, by=c("set", "pair"))) %>%
  mutate(n54 = replace_na(n54, 0)) %>% 
  mutate(frac54 = n54/n1000) %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  left_join(screenpairflip) %>%
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>%
  # filter(frac54>=0.2) %>%
  # filter(frac54>=0.25) %>%
  filter(n1000>=100) %>%
  # filter(set=="background_resamples" & screenpair==1) %>% pull(pair) %>% cat(sep="\n") 
  filter(!(set=="background_resamples" & screenpairflip==1)) %>%   # Filter out background_resample RBP pairs that are screen hits (including inverse pairs)
  # filter(!(set=="screen_hits")) %>%
  arrange(-frac54) %>%
  # filter(frac54>=0.5) %>%
  # filter(ntotal>=100) %>%
  mutate(set=recode(set, "screen_hits"="Screen hits", "background_resamples"="Random pairs")) %>%
  ggplot(aes(x=frac54, y=set, colour=set, fill=set)) + 
  geom_boxplot(notch=T, alpha=0.3) + 
  # geom_violin(draw_quantiles = c(0.5), alpha=0.3) + 
  geom_beeswarm(groupOnX = F, shape=16, alpha=0.7) +
  scale_x_continuous(limits = c(0,1)) +
  scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) +
  guides(colour=F, fill=F) +
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  xlab("Fraction of distances <= 54 nt") +
  ylab(NULL) +
  ggsave(paste0("output-figure-R5.pdf"), width=4, height=2)




# left_join(random1000, random54, by=c("set", "pair")) %>%
# bind_rows(left_join(real1000, real54, by=c("set", "pair"))) %>%
myrandom <- left_join(random1000, random54, by=c("set", "pair"))
myrandom %<>% 
  left_join(screenpairflip) %>% 
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>% 
  filter(screenpairflip==0) %>% select(-screenpairflip)
myreal <- left_join(real1000, real54, by=c("set", "pair"))
# myreal
# nrow(myreal)
# Get an even number of screen_hit and background RBP pairs
bind_rows(slice_sample(myrandom, n=nrow(myreal)), myreal) %>%
  mutate(n54 = replace_na(n54, 0)) %>% 
  mutate(frac54 = n54/n1000) %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  left_join(screenpairflip) %>%
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>%
  # filter(frac54>=0.2) %>%
  # filter(frac54>=0.25) %>%
  # filter(n1000>=100) %>%
  # filter(set=="background_resamples" & screenpair==1) %>% pull(pair) %>% cat(sep="\n")
  filter(!(set=="background_resamples" & screenpairflip==1)) %>%   # Filter out background_resample RBP pairs that are screen hits (including inverse pairs)
  # filter(!(set=="screen_hits")) %>%
  # summary
  # filter(n54>=100) %>%
  arrange(-frac54) %>%
  # filter(frac54>=0.5) %>%
  # filter(ntotal>=100) %>%
  mutate(set=recode(set, "screen_hits"="Screen hits", "background_resamples"="Random pairs")) %>%
  ggplot(aes(x=n54+1, y=set, colour=set, fill=set)) + 
  geom_boxplot(notch=T, alpha=0.3) + 
  # geom_violin(draw_quantiles = c(0.5), alpha=0.3) + 
  geom_beeswarm(shape=16, alpha=0.7) +
  scale_x_log10() +
  # scale_y_continuous(limits = c(0,1)) +
  scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) +
  guides(colour=F, fill=F) +
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  # xlab("Fraction of distances <= 54 nt") +
  # ylab(NULL) +
  ggsave(paste0("output-figure-R7C.pdf"), width=4, height=2)




# left_join(random1000, random54, by=c("set", "pair")) %>%
# bind_rows(left_join(real1000, real54, by=c("set", "pair"))) %>%
myrandom <- left_join(random1000, random54, by=c("set", "pair"))
myrandom %<>% 
  left_join(screenpairflip) %>% 
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>% 
  filter(screenpairflip==0) %>% select(-screenpairflip)
myreal <- left_join(real1000, real54, by=c("set", "pair"))
# myreal
# nrow(myreal)
# Get an even number of screen_hit and background RBP pairs
bind_rows(slice_sample(myrandom, n=nrow(myreal)), myreal) %>%
  mutate(n54 = replace_na(n54, 0)) %>% 
  mutate(frac54 = n54/n1000) %>%
  left_join(screenpairs) %>%
  mutate(screenpair = replace_na(screenpair, 0)) %>%
  left_join(screenpairflip) %>%
  mutate(screenpairflip = replace_na(screenpairflip, 0)) %>%
  # filter(frac54>=0.2) %>%
  # filter(frac54>=0.25) %>%
  # filter(n1000>=100) %>%
  # filter(set=="background_resamples" & screenpair==1) %>% pull(pair) %>% cat(sep="\n")
  filter(!(set=="background_resamples" & screenpairflip==1)) %>%   # Filter out background_resample RBP pairs that are screen hits (including inverse pairs)
  # filter(!(set=="screen_hits")) %>%
  # summary
  arrange(-frac54) %>%
  # filter(n54>=100) %>%
  mutate(n54_100 = ifelse(n54>=100, 1, 0)) %>%
  group_by(set) %>%
  summarise(total=n(), n=sum(n54_100)) %>%
  mutate(n54_frac = n/max(n)) %>%
  # filter(frac54>=0.5) %>%
  # filter(ntotal>=100) %>%
  mutate(set=recode(set, "screen_hits"="Screen hits", "background_resamples"="Random pairs")) %>%
  ggplot(aes(x=fct_rev(set), y=n54_frac, colour=set, fill=set)) + 
  geom_bar(stat="identity") + 
  # geom_violin(draw_quantiles = c(0.5), alpha=0.3) + 
  # geom_beeswarm(shape=16, alpha=0.7) +
  # scale_x_log10() +
  # scale_y_continuous(limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) +
  guides(colour=F, fill=F) +
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  guides(fill=F) +
  xlab("Pairs with distances <= 54 nt") +
  ylab(NULL) +
  ggsave(paste0("output-figure-R7D.pdf"), width=4, height=2)




left_join(random54, random5) %>%
  mutate(n5 = replace_na(n5, 0)) %>% 
  mutate(frac5 = n5/n54) %>%
  left_join(hcpairs) %>%
  mutate(hc = replace_na(hc, 0)) %>%
  # filter(frac5>=0.2) %>%
  # filter(frac5>=0.25) %>%
  # filter(frac5>=0.5) %>%
  filter(frac5>=1) %>%
  left_join(random) %>%
  filter(mindist<1000) %>%
  # print(n=10000)
  group_by(pair) %>%
  filter(mindist<=54) %>% print(n=100)
# summarise(n=n())
