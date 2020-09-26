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


# myset <- "negatives"
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





