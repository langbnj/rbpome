#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table] [eclip type] [screen score type] [mindist threshold for cobinding] [motif search method] [motif 5' extension 0/1] [number of resamples] [RBP symbol 1] [RBP symbol 2] [control to compare to] [-qvalue]\n\n -qvalue: Use motif hits with significant q-values, rather than p-values (default).\n\nExample: $0 rbpome eclip_encode_12 SUM_IS 54 fimo 1 100 FMR1 FXR2\n";
($table, $type, $scoretype, $mindist_threshold, $motif_method, $motif_extend5, $resamples, $symbol1, $symbol2) = args(9);

$tmpextend5 = '';
$tmpextend5 = '_extend5' if ($motif_extend5 == 1);

$tmpqvalue = '';
$tmpqvalue = '-qvalue' if (switch('qvalue'));

$tmpsig = 'psig=1';
$tmpsig = 'qsig=1' if (switch('qvalue'));


$infile = "../output-motif_analysis-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5.txt";
$countfile = "../output/output-motif_counts-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-$symbol1-$symbol2$tmpqvalue.txt";
$ratiofile = "../output/output-motif_ratios-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-$symbol1-$symbol2$tmpqvalue.txt";



state("Reading from '$infile'...");
state("Analysing motif pair occurrence ratios in 'close' vs. 'closest_is_distant' and 'solo' binding regions in R...");
starttime();

startr();

state(runr(qq(

rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

# source("~/include/mysql.R")
library(broom)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(naturalsort)
library(progress)
library(readr)
library(reshape2)
library(scales)
library(stringr)
library(tibble)
library(tidyr)

# setwd("~/Documents/Projects/RBPome/mindist_motifs")

# rbp1 <- "FMR1"; rbp2 <- "FXR2" # FXR2 near FMR1
# rbp1 <- "FXR2"; rbp2 <- "FMR1" # FMR1 near FXR2

q <- read.table(paste0("$infile"), header=T, sep="\t", quote="")
q <- as_tibble(q)
q\$pair <- paste0(q\$symbol1,'|',q\$symbol2)


# Limit everything to the pair of interest
q %<>% filter(symbol1 == "$symbol1" & symbol2 == "$symbol2")


# # Plotting motifcount1 vs. motifcount2:
# # linear hex
# ggplot(q, aes(x=motifcount1, y=motifcount2)) + geom_hex()
# # log hex
# ggplot(q, aes(x=log10(motifcount1 + 1), y=log10(motifcount2 + 1))) + geom_hex()
# # points (too slow)
# # ggplot(q, aes(x=motifcount1, y=motifcount2)) + geom_point()
# # density
# ggplot(q, aes(x=motifcount1, y=motifcount2)) + geom_density2d()
# cor.test(q\$motifcount1, q\$motifcount2)
# # Pearson says: anticorrelated
# cor.test(q\$motifcount1, q\$motifcount2, method="spearman")
# # Spearman says: correlated
# # >> There doesn't seem to be a strong correlation overall. Perhaps for some pairs?

# # Make plots
# mypair <- "FXR2|FMR1"
# mypair <- "HNRNPK|PCBP1"
# # for (mypair in unique(q\$pair)) {
# for (mypair in c("FXR2|FMR1", "HNRNPK|PCBP1", "IGF2BP2|PCBP2", "NONO|EWSR1", "NONO|SFPQ")) {
#   qp %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
#     xlab(paste0(rbp1," motif occurrences\\nwithin its binding region (cobound)")) + ylab(paste0(rbp2," motif occurrences\\nwithin its binding region (cobound)")) + ggtitle(paste0(rbp2," near ",rbp1))
#   # q %>% filter(pair=="FMR1|FXR2") %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA))
#   ggsave(paste0("output-motifcounts-",str_replace(mypair, '\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*1, height=91.5*1, units="mm")
# }






# # Calculate correlation statistics
# # # Get the columns needed...
# # qp <- q %>% filter(pair==mypair)
# # stats <- cor.test(qp\$motifcount1, qp\$motifcount2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1)
# # # ...then delete the row
# # stats %<>% filter(NA)
# stats <- tibble(pair=character(), n=numeric(), pearson_p=double(), pearson_r=double(), spearman_p=double(), spearman_rho=double(), fisher_p=double(), fisher_odds=double())
# # mypair <- "FXR2|FMR1"
# # mypair <- "HNRNPK|PCBP1"
# # mypair <- "SRSF9|TAF15"
# for (mypair in unique(q\$pair)) {
#   rbp1 <- str_split(mypair, '\\|', simplify=T)[1]
#   rbp2 <- str_split(mypair, '\\|', simplify=T)[2]
#
#   qp <- q %>% filter(pair==mypair, bindtype=="close")
#   # if ((nrow(qp) <= 2))
#   if ((nrow(qp) < 3) | sd(qp\$motifcount1) == 0 | sd(qp\$motifcount2) == 0)
#   {
#     next
#   }
#
#   # Correlation tests
#   # print(paste0(" >> ",mypair))
#   # print(paste0("   >> pearson"))
#   # stats %<>% add_row(tibble_row(method="pearson", cor.test(qp\$motifcount1, qp\$motifcount2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate)))
#   # suppressWarnings(stats %<>% add_row(tibble_row(method="spearman", cor.test(qp\$motifcount1, qp\$motifcount2, method="spearman") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate))))
#
#
#   # Do correlation tests
#   suppressWarnings(rm(pear))
#   suppressWarnings(rm(spear))
#   suppressWarnings(rm(fish))
#   pear <- cor.test(qp\$motifcount1, qp\$motifcount2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate)
#   suppressWarnings(spear <- cor.test(qp\$motifcount1, qp\$motifcount2, method="spearman") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate))
#
#
#   # Contingency table for a Fisher test
#   # qp
#   # summary(qp)
#   # summary(qp\$bindtype)
#   fish <- fisher.test(matrix(c(qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2==0) %>% nrow(),
#                                qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2==0) %>% nrow()), ncol=2)) %>% tidy()
#
#
#
#   # Add to stats tibble
#   stats %<>% add_row(pair=mypair, n=nrow(qp), pearson_p=pear\$p.value, pearson_r=pear\$estimate, spearman_p=spear\$p.value, spearman_rho=spear\$estimate, fisher_p=fish\$p.value, fisher_odds=fish\$estimate)
#
#   # stats
#   # str(stats)
#   # >> Strong correlations for FXR2|FMR1 (0.35/0.34), HNRNPK|PCBP1 (0.2/0.23), NONO|SFPQ (0.29/0.14)
#
#   # tidy(cor.test(qp\$motifcount1, qp\$motifcount2))
#   # augment(cor.test(qp\$motifcount1, qp\$motifcount2))
#   # glance(cor.test(qp\$motifcount1, qp\$motifcount2))
#
#
#   # # Make individual RBP pair plots
#   # qp %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
#   #   xlab(paste0(rbp1," motif occurrences\\nwithin its binding region (cobound)")) + ylab(paste0(rbp2," motif occurrences\\nwithin its binding region (cobound)")) + ggtitle(paste0(rbp2," near ",rbp1))
#   # # q %>% filter(pair=="FMR1|FXR2") %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA))
#   # ggsave(paste0("output-motifcounts-",str_replace(mypair, '\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*1, height=91.5*1, units="mm")
# }
#
#
#
#
# summary(stats\$pearson_r)
# table(stats\$pearson_r)
# stats %<>% arrange(-pearson_r)
# stats
# # >> Strong correlations for FXR2|FMR1 (0.35/0.34), HNRNPK|PCBP1 (0.2/0.23), NONO|SFPQ (0.29/0.14)
#
# # Plot
# # stats %>% add_column(label=as.character(NA))
# stats_tmp <- stats
# # as.data.frame(stats_tmp)
# # stats_tmp\$label <- as.character(NA)
# stats_tmp\$label <- ''
# stats_tmp
# # stats_tmp %>% filter(pearson_r >= 0.2 & spearman_rho >= 0.2) %>% mutate(label=pair)
# stats_tmp[stats_tmp\$pearson_r>=0.2 & stats_tmp\$spearman_rho>=0.2 & stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05,]\$label <- stats_tmp[stats_tmp\$pearson_r>=0.2 & stats_tmp\$spearman_rho>=0.2 & stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05,]\$pair
# stats_tmp
# stats_tmp %>% ggplot(aes(x=pearson_r, y=spearman_rho, colour=fisher_p<0.05)) + geom_point() + geom_text_repel(aes(label=label)) + geom_vline(xintercept=0.2, linetype="dotted") + geom_hline(yintercept=0.2, linetype="dotted") + theme_bw() +
#   labs(title='Correlation of motif counts for RBP pairs at "cobound sites"',
#        caption="Cobound sites: cases where both RBPs have binding regions in proximity (≤54 nt).\\n\\nPearson and Spearman correlations of motif counts per RBP pair.\\nThe dotted lines indicate 0.2 as thresholds of strong correlation.\\nEvery pair above these thresholds is labelled.\\n\\nThe Fisher test further indicates whether having a motif at all\\n is correlated between the two RBPs.") +
#   theme(plot.caption = element_text(hjust = 0))
# ggsave(paste0("../output-motif_analysis-correlations-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*2, units="mm")
#




# Motif types
qm <- q
qm\$motif1 <- str_split(qm\$motifs1, ',')
qm\$motif2 <- str_split(qm\$motifs2, ',')
qm
qm %>% unnest(cols=motif1)
head(qm, n=4)
# head(qm, n=5) %>% unnest(cols=motif1)
head(qm, n=4) %>% unnest(cols=motif1)
# head(qm, n=4) %>% unnest(cols=c(motif1, motif2))
head(qm, n=4) %>% unnest(cols=motif1) %>% unnest(cols=motif2)
head(qm, n=1)
head(qm, n=1) %>% unnest(cols=motif1) %>% unnest(cols=motif2)
# Unnest the motif instances (puts each motif instance in its own row) and drop the old text columns listing all motifs
qm <- qm %>% unnest(cols=motif1) %>% unnest(cols=motif2) %>% select(-motifs1, -motifs2)
qm

# Replace '' with NA
# qm %>% mutate(motif1=replace(motif1, "", NA))
# qm %>% mutate(motif2=replace(motif2, "", NA))
# Convert to factors
qm\$motifpair <- as.factor(paste0(qm\$motif1,' & ',qm\$motif2))
qm\$motif1 <- as.factor(qm\$motif1)
qm\$motif2 <- as.factor(qm\$motif2)

# # Plot top n% of motif combinations (that aren't '')
# qm %>%
#   filter(motif1!='' & motif2!='') %>%
#   filter(bindtype=='close') %>%
#   group_by(pair, bindtype, motifpair) %>%
#   summarise(n=n()) %>%
#   arrange(-n) %>%
#   head(n=nrow(.) * 0.005) %>%
#   ggplot(aes(x=fct_reorder(motifpair, n), y=n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
# ggsave(paste0("output-motif-combinations-top0.5percent.pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*2, units="mm")






# # Get ratios & plot them (extremely slow, though - I think the ratio calculation step always scans the entire tibble)
# set.seed(1)
# # qmsam <- sample_n(subset(qm %>% arrange(), pair \%in% c("FMR1|FXR2", "HNRNPK|PCBP1")), 50, replace=F) %>% arrange(symbol1, symbol2, bindtype, motif1, motif2)
# # qmsam %>% print(n=100)
# qmr <- qm %>%
#   filter(motif1!='' & motif2!='') %>%
#   group_by(pair, motifpair, bindtype, .drop=F) %>% 
#   summarise(n=n()) %>%
#   group_by(pair, motifpair) %>%
#   summarise(ratio = n[bindtype=="close"]/n[bindtype=="closest_is_distant"]) %>%
#   filter(!is.na(ratio) & !is.infinite(ratio)) %>%
#   arrange(-ratio)
# # qmr %<>% filter(!is.na(ratio) & !is.infinite(ratio))
# # group_by(bindtype, add=TRUE) %>%
#   # summarise(nmin=min(n), nmax=max(n))
#   # summarise(ratio=n/sum(n))
#   # summarise(nmin=min(n), nmax=max(n)) %>%
#   # mutate(ratio=nmin/nmax)
#   # summarise(nclose=n[bindtype=="close"], ndist=n[bindtype=="closest_is_distant"])
# # %>%
# # arrange(-n) %>%
#   # filter(n>1000) %>%
#   # summarise(n_by_bindtype=sum(n))
#   # mutate(n_by_bindtype=sum(n)) %>% print(n=100)
# summary(qmr\$ratio)
# qmr
# 
# # Filter out ''
# qmr %<>% filter(!grepl("(^ & |^ & \$| & \$)", motifpair))
# 
# qmr %<>% arrange(-ratio)
# qmr
# qmr %>%
#   # head(n=nrow(qmr) * 0.01) %>%     # top 1%
#   head(n=50) %>%     # top 50
# ggplot(aes(x=fct_reorder(motifpair, ratio), y=ratio)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
# ggsave(paste0("output-motif-ratios-top50.pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*2, units="mm")
# 
# qmr %>%
#   # head(n=nrow(qmr) * 0.001) %>%     # top 0.1%
#   head(n=50) %>%     # top 50
#   print(n=100)
#   
# 
# Limit the bootstrapping below to the "top motif pairs":
# qmr %>%
#   arrange(-ratio) %>%
#   filter(ratio>=20)
# resr %>% filter(resample==1)
#
# top_motifpairs <- qmr %>%
#   arrange(-ratio) %>%
#
#   # Get top motif pairs only:
#   # head(n=nrow(qmr) * 0.01) %>%     # top 1% (129 motif pairs)
#   filter(ratio>=20) %>%               # ratio ≥ 20 (11 motif pairs)
#   # filter(ratio>=8) %>%               # ratio ≥ 8
#   # filter(ratio>=4) %>%               # ratio ≥ 4
#
#   ungroup() %>%
#   select(motifpair)
# # top_motifpairs <- sort(unique(as.character(top_motifpairs\$motifpair)))
# top_motifpairs <- unique(as.character(top_motifpairs\$motifpair))
# tibble(top_motifpairs)
# # Number of motifs:
# qm %>% filter(motifpair \%in% top_motifpairs)
# # RBP pairs with "top motif pairs":
# qm %>% filter(motifpair \%in% top_motifpairs) %>% select(symbol1, symbol2) %>% distinct() %>% print(n=100)
# # All RBP pairs:
# qm %>% select(symbol1, symbol2) %>% distinct()
# # >> 74 (as expected)







# Resample to get error bars
resamples <- $resamples
all_motifpairs <- qm %>% filter(motif1!='' & motif2!='') %>% pull(motifpair) %>% unique() %>% as.character()

# # Run for top_motifpairs only (to save time)
# my_motifpairs <- top_motifpairs

# Run for all motif pairs (not just top)
my_motifpairs <- all_motifpairs


# my_motifpair <- "attract|NONO_1 & attract|SFPQ_9"

# res <- tibble(symbol1=character(), symbol2=character(), pair=character(), motifpair=character(), n_close=numeric(), n_closest_is_distant=numeric())
res <- tibble(symbol1=character(), symbol2=character(), pair=character(), motifpair=character(), resample=numeric(), bindtype=character(), n=numeric())
resr <- tibble(symbol1=character(), symbol2=character(), pair=character(), motifpair=character(), resample=numeric(), raw_ratio=double(), ratio=double(), raw_ratio_solo=double(), ratio_solo=double())
p <- progress_bar\$new(format=paste0(" >> resampling [:bar] :current/:total (:percent) :eta"), total=length(my_motifpairs) * resamples)
for (my_motifpair in my_motifpairs) {
  
  # Get pair name for convenience
  qmrstmp <- qm %>% filter(motifpair==my_motifpair) %>% head(n=1)
  qmrstmp
  my_symbol1 <- as.character(qmrstmp\$symbol1)
  my_symbol2 <- as.character(qmrstmp\$symbol2)
  my_pair <- as.character(qmrstmp\$pair)

  # p <- progress_bar\$new(format=paste0(" >> ",my_symbol1," >> ",my_symbol2," >> ",my_motifpair," >> resampling [:bar] :current/:total (:percent)"), total=resamples)

  for (i in 1:resamples) {
    # qmrs <- qm %>%
    #   # filter(motif1!='' & motif2!='') %>%
    #   # filter(motifpair \%in% top_motifpairs) %>%
    #   filter(motifpair==my_motifpair) %>%
    #   sample_frac(size=1, replace=T) %>%
    #   # filter(symbol1=="FMR1" & symbol2=="FXR2") %>%
    #   group_by(motifpair, bindtype, .drop=F) %>% 
    #   summarise(n=n()) %>%
    #   group_by(motifpair) %>%
    #   summarise(ratio = n[bindtype=="close"]/n[bindtype=="closest_is_distant"]) %>%
    #   filter(!is.na(ratio) & !is.infinite(ratio)) %>%
    #   arrange(-ratio)
    # # qmrs
    
    qmrs <-
      qm %>%
      filter(motifpair==my_motifpair) %>%
      sample_frac(size=1, replace=T) %>%
      group_by(bindtype, .drop=F) %>% 
      summarise(n=n())
    qmrs
    n_close <- qmrs %>% filter(bindtype=="close") %>% pull(n)
    n_closest_is_distant <- qmrs %>% filter(bindtype=="closest_is_distant") %>% pull(n)
    n_solo <- qmrs %>% filter(bindtype=="solo") %>% pull(n)
    n_close
    n_closest_is_distant
    n_solo

    my_raw_ratio <- n_close / n_closest_is_distant
    my_raw_ratio
    my_raw_ratio_solo <- n_close / n_solo
    my_raw_ratio_solo

    # Use pseudocounts to get a finite, non-zero ratio
    my_ratio <- (n_close+1) / (n_closest_is_distant+1)
    my_ratio
    my_ratio_solo <- (n_close+1) / (n_solo+1)
    my_ratio_solo
    
    
    
    # %>%
    #   group_by(motifpair) %>%
    #   summarise(ratio = n[bindtype=="close"]/n[bindtype=="closest_is_distant"]) %>%
    #   filter(!is.na(ratio) & !is.infinite(ratio)) %>%
    #   arrange(-ratio)
    
    
    # # str(qmrs\$ratio)
    # if (nrow(qmrs) > 1)
    # {
    #   warning(paste0("Error: More than one ratio value for motifpair ",my_motifpair))
    # }

    # # Get ratio for this resample
    # my_ratio <- NA
    # if (nrow(qmrs) == 1) {
    #   my_ratio <- qmrs\$ratio
    # }
    

    # Add to results tibble
    # res %<>% add_row(tibble_row(symbol1=my_symbol1, symbol2=my_symbol2, pair=my_pair, motifpair=my_motifpair, ratio=my_ratio))
    # res %<>% add_row(tibble_row(symbol1=my_symbol1, symbol2=my_symbol2, pair=my_pair, motifpair=my_motifpair, n_close=n_close, n_closest_is_distant=n_closest_is_distant))
    res %<>% add_row(symbol1=my_symbol1, symbol2=my_symbol2, pair=my_pair, motifpair=my_motifpair, resample=i, bindtype=qmrs\$bindtype, n=qmrs\$n)
    resr %<>% add_row(symbol1=my_symbol1, symbol2=my_symbol2, pair=my_pair, motifpair=my_motifpair, resample=i, raw_ratio=my_raw_ratio, ratio=my_ratio, raw_ratio_solo=my_raw_ratio_solo, ratio_solo=my_ratio_solo)
    
    p\$tick(1)
  }
  # p\$tick(1)
}

res %>% write_tsv("$countfile")
resr %>% write_tsv("$ratiofile")


res
res %>% print(n=100)

res %>%
  ggplot(aes(fct_reorder(motifpair, n, desc=T), y=n, colour=fct_rev(bindtype))) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("../output/output-motif_analysis-boxplots-n-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-$symbol1-$symbol2$tmpqvalue.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")



# Ratio (with pseudocounts)
resr %>%
  ggplot(aes(x=fct_reorder(motifpair, ratio, desc=T), y=ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("../output/output-motif_analysis-boxplots-ratios-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-$symbol1-$symbol2$tmpqvalue.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")



# Raw ratio (no pseudocounts)
resr %>%
  ggplot(aes(x=fct_reorder(motifpair, raw_ratio, desc=T), y=raw_ratio)) + coord_flip() + geom_boxplot(notch=T) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("../output/output-motif_analysis-boxplots-raw_ratios-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-$symbol1-$symbol2$tmpqvalue.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*1.5, units="mm")

)));

stoptime();
stopr();

done();


state("Wrote motif counts to '$countfile'");
state("Wrote pseudocount ratios (+1/+1) and raw ratios to '$ratiofile'");

done();
