#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table] [eclip type] [screen score type] [mindist threshold for cobinding] [motif search method] [motif 5' extension 0/1] [-qvalue]\n\n -qvalue: Use motif hits with significant q-values, rather than p-values (default).\n\nExample: $0 rbpome eclip_encode_12 SUM_IS 54 fimo 1\n";
($table, $type, $scoretype, $mindist_threshold, $motif_method, $motif_extend5) = args(6);

$tmpextend5 = '';
$tmpextend5 = '_extend5' if ($motif_extend5 == 1);

$tmpqvalue = '';
$tmpqvalue = '-qvalue' if (switch('qvalue'));

$tmpsig = 'psig=1';
$tmpsig = 'qsig=1' if (switch('qvalue'));


$infile = "output-motif_analysis-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";


$countfiles = "output/output-motif_counts-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-*$tmpqvalue.txt";
$countfile = "output-motif_counts-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";
$ratiofiles = "output/output-motif_ratios-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-*$tmpqvalue.txt";
$ratiofile = "output-motif_ratios-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";



# start


run("Combine motif count files into one", "head -n 1 `ls -1 $countfiles 2>/dev/null | head -n 1` > $countfile");
run("Combine motif count files into one", "cat $countfiles | grep -iPv '^symbol1\\tsymbol2\\tpair\\t' >> $countfile");

run("Combine motif ratio files into one", "head -n 1 `ls -1 $ratiofiles 2>/dev/null | head -n 1` > $ratiofile");
run("Combine motif ratio files into one", "cat $ratiofiles | grep -iPv '^symbol1\\tsymbol2\\tpair\\t' >> $ratiofile");





state("Making 'output-motif_analysis-motif_counts_for_examples-...' plots for some hand-picked example pairs...");
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



)));

# # Plotting count1 vs. count2:
# # linear hex
# ggplot(q, aes(x=count1, y=count2)) + geom_hex()
# # log hex
# ggplot(q, aes(x=log10(count1 + 1), y=log10(count2 + 1))) + geom_hex()
# # points (too slow)
# # ggplot(q, aes(x=count1, y=count2)) + geom_point()
# # density
# ggplot(q, aes(x=count1, y=count2)) + geom_density2d()
# cor.test(q\$count1, q\$count2)
# # Pearson says: anticorrelated
# cor.test(q\$count1, q\$count2, method="spearman")
# # Spearman says: correlated
# # >> There doesn't seem to be a strong correlation overall. Perhaps for some pairs?

# # Make plots for some pairs of interest
# mypair <- "FXR2|FMR1"
# mypair <- "HNRNPK|PCBP1"
# # for (mypair in unique(q\$pair)) {
# for (mypair in c("FXR2|FMR1", "HNRNPK|PCBP1", "IGF2BP2|PCBP2", "NONO|EWSR1", "NONO|SFPQ")) {
#
#   q %>% filter(pair==mypair, bindtype=="close") %>%
#    ggplot(aes(x=count1, y=count2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
#     xlab(paste0(rbp1," motif occurrences\\nwithin its binding region (cobound)")) + ylab(paste0(rbp2," motif occurrences\\nwithin its binding region (cobound)")) + ggtitle(paste0(rbp2," near ",rbp1))
#   # q %>% filter(pair=="FMR1|FXR2") %>% ggplot(aes(x=count1, y=count2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA))
#   ggsave(paste0("output-motif_analysis-motif_counts_for_examples-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-",str_replace(mypair, '\\\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*1, height=91.5*1, units="mm")
# }




state("Reading '$infile'");
state("Calculating & plotting motif count correlation statistics to 'output-motif_analysis-correlations-...'...");

state(runr(qq(
# Calculate correlation statistics
# # Get the columns needed...
# qp <- q %>% filter(pair==mypair)
# stats <- cor.test(qp\$count1, qp\$count2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1)
# # ...then delete the row
# stats %<>% filter(NA)
stats <- tibble(pair=character(), n=numeric(), pearson_p=double(), pearson_r=double(), spearman_p=double(), spearman_rho=double(), fisher_p=double(), fisher_odds=double())
# mypair <- "FXR2|FMR1"
# mypair <- "HNRNPK|PCBP1"
# mypair <- "SRSF9|TAF15"
for (mypair in unique(q\$pair)) {
  rbp1 <- str_split(mypair, '\\\\|', simplify=T)[1]
  rbp2 <- str_split(mypair, '\\\\|', simplify=T)[2]

  qp <- q %>% filter(pair==mypair, bindtype=="close")
  
  # print(qp)
  
  # if ((nrow(qp) <= 2))
  if ((nrow(qp) < 3) | sd(qp\$motifcount1) == 0 | sd(qp\$motifcount2) == 0)
  {
    next
  }
  
  # Pair order: ensure rbp1 is alphanumerically lower (to avoid double plotting of inverse pairs)
  if (rbp1 >= rbp2)
  {
    next
  }

  # Correlation tests
  # print(paste0(" >> ",mypair))
  # print(paste0("   >> pearson"))
  # stats %<>% add_row(tibble_row(method="pearson", cor.test(qp\$motifcount1, qp\$motifcount2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate)))
  # suppressWarnings(stats %<>% add_row(tibble_row(method="spearman", cor.test(qp\$motifcount1, qp\$motifcount2, method="spearman") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate))))


  # Do correlation tests
  suppressWarnings(rm(pear))
  suppressWarnings(rm(spear))
  suppressWarnings(rm(fish))
  pear <- cor.test(qp\$motifcount1, qp\$motifcount2, method="pearson") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate)
  suppressWarnings(spear <- cor.test(qp\$motifcount1, qp\$motifcount2, method="spearman") %>% tidy() %>% add_column("pair"=mypair, .before=1) %>% select(pair, p.value, estimate))


  # Contingency table for a Fisher test
  # qp
  # summary(qp)
  # summary(qp\$bindtype)
  print("Contingency table:")
  print(matrix(c(qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2==0) %>% nrow(),
                               qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2==0) %>% nrow()), ncol=2))
  fish <- fisher.test(matrix(c(qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1>0 & motifcount2==0) %>% nrow(),
                               qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2>0) %>% nrow(), qp %>% filter(bindtype=="close" & motifcount1==0 & motifcount2==0) %>% nrow()), ncol=2), alternative="greater") %>% tidy()



  # Add to stats tibble
  stats %<>% add_row(pair=mypair, n=nrow(qp), pearson_p=pear\$p.value, pearson_r=pear\$estimate, spearman_p=spear\$p.value, spearman_rho=spear\$estimate, fisher_p=fish\$p.value, fisher_odds=fish\$estimate)

  # stats
  # str(stats)
  # >> Strong correlations for FXR2|FMR1 (0.35/0.34), HNRNPK|PCBP1 (0.2/0.23), NONO|SFPQ (0.29/0.14)

  # tidy(cor.test(qp\$motifcount1, qp\$motifcount2))
  # augment(cor.test(qp\$motifcount1, qp\$motifcount2))
  # glance(cor.test(qp\$motifcount1, qp\$motifcount2))


  # # Make individual RBP pair plots
  # qp %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
  #   xlab(paste0(rbp1," motif occurrences\\nwithin its binding region (cobound)")) + ylab(paste0(rbp2," motif occurrences\\nwithin its binding region (cobound)")) + ggtitle(paste0(rbp2," near ",rbp1))
  # # q %>% filter(pair=="FMR1|FXR2") %>% ggplot(aes(x=motifcount1, y=motifcount2)) + geom_hex() + theme_bw() + coord_cartesian(xlim=c(0, NA), ylim=c(0, NA))
  # ggsave(paste0("output-motifcounts-",str_replace(mypair, '\\\\|', '-'),".pdf"), device=cairo_pdf, width=91.5*1, height=91.5*1, units="mm")
}




summary(stats\$pearson_r)
table(stats\$pearson_r)
stats %<>% arrange(-pearson_r)
stats
# >> Strong correlations for FXR2|FMR1 (0.35/0.34), HNRNPK|PCBP1 (0.2/0.23), NONO|SFPQ (0.29/0.14)

# Plot
# stats %>% add_column(label=as.character(NA))
stats_tmp <- stats
# as.data.frame(stats_tmp)
# stats_tmp\$label <- as.character(NA)
stats_tmp\$label <- ''
stats_tmp
# stats_tmp %>% filter(pearson_r >= 0.2 & spearman_rho >= 0.2) %>% mutate(label=pair)
stats_tmp[stats_tmp\$pearson_r>=0.2 & stats_tmp\$spearman_rho>=0.2 & stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05,]\$label <- stats_tmp[stats_tmp\$pearson_r>=0.2 & stats_tmp\$spearman_rho>=0.2 & stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05,]\$pair

stats_tmp\$col <- "Insignificant motif\\ncount correlation\\n(Pearson or Spearman)"
stats_tmp[stats_tmp\$pearson_p < 0.05 & stats_tmp\$spearman_p < 0.05,]\$col <- "Significant"
stats_tmp[stats_tmp\$fisher_p >= 0.05,]\$col <- "Insignificant motif\\npresence correlation\\n(Fisher's exact test)"
stats_tmp\$col <- as.factor(stats_tmp\$col)

print(stats_tmp)

stats_tmp %>%
ggplot(aes(x=pearson_r, y=spearman_rho, colour=fct_rev(col))) + geom_point() + geom_text_repel(aes(label=label)) + geom_vline(xintercept=0.2, linetype="dotted") + geom_hline(yintercept=0.2, linetype="dotted") +
  scale_colour_manual("Correlations", aesthetics = c("colour", "fill"), values = c("Insignificant motif\\ncount correlation\\n(Pearson or Spearman)" = "#808778", "Significant" = "#00a1d9", "Insignificant motif\\npresence correlation\\n(Fisher's exact test)" = "#bcc3b4")) +
  theme_bw() +
  labs(title='Correlation of motif counts for RBP pairs at cobound sites',
       caption="Cobound sites: cases where both RBPs have binding regions in proximity (≤54 nt).\\n\\nPearson and Spearman correlations of motif counts per RBP pair.\\nThe dotted lines indicate 0.2 as arbitrary thresholds of strong correlation.\\nEvery pair above these thresholds is labelled.\\n\\nThe Fisher test further indicates whether having a motif at all is correlated\\nbetween the two RBPs. The Fisher test takes precedence for the colour coding.") +
	   theme(plot.caption = element_text(hjust = 0)) +
  xlab("Pearson's r") +
  ylab("Spearman's rho") +
# ggsave(paste0("output-motif_analysis-correlations-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5.pdf"), device=cairo_pdf, width=91.5*2, height=91.5*2, units="mm")
ggsave(paste0("output-motif_analysis-correlations-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5.pdf"), width=91.5*2, height=91.5*2, units="mm")



# Make a more compact plot
stats_tmp <- stats
stats_tmp\$label <- ''
stats_tmp
# Label anything significant (according to Pearson, Spearman and Fisher)
stats_tmp[stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05 & stats_tmp\$fisher_p<=0.05,]\$label <- stats_tmp[stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05 & stats_tmp\$fisher_p<=0.05,]\$pair
print("All pairs:")
cat(stats_tmp\$pair, sep="\\n")
print("All tests significant:")
cat(stats_tmp[stats_tmp\$pearson_p<=0.05 & stats_tmp\$spearman_p<=0.05 & stats_tmp\$fisher_p<=0.05,]\$pair, sep="\\n")

stats_tmp\$col <- "Insignificant motif\\ncount correlation\\n(Pearson or Spearman)"
stats_tmp[stats_tmp\$pearson_p < 0.05 & stats_tmp\$spearman_p < 0.05,]\$col <- "Significant"
stats_tmp[stats_tmp\$fisher_p >= 0.05,]\$col <- "Insignificant motif\\npresence correlation\\n(Fisher's exact test)"
stats_tmp\$col <- as.factor(stats_tmp\$col)

stats_tmp %>%
ggplot(aes(x=pearson_r, y=spearman_rho, colour=fct_rev(col))) + geom_point() + geom_text_repel(aes(label=label), force=15, force_pull=0.1, box.padding = 0.5, max.overlaps = 0, max.time = 10, max.iter = 1000000) +
 # geom_vline(xintercept=0.2, linetype="dotted") + geom_hline(yintercept=0.2, linetype="dotted") +
  scale_colour_manual("Correlations", aesthetics = c("colour", "fill"), values = c("Insignificant motif\\ncount correlation\\n(Pearson or Spearman)" = "#808778", "Significant" = "#00a1d9", "Insignificant motif\\npresence correlation\\n(Fisher's exact test)" = "#bcc3b4")) +
  scale_x_continuous(breaks = breaks_pretty(6)) +
  scale_y_continuous(breaks = breaks_pretty(6)) +
  theme_minimal() +
  # labs(title='Correlation of motif counts for RBP pairs at cobound sites',
  #      caption="Cobound sites: cases where both RBPs have binding regions in proximity (≤54 nt).\\n\\nPearson and Spearman correlations of motif counts per RBP pair.\\nThe dotted lines indicate 0.2 as arbitrary thresholds of strong correlation.\\nEvery pair above these thresholds is labelled.\\n\\nThe Fisher test further indicates whether having a motif at all is correlated\\nbetween the two RBPs. The Fisher test takes precedence for the colour coding.") +
  # 	   theme(plot.caption = element_text(hjust = 0)) +
  guides(colour=F, fill=F) +
  xlab("Pearson's r") +
  ylab("Spearman's rho") +
  # ggsave(paste0("output-motif_analysis-correlations-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-compact.pdf"), width=91.5, height=91.5, units="mm")
  ggsave(paste0("output-motif_analysis-correlations-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5-compact.pdf"), width=110, height=110, units="mm")

q %>% filter(bindtype=="close") %>% print

)));






# state("Reading '$countfile' and '$ratiofile' and plotting motif pair ratios in 'close' vs. '$vs' binding contexts...");
# state(runr(qq(
#
#
#
# )));





stoptime();
stopr();

done();
