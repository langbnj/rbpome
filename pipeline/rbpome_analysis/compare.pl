#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $table = 'rbpome';

our $usage = "$0 [table (e.g. rbpome)] [eclip type: eclip_encode/eclip_tom/eclip_tom_local] [interaction screen score type: AVG_IS/SUM_IS/11_IS] [analysis mode for the x-axis] [x-axis threshold for labelling] [analysis mode for the y-axis] [y-axis threshold for labelling] [hc: 0/1] [minlog2fold] [minimum interaction screen score] [mindist_threshold] [resamples]\n\nExample: $0 rbpome eclip_encode AVG_IS probability 0.875 jaccard 0.4 0 0 0 10 100";
($table, $type, $scoretype, $modex, $threshx, $modey, $threshy, $hc, $minlog2fold, $minscore, $mindist_threshold, $resamples) = args(12);
#args(0);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

# cat output/output-txt-rbpome-eclip_encode_12-SUM_IS-jaccard-hc0-minlog2fold0-minscore7.1-mindist_threshold50-resamples100.txt

$infilex = "output/output-txt-$table-$type-$scoretype-$modex-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$infiley = "output/output-txt-$table-$type-$scoretype-$modey-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$outfile = "output-compare-pdf-$table-$type-$scoretype-$modex-$threshx-$modey-$threshy-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.pdf";
$outfiletxt = "output-compare-txt-$table-$type-$scoretype-$modex-$threshx-$modey-$threshy-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

nl();
state("Reading '$infilex' ($modex, x-axis)", 1);
state("Reading '$infiley' ($modey, y-axis)", 1);
nl();

state("Comparing '$type' '$scoretype' '$modex' (x-axis) and '$modey' (y-axis)...");

# state("x: $infilex");
# state("y: $infiley");

startr();
setr('fig2e', switch('2E'));
state(runr(qq(

library(ggplot2)
library(ggrepel)
library(Cairo)
library(dplyr)

# source("~/include/mysql.R")

# setwd("~/Documents/Projects/RBPome/compare")

xtable <- read.table("$infilex", header=T, sep="\t")
ytable <- read.table("$infiley", header=T, sep="\t")

# Melt
# q <- melt(table, id.vars = c("table", "type", "mode", "hc", "minlog2fold", "minscore", "mindist_threshold", "resamples", "set", "rbp1", "rbp2"))
# q <- melt(table, measure.vars = c("cobind"))
# ?melt
# q\$variable <- q\$mode

# Merge
# unique(xtable\$mode)
xtable\$$modex <- xtable\$cobind
xtable\$mode <- NULL
xtable\$cobind <- NULL
# unique(ytable\$mode)
ytable\$$modey <- ytable\$cobind
ytable\$mode <- NULL
ytable\$cobind <- NULL

# head(xtable)
# head(ytable)

# Get mode directionality
mode_directionality <- function(mode) {
  if    (mode == 'fraction')							{ mode_directional <- 0 }
  else if (mode == 'freq') 								{ mode_directional <- 0 }
  else if (mode == 'intersection') 						{ mode_directional <- 0 }
  else if (mode == 'jaccard') 							{ mode_directional <- 0 }
  else if (mode == 'mediandist') 						{ mode_directional <- 0 }
  else if (mode == 'mindist') 							{ mode_directional <- 0 }
  else if (mode == 'mindist_thresh') 					{ mode_directional <- 0 }
  else if (mode == 'mindist_thresh_frac') 				{ mode_directional <- 0 }
  else if (mode == 'oddsratio') 						{ mode_directional <- 0 }
  else if (mode == 'probability') 						{ mode_directional <- 1 }
  else if (mode == 'probability_freqscaled') 			{ mode_directional <- 1 }
  else if (mode == 'probability_mindist') 				{ mode_directional <- 1 }
  else if (mode == 'probability_mindist_freqscaled') 	{ mode_directional <- 1 }
  else {
    warning(paste0("Error: Mode directionality not yet defined for mode '",mode,"'"))
    quit(1)
  }
  return(mode_directional)
}

# If modex isn't directional, but modey is: Add flipped pairs to modex
if ((mode_directionality('$modex') == 0) & (mode_directionality('$modey') == 1)) {
	# Add flipped pairs for modex
	flipxtable <- xtable
	flipxtable\$rbptmp <- flipxtable\$rbp1
	flipxtable\$rbp1 <- flipxtable\$rbp2
	flipxtable\$rbp2 <- flipxtable\$rbptmp
	flipxtable\$rbptmp <- NULL
	xtable <- rbind(xtable, flipxtable)
}
# If modey isn't directional, but modex is: Add flipped pairs to modey
if ((mode_directionality('$modex') == 1) & (mode_directionality('$modey') == 0)) {
	# Add flipped pairs for modey
	flipytable <- ytable
	flipytable\$rbptmp <- flipytable\$rbp1
	flipytable\$rbp1 <- flipytable\$rbp2
	flipytable\$rbp2 <- flipytable\$rbptmp
	flipytable\$rbptmp <- NULL
	ytable <- rbind(ytable, flipytable)
}

q <- merge(xtable, ytable)

# Add pairs column
q\$pairs <- paste0(q\$rbp1, 'ǀ', q\$rbp2)


# # Subset to remove excess background_resamples (only want #1 for plotting)
# Can't do this because every analysis (jaccard, probability...) has its own set of background_resample pairs, i.e. I won't be able to merge them in compare.pl except if they happen to be the same pairs
# q <- q[q\$set \%in% c("positives", "screen_hits", "background_resample_1"),]
# # Remake factor to remove extra levels
# q\$set <- as.factor(as.character(q\$set))

# ...but: subsample background_resamples to be the same size as the pairs
real <- subset(q, set \%in\% c("positives", "screen_hits"))
random <- subset(q, !(set \%in\% c("positives", "screen_hits")))

# Remake factor to remove extra levels
real\$set <- as.factor(as.character(real\$set))
random\$set <- as.factor(as.character(random\$set))

# summary(real)
# summary(random)

# random <- sample_n(random, nrow(real))
# # Reproducible:
# random <- head(random, n=nrow(real))
# Also reproducible:
set.seed(0)
random <- sample_n(random, nrow(real))

# Remove positives from screen_hits (so they're not plotted twice)
# Get sets
positives <- subset(real, set=="positives")
screen_hits <- subset(real, set=="screen_hits")

# Get pairs
positive_pairs <- unique(positives\$pairs)
screen_hits_pairs <- unique(screen_hits\$pairs)
# print(positive_pairs \%in\% screen_hits_pairs)

screen_hits <- subset(screen_hits, !(pairs \%in\% positive_pairs))


# Recombine
# q <- rbind(real, random)
q <- rbind(positives, screen_hits, random)


# Convert 'known' screen_hits to 'positives' (plot them in blue)?



# real

# str(q)
# head(q)

# if ((fig2e == 1) | (("$modex" == "probability_mindist") | ("$modex" == "probability_mindist_freq")) & (("$modey" == "jaccard") | ("$modey" == "freq"))) {
if (fig2e == 1) {

	# Fine-tuned plot style (Figure 2E)

	# Generate point labels
	# Label points above thresholds
	q\$label <- NA
	q\$labelled <- 0
	if (nrow(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]) > 0) {
		q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$label <- paste0(toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp1), "ǀ", toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp2))
		# q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$labelled <- 2
		
		# Label flipped pair, too
		flipped_labelled_pairs <- paste0(toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp2), "ǀ", toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp1))
		q[q\$pair \%in\% flipped_labelled_pairs,]\$label <- q[q\$pair \%in\% flipped_labelled_pairs,]\$pairs
	}
	
	# Label positives as well
	q[q\$set=="positives",]\$label <- paste0(toupper(q[q\$set=="positives",]\$rbp1), "ǀ", toupper(q[q\$set=="positives",]\$rbp2))
	q[q\$set=="positives",]\$labelled <- 1
	
	# Highlight screen_hits
	q[q\$pairs \%in\% screen_hits_pairs,]\$labelled <- 2
	# # Label all screen_hits (including the ones below the labelling thresholds above)
	# q[q\$pairs \%in\% screen_hits_pairs,]\$label <- paste0(toupper(q[q\$pairs \%in\% screen_hits_pairs,]\$rbp1), "ǀ", toupper(q[q\$pairs \%in\% screen_hits_pairs,]\$rbp2))

	# Highlight positives
	q[q\$pairs \%in\% positive_pairs,]\$labelled <- 1
	# # Label all positives (including the ones below the labelling thresholds above)
	# q[q\$pairs \%in\% positive_pairs,]\$label <- paste0(toupper(q[q\$pairs \%in\% positive_pairs,]\$rbp1), "ǀ", toupper(q[q\$pairs \%in\% positive_pairs,]\$rbp2))
	
	
	# Turn into factor
	q\$labelled <- factor(q\$labelled)
	# print(str(q\$labelled))
	# # levels(q\$labelled) <- rev(levels(q\$labelled))
	# print(str(q\$labelled))
	# Sort data frame so labelled points get plotted over the others
	q <- arrange(q, labelled)

	# # Subset to screen hits only
	# q <- subset(q, set=="screen_hits")
	# # Subset to screen hits and positives only
	# q <- subset(q, set=="screen_hits" | set=="positives")
	
	print(str(q))
	print(summary(q))

	# Set plot size
	expand <- 0
	xmin <- 0
	xmax <- 1
	ymin <- max(min(q\$$modey) - expand, 0)
	ymax <- min(max(q\$$modey) + expand, 1)
	# , hjust=0, vjust=0
	
	if 	(("$modex" == "probability_mindist") & ("$modey" == "jaccard")) {
		# 2E Option A
		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) + 
		  xlab("Conditional probability of co-binding within $mindist_threshold nt, p(AǀB)") + ylab("Target set similarity (Jaccard index)") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		  coord_cartesian(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + theme_minimal() + theme(legend.position = "none") +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


	} else if (("$modex" == "probability_mindist_freqscaled") & ("$modey" == "jaccard")) {
		# 2E Option B
		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
		 xlab("Conditional probability of co-binding within $mindist_threshold nt, p(AǀB) (scaled by fraction of transcriptome bound by either protein)") + ylab("Target set similarity (Jaccard index)") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		 coord_cartesian(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + theme_minimal() + theme(legend.position = "none") +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


	} else if (("$modex" == "probability_mindist") & ("$modey" == "freq")) {
		# 2E Option C
		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
		 xlab("Conditional probability of co-binding within $mindist_threshold nt, p(AǀB)") + ylab("Fraction of transcriptome bound by either protein") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		 coord_cartesian(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + theme_minimal() + theme(legend.position = "none") +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


	} else if (("$modex" == "oddsratio") & ("$modey" == "jaccard")) {
		# 2E Option D
		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
		 xlab("Co-binding odds ratio") + ylab("Target set similarity (Jaccard index)") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		 theme_minimal() + theme(legend.position = "none") + scale_x_log10() +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


	} else if (("$modex" == "probability_mindist") & ("$modey" == "oddsratio")) {
		# 2E Option F
		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
		 xlab("Conditional probability of co-binding within $mindist_threshold nt, p(AǀB)") + ylab("Co-binding odds ratio") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		 theme_minimal() + theme(legend.position = "none") + coord_cartesian(xlim=c(0, 1)) + scale_y_log10() +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


 	} else if (("$modex" == "oddsratio") & ("$modey" == "fraction")) {
 		# 2E Option G
 		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
 		 xlab("Co-binding odds ratio") + ylab("Overlap coefficient") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
 		 theme_minimal() + theme(legend.position = "none") + scale_x_log10() +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


 	} else if (("$modex" == "oddsratio") & ("$modey" == "intersection")) {
 		# 2E Option H
 		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
 		 xlab("Co-binding odds ratio") + ylab("Number of co-bound RNAs (intersection)") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
 		 theme_minimal() + theme(legend.position = "none") + scale_x_log10() +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


 	} else if (("$modex" == "oddsratio") & ("$modey" == "mindist_thresh")) {
 		# 2E Option I
 		p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
 		 xlab("Co-binding odds ratio") + ylab("Number of co-bound RNAs (intersection) within $mindist_threshold nt") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
 		 theme_minimal() + theme(legend.position = "none") + scale_x_log10() +
		  annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
		  "\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)

  	} else if (("$modex" == "oddsratio") & ("$modey" == "probability_mindist")) {

  	  	# 2E oddsratio vs probability_mindist (used to be "final")
  	  	p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
  	  	 xlab("Co-binding odds ratio") + ylab("Conditional prob. of co-binding within $mindist_threshold nt") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
  	  	 theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10()
  	  	 # theme_minimal() + theme(legend.position = "none") + scale_x_log10(limits=c(3, 20)) + scale_y_continuous(limits=c(0, 0.4)) # Zoomed
  	 	# ggsave("Figure 2E - Final - $type - $scoretype.pdf", device=cairo_pdf, width=183, height=91.5, units="mm")
  	 	# print("Wrote to Figure 2E - Final - $type - $scoretype.pdf")
  	 	ggsave("Figure 2E - oddsratio vs probability_mindist - $type - $scoretype.pdf", device=cairo_pdf, width=183, height=91.5, units="mm")
  	 	print("Wrote to Figure 2E - oddsratio vs probability_mindist - $type - $scoretype.pdf")

  	 	# 	  	# 2E Linear fit (not used)
  	 	# 	  	p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(method="lm", aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
  	 	# 	  	 xlab("Co-binding odds ratio") + ylab("Conditional prob. of co-binding within $mindist_threshold nt") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
  	 	# 	  	 theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10()
  	 	# ggsave("Figure 2E - Final Linear - $type - $scoretype.pdf", device=cairo_pdf, width=183, height=91.5, units="mm")
  	 	# print("Wrote to Figure 2E - Final Linear - $type - $scoretype.pdf")

  		# # 2E Option E
  		# p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + labs(tag = "E") + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
  		#  xlab("Co-binding odds ratio") + ylab("Conditional probability of co-binding within $mindist_threshold nt, p(AǀB)") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
  		#  theme_minimal() + theme(legend.position = "none") + scale_x_log10() +
  		#   annotate("text", x=Inf, y=Inf, label=paste0("Spearman's rho=",round(cor.test(q\$$modex, q\$$modey, method="spearman")\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey, method="spearman")\$p.value, format="g", digits=0),
  		#   "\\nPearson's r=",round(cor.test(q\$$modex, q\$$modey)\$estimate, 2)," p=",formatC(cor.test(q\$$modex, q\$$modey)\$p.value, format="g", digits=0)), vjust=1, hjust=1)


 	} else if (("$modex" == "jaccard") & ("$modey" == "probability")) {

	  	# 2E jaccard vs probability
	  	p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
	  	 xlab("RNA target set similarity (Jaccard index)") + ylab("Conditional probability of co-binding (p(AǀB))") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
	  	 theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1))
	  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10()
	  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10(limits=c(0.1, NA)) # Zoomed
	 	ggsave("Figure 2E - jaccard vs probability - $type - $scoretype.pdf", device=cairo_pdf, width=183, height=91.5, units="mm")
	 	print("Wrote to Figure 2E - jaccard vs probability - $type - $scoretype.pdf")

 	} else if (("$modex" == "jaccard") & ("$modey" == "probability_mindist")) {

 		# 2E jaccard vs probability_mindist (FINAL)
		# Make a bunch of versions for geom_text_repel
		for (dir in c("y", "both")) {
		  	# p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + stat_smooth(aes(colour=NULL), colour="grey") + geom_point() + geom_text_repel(aes(label=label)) +
		  	p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + stat_smooth(aes(colour=NULL), colour="grey95") + geom_point() + geom_text_repel(aes(label=label), direction=dir, force=10, max.iter=100000) +
		  	 xlab("RNA target set similarity (Jaccard index)") + ylab("Conditional prob. of co-binding within $mindist_threshold nt") + scale_colour_manual(values = c("grey", "#1E3D59", "#FF7F00")) +
		  	 # theme_minimal() + theme(legend.position = "none")
		  	 theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 0.9))
		  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 0.9)) + scale_x_log10() # log
		  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1))
		  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_sqrt() # sqrt
		  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10() # log
		  	 # theme_minimal() + theme(legend.position = "none") + coord_cartesian(ylim=c(0, 1)) + scale_x_log10(limits=c(0.1, NA)) # log Zoomed

			for (i in 1:20) {
			 	ggsave(paste0("Figure 2E - jaccard vs probability_mindist - $type - $scoretype - ",dir," - ",i,".pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
			 	print(paste0("Wrote to Figure 2E - jaccard vs probability_mindist - $type - $scoretype - ",dir," - ",i,".pdf"))
			}
		}

	} else {
		warning(paste0("WARNING: Not sure how to plot Figure 2E with modex '$modex' vs. modey '$modey'"))
	}
	# "Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm)."
	# ggsave("Figure 2E - $type - $scoretype.pdf", device=cairo_pdf, width=183, height=91.5, units="mm")
	# print("Wrote to Figure 2E - $type - $scoretype.pdf")
	
} else {

	# Normal plot style

	q <- subset(q, set=="screen_hits")

	# Generate point labels
	q\$label <- NA
	q\$labelled <- 0
	q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$label <- paste0(toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp1), "/", toupper(q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$rbp2))
	q[q\$$modex>=$threshx | q\$$modey>=$threshy,]\$labelled <- 1
	# q\$label <- paste0(toupper(q\$rbp1), "/", toupper(q\$rbp2))
	# # Apply thresholds
	# q[q\$$modex<$threshx | q\$$modey<$threshy,]\$label <- NA
	# q[q\$$modex<$threshx,]\$label <- NA
	# q[q\$$modey<$threshy,]\$label <- NA
	q\$labelled <- factor(q\$labelled)


	# Set plot size
	expand <- 0.1
	xmin <- max(min(q\$$modex) - expand, 0)
	xmax <- min(max(q\$$modex) + expand, 1)
	ymin <- max(min(q\$$modey) - expand, 0)
	ymax <- min(max(q\$$modey) + expand, 1)
	# , hjust=0, vjust=0

	p <- ggplot(q, aes(x=$modex, y=$modey, colour=labelled)) + geom_point() + stat_smooth(aes(colour=NULL)) + geom_text_repel(aes(label=label), colour="#1E3D59") + xlab("$modex") + ylab("$modey") + theme_bw() + scale_colour_manual(values = c("grey", "#1E3D59")) + coord_cartesian(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + labs(title="Cobinding \\"$modex\\" vs. \\"$modey\\"", subtitle="CLIP data: $type, minimum confidence $hc, minlog2fold $minlog2fold\\nmin_$scoretype $minscore, mindist_threshold $mindist_threshold, resamples $resamples")
	ggsave("$outfile", width=7, height=5)
	
}
write.table(q, "$outfiletxt", sep="\t", quote=F)
)));

# position=position_jitter(width=0.01, height=0.01)

stopr();

state("Wrote to '$outfile'");

done();
