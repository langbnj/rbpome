#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# eCLIP parameters (for getting Jaccard etc. from rbpome_analysis output files)
# $table = 'rbpome';
# $scoretype = 'SUM_IS';
# $threshold = 7.1;
# $type = 'eclip_encode';
# $hc = 0;
# $minlog2fold = 0;
# $mindist_threshold = 100;
# $mindist_threshold = 50;
$mindist_threshold_strict = 54;	# "strict": This is the exact peak of the <=1000 positives curve in Figure 4e
$mindist_threshold_loose = 117;	# "looser": This is 95% of the peak of the <=1000 positives curve in Figure 4e
# $resamples = 100;

@mindist_thresholds = (10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000);

our $usage = "$0 [table (e.g. rbpome)] [type (e.g. eclip_encode_12)] [score type (e.g. SUM_IS)] [-randompairs]\n\n -randompairs: Read reproducible background resample random pairs from e.g. output/output-rbpome-randompairs.txt and use these instead of real RBP pairs\n\nExample: $0 rbpome eclip_encode_12 SUM_IS";
($table, $type, $scoretype) = args(3);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

$tmprandompairs = '';
if (switch('randompairs'))
{
	$tmprandompairs = '-randompairs';
}


$outtable = "output-fit_vs_random-overview-$table-$type-$scoretype-pvalues$tmprandompairs.txt";
$outpdf = "output-fit_vs_random-overview-$table-$type-$scoretype$tmprandompairs.pdf";
$outpng = "output-fit_vs_random-overview-$table-$type-$scoretype$tmprandompairs.png";

open(OUT, ">$outtable") or die("Error: Couldn't open '$outtable'");
print OUT "type\tscoretype\tmindist_threshold\twilcox_pvalue\n";


startr();
state(runr(qq(

library(ggplot2)
library(tibble)
library(forcats)
library(scales)
library(dplyr)

)));

foreach $close (@mindist_thresholds)
{
	$infile = "output-fit_vs_random-$table-$type-$scoretype-$close$tmprandompairs.txt";

	state(runr(qq(

	table$close <- read.table("$infile", header=T, sep="\t", quote="")
	
	# str(table)
	# head(table)

	q$close <- rbind(tibble(mindist=$close, class="Observed", value=table$close\$fraction_close_real), tibble(mindist=$close, class="Randomised", value=table$close\$fraction_close_random))
	# str(q)
	# head(q)
	
	
	p$close <- wilcox.test(table$close\$fraction_close_real, table$close\$fraction_close_random)\$p.value

	)));
}


# Make overview box plot (close binding events in the real data vs. in the random distribution)

$qs = '';
foreach $close (@mindist_thresholds)
{
	$qs .= "q$close, ";
}
$qs =~ s/, $//;

state(runr(qq(

# Combine into one tibble

q <- bind_rows($qs)


# Make overview box plot (close binding events in the real data vs. in the random distribution)
str(q)
head(q)

# ggplot(q, aes(x=class, y=value)) + geom_boxplot(notch=T, outlier.shape = NA) + coord_cartesian(ylim=c(NA, 0.45)) + theme_minimal() + xlab("") + ylab("Fraction of binding events in proximity")
# ggplot(q, aes(x=value, y=fct_rev(as.factor(class)))) + geom_boxplot(notch=T, outlier.shape = NA) + coord_cartesian(xlim=c(NA, 0.45)) + theme_minimal() + ylab("") + xlab("Fraction of binding events in proximity")
# p <- ggplot(q, aes(x=value, y=fct_rev(as.factor(class)), colour=class)) + coord_cartesian(xlim=c(0, 1)) + scale_x_continuous(breaks=pretty_breaks(5)) + geom_boxplot(notch=T) + scale_colour_manual(values=c('black', 'grey'), aesthetics=c("colour", "fill"), guide=F) + theme_minimal() + ylab("") + xlab("Fraction of binding events in proximity")
# p <- ggplot(q, aes(x=as.factor(mindist), y=value, colour=class)) + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=pretty_breaks(5)) + geom_boxplot(notch=T) + scale_colour_manual(values=c('black', 'grey'), aesthetics=c("colour", "fill"), guide=F) + theme_minimal() + xlab("Maximum distance considered [nt]") + ylab("Fraction of binding events in proximity")

# Hack ggpubr::stat_compare_means a little by supplying a fake kruskal.test functiont hat returns the resampling p-value from above:
library(ggpubr)
# # kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"="<0.01"))}
# kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"=\$resampling_p_value_for_plot))}

q\$mindist <- as.factor(q\$mindist)
p <- ggplot(q, aes(x=mindist, y=value, colour=class)) + stat_compare_means(method="wilcox.test", label="p.signif", hide.ns=T) + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=pretty_breaks(5)) + geom_boxplot(notch=T) + scale_colour_manual(values=c('black', 'grey'), aesthetics=c("colour", "fill"), guide=F) + theme_minimal() + xlab("Maximum distance considered [nt]") + ylab("Fraction of binding events in proximity")
ggsave("$outpdf", width=91.5, height=91.5, units="mm", device = cairo_pdf)
ggsave("$outpng", width=91.5, height=91.5, units="mm")

# warnings()

)));


foreach $close (@mindist_thresholds)
{
	$p = getr("p$close");
	
	print " >> $close nt >> $p\n";
	print OUT "$type\t$scoretype\t$close\t$p\n";
}

state("Wrote p-values to '$outtable'");
state("Wrote figure to '$outpdf'");
state("Wrote figure to '$outpng'");

done();
