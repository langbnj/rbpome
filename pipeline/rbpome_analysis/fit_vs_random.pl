#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# NOTE - Negative controls are added manually in fit_vs_random.pl

# eCLIP parameters (for getting Jaccard etc. from rbpome_analysis output files)
# $table = 'rbpome';
# $scoretype = 'SUM_IS';
# $threshold = 7.1;
# $type = 'eclip_encode';
# $hc = 0;
# $minlog2fold = 0;
# $mindist_threshold = 100;
$resamples = 100;

our $usage = "$0 [table (e.g. rbpome)] [type (e.g. eclip_encode_12)] [score type (e.g. SUM_IS)] [proximal binding cutoff] [-distance_signed] [-randompairs]\n\n -distance_signed: Retain sign (+/-) in distance plot instead of using absolute values (to probe whether binding is always upstream or downstream for a given pair)\n -randompairs: Read reproducible background resample random pairs from e.g. output/output-rbpome-randompairs.txt and use these instead of real RBP pairs\n\nExample: $0 rbpome eclip_encode_12 SUM_IS 54";
($table, $type, $scoretype, $close) = args(4);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

$tmprandompairs = '';
if (switch('randompairs'))
{
	$tmprandompairs = '-randompairs';
}

if (!switch('distance_signed'))
{
	$tmpsigned = '';
	$tmpsignedmin = 1;	# x-axis minimum for distance_plot
}
else
{
	$tmpsigned = '-signed';
	$tmpsignedmin = '-5e5';	# x-axis minimum for distance_plot
}



$outfile = "output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs.txt";
$outpdf = "output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs.pdf";
$outpdfdelta = "output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta.pdf";
# $outpdfdeltalim = "output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit.pdf";
# $randomfile = "output/output-txt-$table-$type-$scoretype-jaccard-hc0-minlog2fold0-minscore7.1-mindist_threshold50-resamples100.txt";
$randomfile = "output/output-$table-randompairs.txt";


startr();
runr(qq(

# setwd("~/Documents/Projects/RBPome/Peak distance plot comparison vs random")
library(ggplot2)
library(tibble)
library(forcats)
library(scales)
library(ggbeeswarm)

));


# If output file already exists, jump directly to plotting below
if ((!-s $outfile) or (switch('overwrite')))
{
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

	# print OUT "symbol1\tsymbol2\tks\tpvalue\tks_log\tpvalue_log\n";
	# print OUT "symbol1\tsymbol2\tavg_dist_real_lt54\tavg_dist_real_lt100\tavg_dist_real_lt117\tavg_dist_real_lt1000\tavg_dist_real\tavg_dist_random_lt54\tavg_dist_random_lt100\tavg_dist_random_lt117\tavg_dist_random_lt1000\tavg_dist_random\tmedian_dist_real_lt54\tmedian_dist_real_lt100\tmedian_dist_real_lt117\tmedian_dist_real_lt1000\tmedian_dist_real\tmedian_dist_random_lt54\tmedian_dist_random_lt100\tmedian_dist_random_lt117\tmedian_dist_random_lt1000\tmedian_dist_random\tclose_real\tfraction_close_real\tfraction_close_random\tfraction_close_delta\tfraction_close_ratio\tresampling_p\tresampling_wilcox\tresampling_p_sig\tresampling_wilcox_sig\n";
	print OUT "symbol1\tsymbol2\tavg_dist_real_lt54\tavg_dist_real_lt100\tavg_dist_real_lt117\tavg_dist_real_lt250\tavg_dist_real_lt400\tavg_dist_real_lt1000\tavg_dist_real\tavg_dist_random_lt54\tavg_dist_random_lt100\tavg_dist_random_lt117\tavg_dist_random_lt250\tavg_dist_random_lt400\tavg_dist_random_lt1000\tavg_dist_random\tmedian_dist_real_lt54\tmedian_dist_real_lt100\tmedian_dist_real_lt117\tmedian_dist_real_lt250\tmedian_dist_real_lt400\tmedian_dist_real_lt1000\tmedian_dist_real\tmedian_dist_random_lt54\tmedian_dist_random_lt100\tmedian_dist_random_lt117\tmedian_dist_random_lt250\tmedian_dist_random_lt400\tmedian_dist_random_lt1000\tmedian_dist_random\tclose_real\tfraction_close_real\tfraction_close_random\tfraction_close_delta\tfraction_close_ratio\tresampling_p\tresampling_wilcox\tresampling_p_sig\tresampling_wilcox_sig\n";


	# start

	@pairs = ();
	if (switch('randompairs'))
	{
		# open(RAND, $randomfile) or die("Error: Couldn't open '$randomfile'");
		# <RAND>;	# Skip header
		# startme("Reading random RBP pairs from '$randomfile'");
		# while (<RAND>)
		# {
		# 	chomp;
		#
		# 	# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	FMR1	FXR2	0.353067047075606
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	FXR1	FXR2	0.201612903225806
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	NONO	SFPQ	0.192367199587416
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	PCBP1	PTBP1	0.108222607356056
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	PTBP1	QKI	0.364448188711036
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	positives	U2AF1	U2AF2	0.2545532461967
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	APOBEC3C	HNRNPK	0.0594089264173703
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	APOBEC3C	PCBP1	0.110552763819095
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	APOBEC3C	PCBP2	0.0639899623588457
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	CPEB4	CSTF2T	0.0180102915951973
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	CPSF6	ILF3	0.0643253234750462
		# 	# rbpome	eclip_encode_12	SUM_IS	jaccard	0	0	7.1	50	100	screen_hits	CPSF6	NONO	0.0620272314674735
		#
		# 	@a = split(/\t/);
		# 	$set = $a[9];
		# 	$this_symbol1 = $a[10];
		# 	$this_symbol2 = $a[11];
		#
		# 	# Only use the first background resample
		# 	if ($set eq 'background_resample_1')
		# 	# # Use all background resample pairs
		# 	# if ($set =~ /^background_resample_\d+$/)
		# 	{
		# 		# foreach $pair ("$this_symbol1|$this_symbol2", "$this_symbol2|$this_symbol1")
		# 		push(@pairs, "$this_symbol1|$this_symbol2");
		# 		push(@pairs, "$this_symbol2|$this_symbol1");
		# 	}
		#
		# 	stepme(100);
		# 	stepme(100);
		# }
		# close(RAND);
		# stopme();
	
		open(RAND, "$randomfile") or die("Error: Couldn't open '$randomfile'");
		@pairs = ();
		while (<RAND>)
		{
			chomp;
			push(@pairs, $_);
		}
		close(RAND);
	
		@pairs = addflip(@pairs); # Add inverse pairs
		@pairs = unique(@pairs);

		state("Read ".scalar(@pairs)." unique reproducible random pairs from '$randomfile'");
	}
	else
	{
		startme("Getting real RBP pairs from table '$table'");
		$query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
		while (($this_symbol1, $this_symbol2) = Fetch($query))
		{
			push(@pairs, "$this_symbol1|$this_symbol2");
			push(@pairs, "$this_symbol2|$this_symbol1");

			stepme(100);
		}
		stopme();

		@pairs = unique(@pairs);
		state("Read ".scalar(@pairs)." unique pairs");
	}




	# $query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
	# #DEBUG avoiding RBFOX2 (not run yet, newly fixed its symbol so it is retained - it is an eCLIP protein)
	# warn("WARNING: DEBUG: Skipping RBFOX2 for now!");
	# warn("WARNING: DEBUG: Skipping NONO|SFPQ for now!");	# not run somehow?
	# warn("WARNING: DEBUG: Skipping SF3B4|SFPQ for now!");	# not run somehow?
	# $query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 AND symbol1!='RBFOX2' AND symbol2!='RBFOX2' AND NOT (symbol1='NONO' AND symbol2='SFPQ') AND NOT (symbol1='SF3B4' AND symbol2='SFPQ') GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
	# # $query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 AND symbol1 IN ('FMR1', 'U2AF1') AND symbol2 IN ('FXR2', 'U2AF2') GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
	# #END DEBUG
	# startme("Comparing close binding events (within ≤$close nt) in the real and random distance distributions", 0, Numrows($query) * 2);
	startme("Comparing close binding events (within ≤$close nt) in the real and random distance distributions", 0, scalar(@pairs));
	starttime();
	# while (($this_symbol1, $this_symbol2) = Fetch($query))
	# {
		# foreach $pair ("$this_symbol1|$this_symbol2", "$this_symbol2|$this_symbol1")
		foreach $pair (@pairs)
		{
			($symbol1, $symbol2) = split(/\|/, $pair);

			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-parameters.txt
			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-random-parameters.txt

			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-mindists.txt
			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-mindists-random.txt

			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-mindists.rds
			# output-distance-plot-closest_per_peak-SUM_IS-APOBEC3C-HNRNPK-fit-mindists-random.rds

			runr(qq(real <- readRDS("output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$symbol1-$symbol2-fit-mindists$tmpsigned.rds")));
			runr(qq(random <- readRDS("output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$symbol1-$symbol2-fit-mindists$tmpsigned-random.rds")));

			# state(runr(qq(print(str(real)))));
			# state(runr(qq(print(str(random)))));

			# $ks = getr('ks$statistic');
			# $pvalue = getr('ks$p.value');

			# A log transform (even log10(x + 1)) doesn't affect the KS test at all.
			# runr(qq(ks_log <- ks.test(log10(real + 1), log10(random + 1))));
			#
			# $ks_log = getr('ks_log$statistic');
			# $pvalue_log = getr('ks_log$p.value');
			# print OUT "$symbol1\t$symbol2\t$ks\t$pvalue\t$ks_log\t$pvalue_log\n";

			# Get number of close binding events
			runr(qq(close_real <- length(real[abs(real)<=$close])));
			runr(qq(fraction_close_real <- length(real[abs(real)<=$close]) / length(real)));
			# runr(qq(close_random <- length(random[random<=$close])));
			runr(qq(fraction_close_random <- length(random[abs(random)<=$close]) / length(random)));
		
			runr(qq(

			fails = 0
			wilcox_fails = 0
			for (i in 1:$resamples) {
				my_random <- sample(random, length(real), replace=T)

				close_my_random <- length(my_random[abs(my_random)<=$close])
				fraction_close_my_random <- length(my_random[abs(my_random)<=$close]) / length(my_random)

				if (fraction_close_my_random >= fraction_close_real) {
					fails <- fails + 1
				}

				if (wilcox.test(abs(real), abs(my_random), alternative="less")\$p.value >= 0.05) {
					wilcox_fails <- wilcox_fails + 1
				}
			}

			resampling_p <- fails / $resamples
			resampling_wilcox <- wilcox_fails / $resamples

			));

			$close_real = getr('close_real');
			$fraction_close_real = getr('fraction_close_real');
			# $close_random = getr('close_random');
			$fraction_close_random = getr('fraction_close_random');
			$resampling_p = getr('resampling_p');
			$resampling_wilcox = getr('resampling_wilcox');
			
			$resampling_p_sig = 'FALSE';
			$resampling_wilcox_sig = 'FALSE';
			$resampling_p_sig = 'TRUE' if ($resampling_p < 0.05);
			$resampling_wilcox_sig = 'TRUE' if ($resampling_wilcox < 0.05);

			$fraction_close_delta = $fraction_close_real - $fraction_close_random;	# Delta worked better than fraction in 2020-04-16 Sebastian mindist threshold PDF.pdf
			if ($fraction_close_random != 0)
			{
				$fraction_close_ratio = $fraction_close_real / $fraction_close_random;
			}
			else
			{
				$fraction_close_ratio = 'NA';
			}


			$avg_dist_real = getr('mean(real)');
			$avg_dist_random = getr('mean(random)');

			$avg_dist_real_lt54 = getr('mean(real[abs(real)<=54])');
			$avg_dist_random_lt54 = getr('mean(random[abs(random)<=54])');

			$avg_dist_real_lt100 = getr('mean(real[abs(real)<=100])');
			$avg_dist_random_lt100 = getr('mean(random[abs(random)<=100])');

			$avg_dist_real_lt117 = getr('mean(real[abs(real)<=117])');
			$avg_dist_random_lt117 = getr('mean(random[abs(random)<=117])');

			$avg_dist_real_lt250 = getr('mean(real[abs(real)<=250])');
			$avg_dist_random_lt250 = getr('mean(random[abs(random)<=250])');

			$avg_dist_real_lt400 = getr('mean(real[abs(real)<=400])');
			$avg_dist_random_lt400 = getr('mean(random[abs(random)<=400])');

			$avg_dist_real_lt1000 = getr('mean(real[abs(real)<=1000])');
			$avg_dist_random_lt1000 = getr('mean(random[abs(random)<=1000])');


			$median_dist_real = getr('median(real)');
			$median_dist_random = getr('median(random)');

			$median_dist_real_lt54 = getr('median(real[abs(real)<=54])');
			$median_dist_random_lt54 = getr('median(random[abs(random)<=54])');

			$median_dist_real_lt100 = getr('median(real[abs(real)<=100])');
			$median_dist_random_lt100 = getr('median(random[abs(random)<=100])');

			$median_dist_real_lt117 = getr('median(real[abs(real)<=117])');
			$median_dist_random_lt117 = getr('median(random[abs(random)<=117])');

			$median_dist_real_lt250 = getr('median(real[abs(real)<=250])');
			$median_dist_random_lt250 = getr('median(random[abs(random)<=250])');

			$median_dist_real_lt400 = getr('median(real[abs(real)<=400])');
			$median_dist_random_lt400 = getr('median(random[abs(random)<=400])');

			$median_dist_real_lt1000 = getr('median(real[abs(real)<=1000])');
			$median_dist_random_lt1000 = getr('median(random[abs(random)<=1000])');


			# print OUT "$symbol1\t$symbol2\t$ks\t$pvalue\t$close_real\t$close_random\t$fraction_close_real\t$fraction_close_random\n";
			# print OUT "$symbol1\t$symbol2\t$avg_dist_real_lt54\t$avg_dist_real_lt100\t$avg_dist_real_lt117\t$avg_dist_real_lt1000\t$avg_dist_real\t$avg_dist_random_lt54\t$avg_dist_random_lt100\t$avg_dist_random_lt117\t$avg_dist_random_lt1000\t$avg_dist_random\t$median_dist_real_lt54\t$median_dist_real_lt100\t$median_dist_real_lt117\t$median_dist_real_lt1000\t$median_dist_real\t$median_dist_random_lt54\t$median_dist_random_lt100\t$median_dist_random_lt117\t$median_dist_random_lt1000\t$median_dist_random\t$close_real\t$fraction_close_real\t$fraction_close_random\t$fraction_close_delta\t$fraction_close_ratio\t$resampling_p\t$resampling_wilcox\t$resampling_p_sig\t$resampling_wilcox_sig\n";
			print OUT "$symbol1\t$symbol2\t$avg_dist_real_lt54\t$avg_dist_real_lt100\t$avg_dist_real_lt117\t$avg_dist_real_lt250\t$avg_dist_real_lt400\t$avg_dist_real_lt1000\t$avg_dist_real\t$avg_dist_random_lt54\t$avg_dist_random_lt100\t$avg_dist_random_lt117\t$avg_dist_random_lt250\t$avg_dist_random_lt400\t$avg_dist_random_lt1000\t$avg_dist_random\t$median_dist_real_lt54\t$median_dist_real_lt100\t$median_dist_real_lt117\t$median_dist_real_lt250\t$median_dist_real_lt400\t$median_dist_real_lt1000\t$median_dist_real\t$median_dist_random_lt54\t$median_dist_random_lt100\t$median_dist_random_lt117\t$median_dist_random_lt250\t$median_dist_random_lt400\t$median_dist_random_lt1000\t$median_dist_random\t$close_real\t$fraction_close_real\t$fraction_close_random\t$fraction_close_delta\t$fraction_close_ratio\t$resampling_p\t$resampling_wilcox\t$resampling_p_sig\t$resampling_wilcox_sig\n";


			# Make a plot
			runr(qq(

			# setwd("~/Documents/Projects/RBPome/Peak distance plot comparison vs random")

			table <- "$table"
			type <- "$type"
			scoretype <- "$scoretype"
			symbol1 <- "$symbol1"
			symbol2 <- "$symbol2"

			myo <- "#f3702a"  # My Orange
			mybd <- "#003D51" # My Blue Dark
			# myb <- "#1e3d59"  # My Blue
			apg <- "#bec1c0"

			));



			runr(q(

			prevwd <- getwd()
			# setwd("/Volumes/blang/pipeline/rbpome_analysis/output")
			real <- readRDS(paste0("output/output-distance-plot-closest_per_peak-",table,"-",type,"-",scoretype,"-",symbol1,"-",symbol2,"-fit-mindists.rds"))
			random <- readRDS(paste0("output/output-distance-plot-closest_per_peak-",table,"-",type,"-",scoretype,"-",symbol1,"-",symbol2,"-fit-mindists-random.rds"))
			# setwd(prevwd)

			# head(real)
			# head(random)
			#
			# str(real)
			# str(random)
			#
			# summary(real)
			# summary(random)


			# q <- rbind(tibble(category="Observed", value=real), tibble(category="Control", value=random))
			# str(q)
			# head(q)

			# # q <- rbind(tibble(category="Observed peak distances", value=real), tibble(category="Randomised peak positions", value=random))
			# q <- rbind(tibble(category="Observed\npeak distances", value=real), tibble(category="Randomised\npeak distances", value=random))
			# # ggplot(q, aes(x=category, y=value)) + geom_boxplot(notch=T)
			# ggplot(q, aes(x=value, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance [nt]") + ylab("")
			# ggplot(q, aes(x=value, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance [nt]") + ylab("") + theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1))
			# ggplot(q, aes(x=value, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance [nt]") + ylab("") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
			#
			# q <- rbind(tibble(category="Observed", value=real), tibble(category="Randomised", value=random))
			# ggplot(q, aes(x=value, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

			# Subsample the control to be the same size
			real_tib <- tibble(category="Observed", value=real)
			random_tib <- tibble(category="Randomised", value=sample(random, nrow(real_tib), replace=T))
			q <- rbind(real_tib, random_tib)

			));

			runr(qq(

			# Plot

			# ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2))
			# ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest peak [nt]") + ylab("") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5))
			# ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest peak [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2))
			# p <- ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T) + scale_x_log10(limits=c(1, NA), breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest peak [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5))
			# p <- ggplot(q, aes(x=fct_rev(as.factor(category)), y=value + 1)) + geom_boxplot(notch=T) + scale_y_log10(limits=c(1, NA), breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Peak distances") + ylab("Distance to closest peak [nt]") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
			# p <- ggplot(q, aes(x=fct_rev(as.factor(category)), y=value + 1)) + geom_boxplot(notch=T, outlier.shape=NA) + scale_y_log10(limits=c(1, NA), breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Peak distances") + ylab("Distance to closest peak [nt]") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
			# p <- ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)))) + geom_boxplot(notch=T, outlier.shape=NA) + scale_x_log10(limits=c(1, NA), breaks = trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest peak [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5))
			# p <- ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)), fill=fct_rev(as.factor(category)))) + geom_boxplot(notch=T, outlier.shape=NA) + scale_colour_manual(values=c(mybd, myo), aesthetics=c("colour", "fill"), guide=F) + coord_cartesian(xlim=c(1, 5e5)) + scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest peak [nt]") + ylab("Peak distances") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol1, " vs. ", symbol2)) + theme(plot.title = element_text(hjust = 0.5))
			# Box plots
			p <- ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)), fill=fct_rev(as.factor(category)))) + geom_boxplot(notch=T, outlier.shape=NA) + geom_vline(xintercept=$close, linetype="dashed") + scale_colour_manual(values=c(apg, myo), aesthetics=c("colour", "fill"), guide=F) + coord_cartesian(xlim=c($tmpsignedmin, 5e5)) + scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest binding site [nt]") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol2, " near ", symbol1)) + theme(plot.title = element_text(hjust = 0.5))
			ggsave(paste0("output/output-distance-plot-closest_per_peak-",table,"-",type,"-",scoretype,"-mindist_threshold$close-",symbol1,"-",symbol2,"$tmpsigned-real_vs_random.pdf"), width=91.5, height=45.75, unit="mm")
			# Violin plots
			p <- ggplot(q, aes(x=value + 1, y=fct_rev(as.factor(category)), fill=fct_rev(as.factor(category)))) + geom_violin(draw_quantiles = c(0.5)) + geom_vline(xintercept=$close, linetype="dashed") + scale_colour_manual(values=c(apg, myo), aesthetics=c("colour", "fill"), guide=F) + coord_cartesian(xlim=c($tmpsignedmin, 5e5)) + scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x), labels=comma) + theme_minimal() + xlab("Distance to closest binding site [nt]") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(paste0(symbol2, " near ", symbol1)) + theme(plot.title = element_text(hjust = 0.5))
			ggsave(paste0("output/output-distance-plot-closest_per_peak-",table,"-",type,"-",scoretype,"-mindist_threshold$close-",symbol1,"-",symbol2,"$tmpsigned-real_vs_random.pdf"), width=91.5, height=45.75, unit="mm")

			));


			stepme(1);
		}
	# }
	close(OUT);
	stopme();
	stoptime();

	showmeall(1);

	state("Wrote to '$outfile'");
}
else
{
	state("Output file '$outfile' already existed, reading it (use -overwrite to overwrite it)");
}





# Make overview box plot (close binding events in the real data vs. in the random distribution)
state(runr(qq(

# Make overview box plot (close binding events in the real data vs. in the random distribution)
# library(ggplot2)
# library(tibble)
# library(forcats)
# library(scales)
# library(ggbeeswarm)

# setwd("~/Dropbox/RBPome_paper/FIGURES/Figure4-New RBP interactions predict RBP co-binding in proximity/Panel e (mindist distributions for examples)")

table <- read.table("$outfile", header=T, sep="\t", quote="")

str(table)
head(table)

q <- rbind(tibble(class="Observed", value=table\$fraction_close_real), tibble(class="Randomised", value=table\$fraction_close_random))
str(q)
head(q)

# ggplot(q, aes(x=class, y=value)) + geom_boxplot(notch=T, outlier.shape = NA) + coord_cartesian(ylim=c(NA, 0.45)) + theme_minimal() + xlab("") + ylab("Fraction of close-binding events (≤$close nt)")
# ggplot(q, aes(x=value, y=fct_rev(as.factor(class)))) + geom_boxplot(notch=T, outlier.shape = NA) + coord_cartesian(xlim=c(NA, 0.45)) + theme_minimal() + ylab("") + xlab("Fraction of close-binding events (≤$close nt)")
p <- ggplot(q, aes(x=value, y=fct_rev(as.factor(class)), colour=class)) + coord_cartesian(xlim=c(0, 1)) + scale_x_continuous(breaks=pretty_breaks(5)) + geom_boxplot(notch=T) + scale_colour_manual(values=c('black', 'grey'), aesthetics=c("colour", "fill"), guide=F) + theme_minimal() + ylab("") + xlab("Fraction of close-binding events (≤$close nt)")
ggsave("$outpdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# warnings()

wilcox.test(table\$fraction_close_real, table\$fraction_close_random)

# Plot difference between median_dist_real and median_dist_random
table\$median_dist_delta <- table\$median_dist_random - table\$median_dist_real
q <- tibble(value=table\$median_dist_delta)
str(q)
head(q)

# p <- ggplot(q, aes(x=value)) + coord_cartesian(xlim=c($tmpsignedmin, 5e5)) + geom_boxplot(notch=T) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")

# Boxplot
p <- ggplot(q, aes(x=value)) + geom_boxplot(notch=T) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("$outpdfdelta", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot, filtered to significant pairs
p <- ggplot(subset(table, resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit "close"
p <- ggplot(subset(table, abs(median_dist_real) <= $close), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-$close, $close)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit$close.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit "close", filtered to significant pairs
p <- ggplot(subset(table, abs(median_dist_real) <= $close & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-$close, $close)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit$close-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 100, filtered to significant pairs
p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit100-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 1000
p <- ggplot(subset(table, abs(median_dist_real) <= 1000), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-1000, 1000)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit1000.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)

# Boxplot limit 1000, filtered to significant pairs
p <- ggplot(subset(table, abs(median_dist_real) <= 1000 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_boxplot(notch=T) + coord_cartesian(xlim=c(-1000, 1000)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-limit1000-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)



# Violin limit 100, filtered to significant pairs
# p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=NA)) + geom_violin() + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real)) + geom_density() + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=0)) + geom_violin() + geom_beeswarm(groupOnX=T) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
# p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(x=median_dist_real, y=factor(""))) + geom_violin() + geom_point(position="jitter", alpha=0.3) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
p <- ggplot(subset(table, abs(median_dist_real) <= 100 & resampling_p<0.05 & resampling_wilcox<0.05), aes(y=median_dist_real)) + geom_beeswarm(groupOnX=T) + coord_cartesian(xlim=c(-100, 100)) + theme_minimal() + ylab("") + xlab("Median distance difference [nt]")
ggsave("output-fit_vs_random-$table-$type-$scoretype-$close$tmpsigned$tmprandompairs-delta-violin-limit100-filtered.pdf", width=91.5, height=45.75, units="mm", device = cairo_pdf)



)));




done();











# Define utility function: addflip
sub addflip
{
	my @pos = @_;
	
	# print "before: ".scalar(@pos)."\n";
	# show(@pos);
	my @new = ();
	foreach $pair (@pos)
	{
		my ($rbp1, $rbp2) = split(/\|/, $pair);
		push(@new, "$rbp1|$rbp2");
		push(@new, "$rbp2|$rbp1");
	}
	@pos = @new;
	# print "after: ".scalar(@pos)."\n";
	# show(@pos);
	
	# die("Error: Duplicates found during addflip (shouldn't happen)") if (scalar(@pos) != scalar(unique(@pos)));
	
	return(@pos);
}
