#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $table = 'rbpome';

# eCLIP parameters
$type = 'eclip_encode_12';
$hc = 0;
$minlog2fold = 0;
$mindist_threshold = 50;
$resamples = 100;

our $usage = "$0 [table (e.g. 'rbpome')]\n\nExample: $0 rbpome";
($table) = args(1);

$infile = "input/merged_nanobret_results_07-Feb-2020221237_ben_all.txt";
$outfile = "output-$table-compare-scoretypes.txt";
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "source\tscoretype\tpositives\tnegatives\ttested\tvalidation_ratio\tTP\tFP\tTN\tFN\tF1\tMCC\teclip_TP\teclip_FP\teclip_TN\teclip_FN\teclip_F1\teclip_MCC\teclip_pairs\tavg_jaccard\tavg_oddsratio\tavg_probability_mindist\tmedian_jaccard\tmedian_oddsratio\tmedian_probability_mindist\n";


# start

# Get fullsymbol-symbol mapping
starttime();
$query = Query("SELECT DISTINCT fullsymbol1, symbol1, eclip1 FROM $table");
startme("Getting fullsymbol-symbol mapping from table '$table' (step 1/2)", 0, Numrows($query));
%symbol = ();
%eclip = ();
while (($fullsymbol, $symbol, $eclip) = Fetch($query))
{
	$fullsymbol = uc($fullsymbol);
	$symbol = uc($symbol);
	
	if (exists($symbol{$fullsymbol}))
	{
		if ($symbol{$fullsymbol} ne $symbol)
		{
			die("Error: Multiple symbols for fullsymbol '$fullsymbol': '$symbol{$fullsymbol}' and '$symbol'");
		}
	}

	$symbol{$fullsymbol} = $symbol;
	$eclip{$symbol} = 1 if ($eclip == 1);

	stepme(100);
}
stopme();

$query = Query("SELECT DISTINCT fullsymbol2, symbol2, eclip2 FROM $table");
startme("Getting fullsymbol-symbol mapping from table '$table' (step 2/2)", 0, Numrows($query));
while (($fullsymbol, $symbol, $eclip) = Fetch($query))
{
	$fullsymbol = uc($fullsymbol);
	$symbol = uc($symbol);

	if (exists($symbol{$fullsymbol}))
	{
		if ($symbol{$fullsymbol} ne $symbol)
		{
			die("Error: Multiple symbols for fullsymbol '$fullsymbol': '$symbol{$fullsymbol}' and '$symbol'");
		}
	}

	$symbol{$fullsymbol} = $symbol;
	$eclip{$symbol} = 1 if ($eclip == 1);

	stepme(100);
}
stopme();

state("Got mapping for ".commify(scalar(keys(%symbol)))." fullsymbols");
state("Got ".scalar(keys(%eclip))." eCLIP symbols");




# Getting NanoBRET data
%tested = ();
%positive = ();
%negative = ();
startme("Reading '$infile'");
<IN>;	# Skip header
while (<IN>)
{
	# pair_ID	p1	p2	nb_score	nb_std	merged_pairnames	classification_zscore095	Index
	# R_59	cpeb2	cpeb2	47.78505355	1.537508115	cpeb2_cpeb2	positive	1
	# J103_02	kif16b_cter	ncbp3	0.951010652	0.327433707	kif16b_cter_ncbp3	negative	2
	# R_48	srsf9	rbmx	39.77018046	1.315552926	srsf9_rbmx	positive	3
	# J103_04	tsnax	bicd2	1.247319423	0.282866362	tsnax_bicd2	negative	4
	
	@a = split(/\t/);
	$fullsymbol1 = $a[1];
	$fullsymbol2 = $a[2];
	$verdict = $a[6];

	$fullsymbol1 = uc($fullsymbol1);
	$fullsymbol2 = uc($fullsymbol2);
	
	if (!exists($symbol{$fullsymbol1}))
	{
		# die("Error: No symbol for fullsymbol '$fullsymbol1'");
		addme("nanobret pair skipped because no symbol for fullsymbol (this also indicates the protein wasn't in the screen) (skipped)", $fullsymbol1);
		addme("nanobret pair skipped because no symbol for fullsymbol for pair (this also indicates the protein wasn't in the screen) (skipped)", "$fullsymbol1|$fullsymbol2");
		next;
	}
	if (!exists($symbol{$fullsymbol2}))
	{
		# die("Error: No symbol for fullsymbol '$fullsymbol2'");
		addme("nanobret pair skipped because no symbol for fullsymbol (this also indicates the protein wasn't in the screen) (skipped)", $fullsymbol2);
		addme("nanobret pair skipped because no symbol for fullsymbol for pair (this also indicates the protein wasn't in the screen) (skipped)", "$fullsymbol1|$fullsymbol2");
		next;
	}
	
	# Translate fullsymbol into symbol
	$symbol1 = $symbol{$fullsymbol1};
	$symbol2 = $symbol{$fullsymbol2};
	
	# Tested
	if ($symbol1 lt $symbol2)
	{
		# "Correct pair"
		$tested{"$symbol1|$symbol2"} = 1;
	}
	else
	{
		# "Wrong pair"
		$tested{"$symbol2|$symbol1"} = 1;
	}
	
	# Positives
	if ($verdict eq 'positive')
	{
		if ($symbol1 lt $symbol2)
		{
			# "Correct pair"
			$positive{"$symbol1|$symbol2"} = 1;
		}
		else
		{
			# "Wrong pair"
			$positive{"$symbol2|$symbol1"} = 1;
		}
	}
	elsif ($verdict eq 'negative')
	{
		if ($symbol1 lt $symbol2)
		{
			# "Correct pair"
			$negative{"$symbol1|$symbol2"} = 1;
		}
		else
		{
			# "Wrong pair"
			$negative{"$symbol2|$symbol1"} = 1;
		}
	}
	else
	{
		die;
	}
	
	stepme(100);
}
close(IN);
stopme();





# Cycle through score types and get NanoBRET validation ratios
state("Getting NanoBRET validation ratios:");
$typequery = Query("SELECT DISTINCT scoretype, threshold, source FROM $table ORDER BY scoretype");
%screen_hit = ();
while (($scoretype, $threshold, $source) = Fetch($typequery))
{
	print " >> $scoretype\n";
	
	# $outfile = "output-$scoretype.txt";
	# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	
	$mainquery = Query("SELECT DISTINCT symbol1, symbol2, eclip1, eclip2 FROM $table WHERE scoretype='$scoretype' AND hc=1 AND homodimer=0");
	$total_tested = 0;
	$total_positive = 0;
	$total_negative = 0;
	$TP = 0;
	$FP = 0;
	$TN = 0;
	$FN = 0;
	$eclip_TP = 0;
	$eclip_FP = 0;
	$eclip_TN = 0;
	$eclip_FN = 0;
	while (($symbol1, $symbol2, $eclip1, $eclip2) = Fetch($mainquery))
	{
		# Screen hit
		$screen_hit{"$symbol1|$symbol2"} = 1;
		
		# Tested
		$tested = 0;
		if (exists($tested{"$symbol1|$symbol2"}))
		{
			$tested = 1;
		}
		
		# Negative
		$negative = 0;
		if (exists($negative{"$symbol1|$symbol2"}))
		{
			$negative = 1;
		}
			
		# A positive result overrides negatives
		$positive = 0;
		if (exists($positive{"$symbol1|$symbol2"}))
		{
			$positive = 1;
			$negative = 0;
		}
		
		if ($positive == 1)
		{
			$TP++;
			$eclip_TP++ if (($eclip1 == 1) and ($eclip2 == 1));
		}
		elsif ($negative == 1)
		{
			$FP++;
			$eclip_FP++ if (($eclip1 == 1) and ($eclip2 == 1));
		}
		
		$total_tested += $tested;
		$total_positive += $positive;
		$total_negative += $negative;
	}
	
	# Get TN/FN by cycling through the NanoBRET-tested pairs that weren't screen hits
	foreach $pair (keys(%tested))
	{
		($symbol1, $symbol2) = split(/\|/, $pair);
		
		if (!exists($screen_hit{$pair}))
		{
			# Negative
			$negative = 0;
			if (exists($negative{$pair}))
			{
				$negative = 1;
			}
			
			# A positive result overrides negatives
			$positive = 0;
			if (exists($positive{$pair}))
			{
				$positive = 1;
				$negative = 0;
			}
			
			if ($positive == 1)
			{
				$FN++;
				$eclip_FN++ if (exists($eclip{$symbol1}) and exists($eclip{$symbol2}));
			}
			elsif ($negative == 1)
			{
				$TN++;
				$eclip_TN++ if (exists($eclip{$symbol1}) and exists($eclip{$symbol2}));
			}
		}
	}
	
	$validation = $total_positive / $total_tested;
	
	print "   >> $total_positive / $total_tested\n";
	print "     >> ".sprintf("%.1f", $validation * 100)."%\n\n";
	
	# Calculate F1 score and Matthew's Correlation Coefficient
	$f1 = (2 * $TP) / ((2 * $TP) + $FP + $FN);
	$mcc = (($TP * $TN) - ($FP * $FN)) / sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN));
	
	# Calculate F1 score and Matthew's Correlation Coefficient for eCLIP pairs only
	$eclip_f1 = (2 * $eclip_TP) / ((2 * $eclip_TP) + $eclip_FP + $eclip_FN);
	$eclip_mcc = (($eclip_TP * $eclip_TN) - ($eclip_FP * $eclip_FN)) / sqrt(($eclip_TP + $eclip_FP) * ($eclip_TP + $eclip_FN) * ($eclip_TN + $eclip_FP) * ($eclip_TN + $eclip_FN));
	
	# Round values
	$validation = round($validation, 3);
	$f1 = round($f1, 3);
	$mcc = round($mcc, 3);
	$eclip_f1 = round($eclip_f1, 3);
	$eclip_mcc = round($eclip_mcc, 3);
	
	
	# Get eclip_pairs number (how many screen hits are eCLIP-eCLIP interactions?)
	$query = Query("SELECT COUNT(DISTINCT symbol1, symbol2) FROM $table WHERE scoretype='$scoretype' AND homodimer=0 AND hc=1 AND eclip1=1 AND eclip2=1");
	($eclip_pairs) = FetchOne($query);
	
	
	# Get Jaccard / oddsratio / probability_mindist from ~/pipeline/rbpome_analysis/output/ table files
	# Average
	$mean_jaccard = get_cobind('jaccard', 'mean');
	$mean_oddsratio = get_cobind('oddsratio', 'mean');
	$mean_probability_mindist = get_cobind('probability_mindist', 'mean');
	# Median
	$median_jaccard = get_cobind('jaccard', 'median');
	$median_oddsratio = get_cobind('oddsratio', 'median');
	$median_probability_mindist = get_cobind('probability_mindist', 'median');
	
	# Round values
	# Average
	$mean_jaccard = round($mean_jaccard, 3);
	$mean_oddsratio = round($mean_oddsratio, 3);
	$mean_probability_mindist = round($mean_probability_mindist, 3);
	# Median
	$median_jaccard = round($median_jaccard, 3);
	$median_oddsratio = round($median_oddsratio, 3);
	$median_probability_mindist = round($median_probability_mindist, 3);
	
	
	print OUT "$source\t$scoretype\t$total_positive\t$total_negative\t$total_tested\t$validation\t$TP\t$FP\t$TN\t$FN\t$f1\t$mcc\t$eclip_TP\t$eclip_FP\t$eclip_TN\t$eclip_FN\t$eclip_f1\t$eclip_mcc\t$eclip_pairs\t$mean_jaccard\t$mean_oddsratio\t$mean_probability_mindist\t$median_jaccard\t$median_oddsratio\t$median_probability_mindist\n";
}
nl();
stoptime();



showmeall(1);

done();




sub get_cobind
{
	($mode, $median) = @_;
	
	$tmpfile = "../rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
	open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
	<TMP>;	# Skip header
	@cobind = ();
	while (<TMP>)
	{
		chomp;
		
		# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	PCBP1	0.16151090849886
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	SF3B4	0.135072908672295
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FMR1	FXR2	0.420327027893556
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FXR1	FXR2	0.15188679245283
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	HNRNPK	U2AF2	0.234877126654064
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	NONO	SFPQ	0.203641590800192
	
		@a = split(/\t/);
		$set = $a[9];
		$cobind = $a[12];
		
		# Only use "screen_hits" (not positive controls, not background_resamples)
		next if ($set ne 'screen_hits');
		
		push(@cobind, $cobind);
	}
	close(TMP);
	
	if ($median eq 'median')
	{
		$res = median(@cobind);
	}
	else
	{
		$res = mean(@cobind);
	}

	return($res);
}





