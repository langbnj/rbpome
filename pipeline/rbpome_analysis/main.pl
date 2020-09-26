#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

# $table = 'rbpome';
# $mindist_threshold = 100;	# Minimum distance where we consider it co-binding as a complex
# $mindist_threshold = 50;
$mindist_threshold_strict = 54;	# "strict": This is the exact peak of the <=1000 positives curve in Figure 4e
$mindist_threshold_loose = 117;	# "looser": This is 95% of the peak of the <=1000 positives curve in Figure 4e

our $usage = "$0 [table (e.g. rbpome)] [eclip type: e.g. eclip_encode_12] [interaction screen score type: AVG_IS/SUM_IS/11_IS/F1/MCC] [analysis mode: jaccard/correlation/mindist/mediandist/mindist_thresh/mindist_thresh_frac] [hc: 0/1] [minlog2fold: e.g. 3/4/5] [avgIS: e.g. 1/3.25/9] [mindist threshold: e.g. 10/100/1000] [resamples: e.g. 20/100/1000] [-distance_plot] [-distance_signed] [-nofit] [-overwrite] [-random] [-randompairs]\n\n -distance_plot: Plot a diagnostic peak distance plot and exit (for rationally identifying a good mindist threshold).\n -distance_signed: Retain sign (+/-) in -distance_plot instead of using absolute values (to probe whether binding is always upstream or downstream for a given pair)\n -motif_plot: Plot motif_fraction (the fraction of all peaks that contain motifs) in cobound peaks (closest_per_peak).\n -nofit: With -distance_plot: Do not perform bimodal fit (much faster)\n -overwrite: Overwrite previous output files.\n -random: Use random 'peak' coordinates from table 'clip_cobinding_random'.\n -randompairs: Use random background resample pairs for -distance_plot, instead of real RBP pairs.\n\nExample: $0 rbpome eclip_encode_12 SUM_IS jaccard 0 0 7.1 $mindist_threshold_strict 100\n";
($table, $type, $scoretype, $mode, $hc, $minlog2fold, $minscore, $mindist_threshold, $resamples) = args(9);
#args(0);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

die("Error: -motif_plot requires -distance_plot to be on, too") if (switch('motif_plot') and !switch('distance_plot'));
die("Error: -randompairs requires -distance_plot to be on, too") if (switch('randompairs') and !switch('distance_plot'));

# Add mode axis labels (mode descriptions)
if    ($mode eq "fraction")							{ $axis = "Fraction of genes bound by both RBPs (relative to the smaller set)"; }
elsif ($mode eq "freq")								{ $axis = "Gene binding frequency for either of a pair of RBPs (A binding or B binding, divided by the total number of genes)"; }
elsif ($mode eq "intersection")						{ $axis = "Absolute number of RNA genes bound by both RBPs"; }
elsif ($mode eq "jaccard")							{ $axis = "RNA target set similarity (Jaccard index)"; }
elsif ($mode eq "mediandist")						{ $axis = "Minimum binding site distance for a pair of RBPs (median across their shared target genes)"; }
elsif ($mode eq "mindist")							{ $axis = "Minimum binding site distance for a pair of RBPs (averaged across their shared target genes)"; }
elsif ($mode eq "mindist_thresh")					{ $axis = "Absolute number of RNA genes bound by both RBPs within a minimum binding distance i.e. probably as a complex"; }
elsif ($mode eq "mindist_thresh_frac")				{ $axis = "Absolute number of RNA genes bound by both RBPs within a minimum binding distance as a fraction of their shared target genes"; }
elsif ($mode eq "oddsratio")						{ $axis = "Gene binding odds ratio for a pair of RBPs"; }
elsif ($mode eq "probability")						{ $axis = "Gene co-binding probability for a pair of RBPs p(A|B)"; }
elsif ($mode eq "probability_freqscaled")			{ $axis = "Gene co-binding probability for a pair of RBPs p(A|B) multiplied by the probability of neither RBP binding"; }
elsif ($mode eq "probability_mindist")				{ $axis = "Gene co-binding probability for a pair of RBPs p(A|B)"; }
elsif ($mode eq "probability_mindist_freqscaled")	{ $axis = "Gene co-binding probability for a pair of RBPs p(A|B) multiplied by the probability of neither RBP binding"; }
else
{
	die("Error: Unhandled mode '$mode'");
}


if ($type eq 'clip_postar')
{
	# POSTAR
	$clip_raw_gene = 'clip_postar_raw';
	$clip_gene = 'clip_postar';
}
else
{
	# Any other CLIP
	$clip_raw_gene = 'clip_raw_gene';
	$clip_gene = 'clip_gene';
}

if (!switch('distance_signed'))
{
	$tmpsigned = '';
	$tmpsignedmin = 0;	# x-axis minimum for distance_plot
}
else
{
	$tmpsigned = '-signed';
	$tmpsignedmin = -1000;	# x-axis minimum for distance_plot
}


if (!switch('random'))
{
	# $clip_cobinding = 'clip_cobinding';
	$tmprandom = '';
}
else
{
	# Random "peaks"
	# $clip_cobinding = 'clip_cobinding_random';
	if (!switch('random_balanced'))
	{
		$tmprandom = '-random';
	}
	else
	{
		$tmprandom = '-random_balanced';
	}
}

# Define mode directionality (if 1, add "inverse pairs" via addflip() below)
if    ($mode eq 'fraction')							{ $mode_directional = 0; }
elsif ($mode eq 'freq') 							{ $mode_directional = 0; }
elsif ($mode eq 'intersection') 					{ $mode_directional = 0; }
elsif ($mode eq 'jaccard') 							{ $mode_directional = 0; }
elsif ($mode eq 'mediandist') 						{ $mode_directional = 0; }
elsif ($mode eq 'mindist') 							{ $mode_directional = 0; }
elsif ($mode eq 'mindist_thresh') 					{ $mode_directional = 0; }
elsif ($mode eq 'mindist_thresh_frac') 				{ $mode_directional = 0; }
elsif ($mode eq 'oddsratio') 						{ $mode_directional = 0; }
elsif ($mode eq 'probability') 						{ $mode_directional = 1; }
elsif ($mode eq 'probability_freqscaled') 			{ $mode_directional = 1; }
elsif ($mode eq 'probability_mindist') 				{ $mode_directional = 1; }
elsif ($mode eq 'probability_mindist_freqscaled') 	{ $mode_directional = 1; }
else
{
	die("Error: Mode directionality not yet defined for mode '$mode'");
}
# # Everything needs to be directional, for the comparisons! Implementing this in compare.pl though to avoid reprocessing.
# $mode_directional = 1;

# $infile = "../input/screen_input_proteins_with_eCLIPdata.txt";
$outfile = "../output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$pdffile = "../output/output-pdf-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.pdf";
$pvalfile = "../output/output-pvalues-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$hitsfile = "../output/output-hits-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.txt";

starttime2();
if ((!-s $outfile) or (switch('debug')) or (switch('debug2')) or (switch('distance_plot')) or (switch('motif_plot')) or (switch('overwrite')))
{
	# start



	# Get RBPs that were Y2H-screened
	# startme("Reading RBP list from '$infile'");
	startme("Getting RBPs from table '$table'");
	@rbps_both = ();
	$mainquery = Query("SELECT DISTINCT symbol1 FROM $table WHERE scoretype='$scoretype'");
	while (($symbol) = Fetch($mainquery))
	{
		stepme(1000);

		# Check if it has eCLIP data
		if ($type ne 'clip_postar')
		{
			if (($hc == 0) and ($minlog2fold == 0))
			{
				# No thresholds: use clip_raw_gene directly
				$query = Query("SELECT id FROM $clip_raw_gene WHERE type='$type' AND symbol='$symbol' LIMIT 1");
			}
			else
			{
				$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$symbol' AND hc>=$hc AND log2fold>=$minlog2fold LIMIT 1");
			}
		}
		else
		{
			$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND symbol='$symbol' AND hc>=$hc LIMIT 1");
		}
		if (Numrows($query) == 0)
		{
			# die("Error: '$symbol' is not an eCLIP RBP");
			# This can happen now because of eclip_tom, which doesn't have peaks for some RBPs.
			addme("not an eclip rbp for symbol", $symbol);
			next;
		}
		
		# Check if it has motif data if -motif_plot is active
		if (switch('motif_plot'))
		{
			$query = Query("SELECT id FROM rbp_motifs_extend5_fimo_$type WHERE species='human' AND type='$type' AND symbol='$symbol' LIMIT 1");
			if (Numrows($query) == 0)
			{
				addme("no known motif for symbol (skipped because of -motif_plot)", $symbol);
				next;
			}
		}
	
		# These both have eCLIP data and are in Sebastian Maurer's screen.
		push(@rbps_both, $symbol);
	}
	$mainquery = Query("SELECT DISTINCT symbol2 FROM $table WHERE scoretype='$scoretype'");
	while (($symbol) = Fetch($mainquery))
	{
		stepme(1000);

		# Check if it has eCLIP data
		if ($type ne 'clip_postar')
		{
			if (($hc == 0) and ($minlog2fold == 0))
			{
				# No thresholds: use clip_raw_gene directly
				$query = Query("SELECT id FROM $clip_raw_gene WHERE type='$type' AND symbol='$symbol' LIMIT 1");
			}
			else
			{
				$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$symbol' AND hc>=$hc AND log2fold>=$minlog2fold LIMIT 1");
			}
		}
		else
		{
			$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND symbol='$symbol' AND hc>=$hc LIMIT 1");
		}
		if (Numrows($query) == 0)
		{
			# die("Error: '$symbol' is not an eCLIP RBP");
			# This can happen now because of eclip_tom, which doesn't have peaks for some RBPs.
			addme("not an eclip rbp for symbol", $symbol);
			next;
		}
	
		# Check if it has motif data if -motif_plot is active
		if (switch('motif_plot'))
		{
			$query = Query("SELECT id FROM rbp_motifs_extend5_fimo_$type WHERE species='human' AND type='$type' AND symbol='$symbol' LIMIT 1");
			if (Numrows($query) == 0)
			{
				addme("no known motif for symbol (skipped because of -motif_plot)", $symbol);
				next;
			}
		}
	
		# These both have eCLIP data and are in Sebastian Maurer's screen.
		push(@rbps_both, $symbol);
	}
	@rbps_both = unique(@rbps_both);
	stopme();
	
	if (!switch('motif_plot'))
	{
		state("With eCLIP data: ".scalar(@rbps_both));
		die("Error: No RBPs with eCLIP data of type '$type' are in the screen") if (scalar(@rbps_both) == 0);
	}
	else
	{
		state("With eCLIP and motif data: ".scalar(@rbps_both));
		die("Error: No RBPs with eCLIP and motif data of type '$type' are in the screen") if (scalar(@rbps_both) == 0);
	}


	# #DEBUG
	# @rbps_both = ('FXR2', 'FMR1') if (switch('debug2') or switch('distance_plot'));
	# @rbps_both = ('FMR1', 'FXR1', 'FXR2') if (switch('debug2') or switch('distance_plot'));
	# @rbps_both = ('U2AF1', 'U2AF2') if (switch('debug2') or switch('distance_plot'));
	# @rbps_both = ('CSTF2T', 'EWSR1') if (switch('debug2') or switch('distance_plot'));
	# @rbps_both = ('FMR1', 'U2AF1', 'SLBP', 'FASTKD2', 'FXR1', 'FXR2', 'U2AF1', 'U2AF2', 'CSTF2T', 'EWSR1') if (switch('debug2') or switch('distance_plot'));
	@rbps_both = ('FMR1', 'U2AF1', 'SLBP', 'FASTKD2', 'FXR1', 'FXR2', 'U2AF1', 'U2AF2', 'CSTF2T', 'EWSR1') if (switch('debug2'));
	# #END DEBUG
	
	if (switch('distance_plot'))
	{
		# Filter @rbps_both to contain only RBPs from @rbps_both that interact with another one (only eCLIP-protein complexes, and possibly filtered for RBPs with known motifs as well)
		# $query = Query("SELECT symbol1, symbol2 FROM `$table` WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
		$query = Query("SELECT symbol1, symbol2 FROM `$table` WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 AND symbol1 IN ('".join("', '", @rbps_both)."') AND symbol2 IN ('".join("', '", @rbps_both)."') GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
		@rbps_both = ();
		while (($symbol1, $symbol2) = Fetch($query))
		{
			push(@rbps_both, $symbol1);
			push(@rbps_both, $symbol2);
		}
		@rbps_both = unique(@rbps_both);
	}

	# # Negative control (these two don't interact in the screen, of course)
	# # @rbps_both = ('FASTKD2', 'SLBP') if (switch('debug3'));
	# push(@rbps_both, 'FASTKD2');
	# push(@rbps_both, 'SLBP');
	


	# Get positive control
	@pos = ();
	# push(@pos, "fmr1|fxr1");
	# push(@pos, "fmr1|fxr2");
	# push(@pos, "fxr1|fxr2");
	# push(@pos, "nono|sfpq");
	# push(@pos, "pcbp1|ptbp1");
	# push(@pos, "u2af1|u2af2");
	# Direct query, which gets slightly more (12 pairs):
	# $query = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE eclip1=1 AND eclip2=1 AND bg_alldirect=1 AND homodimer=0 AND scoretype='$scoretype'");
	# our $superloudmysql = 1;
	$mainquery = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND avg_is>=$minscore AND homodimer=0 AND bg_alldirect=1");
	$workaround = 0;
	if (Numrows($mainquery) == 0)
	{
		if ((uc($scoretype) ne '11_IS') and (uc($scoretype) ne 'SUM_IS'))
		{
			die("Error: No positives for score type '$scoretype' in table '$table' (they should have homodimer=0 and bg_alldirect=1)");
		}
		else
		{
			# Jae-Seong's tables don't have bg_alldirect and homodimer information, so I'll get positive pairs from Sebastian's original table (AVG_IS) instead
			$mainquery = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE scoretype='AVG_IS' AND homodimer=0 AND bg_alldirect=1");
			warn("WARNING: Using workaround via AVG_IS to get positive controls (there weren't any with homodimer=0 AND bg_alldirect=1 for scoretype '$scoretype')");
			$workaround = 1;
		}
	}
	while (($rbp1, $rbp2) = Fetch($mainquery))
	{
		# # Skip any pairs that only have been found by one method (require at least 2)
		# $query = Query("SELECT COUNT(DISTINCT system) FROM biogrid WHERE gene1='$rbp1' AND gene2='$rbp2' AND direct=1");
		# ($systemcount) = FetchOne($query);
		# if ($systemcount < 2)
		# {
		# 	state("positive pair skipped because it was only found using a single direct method (requiring 2) for pair (skipped): $rbp1|$rbp2");
		# 	addme("positive pair skipped because it was only found using a single direct method (requiring 2) for pair (skipped)", "$rbp1|$rbp2");
		# 	next;
		# }
		# Skip any pairs that only have been found by one study (require at least 2)
		$query = Query("SELECT COUNT(DISTINCT pmid) FROM biogrid WHERE gene1='$rbp1' AND gene2='$rbp2' AND direct=1");
		($pmidcount) = FetchOne($query);
		if ($pmidcount < 2)
		{
			# state("positive pair skipped because it was only found using a single direct study (requiring 2) for pair (skipped): $rbp1|$rbp2");
			addme("positive pair skipped because it was only found using a single direct study (requiring 2) for pair (skipped)", "$rbp1|$rbp2");
			next;
		}
		
		
		if ($workaround == 1)
		{
			# Check if this pair is present for this score type (else skip)
			$query = Query("SELECT id FROM $table WHERE scoretype='$scoretype' AND symbol1='$rbp1' AND symbol2='$rbp2'");
			if (Numrows($query) == 0)
			{
				addme("WARNING: pair skipped because of AVG_IS workaround", "$rbp1|$rbp2");
				next;
			}
		}
		
		# Check if rbp1 is an eCLIP RBP (of type $type)
		if ($type ne 'clip_postar')
		{
			if (($hc == 0) and ($minlog2fold == 0))
			{
				# No thresholds: use clip_raw_gene directly
				$query = Query("SELECT id FROM $clip_raw_gene WHERE type='$type' AND symbol='$rbp1' LIMIT 1");
			}
			else
			{
				$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$rbp1' AND hc>=$hc AND log2fold>=$minlog2fold LIMIT 1");
			}
		}
		else
		{
			$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND symbol='$rbp1' AND hc>=$hc LIMIT 1");
		}
		if (Numrows($query) == 0)
		{
			addme("not an eclip rbp for symbol", $rbp1);
			next;
		}
		# Check if rbp2 is an eCLIP RBP (of type $type)
		if ($type ne 'clip_postar')
		{
			if (($hc == 0) and ($minlog2fold == 0))
			{
				# No thresholds: use clip_raw_gene directly
				$query = Query("SELECT id FROM $clip_raw_gene WHERE type='$type' AND symbol='$rbp2' LIMIT 1");
			}
			else
			{
				$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$rbp2' AND hc>=$hc AND log2fold>=$minlog2fold LIMIT 1");
			}
		}
		else
		{
			$query = Query("SELECT id FROM $clip_gene WHERE species='human' AND symbol='$rbp2' AND hc>=$hc LIMIT 1");
		}
		if (Numrows($query) == 0)
		{
			addme("not an eclip rbp for symbol", $rbp2);
			next;
		}
	
		# Keep
		push(@pos, "$rbp1|$rbp2");
	}
	@pos = unique(@pos);
	@pos = addflip(@pos) if ($mode_directional == 1); # Add inverse pairs for cobind_probability (which is directional, p(A|B) and p(B|A)) and the other directional analyses

	# DEBUG
	state("POSITIVES:");
	show(@pos);
	# exit;
	# END DEBUG
	

	# Get screen hits (foreground) (skipping homomers)
	# our $superloudmysql = 1;
	$query = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND avg_is>=$minscore AND homodimer=0 AND symbol1 IN ('".join("', '", @rbps_both)."') AND symbol2 IN ('".join("', '", @rbps_both)."')");
	# our $superloudmysql = 0;
	startme("Getting screen hits (foreground) from table '$table'", 0, Numrows($query));
	@fg = ();
	while (($symbol1, $symbol2) = Fetch($query))
	{
		push(@fg, "$symbol1|$symbol2");
	
		stepme(100);
	}
	@fg = unique(@fg);
	@fg = addflip(@fg) if ($mode_directional == 1); # Add inverse pairs for cobind_probability (which is directional, p(A|B) and p(B|A)) and the other directional analyses
	stopme();






	# Get a hash of ENSGs per RBP (incorporating 'hc' and log2fold minimum)
	starttime();
	# $query = Query("SELECT symbol, ensgv FROM clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol IN ('".join("', '", @rbps_both)."') AND hc>=$hc AND log2fold>=$minlog2fold");
	if ($type ne 'clip_postar')
	{
		if (($hc == 0) and ($minlog2fold == 0))
		{
			# No thresholds: use clip_raw_gene directly
			$query = Query("SELECT DISTINCT symbol, ensgv FROM $clip_raw_gene WHERE type='$type' AND symbol IN ('".join("', '", @rbps_both)."')");
		}
		else
		{
			$query = Query("SELECT DISTINCT symbol, ensgv FROM $clip_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol IN ('".join("', '", @rbps_both)."') AND hc>=$hc AND log2fold>=$minlog2fold");
		}
	}
	else
	{
		$query = Query("SELECT DISTINCT c.symbol, t.ensgv FROM $clip_gene c, gencode_gff3_transcript t WHERE c.enstv=t.enstv AND t.species='human' AND c.species='human' AND c.symbol IN ('".join("', '", @rbps_both)."') AND c.hc>=$hc");
	}
	startme("Getting target genes per RBP (incorporating 'hc' and log2fold minimum)", 0, Numrows($query));
	%ensgvs = ();
	%total_ensgvs = ();
	@total_ensgvs = ();
	while (($symbol, $ensgv) = Fetch($query))
	{
		# # Convert symbol to lower case
		# $symbol = lc($symbol);

		# Add to hash
		if (exists($ensgvs{$symbol}))
		{
			$ensgvs{$symbol} .= "|$ensgv";
		}
		else
		{
			$ensgvs{$symbol} = $ensgv;
		}
		$total_ensgvs{$ensgv} = 1;
		# push(@total_ensgvs, $ensgv);
		stepme(100000);
	}
	@total_ensgvs = unique(keys(%total_ensgvs));
	$total_ensgvs = scalar(@total_ensgvs);
	stopme();
	stoptime();
	
	state("Total RNA ENSGVs: ".commify($total_ensgvs));

	if (!switch('nosort'))
	{
		# Sort & unique genes
		startme("Sorting and uniquing target genes", 0, scalar(keys(%ensgvs)));
		starttime();
		$ensgvs = 0;
		foreach $key (keys(%ensgvs))
		{
			$ensgvs{$key} = join("|", unique(split(/\|/, $ensgvs{$key})));
			$ensgvs += scalar(split(/\|/, $ensgvs{$key}));
			stepme(100000);
		}
		stopme();
		state("Target gene count: ".commify($ensgvs), 1);
		stoptime();
	}




	# If in 'mindist' or 'mediandist' mode: Get peak 5' ends per RBP and gene
	# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
	if (($mode eq 'mindist') or ($mode eq 'mediandist') or ($mode eq 'mindist_thresh') or ($mode eq 'mindist_thresh_frac') or ($mode eq 'probability_mindist') or ($mode eq 'probability_mindist_freqscaled'))
	{
		# Get peak 5' ends per RBP and target gene
		starttime();
		if ($type ne 'clip_postar')
		{
			$query = Query("SELECT symbol, ensgv, chr, start, stop, strand FROM $clip_raw_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol IN ('".join("', '", @rbps_both)."') AND log2fold>=$minlog2fold");
		}
		else
		{
			$query = Query("SELECT c.symbol, t.ensgv, c.chr, c.start, c.stop, c.strand FROM $clip_raw_gene c, gencode_gff3_transcript t WHERE c.enstv=t.enstv AND t.species='human' AND c.species='human' AND c.symbol IN ('".join("', '", @rbps_both)."')");
		}
		startme("Getting peaks per RBP and target gene", 0, Numrows($query));
		%peaks = ();
		# Store ENSGVs for %motifpeaks below
		%ensgvs_for_peak = ();
		while (($symbol, $ensgv, $chr, $start, $stop, $strand) = Fetch($query))
		{
			# Convert symbol to lower case
			# $symbol = lc($symbol);
			
			# Store ENSGVs for %motifpeaks below
			if (exists($ensgvs_for_peak{"$chr|$start|$stop|$strand"}))
			{
				$ensgvs_for_peak{"$chr|$start|$stop|$strand"} .= "|$ensgv";
			}
			else
			{
				$ensgvs_for_peak{"$chr|$start|$stop|$strand"} = $ensgv;
			}
			

			# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
			# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
			# $peak_center = $start + (($stop - $start) / 2);

			if ($strand eq '+')
			{
				$peak = $start;
			}
			else
			{
				$peak = $stop;
			}
			
			if (switch('random'))
			{
				# Get gene coordinates
				$genequery = Query("SELECT start, stop FROM gencode_gff3_gene WHERE species='human' AND ensgv='$ensgv'");
				($gene_start, $gene_stop) = FetchOne($genequery);
				
				foreach $i (1..$resamples)
				{
					# Set a random "peak" within the gene
					$peak = $gene_start + int(rand(($gene_stop - $gene_start) + 1));

					# Add to hash
					if (exists($peaks{"$symbol|$ensgv|$i"}))
					{
						$peaks{"$symbol|$ensgv|$i"} .= "|$peak";
					}
					else
					{
						$peaks{"$symbol|$ensgv|$i"} = $peak;
					}
				}
			}
			else
			{
				# Add to hash
				if (exists($peaks{"$symbol|$ensgv"}))
				{
					$peaks{"$symbol|$ensgv"} .= "|$peak";
				}
				else
				{
					$peaks{"$symbol|$ensgv"} = $peak;
				}
			}
			
			stepme(100000);
		}
		stopme();
		stoptime();
	
		if (!switch('nosort'))
		{
			# Sort & unique peaks
			startme("Sorting and uniquing peaks", 0, scalar(keys(%peaks)));
			starttime();
			$peaks = 0;
			foreach $key (keys(%peaks))
			{
				$peaks{$key} = join("|", unique(split(/\|/, $peaks{$key})));
				$peaks += scalar(split(/\|/, $peaks{$key}));
				stepme(100000);
			}
			stopme();
			state("Peak count: ".commify($peaks), 1);
			stoptime();
		}

		
		
		# Get peaks with motifs, if needed (only needed for motif_plot, which is nested within distance_plot (so both switches need to be on))
		# Can't randomise these
		if (switch('motif_plot'))
		{
			# Get peak 5' ends per RBP and target gene
			starttime();
			if ($type ne 'clip_postar')
			{
				# $query = Query("SELECT symbol, ensgv, start, stop, strand FROM $clip_raw_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol IN ('".join("', '", @rbps_both)."') AND log2fold>=$minlog2fold");
				# $query = Query("SELECT symbol, chr, start, stop, strand FROM rbp_motifs_extend5_fimo_$type WHERE species='human' AND type='$type' AND symbol IN ('".join("', '", @rbps_both)."') AND psig=1 GROUP BY start, stop, strand ORDER BY start, stop, strand");
				$query = Query("SELECT DISTINCT symbol, chr, start, stop, strand FROM rbp_motifs_extend5_fimo_$type WHERE species='human' AND type='$type' AND symbol IN ('".join("', '", @rbps_both)."') AND psig=1");
			}
			else
			{
				# $query = Query("SELECT c.symbol, t.ensgv, c.start, c.stop, c.strand FROM $clip_raw_gene c, gencode_gff3_transcript t WHERE c.enstv=t.enstv AND t.species='human' AND c.species='human' AND c.symbol IN ('".join("', '", @rbps_both)."')");
				die("Error: 'clip_postar' not implemented here yet");
			}
			startme("Getting peaks with motifs per RBP and target gene", 0, Numrows($query));
			%motifpeaks = ();
			while (($symbol, $chr, $start, $stop, $strand) = Fetch($query))
			{
				
				# Remove 5' extension like in analyse.R:
			    # q[q$strand=="+",]$start <- q[q$strand=="+",]$start + 50
			    # # - strand
			    # q[q$strand=="-",]$stop <- q[q$strand=="-",]$stop - 50
				
				# Remove 5' extension so the coordinates in rbp_motifs_extend5_... here match up with clip_raw_gene
				if ($strand eq '+')
				{
					$start += 50;
				}
				else
				{
					$stop -= 50;
				}
				
				
				# Convert symbol to lower case
				# $symbol = lc($symbol);

				# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
				# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
				# $peak_center = $start + (($stop - $start) / 2);

				if ($strand eq '+')
				{
					$peak = $start;
				}
				else
				{
					$peak = $stop;
				}
			
				# Can't randomise these
				# if (switch('random'))
				# {
				# 	# Get gene coordinates
				# 	$genequery = Query("SELECT start, stop FROM gencode_gff3_gene WHERE species='human' AND ensgv='$ensgv'");
				# 	($gene_start, $gene_stop) = FetchOne($genequery);
				#
				# 	foreach $i (1..$resamples)
				# 	{
				# 		# Set a random "peak" within the gene
				# 		$peak = $gene_start + int(rand(($gene_stop - $gene_start) + 1));
				#
				# 		# Add to hash
				# 		if (exists($motifpeaks{"$symbol|$ensgv|$i"}))
				# 		{
				# 			$motifpeaks{"$symbol|$ensgv|$i"} .= "|$peak";
				# 		}
				# 		else
				# 		{
				# 			$motifpeaks{"$symbol|$ensgv|$i"} = $peak;
				# 		}
				# 	}
				# }
				# else
				# {


					# Find out ensgvs for this peak
					@ensgvs_for_peak = unique(split(/\|/, $ensgvs_for_peak{"$chr|$start|$stop|$strand"}));
				
					foreach $ensgv (@ensgvs_for_peak)
					{
						# Add to hash
						if (exists($motifpeaks{"$symbol|$ensgv"}))
						{
							$motifpeaks{"$symbol|$ensgv"} .= "|$peak";
						}
						else
						{
							$motifpeaks{"$symbol|$ensgv"} = $peak;
						}
					}


				# }
			
				stepme(100000);
			}
			stopme();
			stoptime();
	
			if (!switch('nosort'))
			{
				# Sort & unique peaks
				startme("Sorting and uniquing peaks with motifs", 0, scalar(keys(%motifpeaks)));
				starttime();
				$motifpeaks = 0;
				foreach $key (keys(%motifpeaks))
				{
					$motifpeaks{$key} = join("|", unique(split(/\|/, $motifpeaks{$key})));
					$motifpeaks += scalar(split(/\|/, $motifpeaks{$key}));
					stepme(100000);
				}
				stopme();
				state("Motif peak count: ".commify($motifpeaks), 1);
				stoptime();
			}
		}


	}
	elsif (($mode eq 'fraction') or ($mode eq 'freq') or ($mode eq 'intersection') or ($mode eq 'jaccard') or ($mode eq 'oddsratio') or ($mode eq 'probability') or ($mode eq 'probability_freqscaled'))
	{
		# Don't need peak coordinates
	}
	else
	{
		die("Error: Not sure if mode '$mode' requires eCLIP peak coordinates (not handled yet)");
	}
	






	# Analysis
	
	if (switch('distance_plot'))
	{
		# distance_plot('CSTF2T|EWSR1');
		# distance_plot('EWSR1|CSTF2T');
		#
		# distance_plot('FMR1|U2AF1');
		# distance_plot('U2AF1|FMR1');
		
		if (switch('debug3'))
		{
			state("Distance plots: FASTKD2 vs. SLBP...", 1);
			distance_plot('FASTKD2|SLBP');
			state("Distance plots: SLBP vs. FASTKD2 (reverse)...", 1);
			distance_plot('SLBP|FASTKD2');
			exit;
		}

		# distance_plot('FXR1|FXR2');
		# distance_plot('FXR2|FXR1');
		#
		# distance_plot('FMR1|FXR2');
		# distance_plot('FXR2|FMR1');
		#
		# distance_plot('FMR1|FXR1');
		# distance_plot('FXR1|FMR1');
		#
		# distance_plot('U2AF1|U2AF2');
		# distance_plot('U2AF2|U2AF1');
		
		# $query = Query("SELECT symbol1, symbol2 FROM `$table` WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
		# our $superloudmysql=1;
		
		if (switch('randompairs'))
		{
			# Use random pairs
			
			# show(@rbps_both);
			# show(@fg);
			# state("PEAKS:");
			# show(%peaks);
			# exit;
			# Get values for randomly resampled backgrounds (from the rbps_both list)
			
			# # Resample $resamples times
			# print " >> Resampling background ($resamples times) (".scalar(@fg)." pairs)\n";
			# @bg_all = ();
			# for ($i = 1; $i<=$resamples; $i++)
			# {
			# 	# Get random list of RBP pairs (as many as in the foreground, scalar(@fg))
			# 	@bg = resample(scalar(@fg));
			# 	@bg = addflip(@bg) if ($mode_directional == 1); # Add inverse pairs for cobind_probability (which is directional, p(A|B) and p(B|A)) and the other directional analyses
			#
			# 	push(@bg_all, @bg);
			# }
			# state("pre-unique: ".scalar(@bg_all));
			# @bg_all = unique(@bg_all);
			# state("post-unique: ".scalar(@bg_all));
			# show(@bg_all);
			
			# Resample once
			print " >> Resampling background (once) (".scalar(@fg)." pairs)\n";
			@bg_all = ();
			# Get random list of RBP pairs (as many as in the foreground, scalar(@fg))
			# @bg = resample(scalar(@fg));
			if ($mode_directional == 1)
			{
				# Only use half since the foreground has already had inverse pairs added to it
				@bg = resample_without_replacement(scalar(@fg) / 2);
			}
			else
			{
				@bg = resample_without_replacement(scalar(@fg));
			}
			# # Add negative controls
			# push(@bg, "FASTKD2|SLBP");
			# # End negative controls
			@bg = addflip(@bg) if ($mode_directional == 1); # Add inverse pairs for cobind_probability (which is directional, p(A|B) and p(B|A)) and the other directional analyses
			# Unique
			@bg = unique(@bg);
			
			
			# Check how many screen hit pairs the random pairs contain
			state("Background (".scalar(@bg)." random pairs) contains ".scalar(intersection(\@bg, \@fg))." pairs that are also in the foreground (".scalar(@fg)." screen hit pairs)");
			# exit;
			
			
			# Write to file for later reproducibility in fit_vs_random.pl
			$randomfile = "../output/output-$table-randompairs.txt";
			if (-s $randomfile)
			{
				# Read it in if it exists
				open(RAND, "$randomfile") or die("Error: Couldn't open '$randomfile'");
				@bg = ();
				while (<RAND>)
				{
					chomp;
					push(@bg, $_);
				}
				close(RAND);
				
				state("Read ".scalar(@bg)." reproducible random pairs from '$randomfile'");
			}
			else
			{
				# Write it out if it doesn't
				open(RAND, ">$randomfile") or die("Error: Couldn't open '$randomfile'");
				foreach $pair (@bg)
				{
					print RAND "$pair\n";
				}
				close(RAND);

				state("Wrote ".scalar(@bg)." reproducible random pairs to '$randomfile'");
			}

			push(@bg_all, @bg);
			state("pre-unique: ".scalar(@bg_all));
			@bg_all = unique(@bg_all);
			state("post-unique: ".scalar(@bg_all));
			# show(@bg_all);

			@pairs = @bg_all;
			
			# # Read random pairs reproducibly from the jaccard mindist_threshold50 resampling file
			# $randomfile = "../output/output-txt-$table-$type-$scoretype-jaccard-hc0-minlog2fold0-minscore7.1-mindist_threshold50-resamples100.txt";
			# @pairs = ();
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

			# # Add negative controls
			# push(@pairs, "FASTKD2|SLBP");
			# push(@pairs, "SLBP|FASTKD2");
			# # End negative controls

			@pairs = unique(@pairs);
			state("Read ".scalar(@pairs)." unique pairs");
		}
		else
		{
			startme("Getting real RBP pairs from table '$table'");
			$query = Query("SELECT symbol1, symbol2 FROM `$table` WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 AND symbol1 IN ('".join("', '", @rbps_both)."') AND symbol2 IN ('".join("', '", @rbps_both)."') GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
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
		
		
		# $query = Query("SELECT symbol1, symbol2 FROM `$table` WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 AND symbol1 IN ('".join("', '", @rbps_both)."') AND symbol2 IN ('".join("', '", @rbps_both)."') GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
		# # $superloudmysql=0;
		# while (($symbol1, $symbol2) = Fetch($query))
		@tmppairs = unique(addflip(@pairs));
		foreach $pair (@tmppairs)
		{
			($symbol1, $symbol2) = split(/\|/, $pair);
			
			addme("distance plot made for symbol", $symbol1);
			addme("distance plot made for symbol", $symbol2);
			
			state("Distance plots: $symbol1 vs. $symbol2...", 1);
			distance_plot("$symbol1|$symbol2");
			# state("Distance plots: $symbol2 vs. $symbol1 (reverse)...", 1);
			# distance_plot("$symbol2|$symbol1");
		}
		
		exit;
	}

	# #DEBUG
	# cobind_mindist('fmr1|fxr1') if (switch('debug'));
	# # cobind_mindist_thresh('fmr1|fxr1') if (switch('debug'));
	# exit if (switch('debug'));
	# #END DEBUG

	state(" >> Creating output file '$outfile'");
	
	# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	open(PVAL, ">$pvalfile") or die("\nError: Couldn't open '$pvalfile'\n\n");
	print OUT "table\ttype\tscoretype\tmode\thc\tminlog2fold\tminscore\tmindist_threshold\tresamples\tset\trbp1\trbp2\tcobind\n";
	print PVAL "table\ttype\tscoretype\tmode\thc\tminlog2fold\tminscore\tmindist_threshold\tresamples\tvalues\tp_value\tresampling_p_value\n";
	if ($mode eq 'mindist_thresh')
	{
		# For 'mindist_thresh' mode only:
		# Write a table listing RBP1-RBP2-ensgv for very close binding sites (i.e. those below $mindist_threshold)
		open(HITS, ">$hitsfile") or die("\nError: Couldn't open '$hitsfile'\n\n");
		print HITS "table\ttype\tscoretype\tmode\thc\tminlog2fold\tminscore\tmindist_threshold\tresamples\tset\trbp1\trbp2\tensgv\tmindist\n";
	}



	# Get values for the positive control
	@posv = ();
	foreach $pair (@pos)
	{
		my ($rbp1, $rbp2) = split(/\|/, $pair);
		# $cobind = &{"cobind_$mode"}($pair);
		$cobind = &{"cobind_$mode"}($pair, 'positives');
		if (defined($cobind))
		{
			push(@posv, $cobind);
			print OUT "$table\t$type\t$scoretype\t$mode\t$hc\t$minlog2fold\t$minscore\t$mindist_threshold\t$resamples\tpositives\t$rbp1\t$rbp2\t$cobind\n";
		}
	}




	# Get values for the screen hits (foreground)
	@fgv = ();
	foreach $pair (@fg)
	{
		my ($rbp1, $rbp2) = split(/\|/, $pair);
		# if ($mode eq 'mindist_thresh')
		# {
		# 	# print "\n\ncobind_$mode(\"$pair|screen_hits\")\n\n";
		# 	$cobind = &{"cobind_$mode"}($pair, 'screen_hits');
		# }
		# else
		# {
		# 	$cobind = &{"cobind_$mode"}($pair);
		# }
		$cobind = &{"cobind_$mode"}($pair, 'screen_hits');
		if (defined($cobind))
		{
			# Keep value in any case for the resampling p-value below
			push(@fgv, $cobind);

			# Avoiding duplicates in the box plots:
			# Remove "pos" from the "fg" set
			# @fg = lonly(\@fg, \@pos);
			# next if (contains($pair, @pos));
			# if (contains($pair, @pos)) { state("duplicate POS pair skipped for FG: $pair", 1) if (switch('debug')); next; }
			print OUT "$table\t$type\t$scoretype\t$mode\t$hc\t$minlog2fold\t$minscore\t$mindist_threshold\t$resamples\tscreen_hits\t$rbp1\t$rbp2\t$cobind\n";
		}
	}

	if (scalar(@fgv) == 0)
	{
		run("touch", "touch $pdffile");
		# die("Error: No 'Screen hits' values (wrote empty output file)");
		die(" >> No 'Screen hits' values (wrote empty output file)");
		exit;
	}

	# Get values for randomly resampled backgrounds (from the rbps_both list)
	@resamples = ();
	print " >> Resampling background ($resamples times) (".scalar(@fg)." pairs)\n";
	@bgv_all = ();
	for ($i = 1; $i<=$resamples; $i++)
	{
		# Get random list of RBP pairs (as many as in the foreground, scalar(@fg))
		if ($mode_directional == 1)
		{
			# Only use half since the foreground has already had inverse pairs added to it
			@bg = resample(scalar(@fg) / 2);
		}
		else
		{
			@bg = resample(scalar(@fg));
		}
		# #DEBUG (adding this since it's already described in the text)
		# push(@bg, "EIF3G|RPS3");
		# push(@bg, "RPS3|EIF3G");
		# #END DEBUG
		@bg = addflip(@bg) if ($mode_directional == 1); # Add inverse pairs for cobind_probability (which is directional, p(A|B) and p(B|A)) and the other directional analyses
		# # Unique
		# @bg = unique(@bg);
		# Check how many screen hit pairs the random pairs contain
		state("     >> Background (".scalar(@bg)." random pairs) contains ".scalar(intersection(\@bg, \@fg))." pairs that are also in the foreground (".scalar(@fg)." screen hit pairs)", 1);
		# exit;


		# state("BG:") if (switch('debug'));
		# show(@bg) if (switch('debug'));
		@bgv = ();
		foreach $pair (@bg)
		{
			my ($rbp1, $rbp2) = split(/\|/, $pair);
			# $cobind = &{"cobind_$mode"}($pair);
			$cobind = &{"cobind_$mode"}($pair, "background_resample_$i");
			if (defined($cobind))
			{
				# Keep values in any case for the resampling p-value below
				push(@bgv, $cobind);
				push(@bgv_all, $cobind);

				# Avoiding duplicates in the box plots:
				# Remove "pos" and "fg" from the background set
				# @bgv_all = lonly(\@bgv_all, \@posv);
				# @bgv_all = lonly(\@bgv_all, \@fgv);
				# next if (contains($pair, @pos));
				# next if (contains($pair, @fg));
				if (contains($pair, @pos)) { state("duplicate POS pair skipped for BG: $pair", 1) if (switch('debug')); next; }
				if (contains($pair, @fg)) { state("duplicate FG pair skipped for BG: $pair", 1) if (switch('debug')); next; }
				print OUT "$table\t$type\t$scoretype\t$mode\t$hc\t$minlog2fold\t$minscore\t$mindist_threshold\t$resamples\tbackground_resample_$i\t$rbp1\t$rbp2\t$cobind\n";
			}
		}
	
		# These should be LOW
		if (($mode eq 'mindist') or ($mode eq 'mediandist'))
		{
			# If the random resample results in an equal or even in a lower average value than for the screen hits, this is a fail:
			if (avg(@bgv) <= avg(@fgv))
			{
				push(@resamples, 1);
				print "   >> $i / $resamples >> ".avg(@fgv)." vs. ".avg(@bgv)." >> 1 (fail)\n";
			}
			else
			{
				# If the random resample results in a lower average value, success!
				push(@resamples, 0);
				print "   >> $i / $resamples >> ".avg(@fgv)." vs. ".avg(@bgv)." >> 0 (SUCCESS!)\n";
			}
		}
		# The rest should be HIGH
		# elsif (($mode eq 'fraction') or ($mode eq 'jaccard') or ($mode eq 'mindist_thresh') or ($mode eq 'mindist_thresh_frac') or ($mode eq 'probability'))
		else
		{
			# If the random resample results in an equal or even in a higher average value than for the screen hits, this is a fail:
			if (avg(@bgv) >= avg(@fgv))
			{
				push(@resamples, 1);
				print "   >> $i / $resamples >> ".avg(@fgv)." vs. ".avg(@bgv)." >> 1 (fail)\n";
			}
			else
			{
				# If the random resample results in a lower average value, success!
				push(@resamples, 0);
				print "   >> $i / $resamples >> ".avg(@fgv)." vs. ".avg(@bgv)." >> 0 (SUCCESS!)\n";
			}
		}
	}

	# Get resampling p-value
	$resampling_p_value = avg(@resamples);

	# Get Wilcoxon p-value (screen hits (foreground) vs. combined background across resamples)
	# These should be LOW
	if (($mode eq 'mindist') or ($mode eq 'mediandist'))
	{
		$p_value = wilcoxon(\@fgv, \@bgv_all, 0, 'less');
	}
	# These should be HIGH
	# if (($mode eq 'fraction') or ($mode eq 'jaccard') or ($mode eq 'mindist_thresh') or ($mode eq 'mindist_thresh_frac'))
	else
	{
		$p_value = wilcoxon(\@fgv, \@bgv_all, 0, 'greater');
	}

	# Show p-values
	nl();
	state("Resampling p-value: $resampling_p_value", 1);
	state("Wilcoxon p-value: $p_value", 1);
	nl();


	# Print both p-values to p-value file
	print PVAL "$table\t$type\t$scoretype\t$mode\t$hc\t$minlog2fold\t$minscore\t$mindist_threshold\t$resamples\t".scalar(@fgv)."\t$p_value\t$resampling_p_value\n";


	# Make plots
	# if (scalar(@posv) == 0)
	# {
	# 	run("touch", "touch $pdffile");
	# 	die("Error: No 'Positive control' values (wrote empty output file)");
	# }
	# if (scalar(@fgv) == 0)
	# {
	# 	run("touch", "touch $pdffile");
	# 	die("Error: No 'Screen hits' values (wrote empty output file)");
	# }
	# if (scalar(@bgv_all) == 0)
	# {
	# 	run("touch", "touch $pdffile");
	# 	die("Error: No 'Random background' values (wrote empty output file)");
	# }

	# setr('pos', @posv);
	# setr('fg', @fgv);
	# setr('bg', @bgv_all);
	close(OUT);
}
else
{
	print " >> Output file '$outfile' already exists, reading it & plotting...\n";
}



$mindist_thresholds = $mindist_threshold;
$minscores = $minscore;
if (switch('2C') or switch('2C_old') or switch('2D') or switch('2D_old'))
{
	# Get all $mindist_thresholds instead of just the $mindist_threshold given as a parameter
	# $mindist_thresholds = 'c(10, 20, 40, 50, 100, 400, 1000, 10000)';
	$mindist_thresholds = "c(10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000)";
}
if (switch('2F') or switch('2G'))
{
	# Get all $mindist_thresholds instead of just the $mindist_threshold given as a parameter
	# $minscores = 'c(0, 1, 2, 3.25, 9, 20)';
	# $minscores = 'c(0, 1.07, 2.39, 3.072, 4.38, 5.4, 7.1, 10)';
	$minscores = 'c(0, 1.07, 2.39, 3.072, 4.38, 5.4, 7.1, 10)';
}



state(runr(qq(

library(ggplot2)
library(ggbeeswarm)
library(Cairo)
library(naturalsort)
library(scales)

# str(pos)
# str(fg)
# str(bg)


# Read p-values
for (outfile in paste0("../output/output-pvalues-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore",$minscores,"-mindist_threshold",$mindist_thresholds,"-resamples$resamples.txt")) {
	if (!exists("pv", mode="list")) {
		pv <- read.table(outfile, header=T, quote="", sep="\\t")
		print(paste0("Reading ", outfile))
	} else {
		pv <- rbind(pv, read.table(outfile, header=T, quote="", sep="\\t"))
		print(paste0("Reading ", outfile))
	}
}
p_value <- pv[1,]\$p_value
resampling_p_value <- pv[1,]\$resampling_p_value


# Read table(s)
for (outfile in paste0("../output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore",$minscores,"-mindist_threshold",$mindist_thresholds,"-resamples$resamples.txt")) {
	if (!exists("q", mode="list")) {
		q <- read.table(outfile, header=T, quote="", sep="\\t")
		print(paste0("Reading ", outfile))
	} else {
		q <- rbind(q, read.table(outfile, header=T, quote="", sep="\\t"))
		print(paste0("Reading ", outfile))
	}
}
# print(str(q))
# print(unique(q\$mindist_threshold))
q\$$mode <- q\$cobind
# print(head(q))
# pos <- subset(tmpq, set=="positives")\$cobind
# fg <- subset(tmpq, set=="screen_hits")\$cobind
# bg <- subset(tmpq, (set != "positives") & (set != "screen_hits"))\$cobind

# # Filter positive and foreground duplicates out from the background for the plots
# # print(length(unique(pos)))
# # print(length(unique(fg)))
# # print(length(unique(bg)))
# # print(length(pos))
# # print(length(fg))
# # print(length(bg))
# # print(str(bg))
# bg_filtered <- bg[!(bg \%in% c(pos, fg))]
# # print(str(bg_filtered))
# # print(length(unique(bg[bg \%in% pos])))
# # print(length(unique(bg[bg \%in% fg])))
# # print(length(unique(bg[bg \%in% c(pos, fg)])))
# # print(length(unique(bg_filtered)))
# bg <- bg_filtered

# q <- data.frame(set=character(), $mode=double())
)));

# Get p-values from output-pvalues-... via R
$r_p_value = getr('p_value');
$r_resampling_p_value = getr('resampling_p_value');

# Verify them (if they've been calculated above)
if (switch('2F'))
{
	# Warn only (because I'm arbitrarily using the p-value from the first score cutoff above)
	warn("Warning: Wilcoxon p-values don't match: original $p_value vs. R's $r_p_value read from the output-pvalues file (going with freshly calculated value) (accepting this for Figure 2F)") if (defined($p_value) and ($p_value != $r_p_value));
	warn("Warning: Resampling p-values don't match: original $resampling_p_value vs. R's $r_resampling_p_value read from the output-pvalues file (going with freshly calculated value) (accepting this for Figure 2F)") if (defined($resampling_p_value) and ($resampling_p_value != $r_resampling_p_value));
}
else
{
	# Error
	die("Error: Wilcoxon p-values don't match: original $p_value vs. R's $r_p_value read from the output-pvalues file (going with freshly calculated value)") if (defined($p_value) and ($p_value != $r_p_value));
	die("Error: Resampling p-values don't match: original $resampling_p_value vs. R's $r_resampling_p_value read from the output-pvalues file (going with freshly calculated value)") if (defined($resampling_p_value) and ($resampling_p_value != $r_resampling_p_value));
}

# Use R's p-values (freshly calculated)
$p_value = $r_p_value;
$resampling_p_value = $r_resampling_p_value;

# Correct the resampling p-value (<0.01 if 100 resamples)
$resampling_p_value_for_plot = $resampling_p_value;
if ($resampling_p_value == 0)
{
	# $resampling_p_value = 1 / $resamples;
	# $resampling_p_value = 0.0001;	# ***
	
	# Set this for the ggpubr::stat_compare_means plot only
	$resampling_p_value_for_plot = 0.0001;	# ***
}

# if (scalar(@posv) > 0)
# {
# 	runr(qq(
# 	q <- rbind(q, data.frame(set="Positive control", $mode=pos))
# 	));
# }
#
# if (scalar(@fgv) > 0)
# {
# 	runr(qq(
# 	q <- rbind(q, data.frame(set="Screen hits", $mode=fg))
# 	));
# }
#
# if (scalar(@bgv_all) > 0)
# {
# 	runr(qq(
# 	q <- rbind(q, data.frame(set="Random background", $mode=bg))
# 	));
# }

# runr(qq(
# 	# q <- rbind(q, data.frame(set="Positive control", $mode=pos))
# 	# q <- rbind(q, data.frame(set="Screen hits", $mode=fg))
# 	# q <- rbind(q, data.frame(set="Random background", $mode=bg))
#
# 	q <- rbind(q, data.frame(set="Positive control", $mode=pos))
# 	q <- rbind(q, data.frame(set="Identified interactions", $mode=fg))
# 	q <- rbind(q, data.frame(set="Random RBP pairs", $mode=bg))
# ));

# state(runr(qq(
# 	print(str(q))
# 	)));
state(runr(qq(
	# print(str(q))
	q\$set <- as.character(q\$set)
	
	# Only use one of the background resamples for plotting (otherwise the background boxplot notches become tiny because the sets are so unbalanced)
	q <- q[q\$set \%in% c("positives", "screen_hits", "background_resample_1"),]
	# q[q\$set=="background_resample_1",]\$set <- "Random RBP pairs"
	
	# q[q\$set!="positives" & q\$set!="screen_hits",]\$set <- "Random RBP pairs"
	# q[q\$set!="positives" & q\$set!="screen_hits",]\$set <- paste0("Random RBP pairs\\nn=",nrow(q[q\$set!="positives" & q\$set!="screen_hits",]))
	q[q\$set!="positives" & q\$set!="screen_hits",]\$set <- paste0("Random RBP pairs")
	# q[q\$set=="screen_hits",]\$set <- "Identified interactions"
	# q[q\$set=="screen_hits",]\$set <- paste0("Identified interactions\\nn=",nrow(q[q\$set=="screen_hits",]))
	q[q\$set=="screen_hits",]\$set <- paste0("Identified interactions")
	if (nrow(q[q\$set=="positives",]) > 0) {
		# q[q\$set=="positives",]\$set <- "Positive control"
		# q[q\$set=="positives",]\$set <- paste0("Positive control\\nn=",nrow(q[q\$set=="positives",]))
		q[q\$set=="positives",]\$set <- paste0("Positive control")
		q\$set <- factor(q\$set)
		# q\$set <- factor(q\$set, levels(q\$set)[c(2, 1, 3)])
		q\$set <- factor(q\$set, levels(q\$set)[c(3, 1, 2)])
	}
)));



# mindist: Plot with log10-scale y-axis
$tmpscale = '';
if (($mode eq 'mindist') or ($mode eq 'mediandist') or ($mode eq 'mindist_thresh') or ($mode eq 'oddsratio'))
{
	$tmpscale = ' + scale_y_log10()';
}


# Plot

if (switch('2A'))
{
	# Make Figure 2A plot
	state(runr(qq(
	# library(ggbeeswarm)
	# Plot
	# p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + labs(tag = "A") + geom_violin() + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + labs(tag = "A") + geom_boxplot(notch=T) + geom_beeswarm(alpha=0.3, colour="transparent") + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=$mode)) + labs(tag = "A") + geom_boxplot(aes(colour=set), notch=T) + geom_point(aes(colour=set), shape=16, alpha=0.3, position="jitter") + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + labs(tag = "A") + geom_boxplot(notch=T) + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + labs(tag = "A") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.3, size=0.5) + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + labs(tag = "A") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm() + xlab("") + ylab("Target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# ggsave("../Figure 2A - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	# p <- ggplot(q, aes(x=$mode, y=set, colour=set)) + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + xlab("") + ylab("") + ggtitle("RNA target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	# p <- ggplot(q, aes(x=$mode, y=set, colour=set)) + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + stat_compare_means(comparisons=list(c("Random RBP pairs", "Identified interactions"))) + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + xlab("") + ylab("") + ggtitle("RNA target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale

	# Hack ggpubr::stat_compare_means a little by supplying a fake kruskal.test function that returns the resampling p-value from above:
	library(ggpubr)
	# kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"="<0.01"))}
	kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"=$resampling_p_value_for_plot))}

	p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + coord_flip() + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=T) + stat_compare_means(method="kruskal.test", comparisons=list(c("Random RBP pairs", "Identified interactions")), label="p.signif") + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + xlab("") + ylab("p=$p_value (resampling: $resampling_p_value)") + ggtitle("RNA target set similarity (Jaccard index)") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	ggsave("../Figure 2A - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=45.75, units="mm")
	)));
}
elsif (switch('2B'))
{
	# Make Figure 2B plot
	state(runr(qq(
	# Plot
	# p <- ggplot(q, aes(x=set, y=cobind, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T) + xlab("") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=cobind, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.3, size=0.5) + xlab("") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=cobind, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm() + xlab("") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=set, y=cobind, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.2, shape=16) + xlab("") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# ggsave("../Figure 2B - $type - $scoretype.pdf", device=cairo_pdf, width=45.75, height=91.5, units="mm")
	# p <- ggplot(q, aes(x=cobind, y=set, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + coord_cartesian(xlim=c(0,1)) + xlab("Conditional probability of co-binding") + ylab("") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "black", "grey")) + scale_fill_manual(values = c("grey", "black", "grey"))$tmpscale
	# p <- ggplot(q, aes(x=cobind, y=set, colour=set)) + labs(tag = "B") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + coord_cartesian(xlim=c(0,1)) + xlab("Conditional probability of co-binding") + ylab("") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill"))$tmpscale
	# p <- ggplot(q, aes(x=cobind, y=set, colour=set)) + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + coord_cartesian(xlim=c(0,1)) + xlab("") + ylab("") + ggtitle("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	p <- ggplot(q, aes(x=cobind, y=set, colour=set)) + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=F) + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + coord_cartesian(xlim=c(0,1)) + xlab("") + ylab("") + ggtitle("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	ggsave("../Figure 2B - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=45.75, units="mm")
	)));
}
elsif (switch('2C_old'))
{
	# Make Figure 2C_old plot
	state(runr(qq(
	# Plot
	p <- ggplot(subset(q, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode, colour=as.factor(mindist_threshold))) + labs(tag = "C") + geom_boxplot(notch=T) + xlab("Maximum distance considered [nt]") + ylab("Fraction of shared targets bound in proximity") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "grey", "grey", "grey", "grey", "black", "grey", "grey"))$tmpscale
	ggsave("../Figure 2C_old - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
elsif (switch('2C'))
{
	# Make Figure 2C plot
	state(runr(qq(
	# Calculate delta: Subtract random pair median for each mindist_threshold
	q1 <- q
	for (i in $mindist_thresholds) {
		q1[q1\$set=="Identified interactions" & q1\$mindist_threshold==i,]\$$mode <- q1[q1\$set=="Identified interactions" & q1\$mindist_threshold==i,]\$$mode - median(q1[q1\$set=="Random RBP pairs" & q1\$mindist_threshold==i,]\$$mode)
	}
	# Plot
	# p <- ggplot(subset(q1, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode, colour=as.factor(mindist_threshold))) + labs(tag = "C") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.3) + xlab("Maximum distance considered [nt]") + ylab("Fraction of shared targets bound in proximity") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "grey", "grey", "black", "grey", "grey", "grey", "grey"))$tmpscale
	p <- ggplot(subset(q1, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode)) + labs(tag = "C") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.3, size=0.5) + xlab("Maximum distance considered [nt]") + ylab("Fraction of shared targets bound in proximity") + theme_minimal() + theme(legend.position = "none")$tmpscale
	ggsave("../Figure 2C - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
elsif (switch('2D_old'))
{
	# Make Figure 2D_old plot
	state(runr(qq(
	# Plot
	p <- ggplot(subset(q, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode, colour=as.factor(mindist_threshold))) + labs(tag = "D") + geom_boxplot(notch=T) + xlab("Maximum distance considered [nt]") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "grey", "grey", "grey", "grey", "black", "grey", "grey"))$tmpscale
	ggsave("../Figure 2D_old - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
elsif (switch('2D'))
{
	# Make Figure 2D plot
	state(runr(qq(
	# Calculate delta: Subtract random pair median for each mindist_threshold
	q1 <- q
	for (i in $mindist_thresholds) {
		q1[q1\$set=="Identified interactions" & q1\$mindist_threshold==i,]\$$mode <- q1[q1\$set=="Identified interactions" & q1\$mindist_threshold==i,]\$$mode - median(q1[q1\$set=="Random RBP pairs" & q1\$mindist_threshold==i,]\$$mode)
	}
	# Plot
	# p <- ggplot(subset(q1, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode, colour=as.factor(mindist_threshold))) + labs(tag = "D") + geom_boxplot(notch=T) + xlab("Maximum distance considered [nt]") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "grey", "grey", "black", "grey", "grey", "grey", "grey"))$tmpscale
	p <- ggplot(subset(q1, set=="Identified interactions"), aes(x=as.factor(mindist_threshold), y=$mode)) + labs(tag = "D") + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.3, size=0.5) + xlab("Maximum distance considered [nt]") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none")$tmpscale
	ggsave("../Figure 2D - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
elsif (switch('2F'))
{
	# Make Figure 2F plot
	state(runr(qq(
	# Plot
	p <- ggplot(subset(q, set=="Identified interactions"), aes(x=naturalfactor(paste0('≥', minscore)), y=$mode, colour=as.factor(minscore))) + labs(tag = "F") + geom_boxplot(notch=T) + xlab("Average protein interaction score") + ylab("Fraction of co-bound RNAs bound within $mindist_threshold nt") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey", "black", "grey"))$tmpscale
	ggsave("../Figure 2F - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
elsif (switch('2G'))
{
	# Make Figure 2G plot
	state(runr(qq(
	# Plot
	p <- ggplot(subset(q, set=="Identified interactions"), aes(x=naturalfactor(paste0('≥', minscore)), y=$mode, colour=as.factor(minscore))) + labs(tag = "G") + geom_boxplot(notch=T) + xlab("Average protein interaction score") + ylab("Conditional probability of co-binding") + theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values = c("black", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))$tmpscale
	ggsave("../Figure 2G - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
	)));
}
# else
# {
# 	# Make normal plot
# 	$p_value = 'NA' if (!defined($p_value));
# 	$resampling_p_value = 'NA' if (!defined($resampling_p_value));
# 	runr(qq(
# 	p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + geom_boxplot(notch=T) + xlab("") + labs(title="Cobinding \\"$mode\\" analysis", subtitle="CLIP data: $type, minimum confidence $hc, minlog2fold $minlog2fold\\nmin_$scoretype $minscore, mindist_threshold $mindist_threshold, resamples $resamples", caption="Wilcoxon p-value: $p_value\\nResampling p-value: $resampling_p_value") + theme_bw() + theme(legend.position = "none")$tmpscale
# 	ggsave("$pdffile", width=5.5, height=4.5)
# 	));
# }
else
{
	# Make normal plot
	$p_value = 'NA' if (!defined($p_value));
	$resampling_p_value = 'NA' if (!defined($resampling_p_value));
	$resampling_p_value_for_plot = 'NA' if (!defined($resampling_p_value_for_plot));
	
	runr(qq(

	# Hacking ggpubr::stat_compare_means a little, by supplying a fake "kruskal.test" function that returns the resampling p-value from above:
	library(ggpubr)
	# kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"="<0.01"))}
	kruskal.test <- function(x,g,formula,data,subset,na.action,paired){return(list("p.value"=$resampling_p_value_for_plot))}

	# Plot horizontal
	p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + coord_flip() + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=T) + stat_compare_means(method="kruskal.test", comparisons=list(c("Random RBP pairs", "Identified interactions")), label="p.signif") + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + xlab("") + ylab("p=$p_value (resampling: $resampling_p_value)") + ggtitle("$mode: $axis") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	ggsave("../Figure 2 - $table - $mode - $type - $scoretype.pdf", device=cairo_pdf, width=91.5, height=45.75, units="mm")
	
	# Plot vertical (old school)
	p <- ggplot(q, aes(x=set, y=$mode, colour=set)) + geom_boxplot(notch=T, outlier.shape=NA) + geom_beeswarm(alpha=0.25, shape=16, groupOnX=T) + stat_compare_means(method="kruskal.test", comparisons=list(c("Random RBP pairs", "Identified interactions")), label="p.signif") + scale_colour_manual(values = c("#BEC1C0", "#FF7F00", "#1E3D59"), aesthetics=c("colour", "fill")) + xlab("") + ylab("p=$p_value (resampling: $resampling_p_value)") + ggtitle("$axis") + theme_minimal() + theme(legend.position = "none") + theme(plot.title=element_text(hjust=0))$tmpscale
	ggsave("$pdffile", width=5.5, height=4.5)
	));
}





# Functions

# Resample background (returns a list of RBP pairs) (sampling with replacement)
sub resample
{
	my ($n) = @_;
	# our @rbps_both;
	
	our $mode_directional;
	
	my @bg = ();
	while (scalar(@bg) < $n)
	{
		# Choose random RBPs
		$rbp1 = $rbps_both[rand(@rbps_both)];
		$rbp2 = $rbps_both[rand(@rbps_both)];
		
		# Skip homomers
		next if ($rbp1 eq $rbp2);
		
		# Skip if rbp1 isn't alphanumerically before rbp2 (to avoid duplicates - not really necessary, but that's how the screen was done)
		if ($mode_directional == 1)
		{
			next if ($rbp1 gt $rbp2);
		}
		
		# Insert into array
		push(@bg, "$rbp1|$rbp2");
		
		# Actually that's fine (sampling with replacement)
		# # Unique it so the same pairs don't show up twice
		# @bg = unique(@bg);
	}
	
	return(@bg);
}


# Resample background (returns a list of RBP pairs) (sampling with replacement)
sub resample_without_replacement
{
	my ($n) = @_;
	# our @rbps_both;
	
	our $mode_directional;
	our @fg;
	
	my @bg = ();
	while (scalar(@bg) < $n)
	{
		# Choose random RBPs
		$rbp1 = $rbps_both[rand(@rbps_both)];
		$rbp2 = $rbps_both[rand(@rbps_both)];
		
		# Skip homomers
		next if ($rbp1 eq $rbp2);
		
		# Skip if rbp1 isn't alphanumerically before rbp2 (to avoid duplicates - not really necessary, but that's how the screen was done)
		if ($mode_directional == 1)
		{
			next if ($rbp1 gt $rbp2);
		}
		
		# Skip foreground (screen hit) pairs (this function is only used by the -randompairs switch, not for resampling p-values)
		next if (contains("$rbp1|$rbp2", @fg));
		
		# Insert into array
		push(@bg, "$rbp1|$rbp2");
		
		# Sampling without replacement:
		# Unique it so the same pairs don't show up twice
		@bg = unique(@bg);
	}
	
	return(@bg);
}


# Define analysis function: cobind_jaccard
# Function to calculate gene co-binding Jaccard index for a pair of RBPs
sub cobind_jaccard
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));

		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		my @union = union(\@a, \@b);
		
		if (scalar(@union) > 0)
		{
			$jaccard = scalar(@ensgvs) / scalar(@union);
		}
		else
		{
			$jaccard = 0;
		}

		return($jaccard);
	}
	else
	{
		return(0);
	}
}

# Define analysis function: cobind_probability
# Function to calculate gene co-binding probability for a pair of RBPs (A binding, given B binding) (rbp1 binding, given rbp2 binding)
sub cobind_probability
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> probability of $rbp1 binding, given $rbp2 binding\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# Calculate conditional probability if rbp2 binds anything
		if (scalar(@b) > 0)
		{
			$prob = scalar(@ensgvs) / scalar(@b);
			print "       >> p($rbp1|$rbp2) = ".scalar(@ensgvs)." / ".scalar(@b)." = $prob\n" if (switch('debug'));
		}
		else
		{
			$prob = undef;
		}

		return($prob);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_oddsratio
# Function to calculate gene binding odds ratio for a pair of RBPs
sub cobind_oddsratio
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> odds ratio of $rbp1 and $rbp2 co-binding" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		# Get intersection
		my @a_b = intersection(\@a, \@b);
		
		# Not bound
		my @nota = lonly(\@total_ensgvs, \@a);
		my @notb = lonly(\@total_ensgvs, \@b);
		
		# Combinations
		my @nota_b = intersection(\@nota, \@b);
		my @a_notb = intersection(\@a, \@notb);
		my @nota_notb = intersection(\@nota, \@notb);
		
		
		print "\n   >> TOTAL >> ".scalar(@total_ensgvs)."\n" if (switch('debug'));
		print "   >> a >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "     >> NOT a >> $rbp1 >> ".scalar(@nota)."\n" if (switch('debug'));
		print "   >> b >> $rbp1 >> ".scalar(@b)."\n" if (switch('debug'));
		print "     >> NOT b >> $rbp1 >> ".scalar(@notb)."\n" if (switch('debug'));
		print "   >> a_b >> ".scalar(@a_b)."\n" if (switch('debug'));
		print "   >> nota_b >> ".scalar(@nota_b)."\n" if (switch('debug'));
		print "   >> a_notb >> ".scalar(@a_notb)."\n" if (switch('debug'));
		print "   >> nota_notb >> ".scalar(@nota_notb)."\n" if (switch('debug'));
		
		# Calculate odds ratio between proteins
		$oddsratio = (scalar(@a_b) / scalar(@a_notb)) / (scalar(@nota_b) / scalar(@nota_notb));
		print "       >> odds ratio ($rbp1, $rbp2) = (".scalar(@a_b)." / ".scalar(@a_notb).") / (".scalar(@nota_b)." / ".scalar(@nota_notb).") = $oddsratio\n" if (switch('debug'));
		# equivalent
		# $oddsratio = (scalar(@a_b) / scalar(@nota_b)) / (scalar(@a_notb) / scalar(@nota_notb));
		# print "       >> odds ratio ($rbp1, $rbp2) = (".scalar(@a_b)." / ".scalar(@nota_b).") / (".scalar(@a_notb)." / ".scalar(@nota_notb).") = $oddsratio\n" if (switch('debug'));
		# exit if (switch('debug'));

		return($oddsratio);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_freq
# Function to calculate gene binding frequency for either of a pair of RBPs (A binding or B binding, divided by the total number of genes)
sub cobind_freq
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> frequency of $rbp1 binding or $rbp2 binding, divided by the total number of genes\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get union
		my @ensgvs = union(\@a, \@b);
		print "     >> union >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# Calculate frequency of either protein binding
		$prob = scalar(@ensgvs) / $total_ensgvs;
		print "       >> p($rbp1 or $rbp2) = ".scalar(@ensgvs)." / ".$total_ensgvs." = $prob\n" if (switch('debug'));

		return($prob);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_probability_freqscaled
# Function to calculate gene co-binding probability for a pair of RBPs (A binding, given B binding) (rbp1 binding, given rbp2 binding), multiplied by the probability of neither binding in order to highlight rare cases (* (1 - p(A u B)))
sub cobind_probability_freqscaled
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> probability of $rbp1 binding, given $rbp2 binding\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# Calculate conditional probability if rbp2 binds anything
		if (scalar(@b) > 0)
		{
			$prob = (scalar(@ensgvs) / scalar(@b)) * (1 - (scalar(union(\@a, \@b)) / $total_ensgvs));
			print "       >> p($rbp1|$rbp2) = ".scalar(@ensgvs)." / ".scalar(@b)." = $prob\n" if (switch('debug'));
		}
		else
		{
			$prob = undef;
		}

		return($prob);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_probability_mindist
# Function to calculate gene co-binding probability for a pair of RBPs (A binding, given B binding) (rbp1 binding, given rbp2 binding)
# ...with a minimum distance filter
sub cobind_probability_mindist
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> probability of $rbp1 binding within $mindist_threshold nt, given $rbp2 binding\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
		
		# Apply mindist filter
		startme("Getting mindist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
		my @mindists = ();
		foreach $ensgv (@ensgvs)
		{
			my $mindist = mindist($rbp1, $rbp2, $ensgv);
			if (defined($mindist) and ($mindist <= $mindist_threshold))
			{
				push(@mindists, $mindist);
				stepme(100) if (switch('debug'));
			}
		}
		stopme() if (switch('debug'));
	
		characterise(@mindists) if (switch('debug') or (switch('debug2')));

		# Calculate conditional probability if rbp2 binds anything
		if (scalar(@b) > 0)
		{
			$prob = scalar(@mindists) / scalar(@b);
			print "       >> p($rbp1|$rbp2) = ".scalar(@mindists)." (out of ".scalar(@ensgvs).") / ".scalar(@b)." = $prob\n" if (switch('debug'));
		}
		else
		{
			$prob = undef;
		}

		return($prob);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_probability_mindist_freqscaled
# Function to calculate gene co-binding probability for a pair of RBPs (A binding, given B binding) (rbp1 binding, given rbp2 binding), multiplied by the probability of neither binding in order to highlight rare cases (* (1 - p(A u B)))
# ...with a minimum distance filter
sub cobind_probability_mindist_freqscaled
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> probability of $rbp1 binding within $mindist_threshold nt, given $rbp2 binding\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
		
		# Apply mindist filter
		startme("Getting mindist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
		my @mindists = ();
		foreach $ensgv (@ensgvs)
		{
			my $mindist = mindist($rbp1, $rbp2, $ensgv);
			if (defined($mindist) and ($mindist <= $mindist_threshold))
			{
				push(@mindists, $mindist);
				stepme(100) if (switch('debug'));
			}
		}
		stopme() if (switch('debug'));
	
		characterise(@mindists) if (switch('debug') or (switch('debug2')));

		# Calculate conditional probability if rbp2 binds anything
		if (scalar(@b) > 0)
		{
			$prob = (scalar(@mindists) / scalar(@b)) * (1 - (scalar(union(\@a, \@b)) / $total_ensgvs));
			print "       >> p($rbp1|$rbp2) = ".scalar(@mindists)." (out of ".scalar(@ensgvs).") / ".scalar(@b)." = $prob\n" if (switch('debug'));
		}
		else
		{
			$prob = undef;
		}

		return($prob);
	}
	else
	{
		return(undef);
	}
}

# Define analysis function: cobind_fraction (the fraction that's in the intersection (for the smaller set of the two, the 'better set'))
sub cobind_fraction
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
		
		# If there is some overlap in genes between these two RBPs...
		my @union = union(\@a, \@b);
		
		if (scalar(@union) > 0)
		{
			$frac_a = scalar(@ensgvs) / scalar(@a);
			print "       >> fraction 1 >> $frac_a\n" if (switch('debug'));
			$frac_b = scalar(@ensgvs) / scalar(@b);
			print "       >> fraction 2 >> $frac_b\n" if (switch('debug'));
			
			# Return the higher fraction
			return(max($frac_a, $frac_b));
		}
		else
		{
			return(0)
		}
	}
	else
	{
		return(0);
	}
}

# Define analysis function: cobind_intersection (the absolute number of RNA genes bound by both RBPs)
sub cobind_intersection
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});
		
		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> intersection >> ".scalar(@ensgvs)."\n" if (switch('debug'));
		
		return(scalar(@ensgvs));
	}
	else
	{
		return(0);
	}
}

# Define analysis function: cobind_mindist_thresh_frac
# Function to get the number of cases below a minimum binding distance for a pair of RBPs (as a fraction of their shared target genes)
# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
sub cobind_mindist_thresh_frac
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	our $mindist_threshold;
	
	print " >> $rbp1 vs. $rbp2 >> threshold $mindist_threshold\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));

		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> overlap >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		if (scalar(@ensgvs) > 0)
		{
			startme("Getting mindist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
			my @mindists = ();
			foreach $ensgv (@ensgvs)
			{
				my $mindist = mindist($rbp1, $rbp2, $ensgv);
				if (defined($mindist) and ($mindist <= $mindist_threshold))
				{
					push(@mindists, $mindist);
					stepme(100) if (switch('debug'));
				}
			}
			stopme() if (switch('debug'));
		
			characterise(@mindists) if (switch('debug') or (switch('debug2')));
			
			state(scalar(@mindists) / scalar(@ensgvs)) if (switch('debug'));
			
			return(scalar(@mindists) / scalar(@ensgvs));
		}
		else
		{
			return(0);
		}
	}
	else
	{
		return(0);
	}
}

# Define analysis function: cobind_mindist_thresh
# Function to get the absolute number of genes in which RBP1 and RBP2 both bind within a minimum binding distance, i.e. probably as a complex
# The distance measure used is the distance between the peaks.
# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
sub cobind_mindist_thresh
{
	my ($pair, $set) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	our $mindist_threshold;
	
	print " >> $rbp1 vs. $rbp2 >> threshold $mindist_threshold\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));

		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> overlap >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		if (scalar(@ensgvs) > 0)
		{
			startme("Getting mindist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
			my @mindists = ();
			foreach $ensgv (@ensgvs)
			{
				my $mindist = mindist($rbp1, $rbp2, $ensgv);
				if (defined($mindist) and ($mindist <= $mindist_threshold))
				{
					push(@mindists, $mindist);
					
					# Write to hits file
					print HITS "$table\t$type\t$scoretype\t$mode\t$hc\t$minlog2fold\t$minscore\t$mindist_threshold\t$resamples\t$set\t$rbp1\t$rbp2\t$ensgv\t$mindist\n";

					stepme(100) if (switch('debug'));
				}
			}
			stopme() if (switch('debug'));
		
			characterise(@mindists) if (switch('debug2'));
	
			return(scalar(@mindists));
		}
		else
		{
			return(0);
		}
	}
	else
	{
		return(0);
	}
}

# Define analysis function: cobind_mindist
# Function to get the minimum binding site distance for a pair of RBPs (averaged across their shared target genes)
# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
sub cobind_mindist
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));

		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> overlap >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		if (scalar(@ensgvs) > 0)
		{
			startme("Getting mindist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
			my @mindists = ();
			foreach $ensgv (@ensgvs)
			{
				my $mindist = mindist($rbp1, $rbp2, $ensgv);
				if (defined($mindist) and ($mindist <= $mindist_threshold))
				{
					push(@mindists, $mindist);
				}
		
				stepme(100) if (switch('debug'));
			}
			stopme() if (switch('debug'));
		
			characterise(@mindists) if (switch('debug2'));
	
			return(avg(@mindists));
		}
		else
		{
			return(undef);
		}
	}
	else
	{
		return(undef);
	}
}

# Using clip_cobinding table, but it's just too problematic (since it doesn't have ENSGV information and it's specific for each cell type and replicate, as opposed to the original mindist function here)
# # Define utility function: mindist
# sub mindist
# {
# 	my ($rbp1, $rbp2, $ensgv) = @_;
# 	# our %peaks;
#
# 	# my @peaks1 = ();
# 	# if (exists($peaks{"$rbp1|$ensgv"}))
# 	# {
# 	# 	@peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv"});
# 	# }
# 	# my @peaks2 = ();
# 	# if (exists($peaks{"$rbp2|$ensgv"}))
# 	# {
# 	# 	@peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv"});
# 	# }
#
# 	my $mindist = undef;
#
# 	# Get peaks within this ENSGV for RBP1 (to look up the closest distance binding site of RBP2 in clip_cobinding)
# 	my $peakquery = Query("SELECT DISTINCT chr, start, stop, strand FROM $clip_raw_gene WHERE species='human' AND type='$type' AND symbol='$rbp1' AND ensgv='$ensgv'");
# 	while (($chr, $start, $stop, $strand) = Fetch($peakquery))
# 	{
# 		# Look up minimum distance in clip_cobinding
# 		my $query = Query("SELECT MIN(ABS(mindist)) FROM $clip_cobinding WHERE type='$type' AND species='human' AND symbol='$rbp1' AND chr='$chr' AND start='$start' AND stop='$stop' AND strand='$strand' AND other_symbol='$rbp2'");
# 		if (Numrows($query) > 0)
# 		{
# 			my ($tmp_mindist) = FetchOne($query);
#
# 			if (!defined($mindist) or (defined($tmp_mindist) and ($tmp_mindist < $mindist)))
# 			{
# 				$mindist = $tmp_mindist;
# 				print "     >> $ensgv >> $mindist\n" if (switch('debug'));
# 			}
# 		}
# 	}
#
#
# 	# foreach my $peak1 (@peaks1)
# 	# {
# 	# 	foreach my $peak2 (@peaks2)
# 	# 	{
# 	# 		if (!defined($mindist) or ($mindist > abs($peak2 - $peak1)))
# 	# 		{
# 	# 			$mindist = abs($peak2 - $peak1);
# 	# 			print "     >> $ensgv >> $peak1 / $peak2 >> $mindist\n" if (switch('debug'));
# 	# 		}
# 	# 	}
# 	# }
#
# 	return($mindist);
# }

# Slower method without using the unfinished pre-calculated clip_cobinding table
# Define utility function: mindist
sub mindist
{
	my ($rbp1, $rbp2, $ensgv) = @_;
	our %peaks;
	
	my $tmpi = '';
	$tmpi = "|1" if (switch('random'));

	my @peaks1 = ();
	if (exists($peaks{"$rbp1|$ensgv$tmpi"}))
	{
		@peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv$tmpi"});
	}
	my @peaks2 = ();
	if (exists($peaks{"$rbp2|$ensgv$tmpi"}))
	{
		@peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv$tmpi"});
	}

	my $mindist = undef;

	foreach my $peak1 (@peaks1)
	{
		foreach my $peak2 (@peaks2)
		{
			if (!switch('distance_signed'))
			{
				# Absolute
				if (!defined($mindist) or ($mindist > abs($peak2 - $peak1)))
				{
					$mindist = abs($peak2 - $peak1);
				}
			}
			else
			{
				# Signed
				if (!defined($mindist) or (abs($mindist) > abs($peak2 - $peak1)))
				{
					$mindist = $peak2 - $peak1;
				}
			}

			print "     >> $ensgv >> $peak1 / $peak2 >> $mindist\n" if (switch('debug'));
		}
	}

	return($mindist);
}

# Define analysis function: cobind_mediandist
# Function to get the minimum binding site distance for a pair of RBPs (averaged across their shared target genes)
# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
sub cobind_mediandist
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));

		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> overlap >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		if (scalar(@ensgvs) > 0)
		{
			startme("Getting mediandist values for $rbp1|$rbp2", 0, scalar(@ensgvs)) if (switch('debug'));
			my @mediandists = ();
			foreach $ensgv (@ensgvs)
			{
				push(@mediandists, mediandist($rbp1, $rbp2, $ensgv));
		
				stepme(100) if (switch('debug'));
			}
			stopme() if (switch('debug'));
		
			characterise(@mediandists) if (switch('debug2'));
	
			return(median(@mediandists));
		}
		else
		{
			return(undef);
		}
	}
	else
	{
		return(undef);
	}
}

# Define utility function: mediandist
sub mediandist
{
	my ($rbp1, $rbp2, $ensgv) = @_;
	our %peaks;
	
	my $tmpi = '';
	$tmpi = "|1" if (switch('random'));

	my @peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv$tmpi"});
	my @peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv$tmpi"});
	
	my @dists = ();
	
	foreach my $peak1 (@peaks1)
	{
		foreach my $peak2 (@peaks2)
		{
			$dist = abs($peak2 - $peak1);
			push(@dists, $dist);
			print "     >> $ensgv >> $peak1 / $peak2 >> $dist\n" if (switch('debug'));
		}
	}
	
	return(median(@dists));
}



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





# Define analysis function: distance_plot
# Function to plot the distribution of binding site distances for a pair of RBPs in all their shared target genes
# Using the 5' end of peak, rather than its full extent (since according to PMID 29883606, this is a good proxy of the binding site)
# # Could alternatively use the center (easiest, and in the case of a narrow peak and a broad peak with almost the same center, it'd be more accurate.)
sub distance_plot
{
	my ($pair) = @_;
	my ($rbp1, $rbp2) = split(/\|/, $pair);
	our %ensgvs;
	
	print " >> $rbp1 vs. $rbp2\n" if (switch('debug'));
	
	# Get shared genes (which both RBPs bind)
	if (exists($ensgvs{$rbp1}) and exists($ensgvs{$rbp2}))
	{
		# Create temporary arrays
		my @a = split(/\|/, $ensgvs{$rbp1});
		my @b = split(/\|/, $ensgvs{$rbp2});

		print "   >> 1 >> $rbp1 >> ".scalar(@a)."\n" if (switch('debug'));
		print "   >> 2 >> $rbp2 >> ".scalar(@b)."\n" if (switch('debug'));
		
		# Get intersection
		my @ensgvs = intersection(\@a, \@b);
		print "     >> overlap >> ".scalar(@ensgvs)."\n" if (switch('debug'));
	
		# If there is some overlap in genes between these two RBPs...
		if (scalar(@ensgvs) > 0)
		{
			# startme("Getting mindist values (one_per_ensgv) for $rbp1|$rbp2", 0, scalar(@ensgvs));
			# my @mindists = ();
			# foreach $ensgv (@ensgvs)
			# {
			# 	# push(@mindists, mindists_for_plot($rbp1, $rbp2, $ensgv));
			# 	push(@mindists, mindist($rbp1, $rbp2, $ensgv));
			#
			# 	# #DEBUG
			# 	# if (mindist($rbp1, $rbp2, $ensgv) != mindist_old($rbp1, $rbp2, $ensgv))
			# 	# {
			# 	# 	die("Error: mindist() and mindist_old() disagree for mindist($rbp1, $rbp2, $ensgv):\n\nmindist says: ".mindist($rbp1, $rbp2, $ensgv)."\nmindist_old says: ".mindist_old($rbp1, $rbp2, $ensgv)."\n\n");
			# 	# }
			# 	# #END DEBUG
			#
			# 	stepme(100);
			# }
			# stopme();
			#
			# state("Mindist values (one_per_ensgv):");
			# characterise(@mindists);
			#
			# startr();
			# setr("mindists", @mindists);
			# state(runr(qq(
			#
			# print(str(mindists))
			#
			# library(ggplot2)
			# library(scales)
			#
			# # qplot(mindists + 1) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw()
			# # qplot(mindists + 1) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = comma) + annotation_logticks(sides="b") + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("Closest peak distance per gene")
			# qplot(mindists, binwidth=10) + coord_cartesian(xlim=c(0, 1000)) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("Closest peak distance per gene")
			# ggsave("../output-distance-plot-one_per_ensgv-$scoretype-$rbp1-$rbp2$tmprandom.pdf", width=7, height=5)
			#
			# )));

			
			# Closest per peak: This is the best, most interesting analysis (how many of A's binding events are in proximity of B's?)
			startme("Getting mindist values (closest_per_peak) for $rbp1|$rbp2", 0, scalar(@ensgvs)) if switch('debug');
			@mindists = ();
			$tmpdistfile = "../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2-fit-mindists-detailed$tmpsigned$tmprandom.txt";
			open(TMPDISTS, ">$tmpdistfile") or die("Error: Couldn't write to '$tmpdistfile'");
			print TMPDISTS "rbp1\trbp2\tensgv\tdist\n";
			foreach $ensgv (@ensgvs)
			{
				@more_mindists = mindists_closest_per_peak($rbp1, $rbp2, $ensgv);
				
				# Store
				push(@mindists, @more_mindists);
				
				# Write to output text table
				foreach $tmpdist (@more_mindists)
				{
					print TMPDISTS "$rbp1\t$rbp2\t$ensgv\t$tmpdist\n";
				}
		
				stepme(100) if switch('debug');
			}
			close(TMPDISTS);
			stopme() if switch('debug');
		
			state("Mindist values (closest_per_peak):") if switch('debug');
			characterise(@mindists) if switch('debug');
			
			startr();
			setr("mindists", @mindists);
			state(runr(qq(
			
			# print(head(mindists, n=1000))

			library(ggplot2)
			library(scales)
			
			# qplot(mindists + 1) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw()
			# qplot(mindists + 1) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = comma) + annotation_logticks(sides="b") + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("Closest peak distance per peak")
			qplot(mindists, binwidth=10) + coord_cartesian(xlim=c($tmpsignedmin, 1000)) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("Closest peak distance per peak")
			ggsave("../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2$tmpsigned$tmprandom.pdf", width=7, height=5)



			saveRDS(mindists, file = "../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2-fit-mindists$tmpsigned$tmprandom.rds")
			write.table(mindists, file = "../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2-fit-mindists$tmpsigned$tmprandom.txt", quote=F, sep="\t", row.names=F, col.names=F)
			
			)));
			
			
			if (!switch('nofit'))
			{
				state(runr(qq(
				
				# Mixture model approach (bimodal Gaussian)

				library(mixtools)
				# library(dplyr)
				# library(scales)

				mindists_shifted <- mindists + 1
				q <- data.frame(mindists=mindists_shifted)
				q\$logmindists <- log10(q\$mindists)
				m <- normalmixEM(log10(mindists_shifted), k=2)

				# Find out which distribution is the smaller one
				smaller <- 1
				bigger <- 2
				if (m\$mu[2] < m\$mu[1]) {
				  smaller <- 2
				  bigger <- 1
				}

				# Define function: dnorm values, scaled by lambda
				dnorm_scaled <- function(x, lambda, mean, sd) {
				  return(dnorm(x, mean=mean, sd=sd) * lambda)
				}

				# Splitting:
				# http://tinyheero.github.io/2015/10/13/mixture-model.html
				mcomps <- as.data.frame(cbind(x = m\$x, m\$posterior))
				# head(mcomps)
				smallercomp <- paste0("comp.", smaller)
				biggercomp <- paste0("comp.", bigger)
				mcomps_smaller <- mcomps[mcomps[[smallercomp]] > mcomps[[biggercomp]],]
				# summary(mcomps_smaller)
				close_max_log10 <- max(mcomps_smaller\$x)
				close_max <- 10^close_max_log10
				# close_max
				# Calculate mean and SD on the non-logged distances
				close_mean <- mean(10^mcomps_smaller\$x)
				close_sd <- sd(10^mcomps_smaller\$x)
				# This fraction differs very slightly from m\$lambda[smaller]:
				close_frac <- nrow(mcomps_smaller) / nrow(mcomps)
				# m\$lambda[smaller]
				# close_frac
				# >> Using close_frac since lambda is a scaling value

				# ggplot(q, aes(x=logmindists)) + geom_histogram(binwidth=0.1)
				# ggplot(q) + geom_histogram(aes(x=log10(mindists), y=..density..), binwidth=0.1)
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + geom_density(linetype="dashed")
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1)
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + stat_function(fun=dnorm, args=list(mean=m\$mu[1], sd=m\$sigma[1])) + stat_function(fun=dnorm, args=list(mean=m\$mu[2], sd=m\$sigma[2]))
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[1], mean=m\$mu[1], sd=m\$sigma[1]), colour="red") + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[2], mean=m\$mu[2], sd=m\$sigma[2]), colour="black")
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[1], mean=m\$mu[1], sd=m\$sigma[1]), colour="red", linetype="dotted", size=0.75) + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[2], mean=m\$mu[2], sd=m\$sigma[2]), colour="grey", linetype="dotted", size=0.75)
				# ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[1], mean=m\$mu[1], sd=m\$sigma[1]), colour="red", linetype="dotted") + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[2], mean=m\$mu[2], sd=m\$sigma[2]), colour="grey", linetype="dotted")
				p <- ggplot(q, aes(x=logmindists)) + geom_histogram(aes(y=..density..), binwidth=0.1) + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[smaller], mean=m\$mu[smaller], sd=m\$sigma[smaller]), colour="red") + stat_function(fun=dnorm_scaled, args=list(lambda=m\$lambda[bigger], mean=m\$mu[bigger], sd=m\$sigma[bigger]), colour="grey", linetype="dotted") + geom_vline(xintercept = close_max_log10, linetype="dotted", colour="red") + ggtitle(paste0(sprintf("%.1f", close_frac * 100), "% of $rbp1 peaks\nhave $rbp2 peaks nearby\n(<",round(close_max, 0)," nt, mean ", round(close_mean, 0), "±",round(close_sd, 0)," nt)")) + scale_x_continuous(breaks=seq(0, max(q\$logmindists), 1), labels=format(10^seq(0, max(q\$logmindists), 1), big.mark=",", trim=T, scientific=F)) + theme_bw() + xlab("Distance from $rbp1 peak to nearest $rbp2 peak") + ylab("Density")
				ggsave("../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2-fit$tmpsigned$tmprandom.pdf", device=cairo_pdf, width=91.5, height=91.5, units="mm")
			
				# Write the fit parameters to a table
				close_df <- data.frame("close_max"=close_max, "close_mean"=close_mean, "close_sd"=close_sd, "close_frac"=close_frac)
				write.table(close_df, "../output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$rbp1-$rbp2-fit$tmpsigned$tmprandom-parameters.txt", quote=F, sep="\t", row.names=F, col.names=T)

				)));
			}





			# startme("Getting mindist values (all_pairs) for $rbp1|$rbp2", 0, scalar(@ensgvs));
			# @mindists = ();
			# foreach $ensgv (@ensgvs)
			# {
			# 	push(@mindists, mindists_all_pairs($rbp1, $rbp2, $ensgv));
			#
			# 	stepme(100);
			# }
			# stopme();
			#
			# state("Mindist values (all_pairs):");
			# characterise(@mindists);
			#
			# startr();
			# setr("mindists", @mindists);
			# state(runr(qq(
			#
			# print(str(mindists))
			#
			# library(ggplot2)
			# library(scales)
			#
			# # qplot(mindists + 1) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw()
			# # qplot(mindists + 1) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = comma) + annotation_logticks(sides="b") + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("All-against-all peak distances")
			# qplot(mindists, binwidth=10) + coord_cartesian(xlim=c(0, 1000)) + xlab("Peak distance") + ylab("Peak comparisons") + theme_bw() + ggtitle("All-against-all peak distances")
			# ggsave("../output-distance-plot-all_pairs-$rbp1-$rbp2.pdf", width=7, height=5)
			#
			# )));
			
			
			
			
			
			
			# Motif fraction when cobound: (how many of A's closest binding events to B's contain motifs?)
			startme("Getting mindist_motifs values (closest per peak) for $rbp1|$rbp2", 0, scalar(@ensgvs)) if switch('debug');
			
			
			# Still got mindist values from above, can also use them here
			# @mindists = ();
			# foreach $ensgv (@ensgvs)
			# {
			# 	push(@mindists, mindists_closest_per_peak($rbp1, $rbp2, $ensgv));
			#
			# 	stepme(100) if switch('debug');
			# }
			# stopme() if switch('debug');
			#
			# state("Mindist values (closest_per_peak):") if switch('debug');
			# characterise(@mindists) if switch('debug');
			#
			# startr();
			# setr("mindists", @mindists);
			# state(runr(qq(
			
			
			
			
			if (switch('motif_plot'))
			{
				# @mindists = ();
				my $tmpmotiffile = "../output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$rbp1-$rbp2-data.txt";
				open(TMPMOTIFS, ">$tmpmotiffile") or die("Error: Couldn't write to '$tmpmotiffile");
				print TMPMOTIFS "rbp1\trbp2\tensgv\tpeak1\tpeak2\tmindist\tbindtype\tmotif1\tmotif2\n";
				# @motif_fractions = ();
				foreach $ensgv (@ensgvs)
				{
					# push(@mindists, mindists_closest_per_peak($rbp1, $rbp2, $ensgv));
					# push(@motif_fractions, mindist_motifs($rbp1, $rbp2, $ensgv));
					mindist_motifs($rbp1, $rbp2, $ensgv);

					stepme(100) if switch('debug');
				}
				stopme() if switch('debug');
				close(TMPMOTIFS);

				# # Get background motif fraction for rbp2 (how many of its own sites have motifs, regardless of cobinding?)
				# print TMPMOTIFS "$rbp1\t$rbp2\t$ensgv\t\t\t\tbackground\t0\t0\n";
		

				state("Wrote values for mindist_motifs analysis to: $tmpmotiffile") if switch('debug');
				# characterise(@mindists) if switch('debug');
				# characterise(@motif_fractions) if switch('debug');

				
				# warn("ADD PLOTTING HERE");
				# startr();
				# setr("motif_fractions", @motif_fractions);
				state(runr(qq(
				
				# library(dplyr)
				library(plyr)
				
				q <- read.table("$tmpmotiffile", header=T, sep="\t", quote="")
				rbp1 <- "$rbp1"; rbp2 <- "$rbp2"

				)));
				state(runr(q(
				qo <- q
				# str(q)

				# Improve labelling
				q$motif1 <- as.factor(q$motif1)
				q$motif1 <- mapvalues(q$motif1, c(1, 0), c(paste0(rbp1," binding site\nwith motif"), paste0(rbp1," binding site\nwithout motif")))
				q$motif1 <- factor(q$motif1, levels = c(paste0(rbp1," binding site\nwith motif"), paste0(rbp1," binding site\nwithout motif")))

				q$motif2 <- as.factor(q$motif2)
				q$motif2 <- mapvalues(q$motif2, c(0, 1), c(paste0(rbp2," binding site\nwithout motif"), paste0(rbp2," binding site\nwith motif")))
				q$motif2 <- factor(q$motif2, levels = c(paste0(rbp2," binding site\nwithout motif"), paste0(rbp2," binding site\nwith motif")))

				q$bindtype <- mapvalues(q$bindtype, c("close", "distant"), c("any_close", "any_distant"))
				# q$bindtype <- mapvalues(q$bindtype, c("closest_is_close", "closest_is_distant"), c("closest_close", "closest_distant"))
				q$bindtype <- factor(q$bindtype, levels = c("closest_is_close", "closest_is_distant", "any_close", "any_distant"))

				)));
				state(runr(qq(

				# Plot
				# ggplot(q, aes(x=bindtype, fill=as.factor(motif2))) + scale_fill_manual(values=c("0"="#BEC1C000", "1"="#F3702A")) + geom_bar(position="fill") + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs"))
				# ggplot(q, aes(x=bindtype, fill=as.factor(motif2))) + scale_fill_manual(values=c("0"="grey", "1"="black")) + geom_bar(position="fill") + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs"))
				# ggplot(q, aes(x=bindtype, fill=motif2)) + scale_fill_manual(values=c("#BEC1C000", "#1E3D59")) + geom_bar(position="fill") + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + xlab(paste0("Relative ",rbp2," binding site position")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs")) + guides(fill=guide_legend(title=NULL))
				ggplot(q, aes(x=bindtype, fill=motif2)) + scale_fill_manual(values=c("#BEC1C000", "#1E3D59")) + geom_bar(position="fill") + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + xlab(paste0("Relative ",rbp2," binding site position")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs")) + guides(fill=F)
				ggsave("../output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$rbp1-$rbp2-plot-stacked-bars.pdf", width=91.5, height=91.5, units="mm")


				# str(qo)
				resamples <- $resamples
				)));
				state(runr(q(
				qr <- data.frame(motif1=character(), bindtype=character(), motif_fraction=numeric())
				for (i in 1:resamples) {
				  print(paste0(" >> Resample ",i,"/",resamples))
				  flush.console()
				  q <- qo[sample(round(nrow(qo)/2), replace=T),]
				  for (my_motif1 in unique(qo$motif1)) {
				    for (my_bindtype in unique(qo$bindtype)) {
				      # print(paste0(" >> ",my_motif1))
				      # print(paste0("   >> ",my_bindtype))
				      my_motif_fraction <- nrow(q[q$bindtype==my_bindtype & q$motif1==my_motif1 & q$motif2==1,]) / nrow(q[q$bindtype==my_bindtype & q$motif1==my_motif1,])
				      # print(paste0("     >> ",my_motif_fraction))
      
				      qr <- rbind(qr, data.frame("resample"=i, "motif1"=my_motif1, "bindtype"=my_bindtype, "motif_fraction"=my_motif_fraction))
				    }
				  }
				}
				# str(qr)
				# ggplot(qr, aes(x=bindtype, y=motif_fraction)) + geom_boxplot(notch=T)
				# Improve labelling
				qr$motif1 <- as.factor(qr$motif1)
				qr$motif1 <- mapvalues(qr$motif1, c(1, 0), c(paste0(rbp1," binding site\nwith motif"), paste0(rbp1," binding site\nwithout motif")))
				qr$motif1 <- factor(qr$motif1, levels = c(paste0(rbp1," binding site\nwith motif"), paste0(rbp1," binding site\nwithout motif")))

				qr$bindtype <- mapvalues(qr$bindtype, c("close", "distant"), c("any_close", "any_distant"))
				# qr$bindtype <- mapvalues(qr$bindtype, c("closest_is_close", "closest_is_distant"), c("closest_close", "closest_distant"))
				qr$bindtype <- factor(qr$bindtype, levels = c("closest_is_close", "closest_is_distant", "any_close", "any_distant"))

				)));
				state(runr(qq(

				# Boxplots
				ggplot(qr, aes(x=bindtype, y=motif_fraction, fill=bindtype)) + scale_fill_manual(values=c("closest_is_close"="#FF7F00", "closest_is_distant"="#95B3D7", "any_close"="#F9C07C", "any_distant"="#D1D8E8")) + coord_cartesian(ylim=c(0,1)) + geom_boxplot(notch=T) + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + xlab(paste0("Relative ",rbp2," binding site position")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs")) + guides(fill=F)
				ggsave("../output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$rbp1-$rbp2-boxplots.pdf", width=91.5, height=91.5, units="mm")


				# Bars with error

				#+++++++++++++++++++++++++
				# Function to calculate the mean and the standard deviation
				# for each group
				#+++++++++++++++++++++++++
				# data : a data frame
				# varname : the name of a column containing the variable
				#to be summariezed
				# groupnames : vector of column names to be used as
				# grouping variables
				data_summary <- function(data, varname, groupnames){
				  require(plyr)
				  summary_func <- function(x, col){
				    c(mean = mean(x[[col]], na.rm=TRUE),
				      sd = sd(x[[col]], na.rm=TRUE))
				  }
				  data_sum<-ddply(data, groupnames, .fun=summary_func,
				                  varname)
				  data_sum <- rename(data_sum, c("mean" = varname))
				  return(data_sum)
				}
				qrs <- data_summary(qr, varname="motif_fraction", groupnames=c("bindtype", "motif1"))
				# str(qrs)
				# qrs
				ggplot(qrs, aes(x=bindtype, y=motif_fraction, fill=bindtype)) + scale_fill_manual(values=c("closest_is_close"="#FF7F00", "closest_is_distant"="#95B3D7", "any_close"="#F9C07C", "any_distant"="#D1D8E8")) + coord_cartesian(ylim=c(0,1)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=motif_fraction-sd, ymax=motif_fraction+sd), width=.2, position=position_dodge(.9)) + facet_grid(cols=vars(motif1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0(rbp2," near ",rbp1,"\\n(<=$mindist_threshold nt)")) + xlab(paste0("Relative ",rbp2," binding site position")) + ylab(paste0("Fraction of ",rbp2," binding sites with motifs")) + guides(fill=F)
				ggsave("../output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$rbp1-$rbp2-errorbars.pdf", width=91.5, height=91.5, units="mm")


				


				)));
				
				
				# Clean up
				# run("Clean up", "rm -f '$tmpmotiffile'");
				
				
			}



			
		}
	}
}


# Define utility function: mindists_all_pairs
sub mindists_all_pairs
{
	my ($rbp1, $rbp2, $ensgv) = @_;
	our %peaks;
	
	my @tmp = ();
	
	if (switch('random'))
	{
		@tmp = (1..$resamples);
	}
	else
	{
		@tmp = ('');
	}
	
	my @my_mindists = ();

	foreach my $i (@tmp)
	{
		my $tmpi = '';
		$tmpi = "|$i" if (switch('random'));
		
		my @peaks1 = ();
		if (exists($peaks{"$rbp1|$ensgv$tmpi"}))
		{
			@peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv$tmpi"});
		}
		my @peaks2 = ();
		if (exists($peaks{"$rbp2|$ensgv$tmpi"}))
		{
			@peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv$tmpi"});
		}
	
		# state("ENSGV: $ensgv");
		# state("PEAKS1 $rbp1");
		# show(@peaks1);
		# state("PEAKS2 $rbp2");
		# show(@peaks2);
	
		# $mindist = undef;
		foreach my $peak1 (@peaks1)
		{
			foreach my $peak2 (@peaks2)
			{
				# if (!defined($mindist) or ($mindist > abs($peak2 - $peak1)))
				# {
					$mindist = abs($peak2 - $peak1);
					# print "     >> $ensgv >> $peak1 / $peak2 >> $mindist\n";
					push(@my_mindists, $mindist);
				# }
			}
		}

		# state("MY_MINDISTS");
		# show(@my_mindists);
	}

	return(@my_mindists);
}


# Define utility function: mindists_closest_per_peak
sub mindists_closest_per_peak
{
	my ($rbp1, $rbp2, $ensgv) = @_;
	our %peaks;
	
	my @tmp = ();
	
	if (switch('random'))
	{
		@tmp = (1..$resamples);
	}
	else
	{
		@tmp = ('');
	}
	
	my @my_mindists = ();
	
	# state("TMP::");
	# show(@tmp);

	foreach my $i (@tmp)
	{
		my $tmpi = '';
		$tmpi = "|$i" if (switch('random'));
		
		# print "PEAKS peaks{\"$rbp1|$ensgv$tmpi\"}\n";
		# print "PEAKS peaks{\"$rbp2|$ensgv$tmpi\"}\n";
		
		my @peaks1 = ();
		if (exists($peaks{"$rbp1|$ensgv$tmpi"}))
		{
			@peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv$tmpi"});
		}
		my @peaks2 = ();
		if (exists($peaks{"$rbp2|$ensgv$tmpi"}))
		{
			@peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv$tmpi"});
		}
	
		if (switch('random_balanced'))
		{
			# Balance the peak sets by subsampling them (without replacement) to the size of the smaller one
			@tmp_peaks1 = @peaks1;
			@tmp_peaks2 = @peaks2;
		
			# state("ENSGV: $ensgv");
			# state("TMP_PEAKS1 $rbp1");
			# show(@tmp_peaks1);
			# state("TMP_PEAKS2 $rbp2");
			# show(@tmp_peaks2);

			# peaks1
			@shuffled_indexes = shuffle(0..(scalar(@tmp_peaks1) - 1));
			@pick_indexes = @shuffled_indexes[0..(min(scalar(@tmp_peaks1), scalar(@tmp_peaks2)) - 1)];  
			@peaks1 = @tmp_peaks1[@pick_indexes];
		
			# peaks2
			@shuffled_indexes = shuffle(0..(scalar(@tmp_peaks2) - 1));
			@pick_indexes = @shuffled_indexes[0..(min(scalar(@tmp_peaks1), scalar(@tmp_peaks2)) - 1)];  
			@peaks2 = @tmp_peaks2[@pick_indexes];
		}
	
		# state("ENSGV: $ensgv");
		# state("PEAKS1 $rbp1");
		# show(@peaks1);
		# state("PEAKS2 $rbp2");
		# show(@peaks2);
	
		foreach my $peak1 (@peaks1)
		{
			$mindist = undef;
			foreach my $peak2 (@peaks2)
			{
				if (!switch('distance_signed'))
				{
					# Absolute
					if (!defined($mindist) or ($mindist > abs($peak2 - $peak1)))
					{
						$mindist = abs($peak2 - $peak1);
					}
				}
				else
				{
					# Signed
					if (!defined($mindist) or (abs($mindist) > abs($peak2 - $peak1)))
					{
						$mindist = $peak2 - $peak1;
					}
				}

				# print "     >> $ensgv >> $peak1 / $peak2 >> $mindist\n" if (switch('debug'));
			}

			push(@my_mindists, $mindist);
		}

		# state("MY_MINDISTS");
		# show(@my_mindists);
	}
	
	return(@my_mindists);
}


# Define utility function: mindist_motifs
sub mindist_motifs
{
	my ($rbp1, $rbp2, $ensgv) = @_;
	our %peaks;
	our %motifpeaks;
	
	# my @tmp = ();
	
	# Can't randomise the motifs
	# if (switch('random'))
	# {
	# 	@tmp = (1..$resamples);
	# }
	# else
	# {
		# @tmp = ('');
	# }
	
	my @my_motif_fractions = ();
	
	# state("TMP::");
	# show(@tmp);

	# foreach my $i (@tmp)
	# {
		my $tmpi = '';
		# Can't randomise the motifs
		# $tmpi = "|$i" if (switch('random'));
		
		# print "PEAKS peaks{\"$rbp1|$ensgv$tmpi\"}\n";
		# print "PEAKS peaks{\"$rbp2|$ensgv$tmpi\"}\n";
		
		my @peaks1 = ();
		if (exists($peaks{"$rbp1|$ensgv$tmpi"}))
		{
			@peaks1 = split(/\|/, $peaks{"$rbp1|$ensgv$tmpi"});
		}
		my @peaks2 = ();
		if (exists($peaks{"$rbp2|$ensgv$tmpi"}))
		{
			@peaks2 = split(/\|/, $peaks{"$rbp2|$ensgv$tmpi"});
		}
	
		my %motifpeaks1 = ();
		if (exists($motifpeaks{"$rbp1|$ensgv$tmpi"}))
		{
			foreach my $motifpeak (split(/\|/, $motifpeaks{"$rbp1|$ensgv$tmpi"}))
			{
				$motifpeaks1{$motifpeak} = 1;
			}
		}
		my %motifpeaks2 = ();
		if (exists($motifpeaks{"$rbp2|$ensgv$tmpi"}))
		{
			foreach my $motifpeak (split(/\|/, $motifpeaks{"$rbp2|$ensgv$tmpi"}))
			{
				$motifpeaks2{$motifpeak} = 1;
			}
		}
		
		# Can't randomise the motifs
		# if (switch('random_balanced'))
		# {
		# 	# Balance the peak sets by subsampling them (without replacement) to the size of the smaller one
		# 	@tmp_peaks1 = @peaks1;
		# 	@tmp_peaks2 = @peaks2;
		#
		# 	# state("ENSGV: $ensgv");
		# 	# state("TMP_PEAKS1 $rbp1");
		# 	# show(@tmp_peaks1);
		# 	# state("TMP_PEAKS2 $rbp2");
		# 	# show(@tmp_peaks2);
		#
		# 	# peaks1
		# 	@shuffled_indexes = shuffle(0..(scalar(@tmp_peaks1) - 1));
		# 	@pick_indexes = @shuffled_indexes[0..(min(scalar(@tmp_peaks1), scalar(@tmp_peaks2)) - 1)];
		# 	@peaks1 = @tmp_peaks1[@pick_indexes];
		#
		# 	# peaks2
		# 	@shuffled_indexes = shuffle(0..(scalar(@tmp_peaks2) - 1));
		# 	@pick_indexes = @shuffled_indexes[0..(min(scalar(@tmp_peaks1), scalar(@tmp_peaks2)) - 1)];
		# 	@peaks2 = @tmp_peaks2[@pick_indexes];
		# }
	
		# state("ENSGV: $ensgv");
		# state("PEAKS1 $rbp1");
		# show(@peaks1);
		# state("PEAKS2 $rbp2");
		# show(@peaks2);
		
		

		# The actual comparison part between @peaks1 and @peaks2 begins here:
	
		foreach my $peak1 (@peaks1)
		{
			my $mindist = undef;
			my $closest_peak = undef;
			
			# Find nearest peak
			foreach my $peak2 (@peaks2)
			{
				$dist = abs($peak2 - $peak1);

				if (!defined($mindist) or ($mindist > $dist))
				{
					$mindist = $dist;
					$closest_peak = $peak2;
				}
				
				print "     >> $ensgv >> $peak1 vs. $peak2 >> $dist\n" if (switch('debug'));;

				my $bindtype = undef;
				if ($dist <= $mindist_threshold)
				{
					# If this is a co-binding event...
					print "       >> CLOSE\n" if (switch('debug'));;
					$bindtype = 'close';
				}
				else
				{
					print "       >> DISTANT\n" if (switch('debug'));;
					$bindtype = 'distant';
				}

				my $motif1 = 0;
				$motif1 = 1 if (exists($motifpeaks1{$peak1}));
				my $motif2 = 0;
				$motif2 = 1 if (exists($motifpeaks2{$peak2}));

				print "         >> $ensgv >> PEAK1 $peak1 >> MOTIF $motif1\n" if (switch('debug'));
				print "         >> $ensgv >> PEAK2 $peak2 >> MOTIF $motif2\n" if (switch('debug'));

				print TMPMOTIFS "$rbp1\t$rbp2\t$ensgv\t$peak1\t$peak2\t$dist\t$bindtype\t$motif1\t$motif2\n";
			}
			
			# Write out the closest peak's motif information
			$peak2 = $closest_peak;

			my $motif1 = 0;
			$motif1 = 1 if (exists($motifpeaks1{$peak1}));
			my $motif2 = 0;
			$motif2 = 1 if (exists($motifpeaks2{$peak2}));

			print "     >> $ensgv >> $peak1 vs. $peak2 >> $dist\n" if (switch('debug'));;

			my $bindtype = undef;
			if ($mindist <= $mindist_threshold)
			{
				# If this is a co-binding event...
				print "       >> CLOSEST: CLOSE\n" if (switch('debug'));;
				$bindtype = 'closest_is_close';
			}
			else
			{
				print "       >> CLOSEST: DISTANT\n" if (switch('debug'));;
				$bindtype = 'closest_is_distant';
			}

			print "         >> $ensgv >> $peak1 >> CLOSEST >> $peak2 >> $mindist\n" if (switch('debug'));
			# $bindtype = 'closest';
			print TMPMOTIFS "$rbp1\t$rbp2\t$ensgv\t$peak1\t$peak2\t$mindist\t$bindtype\t$motif1\t$motif2\n";
			
			
			
			
			
			
			# exit;
			# push(@my_mindists, $mindist);
			# push(@my_motif_fractions, $mindist);
			
			
		}

		# state("MY_MINDISTS");
		# show(@my_mindists);
	# }

	
	# return(@my_mindists);
	# return(@my_motif_fractions);
}




showmeall(1);

print "Total "; stoptime2();

done();
