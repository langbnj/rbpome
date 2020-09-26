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
# $resamples = 100;

# our $usage = "$0 [type (e.g. eclip_encode_12)] [score type (e.g. SUM_IS)] [-randompairs]\n\n -randompairs: Read reproducible background resample random pairs from e.g. output/output-rbpome-randompairs.txt and use these instead of real RBP pairs\n\nExample: $0 eclip_encode_12 SUM_IS";
our $usage = "$0 [table (e.g. rbpome)] [type (e.g. eclip_encode_12)] [score type (e.g. SUM_IS)]\n\nExample: $0 rbpome eclip_encode_12 SUM_IS";
($table, $type, $scoretype) = args(3);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

# $tmprandompairs = '';
# if (switch('randompairs'))
# {
# 	$tmprandompairs = '-randompairs';
# }


# $outfile = "output-fit_vs_random-txt-$table-$type-$scoretype$tmprandompairs.txt";
$outfile = "output-fit_vs_random-txt-$table-$type-$scoretype.txt";
# $outpdf = "output-fit_vs_random-cdf-$table-$type-$scoretype.pdf";
$randomfile = "output/output-$table-randompairs.txt";

open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
print OUT "symbol1\tsymbol2\trandompair\trandom\tmindist\n";


# start

# Get RBP pairs
# if (switch('randompairs'))
# {
	@randompairs = ();
	open(RAND, "$randomfile") or die("Error: Couldn't open '$randomfile'");
	while (<RAND>)
	{
		chomp;
		push(@randompairs, $_);
	}
	close(RAND);
	
	@randompairs = addflip(@randompairs); # Add inverse pairs
	@randompairs = unique(@randompairs);

	state("Read ".scalar(@randompairs)." unique reproducible random pairs from '$randomfile'");
# }
# else
# {
	@pairs = ();
	startme("Getting real RBP pairs from table '$table'");
	$query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
	while (($this_symbol1, $this_symbol2) = Fetch($query))
	{
		push(@pairs, "$this_symbol1|$this_symbol2");
		push(@pairs, "$this_symbol2|$this_symbol1");

		stepme(100);
		stepme(100);
	}
	stopme();
	
	# Combine
	push(@pairs, @randompairs);
	
	@pairs = unique(@pairs);
	state("Got ".scalar(@pairs)." total unique pairs");
# }




# $query = Query("SELECT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
# startme("Consolidating mindists from the output files for individual pairs into one table", 0, Numrows($query) * 2);
startme("Consolidating mindists from the output files for individual pairs into one table", 0, scalar(@pairs));
starttime();
startr();
# while (($this_symbol1, $this_symbol2) = Fetch($query))
# {
	# foreach $pair ("$this_symbol1|$this_symbol2", "$this_symbol2|$this_symbol1")
	foreach $pair (@pairs)
	{
		$randompair = 0;
		if (contains($pair, @randompairs))
		{
			$randompair = 1;
		}
		
		($symbol1, $symbol2) = split(/\|/, $pair);
		
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-HNRNPK-fit-mindists-random.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-HNRNPK-fit-mindists.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-PCBP1-fit-mindists-random.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-PCBP1-fit-mindists.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-PCBP2-fit-mindists-random.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-APOBEC3C-PCBP2-fit-mindists.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-CPEB4-CSTF2T-fit-mindists-random.txt
		# output-distance-plot-closest_per_peak-eclip_encode_12-SUM_IS-CPEB4-CSTF2T-fit-mindists.txt
		
		
		$realfile = "output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$symbol1-$symbol2-fit-mindists.txt";
		$randomfile = "output/output-distance-plot-closest_per_peak-$table-$type-$scoretype-$symbol1-$symbol2-fit-mindists-random.txt";

		# Real data
		open(TMP, $realfile) or die("Error: Couldn't open '$realfile");
		while (<TMP>)
		{
			chomp;
			print OUT "$symbol1\t$symbol2\t$randompair\t0\t$_\n";
		}
		close(TMP);
		
		# Random data
		open(TMP, $randomfile) or die("Error: Couldn't open '$randomfile");
		while (<TMP>)
		{
			chomp;
			print OUT "$symbol1\t$symbol2\t$randompair\t1\t$_\n";
		}
		close(TMP);
		
		stepme(1);
	}
# }
close(OUT);
stopme();
stoptime();

showmeall(1);

state("Wrote to '$outfile'");

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
