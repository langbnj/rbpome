#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table (e.g. rbpome)] [type: eclip_encode/...] [analysis mode 1: jaccard/correlation/mindist/mediandist/mindist_thresh/mindist_thresh_frac] [analysis mode 2] [set: positives/screen_hits/background_resamples] [number of top pairs to get]\n\nExample: $0 rbpome eclip_encode jaccard mindist_thresh screen_hits 20\nExample: $0 rbpome eclip_encode jaccard background_resamples screen_hits 20";
($table, $type, $mode1, $mode2, $set, $num) = args(6);
#args(0);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));

# output-pairs-jaccard-screen_hits.txt
$infile1 = "output-pairs-$table-$type-$mode1-$set.txt";
$infile2 = "output-pairs-$table-$type-$mode2-$set.txt";
$outfile = "output-bestpairs-$table-$type-$mode1-$mode2-$set-top$num.txt";

@skip = ();
if (switch('skip_screen_hits'))
{
	# Read screen_hits pairs and skip these below
	$skipfile = "output-bestpairs-$table-$type-$mode1-$mode2-screen_hits-top$num.txt";
	startme("Reading pairs from '$skipfile' and skipping these below");
	open(SKIP, $skipfile) or die("\nError: Couldn't open '$skipfile'\n\n");
	while (<SKIP>)
	{
		# ppig|sf3b4
		# u2af1|u2af2
		# rps3|sf3b4
		chomp;
		
		$pair = $_;
		
		push(@skip, $pair);
		
		stepme(1);
	}
	close(SKIP);
	stopme();
}
if (switch('skip_positives'))
{
	# Read positives pairs and skip these below
	$skipfile = "output-bestpairs-$table-$type-$mode1-$mode2-positives-top$num.txt";
	startme("Reading pairs from '$skipfile' and skipping these below");
	open(SKIP, $skipfile) or die("\nError: Couldn't open '$skipfile'\n\n");
	while (<SKIP>)
	{
		# ppig|sf3b4
		# u2af1|u2af2
		# rps3|sf3b4
		chomp;
		
		$pair = $_;
		
		push(@skip, $pair);
		
		stepme(1);
	}
	close(SKIP);
	stopme();
}

open(IN1, $infile1) or die("\nError: Couldn't open '$infile1'\n\n");
open(IN2, $infile2) or die("\nError: Couldn't open '$infile2'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

startme("Reading $num top pairs from '$infile1'");
@pairs1 = ();
<IN1>;	# Skip header
while (<IN1>)
{
	chomp;
	
	@a = split(/\t/);
	
	$rbp1 = $a[0];
	$rbp2 = $a[1];
	
	# # flip pairs
	# if ($rbp1 gt $rbp2)
	# {
	# 	($rbp2, $rbp1) = ($rbp1, $rbp2);
	# }

	$pair = $rbp1."|".$rbp2;
	
	# if (!contains($pair, @pairs1))
	if (!contains($pair, @pairs1) and !contains($pair, @skip))
	{
		push(@pairs1, $pair);
		stepme(1);
	}
	
	last if (getme() % $num == 0);
}
stopme();

startme("Reading $num top pairs from '$infile2'");
@pairs2 = ();
<IN2>;	# Skip header
while (<IN2>)
{
	chomp;
	
	@a = split(/\t/);
	
	$rbp1 = $a[0];
	$rbp2 = $a[1];
	
	# # flip pairs
	# if ($rbp1 gt $rbp2)
	# {
	# 	($rbp2, $rbp1) = ($rbp1, $rbp2);
	# }

	$pair = $rbp1."|".$rbp2;
	
	# if (!contains($pair, @pairs2))
	if (!contains($pair, @pairs2) and !contains($pair, @skip))
	{
		push(@pairs2, $pair);
		stepme(1);
	}
	
	last if (getme() % $num == 0);
}
stopme();



# Adjust $num if there aren't enough pairs (based on union)
if ($num > scalar(union(\@pairs1, \@pairs2)))
{
	$num = scalar(union(\@pairs1, \@pairs2));
	state("Readjusted top list length to '$num' (only this many total pairs were available in '$infile1' and '$infile2')");
}

startme("Writing intersection to '$outfile'");
%written = ();
$written1 = 0;
$written2 = 0;
# First, write the top $num intersection of the two lists
foreach $pair (intersection(\@pairs1, \@pairs2))
{
	print OUT "$pair\n";

	$written{$pair} = 1;
	foreach $_ (@pairs1)
	{
		if (defined($_) and ($_ eq $pair))
		{
			$_ = undef;
		}
	}
	foreach $_ (@pairs2)
	{
		if (defined($_) and ($_ eq $pair))
		{
			$_ = undef;
		}
	}
	$written1++;
	$written2++;

	stepme(1);
}
stopme();

# state("PAIRS1");
# show(@pairs1);
# state("PAIRS2");
# show(@pairs2);

startme("Filling alternately from list 1 and list2 '$outfile'");
while (scalar(keys(%written)) < $num)
{
	$i = 0;
	while (!defined($pairs1[$i]))
	{
		$i++;
	}
	if (defined($pairs1[$i]))
	{
		$pair = $pairs1[$i];
	
		print OUT "$pair\n";

		state("Adding from list position $i ($pairs1[$i])...");

		$written{$pair} = 1;
		$pairs1[$i] = undef;
		$written1++;

		stepme(1);
	}



	last if (scalar(keys(%written)) == $num);


	$i = 0;
	while (!defined($pairs2[$i]))
	{
		$i++;
	}
	if (defined($pairs2[$i]))
	{
		$pair = $pairs2[$i];

		print OUT "$pair\n";

		state("Adding from list position $i ($pairs2[$i])...");

		$written{$pair} = 1;
		$pairs2[$i] = undef;
		$written2++;

		stepme(1);
	}
}
stopme();


state("Wrote top $num '$set' pairs for '$type' '$mode1' ($written1) and '$mode2' ($written2) to '$outfile'");

showmeall(1);

run("cat", "cat $outfile");

done();
