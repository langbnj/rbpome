#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [motif source]\n\nExample: $0 attract";
($source) = args(1);
# args(0);

$infile = "input/$source.meme";

# start

startme("Loading list of RBPs with motifs from '$infile'");
open(IN, $infile) or die("Error: Couldn't open '$infile'");
%motifs = ();
while (<IN>)
{
	if (/^MOTIF (\S+)$/)
	{
		$symbol = $1;
		$symbol =~ s/_\d+$//;
		
		if (!exists($motifs{$symbol}))
		{
			$motifs{$symbol} = 1;
			stepme(100);
		}
	}
}
stopme();
close(IN);


# Sort motifs by title
$i = 0;
foreach $symbol (unique(keys(%motifs)))
{
	$i++;
	startme(" >> $i / ".scalar(keys(%motifs))." >> $symbol", 1);
	open(IN, $infile) or die("Error: Couldn't open '$infile'");
	$nl = $/;
	$/ = "\nMOTIF ";
	$head = <IN>;
	chomp($head);
	@motifs = ();
	while (<IN>)
	{
		chomp;
		# warn("['$_']");
		# ['A1CF_1
		#
		# letter-probability matrix: alength= 4 w= 11 nsites= 20 E= 0
		#   0.009615	  0.009615	  0.009615	  0.971154
		#   0.009615	  0.009615	  0.971154	  0.009615
		#   0.971154	  0.009615	  0.009615	  0.009615
		#   0.009615	  0.009615	  0.009615	  0.971154
		#   0.009615	  0.971154	  0.009615	  0.009615
		#   0.971154	  0.009615	  0.009615	  0.009615
		#   0.009615	  0.009615	  0.971154	  0.009615
		#   0.009615	  0.009615	  0.009615	  0.971154
		#   0.971154	  0.009615	  0.009615	  0.009615
		#   0.009615	  0.009615	  0.009615	  0.971154
		#   0.971154	  0.009615	  0.009615	  0.009615
		#
		# URL https://attract.cnic.es
		# ']
	
		/^(\S+)_(\d+)\n/ or die("Error: Couldn't parse motif title in motif \n\n'$_'\n\n");
		$thissymbol = $1;
		$thismotifid = $1.'_'.$2;
	
		addme("total symbols with motifs", $thissymbol);
		addme("total motifids", $thismotifid);
	
		if (/^$symbol\_\d+\n/)
		{
			addme("motifids written to output", $thismotifid);
		
			push(@motifs, $_);
			stepme(1, 1);
		}
	}
	stopme(1);
	close(IN);
	# $foot = pop(@motifs);
	if (scalar(@motifs) > 0)
	{
		# startme("Writing sorted motifs to '$outfile'");
		$outfile = "tmp/tmp-motifs-$source-$symbol.meme";
		open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
		print OUT $head;
		foreach $motif (unique(@motifs))
		{
			print OUT $/.$motif;
			# stepme(100);
		}
		# stopme();
		$/ = $nl;
		# print OUT $foot;
		close(OUT);
		# state("Wrote to '$outfile'");
	}
	else
	{
		# state("No motifs found from source '$source' for RBP '$symbol'");
		addme("no motifs found for rbp", $symbol);
	}
}


showmeall(1);

done();
