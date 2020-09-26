#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [input file]\n\n$0 all_PWMs_human.masta";
($infile) = args(1);
# args(0);

$outfile = $infile;
$outfile =~ s/\.masta$/.meme/ or die("Error: Couldn't generate output filename from '$infile'");

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "MEME version 4

ALPHABET= ACGU

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 U 0.25000 

";


# start

startme("Converting MASTA file '$infile' to MEME format");
starttime();
fastabreak();
%titles = ();
while ($entry = <IN>)
{
	chomp($entry);
	
	# >1004_8676391:1367_NCL
	# -0.97243147 -0.97243147  -2.4521041   1.7173493   1.7873355   1.7873355 -0.97243147  -2.4521041  -2.4521041  -2.4521041  -2.4521041   1.7873355   1.7873355  -2.4521041  -2.4521041   1.7873355  -1.5303424  -2.4521041
	#   1.6437941    1.484377 -0.97243147  -1.5303424  -2.4521041  -2.4521041  -2.4521041   1.7873355   1.7873355   1.7873355  -2.4521041  -2.4521041  -2.4521041  -1.5303424  -2.4521041  -2.4521041 -0.25747839  -1.5303424
	#  -2.4521041  -1.5303424   1.6437941  -2.4521041  -2.4521041  -2.4521041  0.72578896  -2.4521041  -2.4521041  -2.4521041   1.7873355  -2.4521041  -2.4521041   1.7173493  -1.5303424  -2.4521041   1.2063992   1.5662864
	#  -2.4521041  -1.5303424  -2.4521041  -2.4521041  -2.4521041  -2.4521041  0.72578896  -2.4521041  -2.4521041  -2.4521041  -2.4521041  -2.4521041  -2.4521041  -2.4521041   1.7173493  -2.4521041 -0.97243147 -0.97243147
	#
	# >1052_17318228:1720_RBMY1A1
	# -2.8408487 0.29446094 -2.8408487  1.8404688  1.8404688
	#  1.8404688 -2.1193573  1.8404688 -2.8408487 -2.8408487
	# -2.8408487 -2.8408487 -2.8408487 -2.8408487 -2.8408487
	# -2.8408487   1.265344 -2.8408487 -2.8408487 -2.8408487
	
	# print "\nENTRY=['$entry']\n\n";
	
	$entry =~ s/^>?(.+)\n// or die("Error: Couldn't get motif title from entry '$entry'");
	$title = $1;
	
	print " >> Title >> $title\n" if (switch('debug'));
	print "   >> Entry:\n$entry\n" if (switch('debug'));
	
	# Process title:
	# From >1004_8676391:1367_NCL
	# To   >NCL
	# (RBP gene symbol only)
	if ($title =~ /^(\S+):(\S+)$/)
	{
		$title = "$2";
		$title =~ s/^\d+_//;
		
		# Add a suffix (since there can be multiple motifs for one RBP gene symbol)
		$i = 1;
		while (exists($titles{$title.'_'.$i}))
		{
			$i++;
		}
		$title = $title.'_'.$i;
		$titles{$title} = 1;
	}
	
	@lines = split(/\n/, $entry);
	die("Error: Expected 4 lines for motif '$title', but got ".scalar(@lines)." in entry:\n\n$entry\n\n") if (scalar(@lines) != 4);
	
	foreach (@lines)
	{
		chomp;
		s/^\s+//;
		s/\s+$//;
	}
	
	@a = split(/ +/, $lines[0]);
	@c = split(/ +/, $lines[1]);
	@g = split(/ +/, $lines[2]);
	@t = split(/ +/, $lines[3]);
	
	# print "\nA=['".join('|', @a)."']\n\n";
	
	die("Error: Expected a motif length of ".scalar(@a).", but got ".scalar(@a)."/".scalar(@c)."/".scalar(@g)."/".scalar(@t)) if ((scalar(@a) != scalar(@c)) or (scalar(@a) != scalar(@c)) or (scalar(@a) != scalar(@c)) or (scalar(@a) != scalar(@c)));
	
	$motiflength = scalar(@a);
	
	$i = 0;
	$firstline = 1;
	while ($i < $motiflength)
	{
		$a = $a[$i];
		$c = $c[$i];
		$g = $g[$i];
		$t = $t[$i];
		
		# Transform log-odds values (base 2) to probabilities (assuming equiprobability, i.e. 0.25/0.25/0.25/0.25)
		$a = 2 ** $a / 4;
		$c = 2 ** $c / 4;
		$g = 2 ** $g / 4;
		$t = 2 ** $t / 4;
		
		# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($i + 1)) if (round($a + $c + $g + $t, 6) != 1);
		# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($i + 1)) if (($a + $c + $g + $t < 0.9) or ($a + $c + $g + $t > 1.2));
		# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($i + 1)) if (($a + $c + $g + $t < 0.65) or ($a + $c + $g + $t > 1.2));
		# if (round($a + $c + $g + $t, 6) != 1)
		if (($a + $c + $g + $t < 0.9) or ($a + $c + $g + $t > 1.1))
		{
			# addme("frequencies don't sum to 1 for motif|site", $title.'|'.($i+1));
			addme("frequencies don't sum to anything near 1 for motif|site (skipped motif)", $title.'|'.($i+1));
			last;
		}
		
		# MEME version 4
		#
		# ALPHABET= ACGU
		#
		# Background letter frequencies (from uniform background):
		# A 0.25000 C 0.25000 G 0.25000 U 0.25000
		#
		# MOTIF RNCMPT00001 A1CF
		#
		# letter-probability matrix: alength= 4 w= 7 nsites= 20 E= 0
		#   0.395329	  0.105514	  0.105514	  0.393643
		#   0.007705	  0.007705	  0.007705	  0.976886
		#   0.976886	  0.007705	  0.007705	  0.007705
		#   0.976886	  0.007705	  0.007705	  0.007705
		#   0.007705	  0.007705	  0.007705	  0.976886
		#   0.007705	  0.007705	  0.007705	  0.976886
		#   0.321131	  0.143804	  0.478371	  0.056694
		#
		# URL http://hugheslab.ccbr.utoronto.ca/supplementary-data/RNAcompete_eukarya
		#
		# MOTIF RNCMPT00002 ANKHD1
		#
		# letter-probability matrix: alength= 4 w= 7 nsites= 20 E= 0
		#   0.773203	  0.075599	  0.075599	  0.075599
		#   0.004305	  0.004305	  0.987086	  0.004305
		#   0.987086	  0.004305	  0.004305	  0.004305
		#   0.004305	  0.987086	  0.004305	  0.004305
		#   0.004305	  0.004305	  0.987086	  0.004305
		#   0.228810	  0.102089	  0.004305	  0.664797
		#   0.427754	  0.079309	  0.099730	  0.393206
		#
		# URL http://hugheslab.ccbr.utoronto.ca/supplementary-data/RNAcompete_eukarya
		
		
		if ($firstline == 1)
		{
			# Write header
			# nsites 20 is the default (and E 0), according to http://meme-suite.org/doc/meme-format.html
			print OUT "MOTIF $title\n";
			print OUT "\nletter-probability matrix: alength= 4 w= $motiflength nsites= 20 E= 0\n";
			$firstline = 0;
		}
		print OUT "  ".round($a, 6)."\t  ".round($c, 6)."\t  ".round($g, 6)."\t  ".round($t, 6)."\n";

		$i++;
	}
	
	if ($firstline != 1)
	{
		print OUT "\nURL http://rbpdb.ccbr.utoronto.ca\n\n";
	}
	
	stepme(1);
}
close(OUT);
normalbreak();
stopme();
stoptime();

showmesome(50);

state("Wrote to '$outfile'");

done();
