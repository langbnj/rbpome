#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0";
args(0);

$infolder = 'input/dominguez/RBNS_PWMs';
$outfile = "input/dominguez.meme";
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "MEME version 4

ALPHABET= ACGU

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 U 0.25000 

strands: +

";


# start


# From the paper:
# Overlapping Specificities of RBPs
# To visualize and compare the primary sequence specificities, we derived sequence motif logos for each RBP by aligning enriched 5-mers (Z score R3, weighted by enrichment above input using an iterative procedure that avoids overlap issues, Figure 1A, top right; STAR Methods). For roughly half of the RBPs (41/78), this method yielded multiple sequence logos, indicating affinity to multiple distinct motifs that may reflect different binding modes or binding by distinct RBDs (motif 5-mers are listed in Table S3). Clustering proteins based on their top logo, paralogs (e.g., PCBP1/2/4, RBFOX2/3) clustered tightly (Conway et al., 2016; Smith et al., 2013) (Figure 2A). However, unexpectedly, many completely unrelated proteins, often containing distinct types of RBDs, were also grouped together. Fifteen clusters of RBPs with highly similar primary motifs (nine with three or more mem- bers) emerged, leaving 18 RBPs with more distinct motifs unclustered (STAR Methods).
# >> Can use all PWMs (or just the top one). They're all significant.

startme("Converting PWMs in '$infolder' to MEME format and writing to '$outfile'");
starttime();
open(LS, "ls -1 $infolder/*_5mer_logo*.PWM|") or die("Error: Couldn't ls input folder '$infolder");
%titles = ();
while (<LS>)
{
	chomp;
	
	$infile = $_;
	
	# A1CF_5mer_logo0.PWM
	$infile =~ /^$infolder\/([^_]+)_5mer_logo(\d+)\.PWM$/ or die("Error: Couldn't parse input filename '$infile'");
	$symbol = $1;
	$title = $symbol.'_'.($2 + 1);
	
	# # Add a suffix (since there can be multiple motifs for one RBP gene symbol)
	# $i = 1;
	# while (exists($titles{$title.'_'.$i}))
	# {
	# 	$i++;
	# }
	# $title = $title.'_'.$i;
	die("Error: Title '$title' isn't unique'") if (exists($titles{$title}));
	$titles{$title} = 1;
	
	open(IN, $infile) or die("Error: Couldn't open '$infile'");
	
	# Pos	A	C	G	U
	# 1	0.39532879396435	0.105513888686126	0.105513888686126	0.393643427745405
	# 2	0.00770456803068082	0.00770456803068082	0.00770456803068082	0.976886297348457
	# 3	0.976886297348457	0.00770456803068082	0.00770456803068082	0.00770456803068082
	# 4	0.976886297348457	0.00770456803068082	0.00770456803068082	0.00770456803068082
	# 5	0.00770456803068082	0.00770456803068082	0.00770456803068082	0.976886297348457
	# 6	0.00770456803068082	0.00770456803068082	0.00770456803068082	0.976886297348457
	# 7	0.321131137484576	0.143803698114993	0.478370946765367	0.0566942181612342
	
	# ID test
	# PO	A	C	G	U
	# 0	0.200580232093	0.120348139256	0.120348139256	0.558723489396
	# 1	0.402160864346	0.077430972389	0.208483393357	0.311924769908
	# 2	0.0	0.0	0.0	1.0
	# 3	0.0	0.0	1.0	0.0
	# 4	0.0	0.0	0.0	1.0
	# 5	0.12975190076	0.33793517407	0.129551820728	0.402761104442
	# 6	0.172669067627	0.262705082033	0.172468987595	0.392156862745

	<IN>; <IN>;	# Skip header
	
	# Write header
	# nsites 20 is the default (and E 0), according to http://meme-suite.org/doc/meme-format.html
	print OUT "MOTIF $title\n";
	print OUT "\nletter-probability matrix: alength= 4 w= ".(chompme(`cat $infile | wc -l`) - 2)." nsites= 20 E= 0\n";
	
	$i = 0;
	while (<IN>)
	{
		chomp;
		
		@a = split(/\t/);
		
		die("Error: Expected 5 columns for motif '$title', but got ".scalar(@a)) if (scalar(@a) != 5);
	
		$a = $a[1];
		$c = $a[2];
		$g = $a[3];
		$u = $a[4];
		
		$i++;
		
		die("Error: Sum is ".($a + $c + $g + $u)." instead of 1 for motif '$title' position $i") if (($a + $c + $g + $u < 0.98) or ($a + $c + $g + $u > 1.02));

		print OUT "  ".round($a, 6)."\t  ".round($c, 6)."\t  ".round($g, 6)."\t  ".round($u, 6)."\n";
	}
	close(IN);

	# Write footer
	print OUT "\nURL https://www.encodeproject.org\n\n";
	
	addme("total RBP gene symbols written", $symbol);
	addme("total motif titles written", $title);
	
	stepme(1);
}
close(OUT);
stopme();
stoptime();

showmeall(1);

state("Wrote to '$outfile'");

done();
