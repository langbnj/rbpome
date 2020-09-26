#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0";
args(0);

$infolder = 'input/rbpmap/RBP_PSSMs';
$outfile = "input/rbpmap.meme";
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "MEME version 4

ALPHABET= ACGU

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 U 0.25000 

strands: +

";


# start


startme("Converting PSSMs in '$infolder' to MEME format and writing to '$outfile'");
starttime();
open(LS, "ls -1 $infolder/*_human_PSSM.txt|") or die("Error: Couldn't ls input folder '$infolder");
%titles = ();
while (<LS>)
{
	chomp;
	
	$infile = $_;
	
	# A1CF_wuaauur_human_PSSM.txt
	$infile =~ /^$infolder\/([^_]+)_([a-z]+)_human_PSSM\.txt$/ or die("Error: Couldn't parse input filename '$infile'");
	$symbol = $1;
	$title = $symbol;
	
	# Add a suffix (since there can be multiple motifs for one RBP gene symbol)
	$i = 1;
	while (exists($titles{$title.'_'.$i}))
	{
		$i++;
	}
	$title = $title.'_'.$i;
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
	
	<IN>;	# Skip header
	
	# Write header
	# nsites 20 is the default (and E 0), according to http://meme-suite.org/doc/meme-format.html
	print OUT "MOTIF $title\n";
	print OUT "\nletter-probability matrix: alength= 4 w= ".(chompme(`cat $infile | wc -l`) - 1)." nsites= 20 E= 0\n";
	
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
	print OUT "\nURL http://rbpmap.technion.ac.il\n\n";
	
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
