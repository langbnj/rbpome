#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0";
args(0);

$annotfile = "input/attract/attract_db.txt";
$infile = "input/attract/pwm.txt";
$tmpfile = "input/attract_unsorted.meme";
$outfile = "input/attract.meme";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(ANNOT, $annotfile) or die("\nError: Couldn't open '$annotfile'\n\n");
open(TMP, ">$tmpfile") or die("\nError: Couldn't open '$tmpfile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print TMP "MEME version 4

ALPHABET= ACGU

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 U 0.25000 

strands: +

";


# start


startme("Reading annotation from '$annotfile'");
%symbol = ();
# %symbollen = ();
%species = ();
%mutated = ();
<ANNOT>;	# Skip header
while (<ANNOT>)
{
	chomp;

	# Gene_name	Gene_id	Mutated	Organism	Motif	Len	Experiment_description	Database	Pubmed	Experiment_description	Family	Matrix_id	Score
	# 3IVK	3IVK	no	Mus_musculus	GAAACA	6	X-RAY DIFFRACTION	PDB	19965478	X-RAY DIFFRACTION	N/A	519	1.000000**
	# 3IVK	3IVK	no	Mus_musculus	UGGG	4	X-RAY DIFFRACTION	PDB	19965478	X-RAY DIFFRACTION	N/A	574	1.000000**
	# 4KZD	4KZD	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	464	1.000000**
	# 4KZE	4KZE	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	437	1.000000**
	# 4Q9Q	4Q9Q	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	423	1.000000**
	# 4Q9R	4Q9R	no	Mus_musculus	CGAAAC	6	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	433	1.000000**
	# A1CF	ENSG00000148584	no	Homo_sapiens	UGAUCAGUAUA	11	UV cross-linking	R	10669759	UV cross-linking	RRM	110	1.000000**
	# A1CF	ENSG00000148584	no	Homo_sapiens	AUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.126952
	# A1CF	ENSG00000148584	no	Homo_sapiens	UUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.126411
	# A1CF	ENSG00000148584	no	Homo_sapiens	AUAAUUG	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.189114**
	# A1CF	ENSG00000148584	no	Homo_sapiens	UUAAUUG	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.188308
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	AUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.029649
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	UUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.04527
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	GUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.075703
	
	# Gene_name	( no need to explain :-) right?)
	# Gene_id	( no need to explain :-) right?)
	# Mutated	(if the target gene is mutated)
	# Organism	( no need to explain :-) right?)
	# Motif	( no need to explain :-) right?)
	# Len	(lenght of the motif)
	# Experiment_description(when available)
	# Database (Database from where the motifs were extracted PDB: Protein data bank, C: Cisbp-RNA, R:RBPDB, S: Spliceaid-F, AEDB:ASD)
	# Pubmed (pubmed ID)
	# Experiment (type of experiment; short description)
	# Family (domain)
	# Matrix_id (linked to the file PWM.txt)
	# Score (Qscore refer to the paper)
	#
	# The field Matrix_id refers to the pwm id that you can find in the pwm.txt file.
	# The position weight matrices are annotated in fasta format.

	@a = split(/\t/);
	
	$symbol = $a[0];
	# $ensg = $a[1];
	$mutated = $a[2];
	$fullspecies = $a[3];
	# $motif = $a[4];
	$motiflen = $a[5];
	# $exptype = $a[6];
	$motifid = $a[11];
	$score = $a[12];
	
	# Remove whitespace from symbol (necessary for NOVA1)
	$symbol =~ s/^\s+//;
	$symbol =~ s/\s+$//;
	
	# From the paper: "Score" is a quality score for SELEX experiments etc. (i.e. where there is data on multiple motifs). Where there is only one motif, the score is 1. Also, the top score seems to be highlighted with two asterisks: **.
	# I'll keep only lines with the asterisks (since I only care about the mapping).
	
	# # Skip mutated motifs
	# next if ($mutated ne 'no');
	
	# Skip non-"top" motifs (they'd be represented as one PWM anyway)
	if ($score !~ /\*\*$/)
	{
		next;
	}
	
	# # Skip non-human motifs
	# if ($fullspecies ne 'Homo_sapiens')
	# {
	# 	addme("skipped non-human motif for fullspecies", $fullspecies);
	# 	next;
	# }
	
	$symbol{$motifid} = $symbol;
	$motiflen{$motifid} = $motiflen;
	$species{$motifid} = $fullspecies;
	$mutated{$motifid} = $mutated;
	
	stepme(1000);
}
stopme();



startme("Converting '$infile' to MEME format and writing to '$tmpfile'");
starttime();
%titles = ();
$firstline = 1;
$out = '';
while (<IN>)
{
	chomp;
	
	# print " >> IN  >> $_\n";
	
	# >904	5
	# 0.00961538461538	0.00961538461538	0.00961538461538	0.971153846154
	# 0.00961538461538	0.00961538461538	0.971153846154	0.00961538461538
	# 0.00961538461538	0.00961538461538	0.971153846154	0.00961538461538
	# 0.00961538461538	0.00961538461538	0.971153846154	0.00961538461538
	# 0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
	# >s36	7
	# 0.844325153374	0.000766871165644	0.154141104294	0.000766871165644
	# 0.0774539877301	0.690950920245	0.0774539877301	0.154141104294
	# 0.000766871165644	0.0774539877301	0.154141104294	0.76763803681
	# 0.76763803681	0.0774539877301	0.000766871165644	0.154141104294
	# 0.0774539877301	0.0774539877301	0.76763803681	0.0774539877301
	# 0.230828220859	0.614263803681	0.0774539877301	0.0774539877301
	# 0.460889570552	0.230828220859	0.154141104294	0.154141104294
	
	# New title
	if (/^>(.+)/)
	{
		$title = $1;

		print " >> Title >> $title\n" if (switch('debug'));
	
		# Process title:
		# From >M120_0.6	7
		# To   >M120_0.6
		# (the "7" is just the length of the motif)
		if ($title =~ /^(\S+)\t(\d+)$/)
		{
			if (defined($motiflength))
			{
				if ($motiflength != $tmplen)
				{
					die("Error: Expected a motif length of $tmplen for title '$prevtitle', but got $motiflength");
				}
			}

			# New motif
			$motiflength = 0;

			$prevtitle = $title;
			$title = $1;
			$tmplen = $2;
			
			if (!exists($symbol{$title}))
			{
				die("Error: No symbol for title '$title'");
			}
			else
			{
				if ($mutated{$title} ne 'no')
				{
					addme("skipped mutated motif for motifid", $title);
					next;
				}
				if ($species{$title} ne 'Homo_sapiens')
				{
					addme("skipped non-human motif for motifid", $title);
					addme("skipped non-human motif for fullspecies", $fullspecies);
					next;
				}

				if ($tmplen ne $motiflen{$title})
				{
					# die("Error: Motif length doesn't match for title '$title' (should be ".$motiflen{$title}.", but is $tmplen)");
					addme("motif length mismatch between attract_db.txt and pwm.txt for title (kept)", $title);
				}
				
				$title = $symbol{$title};

				addme("RBP symbols with motifs", $title);
			}
		
			# Add a suffix (since there can be multiple motifs for one RBP gene symbol)
			$i = 1;
			while (exists($titles{$title.'_'.$i}))
			{
				$i++;
			}
			$title = $title.'_'.$i;
			die("Error: Title '$title' isn't unique in '$infile'") if (exists($titles{$title}));
			$titles{$title} = 1;
		}
		else
		{
			die("Error: Couldn't parse title '$title'");
		}
		
		if ($firstline != 1)
		{
			addme("successfully wrote output for title", $title);
			stepme(100);

			# Write motif
			print TMP $out;
			
			# Write footer
			print TMP "\nURL https://attract.cnic.es\n\n";
		}
		$firstline = 0;

		# Write header
		# nsites 20 is the default (and E 0), according to http://meme-suite.org/doc/meme-format.html
		$out = '';
		print TMP "MOTIF $title\n";
		print TMP "\nletter-probability matrix: alength= 4 w= $tmplen nsites= 20 E= 0\n";
		
		# New motif: Skip to next line
		next;
	}
	
	@a = split(/\t/);
	die("Error: Expected 4 columns for motif '$title', but got ".scalar(@a)) if (scalar(@a) != 4);
	
	($a, $c, $g, $t) = @a;
	# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($motiflength + 1)) if (round($a + $c + $g + $t, 6) != 1);
	die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($motiflength + 1)) if (($a + $c + $g + $t < 0.98) or ($a + $c + $g + $t > 1.02));
	# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($motiflength + 1)) if (($a + $c + $g + $t < 0.9) or ($a + $c + $g + $t > 1.2));
	# die("Error: Sum is ".($a + $c + $g + $t)." instead of 1 for motif '$title' position ".($motiflength + 1)) if (($a + $c + $g + $t < 0.65) or ($a + $c + $g + $t > 1.2));
	# if (round($a + $c + $g + $t, 6) != 1)
	# if (($a + $c + $g + $t < 0.9) or ($a + $c + $g + $t > 1.1))
	# {
	# 	# addme("frequencies don't sum to 1 for motif|site", $title.'|'.($i+1));
	# 	addme("frequencies don't sum to anything near 1 for motif|site (skipped motif)", $title.'|'.($i+1));
	# 	last;
	# }
		
	# print " >> IN  >> $a\t$c\t$g\t$t\n";

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
	
	$out .= "  ".round($a, 6)."\t  ".round($c, 6)."\t  ".round($g, 6)."\t  ".round($t, 6)."\n";
	$motiflength++;

	# print TMP "  ".round($a, 6)."\t  ".round($c, 6)."\t  ".round($g, 6)."\t  ".round($t, 6)."\n";
	# print " >> TMP >> "."  ".round($a, 6)."\t  ".round($c, 6)."\t  ".round($g, 6)."\t  ".round($t, 6)."\n";
	# print "\n\n >> TMP >> $out\n\n\n";
	# print " >> MOTIFLENGTH >> $motiflength\n";
}
if ($firstline != 1)
{
	addme("successfully wrote output for title", $title);
	stepme(100);

	# Write motif
	print TMP $out;
	
	# Write footer
	print TMP "\nURL https://attract.cnic.es\n\n";
}
close(TMP);
normalbreak();
stopme();
stoptime();

showmesome(50);
# showmeall();

state("Wrote to '$tmpfile'");

# Sort motifs by title
startme("Reading motifs in '$tmpfile'");
open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
$nl = $/;
$/ = "\nMOTIF ";
@motifs = ();
$head = <TMP>;
chomp($head);
print OUT $head;
while (<TMP>)
{
	chomp;
	push(@motifs, $_);
	stepme(100);
}
stopme();
close(TMP);
# $foot = pop(@motifs);
startme("Writing sorted motifs to '$outfile'");
foreach $motif (unique(@motifs))
{
	print OUT $/.$motif;
	stepme(100);
}
stopme();
$/ = $nl;
# print OUT $foot;
close(OUT);


state("Wrote to '$outfile'");

done();
