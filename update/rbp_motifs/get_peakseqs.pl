#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
# $infolder = '../../pipeline/tom_functional_enrichment_analysis/input/results_txt/';
# $outtable = 'rbp_motif_instances';
$species = 'human';
$motiftable = 'rbp_motif_list';

our $usage = "$0 [method: grep/fimo] [type: eclip_encode/eclip_tom/...] [motif source] [-extend]\n\n -extend: Extend peaks 5' and 3' by 150 nucleotides, as \"the 5' start of the peak is predicted to correspond to the site of crosslink between the RBP and the RNA.\" (Dominguez et al. 2018). eclip_encode and eclip_encode_12 peaks are always extended 50nt 5' even when this is off, or 150nt in both directions when this is on.\n\nExample: $0 fimo eclip_encode attract\nExample: $0 grep eclip_tom attract";
($method, $type, $source) = args(3);

$tmpextend = '';
if (switch('extend'))
{
	$tmpextend = '-extend';
}

$source = lc($source);
$motiffile = "input/$source.meme";

$infile = "input/GRCh38.p10.genome.fa";
# $outfile = "tmp/tmp-peakseqs$tmpextend-$method-$type-$source-combined.txt";
# $modelfile = "tmp/tmp-peakseqs$tmpextend-$method-$type-$source-combined-model.txt";
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

$map = 'gene';
$intable = "clip_raw_$map";



# start

# Get RBP accessions
startme("Getting RBP symbol-to-acc mapping from table '$intable'");
$query = Query("SELECT symbol, acc FROM `$intable` WHERE species='$species' AND map='$map' GROUP BY symbol, acc ORDER BY symbol, acc");
%acc = ();
while (($symbol, $acc) = Fetch($query))
{
	$acc{$symbol} = $acc;

	stepme(10000);
}
stopme();

if ($method eq 'fimo')
{
	startme("Loading list of RBPs with motifs from '$motiffile'");
	open(MOTIFS, $motiffile) or die("Error: Couldn't open '$motiffile'");
	%motifs = ();
	while (<MOTIFS>)
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
}
elsif ($method eq 'grep')
{
	$query = Query("SELECT DISTINCT symbol FROM `$motiftable` WHERE species='$species' AND source='$source'");
	startme("Loading list of RBPs with '$source' motifs from table '$motiftable'", 0, Numrows($query));
	%motifs = ();
	while (($symbol) = Fetch($query))
	{
		$motifs{$symbol} = 1;
		stepme(100);
	}
	stopme();
}
else
{
	die("Unhandled method '$method'");
}

startme("Loading main chromosomes (plus X/Y/M) from '$infile' into memory");
starttime();
fastabreak();
%chr = ();
while (<IN>)
{
	($title, $seq) = getfasta($_);
	
	$title =~ /^(\S+)/ or die("Error: Couldn't match title '$title'");
	$chr = $1;
	
	if ($chr !~ /^chr(\d+|M|X|Y)$/)
	{
		addme("scaffold skipped during loading (only need main chromosomes)", $chr);
		next;
	}
	addme("chromosome loaded", $chr);
	
	$chr{$chr} = $seq;
	
	stepme(1);
}
normalbreak();
stopme();
stoptime();



state("Getting '$type' peaks from table '$intable' and writing their sequences to output files:");

# Cycle through all RBPs that have motifs and eCLIP data
$rbpquery = Query("SELECT symbol, celltype, rep FROM `$intable` WHERE type='$type' AND species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol, celltype, rep ORDER BY symbol, celltype, rep");
starttime();
$i = 0;
# $totalhits = 0;
%titles = ();
while (($symbol, $celltype, $rep) = Fetch($rbpquery))
{
	# #DEBUG
	# next unless (($symbol eq 'HNRNPC') and ($celltype eq 'HEK293'));
	# #END DEBUG

	$tmpfile = "tmp/tmp-peakseqs$tmpextend-$method-$type-$source-$symbol-$celltype-$rep.txt";
	# $tmpmodelfile = "tmp/tmp-peakseqs$tmpextend-$method-$type-$source-$symbol-$celltype-$rep-model.txt";
	open(TMP, ">$tmpfile") or die("\nError: Couldn't open '$tmpfile'\n\n");

	# Does the RBP have at least one known motif?
	die("Error: Unexpected RBP '$symbol'") if (!exists($motifs{$symbol}));
	
	# Get UniProt accession
	die("Error: No UniProt accession for symbol '$symbol'") if (!exists($acc{$symbol}));
	$acc = $acc{$symbol};

	addme("$type: total rbp symbols", $symbol);
	addme("$type: total uniprot accessions", $acc);
	addme("$type: total rbp symbols|celltype|reps", "$symbol|$celltype|$rep");
	# addme("$type: total peak regions", $title);

	$mainquery = Query("SELECT DISTINCT chr, start, stop, strand FROM `$intable` WHERE type='$type' AND species='$species' AND map='$map' AND symbol='$symbol' AND celltype='$celltype' AND rep='$rep' ORDER BY chr, start, stop, strand");
	$i++;
	startme(" >> $i / ".Numrows($rbpquery)." >> $symbol >> $celltype >> rep$rep", 1, Numrows($mainquery));
	# while (($ensgv, $celltype, $rep, $chr, $start, $stop, $strand) = Fetch($mainquery))
	while (($chr, $start, $stop, $strand) = Fetch($mainquery))
	{
		# $rep = '';
	
		stepme(1000, 1);

		# Got at least one motif for this RBP
		next if (!exists($motifs{$symbol}));
		
		
		
		# 5' and 3' extensions
		$tmpstart = $start;
		$tmpstop = $stop;
		
		# Even with -extend off:
		# Always extend ENCODE CLIPper peaks 5' by 50 nt, since this is what they seem to consider the true binding region.
		# From Dominguez et al. 2018 RNA Bind-N-Seq:
		# "Peaks were also extended 50 nucleotides in the 5' direction as the 5' start of the peak is predicted to correspond to the site of crosslink between the RBP and the RNA."
		if ($type =~ /^eclip_encode/)	# eclip_encode and eclip_encode_12
		{
			if ($strand eq '+')
			{
				$tmpstart -= 50;
			}
			elsif ($strand eq '-')
			{
				$tmpstop += 50;
			}
		}

		# Extend 5' and 3' by 150 nt
		# (For eclip_encode and eclip_encode_12, this is in addition to their 5' 50 nt extension)
		if (switch('extend'))
		{
			# 
			if ($strand eq '+')
			{
				$tmpstart -= 150;
				$tmpstop += 150;
			}
			elsif ($strand eq '-')
			{
				$tmpstart -= 150;
				$tmpstop += 150;
			}
		}
		
		# Get genomic sequence of the peak
		$peakseq = substr($chr{$chr}, $tmpstart - 1, ($tmpstop - $tmpstart) + 1);

		if ($strand eq '-')
		{
			# $peakseq = strandflip($peakseq);
			$peakseq = strandflip_loosely($peakseq);	# Allow N and - (there's one N case that came up)
		}
	
		# Translate to RNA
		$peakseq =~ s/T/U/g;

		# Construct sequence title
		# $title = "$chr:$start-$stop";
		# $title = "$symbol|$celltype|$chr|$start|$stop";
		$title = "$symbol|$celltype|$rep|$chr|$start|$stop|$strand";
		die("Error: Title '$title' occurs multiple times in '$tmpfile' (apparently this region exists as both '+' and '-' strands - might lead to problems with MEME FIMO)") if (exists($titles{$title}));
		$titles{$title} = 1;

		# Write to sequence file
		print TMP ">$title\n".split60($peakseq)."\n";
		# print OUT ">$title\n".split60($peakseq)."\n";

		addme("$type: total peak regions", $title);
		addme("$type: distinct peak regions", "$chr|$start|$stop|$strand");
	}
	stopme(1);
	close(TMP);

	# run("Calculating background model", "./bin/meme-5.0.5/src/fasta-get-markov -norc $tmpfile $tmpmodelfile > /dev/null 2>&1", 1);
	# print "   >> ";
	# run("cat", q(cat ).$tmpmodelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
	# print "\n";
	# # state("Wrote to '$tmpmodelfile'");
}
stoptime();



showmeall(1);

# state("Calculating combined background model:");
# run("Calculating combined background model", "./bin/meme-5.0.5/src/fasta-get-markov -norc $outfile $modelfile > /dev/null 2>&1", 1);
# run("cat", q(cat ).$modelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
# nl();
# state("Wrote to '$modelfile'");

done();
