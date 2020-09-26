#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$species = 'human';
$motiftable = 'rbp_motif_list';

our $usage = "$0";
args(0);

# @grepsources = ('attract', 'dominguez', 'custom');
@fimosources = ('attract', 'cisbp', 'dominguez', 'rbpdb', 'rbpmap', 'rnacompete');

$infile = "input/GRCh38.p10.genome.fa";
$outfile = "tmp/tmp-geneseqs.txt";
$modelfile = "tmp/tmp-geneseqs-model.txt";
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

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

%motifs = ();
foreach $source (@fimosources)
{
	$source = lc($source);
	$motiffile = "input/$source.meme";

	startme("Loading list of RBPs with motifs from '$motiffile'", 1);
	open(MOTIFS, $motiffile) or die("Error: Couldn't open '$motiffile'");
	while (<MOTIFS>)
	{
		if (/^MOTIF (\S+)$/)
		{
			$symbol = $1;
			$symbol =~ s/_\d+$//;
			
			# If it's also an eCLIP protein...
			if (exists($acc{$symbol}))
			{
				$motifs{$symbol} = 1;
				addme("total eCLIP RBPs with motifs", $symbol);
				stepme(100, 1);
			}
		}
	}
	stopme(1);
}


$query = Query("SELECT DISTINCT symbol FROM `$motiftable` WHERE species='$species'");
startme("Loading list of RBPs with motifs from any source from table '$motiftable'", 1, Numrows($query));
while (($symbol) = Fetch($query))
{
	# If it's also an eCLIP protein...
	if (exists($acc{$symbol}))
	{
		$motifs{$symbol} = 1;
		addme("total eCLIP RBPs with motifs", $symbol);
		stepme(100, 1);
	}
}
stopme(1);


state("Total RBPs with motifs: ".scalar(keys(%motifs)));

# show(unique(keys(%motifs)));
# exit;





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



state("Getting genes from table '$intable' and writing their sequences to output files:");

# Cycle through all RBPs that have motifs and eCLIP data
$rbpquery = Query("SELECT symbol FROM `$intable` WHERE species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol ORDER BY symbol");
starttime();
$i = 0;
# $totalhits = 0;
%titles = ();
while (($symbol) = Fetch($rbpquery))
{
	$tmpfile = "tmp/tmp-geneseqs-$symbol.txt";
	$tmpmodelfile = "tmp/tmp-geneseqs-$symbol-model.txt";
	open(TMP, ">$tmpfile") or die("\nError: Couldn't open '$tmpfile'\n\n");

	# Does the RBP have at least one known motif?
	die("Error: Unexpected RBP '$symbol'") if (!exists($motifs{$symbol}));
	
	# Get UniProt accession
	die("Error: No UniProt accession for symbol '$symbol'") if (!exists($acc{$symbol}));
	$acc = $acc{$symbol};

	addme("total rbp symbols", $symbol);
	addme("total uniprot accessions", $acc);

	$mainquery = Query("SELECT DISTINCT ensgv FROM `$intable` WHERE species='$species' AND map='$map' AND symbol='$symbol'");
	$i++;
	startme(" >> $i / ".Numrows($rbpquery)." >> $symbol", 1, Numrows($mainquery));
	# while (($ensgv, $celltype, $rep, $chr, $start, $stop, $strand) = Fetch($mainquery))
	while (($ensgv) = Fetch($mainquery))
	{
		stepme(1000, 1);

		# Got at least one motif for this RBP
		next if (!exists($motifs{$symbol}));
		
		# Get $chr, $start, $stop, $strand
		$query = Query("SELECT chr, start, stop, strand FROM gencode_gff3_gene WHERE ensgv='$ensgv' AND species='$species'");
		($chr, $start, $stop, $strand) = FetchOne($query);
		
		# Get genomic sequence of the gene
		$geneseq = substr($chr{$chr}, $start - 1, ($stop - $start) + 1);

		if ($strand eq '-')
		{
			# $geneseq = strandflip($geneseq);
			$geneseq = strandflip_loosely($geneseq);	# Allow N and - (there's one N case that came up)
		}
	
		# Translate to RNA
		$geneseq =~ s/T/U/g;

		# Construct sequence title
		# $title = "$chr:$start-$stop";
		# $title = "$symbol|$celltype|$chr|$start|$stop";
		$title = "$symbol|$chr|$start|$stop|$strand";
		die("Error: Title '$title' occurs multiple times in '$tmpfile' (apparently this region exists as both '+' and '-' strands - might lead to problems with MEME FIMO)") if (exists($titles{$title}));
		$titles{$title} = 1;

		# Write to sequence file
		print TMP ">$title\n".split60($geneseq)."\n";
		print OUT ">$title\n".split60($geneseq)."\n";

		addme("total gene regions", $title);
		addme("distinct gene regions", "$chr|$start|$stop|$strand");
	}
	stopme(1);
	close(TMP);

	run("Calculating background model", "./bin/meme-5.0.5/src/fasta-get-markov -norc $tmpfile $tmpmodelfile > /dev/null 2>&1", 1);
	print "   >> ";
	run("cat", q(cat ).$tmpmodelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
	print "\n";
	# state("Wrote to '$tmpmodelfile'");
}
stoptime();



showmeall(1);

state("Calculating combined background model:");
run("Calculating combined background model", "./bin/meme-5.0.5/src/fasta-get-markov -norc $outfile $modelfile > /dev/null 2>&1", 1);
run("cat", q(cat ).$modelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
nl();
state("Wrote to '$modelfile'");

done();
