#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$species = 'human';
$type = 'geneseqs';
$map = 'gene';
$intable = "clip_raw_$map";

our $usage = "$0 [p/q-value threshold] [-qvalue]\n\n -qvalue: Use a q-value threshold (instead of a p-value threshold)\n\nExample: $0 0.001\nExample: $0 0.05 -qvalue";
($thresh) = args(1);

$tmpthresh = $thresh;
$tmpthresh =~ s/\./_/g;
$tmp_qval = '';
$allmotiffile = "input/all.meme";

$combinedmodelfile = "tmp-geneseqs-model.txt";
state("Combined background:");
run("cat", q(cat tmp/).$combinedmodelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
nl();

# start

startme("Loading list of RBPs with motifs from '$allmotiffile'");
open(MOTIFS, $allmotiffile) or die("Error: Couldn't open '$allmotiffile'");
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

# Cycle through all RBPs that have motifs and eCLIP data
state("Running MEME FIMO on RBPs with motifs and eCLIP data:");
$rbpquery = Query("SELECT symbol FROM `$intable` WHERE species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol ORDER BY symbol");
starttime();
$i = 0;
while (($symbol) = Fetch($rbpquery))
{
	# #DEBUG
	# next unless (($symbol eq 'HNRNPC') and ($celltype eq 'HEK293'));
	# #END DEBUG
	
	addme("total symbols", $symbol);

	# Set filenames
	$seqfile = "tmp-geneseqs-$symbol.txt";
	if (switch('qvalue'))
	{
		$tmp_qval = '_qval_'.$tmpthresh;
		$outfile = "../output/output-fimo-$type-$symbol-$tmpthresh-qvalue.txt";
		$outfilebg = "../output/output-fimobg-$type-$symbol-$tmpthresh-qvalue.txt";
		$outfilebgi = "../output/output-fimobgi-$type-$symbol-$tmpthresh-qvalue.txt";
		$tmpoutfile = "output/output-fimo-$type-$symbol-$tmpthresh-qvalue.txt";
		$tmpoutfilebg = "output/output-fimobg-$type-$symbol-$tmpthresh-qvalue.txt";
		$tmpoutfilebgi = "output/output-fimobgi-$type-$symbol-$tmpthresh-qvalue.txt";
	}
	else
	{
		$outfile = "../output/output-fimo-$type-$symbol-$tmpthresh.txt";
		$outfilebg = "../output/output-fimobg-$type-$symbol-$tmpthresh.txt";
		$outfilebgi = "../output/output-fimobgi-$type-$symbol-$tmpthresh.txt";
		$tmpoutfile = "output/output-fimo-$type-$symbol-$tmpthresh.txt";
		$tmpoutfilebg = "output/output-fimobg-$type-$symbol-$tmpthresh.txt";
		$tmpoutfilebgi = "output/output-fimobgi-$type-$symbol-$tmpthresh.txt";
	}
	$fimodir = "fimo$type\_$symbol$tmp_qval";
	$fimodirbg = "fimobg$type\_$symbol$tmp_qval";
	$fimodirbgi = "fimobgi$type\_$symbol$tmp_qval";
	$fimofile = "$fimodir/fimo.tsv";
	$fimofilebg = "$fimodirbg/fimo.tsv";
	$fimofilebgi = "$fimodirbgi/fimo.tsv";
	$motiffile = "tmp-motifs-all-$symbol.meme";
	# if (!switch('combinedmodel'))
	# {
		$modelfile = "tmp-geneseqs-$symbol-model.txt";
	# }
	
	# Skip if there were no significant peaks ("eclip_tom..." only)
	if (!-s 'tmp/'.$seqfile)
	{
		# addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol", $symbol);
		# next;
		die("Error: Sequence file '$seqfile' doesn't exist");
	}
	
	# Start
	$i++;
	state(" >> $i / ".Numrows($rbpquery)." >> $symbol", 1);
	print "   >> ";
	# Show ACGU fractions (the 'background model')
	run("cat", q(cat tmp/).$modelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
	print "\n";

	# Run MEME FIMO if necessary
	if (!-e $tmpoutfile)
	{
		cd("tmp", 1);
		if (!switch('qvalue'))
		{
			# p-value threshold
			state("     >> Running FIMO >> '$tmpoutfile'", 1);
			run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfile");
			state("     >> Running FIMO >> '$tmpoutfilebg'", 1);
			run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --bfile $combinedmodelfile --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebg");
			state("     >> Running FIMO >> '$tmpoutfilebgi'", 1);
			run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --bfile $modelfile --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebgi");
		}
		else
		{
			# q-value threshold
			# run("MEME FIMO (q-value threshold of $thresh)", "fimo --verbosity 1 --qv-thresh --thresh $thresh --norc --max-stored-scores 1000000000 --oc $fimodir $motiffile $seqfile 2>&1", 1);
			state("     >> Running FIMO >> '$tmpoutfile'", 1);
			$cmd = "fimo --verbosity 1 --qv-thresh --thresh $thresh --norc --max-stored-scores 1000000000 --oc $fimodir $motiffile $seqfile";
			run("MEME FIMO (q-value threshold of $thresh)", $cmd);
			if (-s $fimofile)
			{
				# Silly experiment because sometimes cat claims the file doesn't exist. Keeping the output directory now for diagnosis...and found out that FIMO simply doesn't produce an output directory in these cases. There are no results. All good
				run("Move output file to final location", "cat $fimofile | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfile", 1);
				run("Remove output directory", "rm -rf $fimodir", 1);
				print " >> ".chompme(`cat $outfile | wc -l`)." hits\n";
			}
			else
			{
				warn("Warning: FIMO doesn't seem to have produced any output from '$cmd' (expecting directory 'tmp/$fimodir')");
			}

			state("     >> Running FIMO >> '$tmpoutfilebg'", 1);
			$cmd = "fimo --verbosity 1 --qv-thresh --thresh $thresh --bfile $combinedmodelfile --norc --max-stored-scores 1000000000 --oc $fimodirbg $motiffile $seqfile";
			run("MEME FIMO (q-value threshold of $thresh)", $cmd);
			if (-s $fimofilebg)
			{
				# Silly experiment because sometimes cat claims the file doesn't exist. Keeping the output directory now for diagnosis...and found out that FIMO simply doesn't produce an output directory in these cases. There are no results. All good
				run("Move output file to final location", "cat $fimofilebg | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebg", 1);
				run("Remove output directory", "rm -rf $fimodirbg", 1);
				print " >> ".chompme(`cat $outfilebg | wc -l`)." hits\n";
			}
			else
			{
				warn("Warning: FIMO doesn't seem to have produced any output from '$cmd' (expecting directory 'tmp/$fimodirbg')");
			}

			state("     >> Running FIMO >> '$tmpoutfilebgi'", 1);
			$cmd = "fimo --verbosity 1 --qv-thresh --thresh $thresh --bfile $modelfile --norc --max-stored-scores 1000000000 --oc $fimodirbgi $motiffile $seqfile";
			run("MEME FIMO (q-value threshold of $thresh)", $cmd);
			if (-s $fimofilebgi)
			{
				# Silly experiment because sometimes cat claims the file doesn't exist. Keeping the output directory now for diagnosis...and found out that FIMO simply doesn't produce an output directory in these cases. There are no results. All good
				run("Move output file to final location", "cat $fimofilebgi | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebgi", 1);
				run("Remove output directory", "rm -rf $fimodirbgi", 1);
				print " >> ".chompme(`cat $outfilebgi | wc -l`)." hits\n";
			}
			else
			{
				warn("Warning: FIMO doesn't seem to have produced any output from '$cmd' (expecting directory 'tmp/$fimodirbgi')");
			}
		}
		# Keeping these now (to stop unnecessary re-runs that won't return any motifs)
		# run("Remove empty output file", "rm -f $outfile", 1) if (-z $outfile);
		# run("Remove empty output file", "rm -f $outfilebg", 1) if (-z $outfilebg);
		# run("Remove empty output file", "rm -f $outfilebgi", 1) if (-z $outfilebgi);

		cd("..", 1);

		addme("ran FIMO for symbol", $symbol);
	}
	else
	{
		state("     >> Skipping run: FIMO output file '$tmpoutfile' already exists", 1);
		addme("already run for symbol", $symbol);
	}
}
stoptime();

showmeall(1);

done();
