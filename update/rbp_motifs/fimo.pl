#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$species = 'human';
$method = 'fimo';

our $usage = "$0 [type: eclip_encode/eclip_tom/...] [motif source] [p/q-value threshold] [-qvalue]\n\n -qvalue: Use a q-value threshold (instead of a p-value threshold)\n\nExample: $0 eclip_encode attract 0.001\nExample: $0 eclip_tom attract 0.05 -qvalue";
# our $usage = "$0 [type: eclip_encode/eclip_tom/...] [motif source] [p/q-value threshold] [-bgs] [-qvalue]\n\n -bgs: Run 'fimobg' and 'fimobgi' as well (with background models).\n -qvalue: Use a q-value threshold (instead of a p-value threshold)\n\nExample: $0 eclip_encode attract 0.001\nExample: $0 eclip_tom attract 0.05 -qvalue";
($type, $source, $thresh) = args(3);

$tmpextend = '';
$tmpextend_2 = '';
if (switch('extend'))
{
	$tmpextend = '-extend';
	$tmpextend_2 = '_extend';
}

$tmpthresh = $thresh;
$tmpthresh =~ s/\./_/g;
$tmp_qval = '';
$allmotiffile = "input/$source.meme";

$map = 'gene';
$intable = "clip_raw_$map";

# # $modelfile = '../input/gencode.human.v27.transcripts.model.txt';
# # $modelfile = '../input/tmp-peakseqs-eclip_encode-attract-model.txt';
# # $modelfile = '../input/tmp-peakseqs-eclip_tom-attract-model.txt';
# # if (switch('combinedmodel'))
# # {
# 	$combinedmodelfile = "tmp-peakseqs$tmpextend-$method-$type-$source-combined-model.txt";
# 	state("Combined background:");
# 	run("cat", q(cat tmp/).$combinedmodelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
# 	nl();
# # }
# # ~/update/rbp_motifs/bin/meme-5.0.5/src >> fasta-get-markov -norc ~/update/gencode_gff3/input/gencode.human.v27.transcripts.fa ../../../tmp2_HUGE_DELETE_ME_SOON/gencode.human.v27.transcripts.model.txt
# # ~/update/rbp_motifs/bin/meme-5.0.5/src >> fasta-get-markov -norc ../../../tmp2_HUGE_DELETE_ME_SOON/tmp-peakseqs-eclip_encode-attract.txt ../../../tmp2_HUGE_DELETE_ME_SOON/tmp-peakseqs-eclip_encode-attract-model.txt
# # ~/update/rbp_motifs/bin/meme-5.0.5/src >> fasta-get-markov -norc ../../../tmp2_HUGE_DELETE_ME_SOON/tmp-peakseqs-eclip_tom-attract.txt ../../../tmp2_HUGE_DELETE_ME_SOON/tmp-peakseqs-eclip_tom-attract-model.txt

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
state("Running MEME FIMO on RBPs with '$source' motifs and eCLIP data:");
$rbpquery = Query("SELECT symbol, celltype, rep FROM `$intable` WHERE type='$type' AND species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol, celltype, rep ORDER BY symbol, celltype, rep");
starttime();
$i = 0;
while (($symbol, $celltype, $rep) = Fetch($rbpquery))
{
	#DEBUG
	# next unless (($symbol eq 'HNRNPC') and ($celltype eq 'HEK293'));
	#END DEBUG
	
	addme("total symbols", $symbol);
	addme("total symbol|celltype|reps", "$symbol|$celltype|$rep");

	# Set filenames
	$seqfile = "tmp-peakseqs$tmpextend-$method-$type-$source-$symbol-$celltype-$rep.txt";
	if (switch('qvalue'))
	{
		$tmp_qval = '_qval_'.$tmpthresh;
		$outfile = "../output/output$tmpextend-fimo-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
		# $outfilebg = "../output/output$tmpextend-fimobg-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
		# $outfilebgi = "../output/output$tmpextend-fimobgi-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
		$tmpoutfile = "output/output$tmpextend-fimo-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
		# $tmpoutfilebg = "output/output$tmpextend-fimobg-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
		# $tmpoutfilebgi = "output/output$tmpextend-fimobgi-$type-$source-$symbol-$celltype-$rep-$tmpthresh-qvalue.txt";
	}
	else
	{
		$outfile = "../output/output$tmpextend-fimo-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
		# $outfilebg = "../output/output$tmpextend-fimobg-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
		# $outfilebgi = "../output/output$tmpextend-fimobgi-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
		$tmpoutfile = "output/output$tmpextend-fimo-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
		# $tmpoutfilebg = "output/output$tmpextend-fimobg-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
		# $tmpoutfilebgi = "output/output$tmpextend-fimobgi-$type-$source-$symbol-$celltype-$rep-$tmpthresh.txt";
	}
	$fimodir = "fimo$tmpextend_2\_$type\_$source\_$symbol\_$celltype\_$rep$tmp_qval";
	# $fimodirbg = "fimobg$tmpextend_2\_$type\_$source\_$symbol\_$celltype\_$rep$tmp_qval";
	# $fimodirbgi = "fimobgi$tmpextend_2\_$type\_$source\_$symbol\_$celltype\_$rep$tmp_qval";
	$fimofile = "$fimodir/fimo.tsv";
	# $fimofilebg = "$fimodirbg/fimo.tsv";
	# $fimofilebgi = "$fimodirbgi/fimo.tsv";
	$motiffile = "tmp-motifs-$source-$symbol.meme";
	# # if (!switch('combinedmodel'))
	# # {
	# 	$modelfile = "tmp-peakseqs$tmpextend-$method-$type-$source-$symbol-$celltype-$rep-model.txt";
	# # }
	
	# Skip if there were no significant peaks ("eclip_tom..." only)
	if (!-s 'tmp/'.$seqfile)
	{
		addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol|celltype|rep", "$symbol|$celltype|$rep");
		addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol", $symbol);
		next;
	}
	
	# Start
	$i++;
	state(" >> $i / ".Numrows($rbpquery)." >> $symbol >> $celltype >> rep$rep", 1);
	# print "   >> ";
	# # Show ACGU fractions (the 'background model')
	# run("cat", q(cat tmp/).$modelfile.q( | perl -ne 'chomp; next if /^#/; ($a, $b) = split(/ /); $b+=0; print "$a $b "'), 1);
	# print "\n";

	# Run MEME FIMO if necessary
	if (!-e $tmpoutfile)
	{
		cd("tmp", 1);
		if (!switch('qvalue'))
		{
			# p-value threshold
			state("     >> Running FIMO >> '$tmpoutfile'", 1);
			run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfile");
			# state("     >> Running FIMO >> '$tmpoutfilebg'", 1);
			# run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --bfile $combinedmodelfile --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebg");
			# state("     >> Running FIMO >> '$tmpoutfilebgi'", 1);
			# run("MEME FIMO (p-value threshold of $thresh)", "fimo --verbosity 1 --thresh $thresh --bfile $modelfile --text --skip-matched-sequence --norc $motiffile $seqfile 2> /dev/null | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebgi");
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

		# 	state("     >> Running FIMO >> '$tmpoutfilebg'", 1);
		# 	$cmd = "fimo --verbosity 1 --qv-thresh --thresh $thresh --bfile $combinedmodelfile --norc --max-stored-scores 1000000000 --oc $fimodirbg $motiffile $seqfile";
		# 	run("MEME FIMO (q-value threshold of $thresh)", $cmd);
		# 	if (-s $fimofilebg)
		# 	{
		# 		# Silly experiment because sometimes cat claims the file doesn't exist. Keeping the output directory now for diagnosis...and found out that FIMO simply doesn't produce an output directory in these cases. There are no results. All good
		# 		run("Move output file to final location", "cat $fimofilebg | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebg", 1);
		# 		run("Remove output directory", "rm -rf $fimodirbg", 1);
		# 		print " >> ".chompme(`cat $outfilebg | wc -l`)." hits\n";
		# 	}
		# 	else
		# 	{
		# 		warn("Warning: FIMO doesn't seem to have produced any output from '$cmd' (expecting directory 'tmp/$fimodirbg')");
		# 	}
		#
		# 	state("     >> Running FIMO >> '$tmpoutfilebgi'", 1);
		# 	$cmd = "fimo --verbosity 1 --qv-thresh --thresh $thresh --bfile $modelfile --norc --max-stored-scores 1000000000 --oc $fimodirbgi $motiffile $seqfile";
		# 	run("MEME FIMO (q-value threshold of $thresh)", $cmd);
		# 	if (-s $fimofilebgi)
		# 	{
		# 		# Silly experiment because sometimes cat claims the file doesn't exist. Keeping the output directory now for diagnosis...and found out that FIMO simply doesn't produce an output directory in these cases. There are no results. All good
		# 		run("Move output file to final location", "cat $fimofilebgi | perl -ne 'if (/^\(\\S+\)_\\d+\\t\\t\([^\\|]+\)\\|/) { print if (\$1 eq \$2); }' > $outfilebgi", 1);
		# 		run("Remove output directory", "rm -rf $fimodirbgi", 1);
		# 		print " >> ".chompme(`cat $outfilebgi | wc -l`)." hits\n";
		# 	}
		# 	else
		# 	{
		# 		warn("Warning: FIMO doesn't seem to have produced any output from '$cmd' (expecting directory 'tmp/$fimodirbgi')");
		# 	}
		}
		# Keeping these now (to stop unnecessary re-runs that won't return any motifs)
		# run("Remove empty output file", "rm -f $outfile", 1) if (-z $outfile);
		# run("Remove empty output file", "rm -f $outfilebg", 1) if (-z $outfilebg);
		# run("Remove empty output file", "rm -f $outfilebgi", 1) if (-z $outfilebgi);

		cd("..", 1);

		addme("ran FIMO for symbol", $symbol);
		addme("ran FIMO for symbol|celltype|rep", "$symbol|$celltype|$rep");
	}
	else
	{
		state("     >> Skipping run: FIMO output file '$tmpoutfile' already exists", 1);
		addme("already run for symbol", $symbol);
		addme("already run for symbol|celltype|rep", "$symbol|$celltype|$rep");
	}
}
stoptime();

showmeall(1);

done();
