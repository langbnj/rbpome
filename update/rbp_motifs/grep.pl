#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$species = 'human';
$motiftable = 'rbp_motif_list';
$table = 'rbp_motifs';	# to be extended by _extend5_grep_$type below
$method = 'grep';

our $usage = "$0 [type: eclip_encode/eclip_tom/...] [motif source] [-extend5] [-firstonly]\n\n -extend5: Set extend5=1 in table rbp_motifs (to indicate that get_peakseqs.pl extended regions 5' by 50 nt, as done in the RNA Bind-N-Seq paper by Chris Burge's lab)\n-firstonly: Use only the 'best' motifs (and store them with a 'significant' q-value, as opposed to a significant p-value)\n\nExample: $0 eclip_encode attract\nExample: $0 eclip_tom dominguez -firstonly";
($type, $source) = args(2);

$extend5 = 0;
$tmpextend5 = '';
if (switch('extend5'))
{
	$table = 'rbp_motifs_extend5';
	$extend5 = 1;
	$tmpextend5 = '-extend5';
}
$table .= '_grep_'.$type;



$map = 'gene';
$intable = "clip_raw_$map";


# Clear table
if (switch('firstonly'))
{
	$query = Query("DELETE FROM `$table` WHERE method='$method' AND type='$type' AND source='$source' AND extend5='$extend5' AND qvalue=0");
	state("Cleared ".commify(Numrows($query))." '$method' '$type' '$source' extend5 '$extend5' entries with q-values from table '$table'");
}
else
{
	$query = Query("DELETE FROM `$table` WHERE method='$method' AND type='$type' AND source='$source' AND extend5='$extend5' AND qvalue=1");
	state("Cleared ".commify(Numrows($query))." '$method' '$type' '$source' extend5 '$extend5' entries without q-values from table '$table'");
}


# start

$query = Query("SELECT DISTINCT symbol FROM `$motiftable` WHERE species='$species' AND source='$source'");
startme("Loading list of RBPs with '$source' motifs from table '$motiftable'", 0, Numrows($query));
%motifs = ();
while (($symbol) = Fetch($query))
{
	$motifs{$symbol} = 1;
	stepme(100);
}
stopme();

# Get RBP accessions
startme("Getting RBP symbol-to-acc mapping from table '$intable'");
$query = Query("SELECT symbol, acc FROM `$intable` WHERE species='$species' AND map='$map' GROUP BY symbol, acc ORDER BY symbol, acc");
%acc = ();
while (($symbol, $acc) = Fetch($query))
{
	$acc{$symbol} = $acc;

	stepme(100);
}
stopme();

# Cycle through all RBPs that have motifs and eCLIP data
state("Running grep on RBPs with '$source' motifs and eCLIP data:");
$rbpquery = Query("SELECT symbol, celltype, rep FROM `$intable` WHERE type='$type' AND species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol, celltype, rep ORDER BY symbol, celltype, rep");
starttime();
$i = 0;
$inserted = 0;
while (($symbol, $celltype, $rep) = Fetch($rbpquery))
{
	addme("total symbols", $symbol);
	addme("total symbol|celltype|reps", "$symbol|$celltype|$rep");

	# Set filenames
	$seqfile = "tmp/tmp-peakseqs$tmpextend5-$method-$type-$source-$symbol-$celltype-$rep.txt";
	# $outfile = "output/output-grep-$type-$source-$symbol-$celltype-$rep.txt";
	
	# Skip if there were no significant peaks ("eclip_tom..." only)
	if (!-e $seqfile)
	{
		addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol|celltype|rep", "$symbol|$celltype|$rep");
		addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol", $symbol);
		next;
	}
	
	# Start
	$i++;

	# %hit = ();
	# if (!switch('nozero'))
	# {
	# 	# Get input peak regions (so I can add hit=0 rows below)
	# 	open(SEQ, $seqfile) or die("Error: Couldn't open sequence file '$seqfile'");
	# 	%hit = ();
	# 	fastabreak();
	# 	while (<SEQ>)
	# 	{
	# 		($title, $seq) = getfasta();
	# 		$seq = '';
	#
	# 		# >HNRNPA1|HepG2|chr20|49284067|49284173|+
	# 		# AUAAGGCAGUCCACAUAACAGAGAGAUAGGCCAGUAUCUUUCUGAGAGGCAGUCCUGACU
	# 		# GAAUUAGGUGAAUCUGUAAACCCUUGGGGUAGGACUCUCCAUUGGUA
	# 		# >HNRNPA1|HepG2|chr20|43459617|43459692|+
	#
	# 		@a = split(/\|/, $title);
	# 		die("Error: Expected 7 fields in title '$title' in '$seqfile'") if (scalar(@a) != 7);
	#
	# 		$chr = $a[3];
	# 		$start = $a[4];
	# 		$stop = $a[5];
	# 		$strand = $a[6];
	#
	# 		$hit{"$chr|$start|$stop|$strand"} = 0;
	# 	}
	# 	normalbreak();
	# 	close(SEQ);
	# }

	# Run grep
	startme(" >> $i / ".Numrows($rbpquery)." >> $symbol >> $celltype >> rep$rep", 1);
	open(SEQ, $seqfile) or die("Error: Couldn't open '$seqfile'");
	fastabreak();
	while (<SEQ>)
	{
		($title, $seq) = getfasta();
		
		@a = split(/\|/, $title);
		die("Error: Expected 7 fields in title '$title' in '$seqfile'") if (scalar(@a) != 7);

		$chr = $a[3];
		$start = $a[4];
		$stop = $a[5];
		$strand = $a[6];

		# Grep for motifs within the region
		if (!switch('firstonly'))
		{
			$motifquery = Query("SELECT DISTINCT motif FROM rbp_motif_list WHERE symbol='$symbol' AND species='$species' AND source='$source'");
		}
		else
		{
			$motifquery = Query("SELECT DISTINCT motif FROM rbp_motif_list WHERE symbol='$symbol' AND species='$species' AND source='$source' AND first=1");
		}
		while ($motif = Fetch($motifquery))
		{
			# $hits = 0;
			while ($seq =~ m/$motif/g)
			{
			# 	$hits++;
			# }
			# if ($hits > 0)
			# {
				
				$tmp_motifstart = (pos($seq) - length($motif)) + 1; # e.g. 1
				$tmp_motifstop = pos($seq); # e.g. 6 (pos() is already one position after the end of the occurrence)
				
				# state("\n\nMOTIF '$motif' IN $chr:$start-$stop ($strand strand) SEQ '$seq':\n\nTMP $tmp_motifstart..$tmp_motifstop\n\n\n");
				
				# Assign motif genomic coordinates
				if ($strand eq '+')
				{
					$motifstart = $start + ($tmp_motifstart - 1);
					$motifstop = $start + ($tmp_motifstop - 1);
				}
				elsif ($strand eq '-')
				{
					$motifstart = $stop - ($tmp_motifstop - 1);
					$motifstop = $stop - ($tmp_motifstart - 1);
				}
				else { die; }

				# state("\n\nMOTIF POSITION $chr:$motifstart-$motifstop\n\n\n");
				# exit;

				die("Error: No UniProt accession for symbol '$symbol'") if (!exists($acc{$symbol}));
				$acc = $acc{$symbol};
				
				if (!switch('firstonly'))
				{
					# Default: significant p-value
					$pvalue = '0';
					$qvalue = '1';
					$psig = 1;
					$qsig = 0;
				}
				else
				{
					# -firstonly: significant q-value
					$pvalue = '1';
					$qvalue = '0';
					$psig = 0;
					$qsig = 1;
				}
				
				# Insert hit into rbp_motifs_...
				# $motifid = 'grep_'.$motif;
				# $motifid = $motif;
				# $s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', celltype='$celltype', rep='$rep', method='$method', type='$type', chr='$chr', start='$start', stop='$stop', strand='$strand', source='$source', motif='$motif', hit='$hits', pvalue='$pvalue', qvalue='$qvalue', psig='$psig', qsig='$qsig', extend5='$extend5'";
				$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', celltype='$celltype', rep='$rep', method='$method', type='$type', chr='$chr', start='$start', stop='$stop', strand='$strand', source='$source', motif='$motif', motifstart='$motifstart', motifstop='$motifstop', pvalue='$pvalue', qvalue='$qvalue', psig='$psig', qsig='$qsig', extend5='$extend5'";
				$s =~ s/=''/=NULL/g;
				Query($s);
				$inserted++;
				
				addme("found motifs", $motif);
				addme("found symbol|motifs", "$symbol|$motif");
				addme("found symbol|celltype|rep|motif|chr|start|stop|strands", "$symbol|$celltype|$rep|$motif|$chr|$start|$stop|$strand");
				
				# # Store hit regions
				# if (!switch('nozero'))
				# {
				# 	$hit{"$chr|$start|$stop|$strand"} = 1;
				# }
			}
		}
		
		stepme(1000, 1);
	}
	normalbreak();
	close(SEQ);
	stopme(1);

	addme("ran grep for symbol", $symbol);
	addme("ran grep for symbol|celltype|rep", "$symbol|$celltype|$rep");
	
	# if (!switch('nozero'))
	# {
	# 	startme("   >> adding hit=0 rows", 1);
	# 	$sites_hit = 0;
	# 	foreach $site (keys(%hit))
	# 	{
	# 		if ($hit{$site} == 0)
	# 		{
	# 			@a = split(/\|/, $site);
	# 			die("Error: Expected 4 fields in '$site'") if (scalar(@a) != 4);
	# 			($chr, $start, $stop, $strand) = @a;
	#
	# 			if (!switch('firstonly'))
	# 			{
	# 				# Default: non-NULL p-value
	# 				$pvalue = '1';
	# 				$qvalue = '';
	# 			}
	# 			else
	# 			{
	# 				# -firstonly: non-NULL q-value
	# 				$pvalue = '';
	# 				$qvalue = '1';
	# 			}
	#
	# 			$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', celltype='$celltype', rep='$rep', method='$method', type='$type', chr='$chr', start='$start', stop='$stop', strand='$strand', source='$source', hit=0, pvalue='$pvalue', qvalue='$qvalue'";
	# 			$s =~ s/=''/=NULL/g;
	# 			Query($s);
	# 			$inserted++;
	# 			stepme(1000, 1);
	# 		}
	# 		else
	# 		{
	# 			$sites_hit++;
	# 		}
	# 	}
	# 	stopme(1);
	# 	state("   >> sites with hits >> ".commify($sites_hit), 1);
	# 	state("   >> total sites >> ".commify(scalar(keys(%hit))), 1);
	# 	state("     >> fraction >> ".sprintf('%.1f', ($sites_hit / scalar(keys(%hit))) * 100)."%", 1);
	# }
}
stoptime();

state("Rows inserted: ".commify($inserted));
showmeall(1);

Optimize($table);

done();
