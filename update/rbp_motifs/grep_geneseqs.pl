#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$species = 'human';
$motiftable = 'rbp_motif_list';
$method = 'grep';

$type = 'geneseqs';
$source = 'geneseqs';
$table = 'rbp_motifs_grep_'.$type;
$map = 'gene';
$intable = "clip_raw_$map";

$extend5 = 0;

our $usage = "$0 [-firstonly]\n\n-firstonly: Use only the 'best' motifs (and store them with a 'significant' q-value, as opposed to a significant p-value)\n\nExample: $0";
args(0);



# Clear table
if (switch('firstonly'))
{
	$query = Query("DELETE FROM `$table` WHERE method='$method' AND type='$type' AND qvalue=0");
	state("Cleared ".commify(Numrows($query))." '$method' '$type' entries with q-values from table '$table'");
}
else
{
	$query = Query("DELETE FROM `$table` WHERE method='$method' AND type='$type' AND qvalue=1");
	state("Cleared ".commify(Numrows($query))." '$method' '$type' entries without q-values from table '$table'");
}


# start

$query = Query("SELECT DISTINCT symbol FROM `$motiftable` WHERE species='$species'");
startme("Loading list of RBPs with motifs from table '$motiftable'", 0, Numrows($query));
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
state("Running grep on RBPs with motifs and eCLIP data:");
$rbpquery = Query("SELECT symbol FROM `$intable` WHERE species='$species' AND map='$map' AND symbol IN ('".join("', '", keys(%motifs))."') GROUP BY symbol ORDER BY symbol");
starttime();
$i = 0;
$inserted = 0;
while (($symbol) = Fetch($rbpquery))
{
	addme("total symbols", $symbol);

	# Set filenames
	$seqfile = "tmp/tmp-geneseqs-$symbol.txt";
	
	# Skip if there were no significant peaks ("eclip_tom..." only)
	if (!-e $seqfile)
	{
		# addme("skipped because there were no significant peaks (empty/nonexistent peak file) for symbol", $symbol);
		# next;
		die("Error: Sequence file '$seqfile' doesn't exist");
	}
	
	# Start
	$i++;

	# Run grep
	startme(" >> $i / ".Numrows($rbpquery)." >> $symbol", 1);
	open(SEQ, $seqfile) or die("Error: Couldn't open '$seqfile'");
	fastabreak();
	while (<SEQ>)
	{
		($title, $seq) = getfasta();
		
		@a = split(/\|/, $title);
		die("Error: Expected 5 fields in title '$title' in '$seqfile'") if (scalar(@a) != 5);
		
		$chr = $a[1];
		$start = $a[2];
		$stop = $a[3];
		$strand = $a[4];

		#DEBUG
		# next if ($strand eq '+');
		#END DEBUG
		
		# Grep for motifs within the region
		if (!switch('firstonly'))
		{
			$motifquery = Query("SELECT DISTINCT source, motif FROM rbp_motif_list WHERE symbol='$symbol' AND species='$species'");
		}
		else
		{
			$motifquery = Query("SELECT DISTINCT source, motif FROM rbp_motif_list WHERE symbol='$symbol' AND species='$species' AND first=1");
		}
		while (($source, $motif) = Fetch($motifquery))
		{
			# $hits = 0;
			foreach $pos (positions($motif, $seq))
			# while ($seq =~ m/$motif/g)
			# {
			# 	$hits++;
			# }
			# if ($hits > 0)
			{
				if ($strand eq '+')
				{
					$thisstart = $start + ($pos - 1);
					$thisstop = $thisstart + (length($motif) - 1);
				}
				else
				{
					$thisstart = $stop - ($pos - 1);
					$thisstop = $thisstart - (length($motif) - 1);
				}
				
				# print " >> Title: $title\n";
				# print "   >> Gene Pos: $chr:$start-$stop:$strand\n";
				# print "     >> Motif Pos: $pos\n";
				# print "       >> Motif: $thisstart|$thisstop\n";
				# print "       >> Motif Seq: $motif\n";
				# print "       >> Gene Seq: $seq\n";
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
				# $s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', method='$method', type='$type', chr='$chr', start='$thisstart', stop='$thisstop', strand='$strand', source='$source', motif='$motif', hit='$hits', pvalue='$pvalue', qvalue='$qvalue', psig='$psig', qsig='$qsig', extend5='$extend5'";
				$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', method='$method', type='$type', chr='$chr', start='$thisstart', stop='$thisstop', strand='$strand', source='$source', motif='$motif', pvalue='$pvalue', qvalue='$qvalue', psig='$psig', qsig='$qsig', extend5='$extend5'";
				$s =~ s/=''/=NULL/g;
				Query($s);
				$inserted++;
				
				addme("found motifs", $motif);
				addme("found symbol|motifs", "$symbol|$motif");
				addme("found symbol|motif|chr|start|stop|strands", "$symbol|$motif|$chr|$start|$stop|$strand");
				
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
}
stoptime();

state("Rows inserted: ".commify($inserted));
showmeall(1);

Optimize($table);

done();
