#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table] [eclip type] [screen score type] [mindist threshold for cobinding] [motif search method] [motif 5' extension 0/1] [-qvalue]\n\n -qvalue: Use motif hits with significant q-values, rather than p-values (default).\n\nExample: $0 rbpome eclip_encode_12 SUM_IS 54 fimo 1\n";
($table, $type, $scoretype, $mindist_threshold, $motif_method, $motif_extend5) = args(6);

$tmpextend5 = '';
$tmpextend5 = '_extend5' if ($motif_extend5 == 1);

$tmpqvalue = '';
$tmpqvalue = '-qvalue' if (switch('qvalue'));

$tmpsig = 'psig=1';
$tmpsig = 'qsig=1' if (switch('qvalue'));


if (!switch('debug'))
{
	$outfile = "output-motif_analysis-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";
}
else
{
	$outfile = "output-motif_analysis-rbpome-eclip_encode_12-SUM_IS-mindist_threshold54-fimo-extend5_1-debug.txt";
}

die("Error: '$outfile' didn't exist yet (please run motif_analysis.pl first)") if (!-s $outfile);

# First, read existing output to get the number of random sites needed for each RBP pair
startme("Getting number of close binding cases from '$outfile' (in order to randomly pick the same number of distant cases as a control)");
open(IN, $outfile) or die("Error: Couldn't read '$outfile' (please run motif_analysis.pl first)");
<IN>;	# Skip header
%closecount = ();
while (<IN>)
{
	chomp;
	
	# symbol1	symbol2	bindtype	mindist	motifcount1	motifcount2	motifs1	motifs2
	# EWSR1	HNRNPK	close	10	5	7	dominguez|EWSR1_1	attract|HNRNPK_2,attract|HNRNPK_4,attract|HNRNPK_6,attract|HNRNPK_9,cisbp|HNRNPK_1,rbpmap|HNRNPK_1,rnacompete|HNRNPK_1
	# EWSR1	HNRNPK	close	17	5	6	dominguez|EWSR1_1	attract|HNRNPK_1,attract|HNRNPK_13,attract|HNRNPK_2,attract|HNRNPK_21,attract|HNRNPK_5,attract|HNRNPK_8
	# EWSR1	HNRNPK	close	11	1	6	dominguez|EWSR1_1	attract|HNRNPK_19,attract|HNRNPK_2,attract|HNRNPK_6,cisbp|HNRNPK_1,rbpmap|HNRNPK_1,rnacompete|HNRNPK_1
	# EWSR1	HNRNPK	close	52	2	16	dominguez|EWSR1_1,dominguez|EWSR1_2	attract|HNRNPK_10,attract|HNRNPK_12,attract|HNRNPK_14,attract|HNRNPK_15,attract|HNRNPK_16,attract|HNRNPK_17,attract|HNRNPK_20,attract|HNRNPK_22,attract|HNRNPK_3,attract|HNRNPK_5,attract|HNRNPK_6,attract|HNRNPK_7,cisbp|HNRNPK_1,rbpmap|HNRNPK_1,rnacompete|HNRNPK_1
	
	@a = split(/\t/);
	
	$symbol1 = $a[0];
	$symbol2 = $a[1];
	$bindtype = $a[2];
	
	next if ($bindtype ne 'close');
	
	if (exists($closecount{$symbol1.'|'.$symbol2}))
	{
		$closecount{$symbol1.'|'.$symbol2}++;
	}
	else
	{
		$closecount{$symbol1.'|'.$symbol2} = 1;
	}
	
	addme("total pairs with close binding sites", $symbol1.'|'.$symbol2);
	
	stepme(10000);
}
close(IN);
stopme();





if (!switch('debug'))
{
	open(OUT, ">>$outfile") or die("\nError: Couldn't append to '$outfile'\n\n");
	# print OUT "symbol1\tsymbol2\tbindtype\tmindist\tmotifcount1\tmotifcount2\tmotifs1\tmotifs2\n";
}



# start

state("Getting 'output/output-mindist_motifs_plot-*-data.txt' files and analysing motif types for these peaks...");
starttime();

$totalfiles = chompme(`ls -1 output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-*-data.txt | wc -l`);

open(LS, "ls -1 output/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-*-data.txt|") or die("Error: Couldn't ls");
$i = 0;
while (<LS>)
{
	$i++;
	chomp;
	
	$infile = $_;
	$infile =~ /^output\/output-mindist_motifs_plot-$table-$type-$scoretype-mindist_threshold$mindist_threshold-(\S+)-(\S+)-data\.txt/ or die("Error: Couldn't parse file name '$infile'");
	$symbol1 = $1;
	$symbol2 = $2;
	
	# $tmpfile = $infile;
	# $tmpfile =~ s/-data\.txt$/-data-tmp.txt/;
	#
	# # Process input file so it no longer contains gene information (drop the ensgv column) (unless it exists already)
	# if (!-s $tmpfile or switch('debug'))
	# {
	# 	run("Process input file", "cat $infile | cut -f1-2,4- | natsort | uniq > $tmpfile");
	# }
	
	startme(" >> $i / $totalfiles >> $symbol1|$symbol2", 1, chompme(`cat $infile | wc -l`));
	print " >> $symbol1|$symbol2 ($symbol2 near $symbol1)\n" if (switch('debug'));

	
	
	# Get ENSGVs for symbol1
	$query = Query("SELECT DISTINCT ensgv FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol1'");
	print "   >> Getting ENSGVs bound by $symbol1\n" if (switch('debug'));
	@ensgvs1 = ();
	while (($ensgv) = Fetch($query))
	{
		push(@ensgvs1, $ensgv);
	}
	
	# Get ENSGVs for symbol2
	$query = Query("SELECT DISTINCT ensgv FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol2'");
	print "   >> Getting ENSGVs bound by $symbol2\n" if (switch('debug'));
	@ensgvs2 = ();
	while (($ensgv) = Fetch($query))
	{
		push(@ensgvs2, $ensgv);
	}

	# Get exclusive genes
	@ensgvs_1_only = lonly(\@ensgvs1, \@ensgvs2);
	@ensgvs_2_only = ronly(\@ensgvs1, \@ensgvs2);
	
	print "     >> Exclusively bound by $symbol1: ".scalar(@ensgvs_1_only)." genes\n" if (switch('debug'));
	print "     >> Exclusively bound by $symbol2: ".scalar(@ensgvs_2_only)." genes\n" if (switch('debug'));
	
	# @exclusives = (@ensgvs_1_only, @ensgvs_2_only);
	
	# Get all solo peaks into arrays
	$query = Query("SELECT chr, start, stop, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol1' AND ensgv IN ('".join("', '", @ensgvs_1_only)."')");
	@solopeaks1 = ();
	while (($chr, $start, $stop, $strand) = Fetch($query))
	{
		push(@solopeaks1, "$chr|$start|$stop|$strand");
	}
	
	$query = Query("SELECT chr, start, stop, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol2' AND ensgv IN ('".join("', '", @ensgvs_2_only)."')");
	@solopeaks2 = ();
	while (($chr, $start, $stop, $strand) = Fetch($query))
	{
		push(@solopeaks2, "$chr|$start|$stop|$strand");
	}
	
	
	# if ((Numrows($mainquery1) < $closecount{$symbol1.'|'.$symbol2}) or (Numrows($mainquery1) < $closecount{$symbol1.'|'.$symbol2}))
	# {
	# 	die("Error: Not enough solo-bound genes for $symbol1|$symbol2 (needed ".$closecount{$symbol1.'|'.$symbol2}.", but got $symbol1: ".Numrows($mainquery1)." and $symbol2: ".Numrows($mainquery2).")");
	# }
	
	
	$picked = 0;
	next if (!exists($closecount{$symbol1.'|'.$symbol2}));
	while ($picked < $closecount{$symbol1.'|'.$symbol2})
	{
		# # Pick a random line (@file doesn't contain the header, so starting from 0 is fine)
		# $_ =  @file[int(rand(scalar(@file)))];
		# $ensgv1 = @ensgvs_1_only[int(rand(scalar(@ensgvs_1_only)))];
		# $ensgv2 = @ensgvs_2_only[int(rand(scalar(@ensgvs_2_only)))];
		
		# Pick random peaks
		$solopeak1 = @solopeaks1[int(rand(scalar(@solopeaks1)))];
		$solopeak2 = @solopeaks2[int(rand(scalar(@solopeaks2)))];
		
		($chr1, $start1, $stop1, $strand1) = split(/\|/, $solopeak1);
		($chr2, $start2, $stop2, $strand2) = split(/\|/, $solopeak2);
		
		# Get a random peak for protein 1
		# $query = Query("SELECT chr, start, stop, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol1' AND ensgv='$ensgv1' ORDER BY RAND() LIMIT 1");
		# ($chr1, $start1, $stop1, $strand1) = FetchOne($query);
		# ($chr1, $start1, $stop1, $strand1) = Fetch($mainquery1);
		if ($strand1 eq '+')
		{
			$peak1 = $start1;
		}
		else
		{
			$peak1 = $stop1;
		}

		# Get a random peak for protein 2
		# $query = Query("SELECT chr, start, stop, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol2' AND ensgv='$ensgv2' ORDER BY RAND() LIMIT 1");
		# ($chr2, $start2, $stop2, $strand2) = FetchOne($query);
		# ($chr2, $start2, $stop2, $strand2) = Fetch($mainquery2);
		if ($strand2 eq '+')
		{
			$peak2 = $start2;
		}
		else
		{
			$peak2 = $stop2;
		}
		
		# chomp;

		# symbol1	symbol2	ensgv	peak1	peak2	mindist	bindtype	motif1	motif2
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641752	8	close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641806	46	close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7642477	717	distant	1	1
		# FMR1	FXR2	ENSG00000076944.15	7641760	7644614	2854	distant	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641752	8	closest_is_close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7642434	7641752	682	distant	1	0
		
		# @a = split(/\t/);
		
		# # Verify symbols (should never happen)
		# die if ($a[0] ne $symbol1);
		# die if ($a[1] ne $symbol2);
		
		# $ensgv = $a[2];
		# $peak1 = $a[3];
		# $peak2 = $a[4];
		# $mindist = $a[5];
		# die("Error: Unexpected distance result '$mindist' for peaks '$peak1' and '$peak2'") if ($mindist ne abs($peak2 - $peak1));
		
		$mindist = 'NA';
		
		# $bindtype = $a[6];
		# $peak1_has_motifs = $a[7];
		# $peak2_has_motifs = $a[8];
		
		
		$bindtype = 'solo';
		


		# # Get ENSG chromosome and strand for symbol1 peak
		# # $query = Query("SELECT chr, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol1' AND ensgv='$ensgv' AND ((start='$peak1' AND strand='+') OR (stop='$peak1' AND strand='-'))");
		# $query = Query("SELECT chr, strand FROM gencode_gff3_gene WHERE species='human' AND ensgv='$ensgv'");
		# ($chr, $strand) = FetchOne($query);
		
		
		
		
		# Apply 50 nt 5' extension
		if ($motif_extend5 == 1)
		{
			# Peak 1
			if ($strand1 eq '+')
			{
				$peak1 -= 50;
			}
			elsif ($strand1 eq '-')
			{
				$peak1 += 50;
			}

			# Peak 2
			if ($strand2 eq '+')
			{
				$peak2 -= 50;
			}
			elsif ($strand2 eq '-')
			{
				$peak2 += 50;
			}
		}

		
		
		# print "   >> ENSGV $ensgv >> PEAK1 $peak1 >> PEAK2 $peak2\n" if (switch('debug'));
		print "   >> $bindtype >> mindist $mindist nt\n" if (switch('debug'));
		print "     >> $symbol1 motifs at $peak1:\n" if (switch('debug'));
		
		
		# Get motifs for symbol1 peak
		# Not caring about cell type or replicate here (using DISTINCT)
		# # Not using chr and strand here to query, hopefully no ambiguity
		# $query = Query("SELECT motif, COUNT(id) AS hits, COUNT(DISTINCT chr) AS chrs, COUNT(DISTINCT strand) AS strands FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol1' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND ((start='$peak1' AND strand='+') OR (stop='$peak1' AND strand='-')) GROUP BY motif");
		# our $superloudmysql = 1;
		$q = "SELECT source, motif, COUNT(id) AS hits FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol1' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND chr='$chr1' AND strand='$strand1' AND ((start='$peak1' AND strand='+') OR (stop='$peak1' AND strand='-')) GROUP BY source, motif";
		$query = Query($q);
		@motifs1 = ();
		$motifcount1 = 0;
		# while (($motif, $hits, $chrs, $strands) = Fetch($query))
		# if ((Numrows($query) == 0) and ($peak1_has_motifs == 1))
		# {
		# 	die("\n\nError: No peaks found for query, even though there should be some:\n\n$q\n\n");
		# }
		while (($source, $motif, $hits) = Fetch($query))
		{
			# # Verification that just querying by the peak start/stop is a unique identifier
			# die("Error: Multiple chromosomes found ($chrs) for symbol '$symbol1' and peak '$peak1'") if ($chrs != 1);
			# die("Error: Multiple strand orientations found for symbol '$symbol1' and peak '$peak1'") if ($strands != 1);
			
			# # Shouldn't happen
			# die if ($peak1_has_motifs == 0);
			
			push(@motifs1, "$source|$motif");
			$motifcount1 += $hits;

			print "       >> $source >> $motif ($hits)\n" if (switch('debug'));
		}
		
		



		print "     >> $symbol2 motifs at $peak2:\n" if (switch('debug'));
		
		# Get motifs for symbol2 peak
		# Not caring about cell type or replicate here (using DISTINCT)
		# # Not using chr and strand here to query, hopefully no ambiguity
		# $query = Query("SELECT motif, COUNT(id) AS hits, COUNT(DISTINCT chr) AS chrs, COUNT(DISTINCT strand) AS strands FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol2' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND ((start='$peak2' AND strand='+') OR (stop='$peak2' AND strand='-')) GROUP BY motif");
		# our $superloudmysql = 1;
		$q = "SELECT source, motif, COUNT(id) AS hits FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol2' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND chr='$chr2' AND strand='$strand2' AND ((start='$peak2' AND strand='+') OR (stop='$peak2' AND strand='-')) GROUP BY source, motif";
		$query = Query($q);
		@motifs2 = ();
		$motifcount2 = 0;
		# while (($motif, $hits, $chrs, $strands) = Fetch($query))
		# if ((Numrows($query) == 0) and ($peak2_has_motifs == 1))
		# {
		# 	die("\n\nError: No peaks found for query, even though there should be some:\n\n$q\n\n");
		# }
		while (($source, $motif, $hits) = Fetch($query))
		{
			# # Verification that just querying by the peak start/stop is a unique identifier
			# die("Error: Multiple chromosomes found ($chrs) for symbol '$symbol2' and peak '$peak2'") if ($chrs != 1);
			# die("Error: Multiple strand orientations found for symbol '$symbol2' and peak '$peak2'") if ($strands != 1);
			
			# # Shouldn't happen
			# die if ($peak2_has_motifs == 0);
			
			push(@motifs2, "$source|$motif");
			$motifcount2 += $hits;

			print "       >> $source >> $motif ($hits)\n" if (switch('debug'));
		}
		
		
		
		
		# Print
		if (!switch('debug'))
		{
			print OUT "$symbol1\t$symbol2\t$bindtype\t$mindist\t$motifcount1\t$motifcount2\t".join(",", @motifs1)."\t".join(",", @motifs2)."\n";
		}
		
		$picked++;
		
		stepme(1, 1);
	}
	stopme(1);
	
	
	
	
	stepme(100);
}
stopme();
stoptime();

showmeall(1);

state("Wrote to '$outfile'");

done();
