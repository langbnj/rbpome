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


$outfile = "output-motif_analysis-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "symbol1\tsymbol2\tbindtype\tmindist\tmotifcount1\tmotifcount2\tmotifs1\tmotifs2\n";



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
	
	startme(" >> $i / $totalfiles >> Reading '$infile'", 1, chompme(`cat $infile | wc -l`));
	print " >> $symbol1|$symbol2 ($symbol2 near $symbol1)\n" if (switch('debug'));
	open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
	# startme(" >> $i / $totalfiles >> Reading '$tmpfile'", 1);
	# open(IN, $tmpfile) or die("\nError: Couldn't open '$tmpfile'\n\n");
	<IN>;	# Skip header
	while (<IN>)
	{
		chomp;

		# symbol1	symbol2	ensgv	peak1	peak2	mindist	bindtype	motif1	motif2
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641752	8	close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641806	46	close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7642477	717	distant	1	1
		# FMR1	FXR2	ENSG00000076944.15	7641760	7644614	2854	distant	1	0
		# FMR1	FXR2	ENSG00000076944.15	7641760	7641752	8	closest_is_close	1	0
		# FMR1	FXR2	ENSG00000076944.15	7642434	7641752	682	distant	1	0
		
		@a = split(/\t/);
		
		# Verify symbols (should never happen)
		die if ($a[0] ne $symbol1);
		die if ($a[1] ne $symbol2);
		
		$ensgv = $a[2];
		$peak1 = $a[3];
		$peak2 = $a[4];
		$mindist = $a[5];
		die("Error: Unexpected distance result '$mindist' for peaks '$peak1' and '$peak2'") if ($mindist ne abs($peak2 - $peak1));
		
		$bindtype = $a[6];
		$peak1_has_motifs = $a[7];
		$peak2_has_motifs = $a[8];
		
		
		# Skip uninteresting lines (these are duplicates of close & distant anyway)
		next if ($bindtype eq 'closest_is_close');
		next if ($bindtype eq 'closest_is_distant');
		
		
		# # For debugging
		# next if (($bindtype eq 'distant') and switch('debug'));
		# next if (($peak1_has_motifs != 1) and switch('debug'));
		# next if (($peak2_has_motifs != 1) and switch('debug'));

		# For speed
		next if ($bindtype eq 'distant');
		# next if ($peak1_has_motifs != 1);
		# next if ($peak2_has_motifs != 1);


		# Get ENSG chromosome and strand for symbol1 peak
		# $query = Query("SELECT chr, strand FROM clip_raw_gene WHERE type='$type' AND symbol='$symbol1' AND ensgv='$ensgv' AND ((start='$peak1' AND strand='+') OR (stop='$peak1' AND strand='-'))");
		$query = Query("SELECT chr, strand FROM gencode_gff3_gene WHERE species='human' AND ensgv='$ensgv'");
		($chr, $strand) = FetchOne($query);
		
		
		
		# Apply 50 nt 5' extension
		if ($motif_extend5 == 1)
		{
			if ($strand eq '+')
			{
				$peak1 -= 50;
				$peak2 -= 50;
			}
			elsif ($strand eq '-')
			{
				$peak1 += 50;
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
		$q = "SELECT source, motif, COUNT(id) AS hits FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol1' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND chr='$chr' AND strand='$strand' AND ((start='$peak1' AND strand='+') OR (stop='$peak1' AND strand='-')) GROUP BY source, motif";
		$query = Query($q);
		@motifs1 = ();
		$motifcount1 = 0;
		# while (($motif, $hits, $chrs, $strands) = Fetch($query))
		if ((Numrows($query) == 0) and ($peak1_has_motifs == 1))
		{
			die("\n\nError: No peaks found for query, even though there should be some:\n\n$q\n\n");
		}
		while (($source, $motif, $hits) = Fetch($query))
		{
			# # Verification that just querying by the peak start/stop is a unique identifier
			# die("Error: Multiple chromosomes found ($chrs) for symbol '$symbol1' and peak '$peak1'") if ($chrs != 1);
			# die("Error: Multiple strand orientations found for symbol '$symbol1' and peak '$peak1'") if ($strands != 1);
			
			# Shouldn't happen
			die("\n\nError: This query gives motifs for '$symbol1' peak1 '$peak1', but there shouldn't be any according to '$infile':\n\n$q\n\n") if ($peak1_has_motifs == 0);
			
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
		$q = "SELECT source, motif, COUNT(id) AS hits FROM rbp_motifs$tmpextend5\_$motif_method\_$type WHERE symbol='$symbol2' AND species='human' AND method='$motif_method' AND extend5='$motif_extend5' AND type='$type' AND $tmpsig AND chr='$chr' AND strand='$strand' AND ((start='$peak2' AND strand='+') OR (stop='$peak2' AND strand='-')) GROUP BY source, motif";
		$query = Query($q);
		@motifs2 = ();
		$motifcount2 = 0;
		# while (($motif, $hits, $chrs, $strands) = Fetch($query))
		if ((Numrows($query) == 0) and ($peak2_has_motifs == 1))
		{
			die("\n\nError: No peaks found for query, even though there should be some:\n\n$q\n\n");
		}
		while (($source, $motif, $hits) = Fetch($query))
		{
			# # Verification that just querying by the peak start/stop is a unique identifier
			# die("Error: Multiple chromosomes found ($chrs) for symbol '$symbol2' and peak '$peak2'") if ($chrs != 1);
			# die("Error: Multiple strand orientations found for symbol '$symbol2' and peak '$peak2'") if ($strands != 1);
			
			# Shouldn't happen
			die("\n\nError: This query gives motifs for '$symbol2' peak2 '$peak2', but there shouldn't be any according to '$infile':\n\n$q\n\n") if ($peak2_has_motifs == 0);
			
			push(@motifs2, "$source|$motif");
			$motifcount2 += $hits;

			print "       >> $source >> $motif ($hits)\n" if (switch('debug'));
		}
		
		
		
		
		# Print
		print OUT "$symbol1\t$symbol2\t$bindtype\t$mindist\t$motifcount1\t$motifcount2\t".join(",", @motifs1)."\t".join(",", @motifs2)."\n";


		
		stepme(1, 1);
	}
	close(IN);
	stopme(1);
	
	
	
	
}
stoptime();

showmeall(1);

done();
