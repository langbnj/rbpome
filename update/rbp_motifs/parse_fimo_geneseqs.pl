#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
# $table = 'rbp_motifs';
$species = 'human';
$map = 'gene';
$intable = "clip_raw_$map";
$table = 'rbp_motifs';	# to be extended by _extend5_grep_$type below
$type = "geneseqs";

our $usage = "$0 [method: fimo/fimobg/fimobgi] [type: eclip_encode/eclip_tom/...] [p/q-value threshold] [-qvalue]\n\n-qvalue: Use a q-value threshold of 0.05 (instead of a p-value of 0.05)\n\nExample: $0 fimo 0.001\nExample: $0 fimobgi 0.05 -qvalue";
($method, $thresh) = args(2);

$tmpthresh = $thresh;
$tmpthresh =~ s/\./_/g;
$tmpqvalue = '';
$psig = 1; $qsig = 0;
if (switch('qvalue'))
{
	$tmpqvalue = '-qvalue';
	$psig = 0; $qsig = 1;
}

$table = 'rbp_motifs_'.$method.'_geneseqs';
$outfile = "tmp-$table-$thresh$tmpqvalue.txt";
# if (!-e $outfile)
# {
	# No header needed for LOAD DATA INFILE
	# print OUT "symbol\tacc\tspecies\tcelltype\trep\tmethod\ttype\tchr\tstart\tstop\tstrand\tsource\tmotif\thit\tpvalue\tqvalue\n";
# }
# else
# {
# 	open(OUT, ">>$outfile") or die("Error: Couldn't append to '$outfile'");
# }




# Check if the table already contains results
if (switch('qvalue'))
{
	$query = Query("SELECT id FROM $table WHERE qvalue!=1 LIMIT 1");
}
else
{
	$query = Query("SELECT id FROM $table WHERE qvalue=1 LIMIT 1");
}
if (Numrows($query) > 0)
{
	state("Table '$table' already contains '$tmpqvalue' entries, skipping run!");
	exit;
}



# Only run if there isn't already a temporary file we can import
if (!-s $outfile)
{
	open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");
	# Get RBP accessions
	startme("Getting RBP symbol-to-acc mapping from table '$intable'");
	$query = Query("SELECT symbol, acc FROM `$intable` WHERE species='$species' AND map='$map' GROUP BY symbol, acc ORDER BY symbol, acc");
	%acc = ();
	while (($symbol, $acc) = Fetch($query))
	{
		$acc{$symbol} = $acc;

		stepme(1000);
	}
	stopme();



	# Start
	nl();
	state("Parsing motif hits from FIMO output files and writing to temporary file '$outfile':", 1);
	open(LS, "ls -1 output/output-$method-$type-*-$tmpthresh$tmpqvalue.txt|") or die("Error: Couldn't ls input files at 'output/output-$method-$type-*-$tmpthresh$tmpqvalue.txt'");
	starttime();
	$inserted = 0;
	while (<LS>)
	{
		chomp;
		# output/output-fimo-eclip_tom-rnacompete-TARDBP-K562-0_05-qvalue.txt
		$infile = $_;
	
		$infile =~ /^output\/output-$method-$type-([^\-]+)-$tmpthresh$tmpqvalue\.txt/ or die("Error: Couldn't match input file name '$infile");
		$symbol = $1;
	
		die("Error: Unexpected input filename from ls: $infile") if ($infile ne "output/output-$method-$type-$symbol-$tmpthresh$tmpqvalue.txt");
	
		# Read results
		open(IN, $infile) or die("Error: Couldn't open '$infile'");

		# startme("Reading '$infile'", 0, chompme(`cat $infile | wc -l`));
		# startme("Reading '$infile'");
		startme(" >> $symbol", 1);
		while (<IN>)
		{
			stepme(100000, 1);

			chomp;
	
			# Skip header
			next if ($_ eq "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence");
	
			# Stop parsing on blank or comment line
			last if (/^$/);
			last if (/^#/);

			# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
			# HNRNPH2_1		AKAP1|HepG2|chr10|100150213|100150263	11	15	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100150373|100150423	9	13	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100150393|100150443	29	33	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151373|100151423	32	36	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151513|100151563	32	36	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151633|100151683	10	14	+	9.75524	0.000977
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151633|100151683	36	40	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151673|100151723	2	6	+	3.27273	0.0156
			# HNRNPH2_1		AKAP1|HepG2|chr10|100151673|100151723	30	34	+	9.75524	0.000977

			@a = split(/\t/);
			$motifid = $a[0];
			$site = $a[2];
			$tmp_motifstart = $a[3];
			$tmp_motifstop = $a[4];
			$motifstrand = $a[5];
			$score = $a[6];
			$pvalue = $a[7];
			if (switch('qvalue'))
			{
				# q-value threshold used
				$qvalue = $a[8];
				$psig = 0;
				$qsig = 1;
			}
			else
			{
				# p-value threshold used
				$qvalue = '1';
				$psig = 1;
				$qsig = 0;
			}
			
			# Verify
			die("Error: FIMO reported - strand hits") if ($motifstrand ne '+');
			die("Error: FIMO motif positions are actually 0-based") if ($motifstart == 0);
			
	
			# Get source from the motifid's number (and remove the code)
			if    ($motifid =~ /^(\w+)_1(\d+)$/) { $motifid = $1.'_'.$2; $source = 'attract'; }
			elsif ($motifid =~ /^(\w+)_2(\d+)$/) { $motifid = $1.'_'.$2; $source = 'cisbp'; }
			elsif ($motifid =~ /^(\w+)_3(\d+)$/) { $motifid = $1.'_'.$2; $source = 'dominguez'; }
			elsif ($motifid =~ /^(\w+)_4(\d+)$/) { $motifid = $1.'_'.$2; $source = 'rbpdb'; }
			elsif ($motifid =~ /^(\w+)_5(\d+)$/) { $motifid = $1.'_'.$2; $source = 'rbpmap'; }
			elsif ($motifid =~ /^(\w+)_6(\d+)$/) { $motifid = $1.'_'.$2; $source = 'rnacompete'; }
			else
			{
				die("Error: Can't determine source from motifid '$motifid'");
			}
			
	
			# Parse site
			@site = split(/\|/, $site);
			die("Error: Expected 7 fields in '$site'") if (scalar(@site) != 7);
			$symbol = $site[0];
			$celltype = $site[1];
			$rep = $site[2];
			$chr = $site[3];
			$start = $site[4];
			$stop = $site[5];
			$strand = $site[6];

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


			# Verify RBP symbol & skip irrelevant motif hits
			$this_rbp = $motifid;
			$this_rbp =~ s/_\d+$//;
			# addme("total RBPs in whose peaks motifs occur", $symbol);
			# addme("total RBPs whose motifs occur", $this_rbp);
			# addme("total motif IDs which occur", $motifid);
			if ($this_rbp ne $symbol)
			{
				# addme("motif '$motifid' occurs in another RBP's peaks for motifid|chr|start|stop|strand (skipped)", "$motifid|$chr|$start|$stop|$strand");
				next;
			}
	
			# addme("total RBPs in whose peaks matching motifs occur", $symbol);
	
			# Get UniProt accession
			if (!exists($acc{$symbol}))
			{
				die("Error: No UniProt acc for symbol '$symbol'");
				# warn("Error: No UniProt acc for symbol '$symbol'");
				# $acc = '';
			}
			else
			{
				$acc = $acc{$symbol};
			}
		
			# # Hit
			# $hit{"$chr|$start|$stop|$strand"} = 1;

			# Insert
			# print OUT "$symbol\t$acc\t$species\t$celltype\t$rep\t$method\t$type\t$chr\t$start\t$stop\t$strand\t$source\t$motifid\t1\t$pvalue\t$qvalue\n";
			# print OUT "$symbol\t$acc\t$species\t$celltype\t$rep\t$method\t$type\t$chr\t$start\t$stop\t$strand\t$source\t$motifid\t1\t$pvalue\t$qvalue\t$psig\t$qsig\n";
			print OUT "$symbol\t$acc\t$species\t$celltype\t$rep\t$method\t$type\t$chr\t$start\t$stop\t$strand\t$source\t$motifid\t$motifstart\t$motifstop\t$score\t$pvalue\t$qvalue\t$psig\t$qsig\n";
			# $s =~ s/=''/=NULL/g;
			# Query($s);
			$inserted++;
	
			# stepme(1000);
		}
		close(IN);
		stopme(1);

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
		# 			if (!switch('qvalue'))
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
		# 			# $s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', celltype='$celltype', rep='$rep', method='$method', type='$type', chr='$chr', start='$start', stop='$stop', strand='$strand', source='$source', hit=0, pvalue='$pvalue', qvalue='$qvalue'";
		# 			print OUT "$symbol\t$acc\t$species\t$celltype\t$rep\t$method\t$type\t$chr\t$start\t$stop\t$strand\t$source\t\t0\t$pvalue\t$qvalue\n";
		# 			# $s =~ s/=''/=NULL/g;
		# 			# Query($s);
		# 			$inserted++;
		# 			stepme(100000, 1);
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
	close(OUT);
	nl();
	stoptime();

	showmeall(1);

	state("Lines written to '$outfile': ".commify($inserted));
}



# Debug workaround: Some of these files have 18 rows (they've already got psig and qsig fields), others still only have 16.
if (-z $outfile)
{
	state("Output file '$outfile' is empty, skipping parsing!");
	exit;
}
$psig = 1; $qsig = 0;
if (switch('qvalue'))
{
	$psig = 0; $qsig = 1;
}
$columns = chompme(`head -n1 $outfile | perl -ne 'chomp; \@a=split(/\t/); print scalar(\@a)'`);
die("Error: Column count is '$columns' in '$outfile' (should be 20)") if ($columns != 20);

state("Loading data into table '$table'...", 1);
starttime();
# Query("ALTER TABLE $table DISABLE KEYS");
# stoptime();

# if ($columns == 16)
# {
# 	$query = Query("LOAD DATA LOCAL INFILE '/users/gt/blang/update/rbp_motifs/$outfile' INTO TABLE $table (symbol, acc, species, celltype, rep, method, type, chr, start, stop, strand, source, motif, hit, pvalue, qvalue) SET id=NULL, extend5=0, psig=$psig, qsig=$qsig");
# }
if ($columns == 20)
{
	$query = Query("LOAD DATA LOCAL INFILE '/users/gt/blang/update/rbp_motifs/$outfile' INTO TABLE $table (symbol, acc, species, celltype, rep, method, type, chr, start, stop, strand, source, motif, hit, pvalue, qvalue, psig, qsig) SET id=NULL, extend5=0");
}
else
{
	die("Error: Column count is '$columns' in '$outfile' (should be 20)");
}
# state("Re-enabling keys...", 1);
# starttime();
# Query("ALTER TABLE $table ENABLE KEYS");
stoptime();
state("Rows inserted: ".commify(Numrows($query)));

# Optimize($table);

done();
