#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $table = "clip_raw";
$metatable = 'clip_encode_metadata';

our $usage = "$0 [species] [assembly] [mapping: gene/transcript/exon] [log2fold threshold] [-log10p threshold] [-12]\n\n -12: Use combined '1, 2' replicate files (only) and insert them as eclip_encode_12\n\nExample: $0 human GRCh38 transcript 3 5";
($species, $assembly, $map, $thresh_log2fold, $thresh_log10p) = args(5);
# args(0);

$type = 'eclip_encode';
if (switch('12'))
{
	$type = 'eclip_encode_12';
}

$table = "clip_raw_$map";

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

if (!switch('debug'))
{
	# Clear($table);
	state("Clearing '$type' '$map' rows from table '$table'...");
	$query = Query("DELETE FROM `$table` WHERE type='$type' AND map='$map'");
}


# Get RBP accessions
startme("Getting RBP symbol-to-acc mapping from UniProt");
$query = Query("SELECT a.primary_acc, p.gene, COUNT(DISTINCT a.primary_acc) FROM uniprot p, uniacc a WHERE p.name=a.name AND p.species='$species' AND a.species='$species' GROUP BY p.gene");
%acc = ();
while (($acc, $symbol, $count) = Fetch($query))
{
	# Skip ambiguous ones
	next if ($count != 1);
	
	$acc{uc($symbol)} = $acc;
	
	stepme(10000);
}
stopme();

# startme("Getting RBP ORF number-to-acc mapping from UniProt");
# $query = Query("SELECT a.primary_acc, p.orf, COUNT(DISTINCT a.primary_acc) FROM uniprot p, uniacc a WHERE p.name=a.name AND p.species='$species' AND a.species='$species' GROUP BY p.orf");
# while (($acc, $orf, $count) = Fetch($query))
# {
# 	# Skip ambiguous ones
# 	next if ($count != 1);
#
# 	die("Error: Overlap between ORFs and symbols for orf '$orf'") if (exists($acc{$orf}));
# 	$acc{uc($orf)} = $acc;
#
# 	stepme(10000);
# }
# stopme();


if ($species eq 'human')
{
	# Use GENCODE
	# $gff3file = "../gencode_gff3/input/sorted/gencode.human.v27.chr_patch_hapl_scaff.annotation.$map.sorted.gff3";
	$gff3file = "../gencode_gff3/input/gencode.human.v27.chr_patch_hapl_scaff.annotation.gff3";
}
elsif ($species eq 'mouse')
{
	# Use GENCODE
	# $gff3file = "../gencode_gff3/input/sorted/gencode.mouse.vM16.chr_patch_hapl_scaff.annotation.$map.sorted.gff3";
	$gff3file = "../gencode_gff3/input/gencode.mouse.vM16.chr_patch_hapl_scaff.annotation.gff3";
}
elsif ($species eq 'yeast')
{
	# Use Ensembl
	# $gff3file = "../ensembl_gff3/input/sorted/Saccharomyces_cerevisiae.R64-1-1.92.$map.sorted.gff3";
	$gff3file = "../ensembl_gff3/input/Saccharomyces_cerevisiae.R64-1-1.92.gff3";
}
elsif ($species eq 'drome')
{
	# Use Ensembl
	# $gff3file = "../ensembl_gff3/input/sorted/Drosophila_melanogaster.BDGP6.92.$map.sorted.gff3";
	$gff3file = "../ensembl_gff3/input/Drosophila_melanogaster.BDGP6.92.gff3";
}
else
{
	die("Error: Unhandled species '$species'");
}




$query = Query("SELECT DISTINCT fullspecies FROM ensembl WHERE species='$species'");
($fullspec) = FetchOne($query);
$fullspec =~ s/_/ /g;

# $mainquery = Query("SELECT DISTINCT file_accession, experiment_target, biosample_term_name, biological_replicate_s, technical_replicate FROM `$metatable` WHERE biosample_organism='$fullspec' AND assembly='$assembly' ORDER BY experiment_target, biological_replicate_s");
# HepG2 and K562 only
$mainquery = Query("SELECT DISTINCT file_accession, experiment_target, biosample_term_name, biological_replicate_s, technical_replicate FROM `$metatable` WHERE biosample_organism='$fullspec' AND assembly='$assembly' AND biosample_term_name IN ('HepG2', 'K562') ORDER BY experiment_target, biological_replicate_s");
state("Getting '$species' '$assembly' datasets from '$metatable', parsing their files and mapping to $map"."s [".Numrows($mainquery)."]:");
starttime();
$inserted = 0;
$threshfail_log2fold = 0;
$threshfail_log10p = 0;
$i = 0;
while (($fileacc, $exptarget, $celltype, $rep, $techrep) = Fetch($mainquery))
{
	$i++;
	
	if (switch('12'))
	{
		# Combined reps from ENCODE (simple intersection with exact coordinate match, which they added in 03/2019)
		next if ($rep ne '1, 2');
		$rep = '1';
	}
	else
	{
		# Individual reps (normal)
		next if ($rep eq '1, 2');
		die("Error: Technical replicate isn't '1'") if ($techrep ne '1');
	}
	
	$bedfile = "input/$fileacc.bed";
	# $infile = "input/$fileacc.intersect.bed";
	# $infile = "input/intersect.$map.$fileacc.bed";
	$infile = "tmp/intersect.$fileacc.bed";

	# Run bedtools intersect, unless the output file exists already
	if (!-s $infile)
	{
		# run("Intersect", "bedtools intersect -s -wb -nonamecheck -a $bedfile -b $gff3file > $infile");
		run("Intersect", "bedtools intersect -s -wo -nonamecheck -a $bedfile -b $gff3file > $infile");
		# run("Intersect", "~/scripts/qsub.sh bedtools intersect -s -wo -nonamecheck -a $bedfile -b $gff3file \\> $infile");
	}
	# next;

	$exptarget =~ /^(\w+)-$species$/ or die("Error: experiment_target '$exptarget' doesn't end in '-$species'");
	$symbol = $1;

	die("Error: No acc for RBP symbol '$symbol'") if (!exists($acc{uc($symbol)}));
	$acc = $acc{uc($symbol)};
	
	open(IN, $infile) or die("Error: Couldn't open '$infile'");
	# startme(" >> Importing '$infile' ($symbol rep $rep) ($i / ".Numrows($mainquery).")", 1);
	startme(" >> $i / ".Numrows($mainquery)." >> Importing '$infile' ($symbol $celltype rep $rep)", 1, chompme(`cat '$infile' | wc -l`));
	starttime();
	while (<IN>)
	{
		stepme(10000, 1);
		# last if (getme() % 10 == 0);

		chomp();
		# # chr14_GL000194v1_random	72393	72489	GEMIN5_K562_rep01	1000	-	3.77792607096336	15.360722739658	-1	-1
		# # chr14_GL000194v1_random	71631	71668	GEMIN5_K562_rep01	1000	-	3.66940161418519	4.09004729666783	-1	-1
		# # chr14_GL000194v1_random	73192	73240	GEMIN5_K562_rep01	200	-	2.51739852074014	4.07700356868956	-1	-1
		# chr7	98881814	98881876	AARS_K562_rep01	200	+	1.23600949878672	12.1431261676372	-1	-1	chr7	HAVANA	transcript	98877933	98890399	.	+	.	ID=ENST00000417523.5;Parent=ENSG00000196367.12;gene_id=ENSG00000196367.12;transcript_id=ENST00000417523.5;gene_type=protein_coding;gene_name=TRRAP;transcript_type=protein_coding;transcript_name=TRRAP-203;level=2;protein_id=ENSP00000401107.1;transcript_support_level=4;tag=alternative_5_UTR,mRNA_end_NF,cds_end_NF;havana_gene=OTTHUMG00000150403.4;havana_transcript=OTTHUMT00000317982.2
		# 																										chr1	ENSEMBL	exon	925922	926013	.	+	.	ID=exon:ENST00000617307.4:2;Parent=ENST00000617307.4;gene_id=ENSG00000187634.11;transcript_id=ENST00000617307.4;gene_type=protein_coding;gene_name=SAMD11;transcript_type=protein_coding;transcript_name=SAMD11-212;exon_number=2;exon_id=ENSE00001763717.1;level=3;protein_id=ENSP00000482090.1;transcript_support_level=5;tag=basic,appris_alternative_2;havana_gene=OTTHUMG00000040719.10
		# chr7	98881814	98881876	AARS_K562_rep01	200	+	1.23600949878672	12.1431261676372	-1	-1	chr7	HAVANA	transcript	98877933	98890399	.	+	.	ID=ENST00000417523.5;Parent=ENSG00000196367.12;gene_id=ENSG00000196367.12;transcript_id=ENST00000417523.5;gene_type=protein_coding;gene_name=TRRAP;transcript_type=protein_coding;transcript_name=TRRAP-203;level=2;protein_id=ENSP00000401107.1;transcript_support_level=4;tag=alternative_5_UTR,mRNA_end_NF,cds_end_NF;havana_gene=OTTHUMG00000150403.4;havana_transcript=OTTHUMT00000317982.2	62
		
		@a = split(/\t/);
	
		# die("Error: Expected 19 columns in line '$_', but found ".scalar(@a)) if (scalar(@a) != 19);
		die("Error: Expected 20 columns in line '$_', but found ".scalar(@a)) if (scalar(@a) != 20);

		$chromo = $a[0];
		$start = $a[1];
		$stop = $a[2];
		$sample = $a[3];
		$score = $a[4];
		$strand = $a[5];
		die if (($strand ne '+') and ($strand ne '-'));
		$log2fold = $a[6];
		die("Error: log2fold is $log2fold") if ($log2fold !~ /^-?\d+\.\d+(e-\d+)?$/);
		$log10p = $a[7];
		die("Error: Couldn't parse log10p '$log10p' in line '$_'") if ($log10p !~ /^\d+(\.\d+)?(e-\d+)?$/);
		die if ($a[8] ne '-1');
		die if ($a[9] ne '-1');
		
		# BED files are zero-based (https://www.biostars.org/p/84686/), hence:
		$start += 1;
		
		# Apply thresholds
		$threshfail_log2fold++ if ($log2fold < $thresh_log2fold);
		$threshfail_log10p++ if ($log10p < $thresh_log10p);

		next if ($log2fold < $thresh_log2fold);
		next if ($log10p < $thresh_log10p);

		
		
		# Parse the mapped fields from the GENCODE/Ensembl GFF3 file (the transcript annotation etc.)

		# Verify that it's the desired feature type (e.g. 'transcript')
		$thistype = $a[12];
		next if ($thistype ne $map);
		
		# Get feature ID (e.g. ENSTV)
		$att = $a[18];
		if ($map eq 'gene')
		{
			$att =~ /^ID=([^;]+);/ or die("Error: Couldn't parse ID from attributes '$att'");
			$ensgv = $1;
			$enstv = '';
			$ensev = '';
		}
		elsif ($map eq 'transcript')
		{
			$att =~ /;gene_id=([^;]+);/ or die("Error: Couldn't parse gene_id from attributes '$att'");
			$ensgv = $1;
			$att =~ /^ID=([^;]+);/ or die("Error: Couldn't parse ID from attributes '$att'");
			$enstv = $1;
			$ensev = '';
		}
		elsif ($map eq 'exon')
		{
			$att =~ /;gene_id=([^;]+);/ or die("Error: Couldn't parse gene_id from attributes '$att'");
			$ensgv = $1;
			$att =~ /;transcript_id=([^;]+);/ or die("Error: Couldn't parse transcript_id from attributes '$att'");
			$enstv = $1;
			$att =~ /;exon_id=([^;]+);/ or die("Error: Couldn't parse exon_id from attributes '$att'");
			$ensev = $1;
		}
		else
		{
			die("Error: Unhandled mapping type '$map'");
		}
		
		addme("total chromosomes", $chromo);
		die("Error: Couldn't parse start coordinate '$start'") if ($start !~ /^\d+$/);
		die("Error: Couldn't parse stop coordinate '$stop'") if ($stop !~ /^\d+$/);
		die("Error: Couldn't parse strand '$strand'") if ($strand !~ /^(\+|\-)$/);
		die if ($score !~ /^\d+$/);
		
		if (switch('12'))
		{
			$sample =~ /^$symbol\_(HepG2|K562)_IDR$/ or die("Error: Couldn't parse sample '$sample'");
		}
		else
		{
			# $sample =~ /^$symbol\_(HepG2|K562|SM-9MVZL_hAdrenalGland)_rep0$rep$/ or die("Error: Couldn't parse sample '$sample'");
			$sample =~ /^$symbol\_(HepG2|K562)_rep0$rep$/ or die("Error: Couldn't parse sample '$sample'");
			# $celltype = $1;
		}
		
		$q = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', ensgv='$ensgv', enstv='$enstv', ensev='$ensev', species='$species', type='$type', map='$map', celltype='$celltype', rep='$rep', score='$score', log2fold='$log2fold', log10p='$log10p', chr='$chromo', start='$start', stop='$stop', strand='$strand'";
		$q =~ s/=''/=NULL/g;
		Query($q) if (!switch('debug'));

		$inserted++;
	}
	close(IN);
	stopme(1);

	# exit;
}
done();
stoptime();

showmesome(50);

nl();
state("Rows inserted: ".commify($inserted), 1);
nl();
state("Peaks failed log2fold threshold: ".commify($threshfail_log2fold), 1);
state("Peaks failed log10p threshold: ".commify($threshfail_log10p), 1);
nl();

Optimize($table);

done();
