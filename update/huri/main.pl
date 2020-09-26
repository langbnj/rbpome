#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'huri';

our $usage = "$0 [type: huri/hi-union]\n\nExample: $0 huri\nExample: $0 hi-union\n";
($type) = args(1);

$infile = "input/$type.tsv";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

# Clear($table);
Query("DELETE FROM `$table` WHERE type='$type'");
state("Deleted rows of type '$type' from table '$table'");


# start

startme("Reading '$infile' and inserting into '$table'", 0, chompme(`cat $infile | wc -l`));
starttime();
<IN>;	# Skip header
while (<IN>)
{
	chomp;
	
	# HuRI:
	# Ensembl_gene_id_a	Ensembl_gene_id_b	in_screen_1	in_screen_2	in_screen_3	in_screen_4	in_screen_5	in_screen_6	in_screen_7	in_screen_8	in_screen_9	in_assay_v1	in_assay_v2	in_assay_v3
	# ENSG00000000005	ENSG00000061656	0	0	0	0	0	1	0	0	0	0	1	0
	# ENSG00000000005	ENSG00000099968	0	0	0	0	1	0	0	0	0	0	1	0
	# ENSG00000000005	ENSG00000104765	0	0	1	0	0	0	0	0	0	1	0	0
	# ENSG00000000005	ENSG00000105383	0	0	0	0	0	0	0	0	1	0	0	1
	# ENSG00000000005	ENSG00000114455	0	0	0	0	0	0	0	0	1	0	0	1
	# ENSG00000000005	ENSG00000124103	0	0	0	0	0	0	1	0	0	0	0	1
	# ENSG00000000005	ENSG00000139637	0	0	0	0	0	0	1	0	0	0	0	1
	# ENSG00000000005	ENSG00000150337	0	0	0	0	0	0	0	0	1	0	0	1
	# ENSG00000000005	ENSG00000157613	0	0	0	0	1	0	0	0	0	0	1	0
	
	# HI-Union:
	# Ensembl_gene_id_a	Ensembl_gene_id_b	in_HI_III	in_tested_space	in_Rual	in_Venkatesan	in_Yu	in_Rolland	in_Yang
	# ENSG00000000005	ENSG00000061656	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000099968	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000104765	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000105383	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000114455	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000124103	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000139637	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000150337	1	0	0	0	0	0	0
	# ENSG00000000005	ENSG00000157613	1	0	0	0	0	0	0

	
	
	@a = split(/\t/);
	
	$ensg1 = $a[0];
	$ensg2 = $a[1];
	
	# Get gene symbols from gencode_gff3_gene (GENCODE v27)
	$query = Query("SELECT symbol FROM gencode_gff3_gene WHERE ensg='$ensg1' AND species='human'");
	($symbol1) = FetchOne($query);

	$query = Query("SELECT symbol FROM gencode_gff3_gene WHERE ensg='$ensg2' AND species='human'");
	($symbol2) = FetchOne($query);
	
	# This extra information isn't there anymore
	# if ($type eq 'huri')
	# {
	# 	die("Error: Expected 14 columns, but got ".scalar(@a)) if (scalar(@a) != 14);
	# 	$screens = $a[2] + $a[3] + $a[4] + $a[5] + $a[6] + $a[7] + $a[8] + $a[9] + $a[10];
	# 	$assays = $a[11] + $a[12] + $a[13];
	# }
	# else
	# {
	# 	die("Error: Expected 9 columns, but got ".scalar(@a)) if (scalar(@a) != 9);
	# 	$screens = $a[2] + $a[3] + $a[4] + $a[5] + $a[6] + $a[7] + $a[8];
	# 	$assays = 0;
	# }
	#
	# $total = $screens + $assays;
	
	addme("total ensgs", $ensg1);
	addme("total ensgs", $ensg2);
	addme("total interactions", "$ensg1|$ensg2");
	
	# addme("total interactions observed in total screens and assays: $total total", "$ensg1|$ensg2");
	# addme("total interactions observed in screens: $screens screens", "$ensg1|$ensg2");
	# addme("total interactions observed in assays: $assays assays", "$ensg1|$ensg2");
	
	# Query("INSERT INTO `$table` SET ensg1='$ensg1', ensg2='$ensg2', type='$type', species='human', screens='$screens', assays='$assays', total='$total'");
	Query("INSERT INTO `$table` SET ensg1='$ensg1', ensg2='$ensg2', symbol1='$symbol1', symbol2='$symbol2', type='$type', species='human'");
	
	stepme(1000);
}
stopme();
stoptime();

showmeallsorted(1);

Optimize($table);

done();
