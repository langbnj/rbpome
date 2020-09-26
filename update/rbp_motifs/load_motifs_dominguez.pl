#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'rbp_motif_list';
$source = 'dominguez';
$species = 'human';

our $usage = "$0 [-firstonly]\n\n -firstonly: Only use the first (top) motif per RBP (skip all lower-ranked significant motifs)\n\nExample: $0 -firstonly";
# ($var) = args(1);
args(0);

$infile = "input/dominguez-s3.txt";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

Query("DELETE FROM `$table` WHERE source='$source'");
state("Cleared '$source' rows from table '$table'");



# start

startme("Reading '$infile' and inserting motifs into table '$table'");
starttime();
<IN>;	# Skip header
$inserted = 0;
while (<IN>)
{
	chomp;
	
	s/\s+$//;
	
	# RBP	ReadLen	Temp	ENCODE Accession ID	MostRConc	Motif 5mer_logonum_stepwiseRminus1
	# A1CF	20	4	ENCSR934TDK	5 nM	AAUUA_1_5.28	AAUCA_1_3.98	UAAUU_1_3.51	AAUAA_1_2.50	AAUGA_1_2.07	CUAAU_1_1.69	UAAUA_1_1.66	UUAAU_1_1.62	AUAAU_1_1.53	AAAUU_1_1.39
	# BOLL	20	4	ENCSR497LIF	20 nM	UUUUU_1_8.79	UGUUU_1_3.86	UUUUA_1_3.50	UUAUU_1_2.70	UUUUC_1_2.49	UUGUU_1_2.47	GUUUU_1_2.40	AUUUU_1_2.38	CUUUU_1_2.30	UUUUG_1_2.28	UUUAU_1_2.06	UGUUA_2_3.25	GUGUU_2_2.04
	# CELF1	40	21	ENCSR992NHR	64 nM	UAUGU_1_2.55	UUUGU_1_2.55	UGUGU_1_2.18	UGUCU_1_1.96	UGUUU_1_1.69	AUGUU_1_1.51	UGUCC_1_1.50	UUGUU_1_1.34	AAUGU_1_1.33
	
	@a = split(/\t/);
	$symbol = $a[0];
	@motifs = splice(@a, 5);

	# Remove whitespace from symbol (necessary for NOVA1)
	$symbol =~ s/^\s+//;
	$symbol =~ s/\s+$//;
	
	addme("total symbols", $symbol);
	
	# Get UniProt accession
	# $query = Query("SELECT DISTINCT name FROM uniprot WHERE gene='$symbol' AND species='$species'");
	# ($name) = FetchOne($query);
	# $query = Query("SELECT DISTINCT acc FROM uniacc WHERE name='$name' AND species='$species'");
	# ($acc) = FetchOne($query);
	$query = Query("SELECT DISTINCT acc FROM clip_gene WHERE symbol='$symbol' AND species='$species'");
	$acc = '';
	if (Numrows($query) > 0)
	{
		($acc) = FetchOne($query);
	}
	
	$first = 1;
	foreach $s (@motifs)
	{
		# AAUUA_1_5.28
		$s =~ /^(\w{5})_(\d+)_([\d\.]+)$/ or die("Error: Couldn't match motif string '$s'");
		$motif = $1;
		$logonum = $2;
		$score = $3;
		
		addme("total symbol|motifs", "$symbol|$motif");
		addme("total symbol|motif|logonums", "$symbol|$motif|$logonum");
		addme("total symbol|motifstrings", "$symbol|$s");
		
		$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', source='$source', motif='$motif', first='$first', logonum='$logonum', score='$score'";
		$s =~ s/=''/=NULL/g;
		Query($s);
		$inserted++;

		if ($first == 1)
		{
			$first = 0;
		}
	}
	
	stepme(100);
}
stopme();
stoptime();

showmeall(1);

state("Rows inserted: ".commify($inserted));

Optimize($table);

done();
