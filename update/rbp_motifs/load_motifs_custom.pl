#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'rbp_motif_list';
$source = 'custom';
$species = 'human';

our $usage = "$0";
# ($var) = args(1);
args(0);

Query("DELETE FROM `$table` WHERE source='$source'");
state("Cleared '$source' rows from table '$table'");


@motifs = ();
# Insert motifs here
push(@motifs, "YBX3|UCCAUCA");	# From https://www.uniprot.org/uniprot/P16989
# Insert motifs here



# start

startme("Inserting custom motifs into table '$table'");
$inserted = 0;
foreach $motif (@motifs)
{
	@a = split(/\|/, $motif);
	$symbol = $a[0];
	$motif = $a[1];

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
	
	addme("total symbol|motifs", "$symbol|$motif");
	
	$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', source='$source', motif='$motif', first=1";
	$s =~ s/=''/=NULL/g;
	Query($s);
	$inserted++;
	
	stepme(100);
}
stopme();

showmeall(1);

state("Rows inserted: ".commify($inserted));

Optimize($table);

done();
