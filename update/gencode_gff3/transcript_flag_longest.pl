#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_transcript';

our $usage = "$0 [type: ensgv/symbol\n\nExample: $0 ensgv";
($type) = args(1);
# args(0);

if ($type eq 'ensgv')
{
	$field = 'longest';
}
elsif ($type eq 'symbol')
{
	$field = 'longest_for_symbol';
}

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

state("Clearing '$field' field in table '$table'...");
Query("UPDATE `$table` SET `$field`=0");
done();


state("Flagging '$field' transcripts in table '$table':");
starttime();

$speciesquery = Query("SELECT DISTINCT species FROM `$table`");
while (($species) = Fetch($speciesquery))
{
	$affected = 0;
	$ensgvquery = Query("SELECT DISTINCT `$type` FROM `$table` WHERE species='$species' AND appris IS NOT NULL");
	startme(" >> Flagging '$field' APPRIS Principal transcript for '$species' '".$type."'s...", 1, Numrows($ensgvquery));
	while (($ensgv) = Fetch($ensgvquery))
	{
		$query = Query("SELECT enstv FROM `$table` WHERE species='$species' AND `$type`='$ensgv' AND appris IS NOT NULL ORDER BY appris LIKE 'P%' DESC, appris, LENGTH(seq) DESC LIMIT 1");
		($enstv) = FetchOne($query);
		
		$query = Query("UPDATE `$table` SET `$field`=1 WHERE species='$species' AND `$type`='$ensgv' AND appris IS NOT NULL AND enstv='$enstv'");
		$affected += Numrows($query);
		
		stepme(1000, 1);
	}
	stopme(1);
	
	state("Rows affected: ".commify($affected));
}
stoptime();

showmeall(1);

Optimize($table);

done();
