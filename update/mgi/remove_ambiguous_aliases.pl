#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0";
# ($var) = args(1);
args(0);

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

$mainquery = Query("SELECT alias, GROUP_CONCAT(DISTINCT symbol ORDER BY symbol SEPARATOR '|'), COUNT(DISTINCT symbol) AS c FROM mgi_aliases GROUP BY alias ORDER BY c DESC;");
startme("Removing ambiguous aliases for all symbols in table 'mgi_aliases'", 0, Numrows($mainquery));
starttime();
$affected = 0;
while (($alias, $symbols, $symbol_count) = Fetch($mainquery))
{
	addme("got $symbol_count symbols for alias", $alias);

	if ($symbol_count > 1)
	{
		if (!switch('debug'))
		{
			$query = Query("DELETE FROM mgi_aliases WHERE alias='".esc($alias)."'");
		}
		else
		{
			$query = Query("SELECT id FROM mgi_aliases WHERE alias='".esc($alias)."'");
		}
		$affected += Numrows($query);
		
		addme("deleted ambiguous alias for alias", $alias);
		foreach $symbol (split(/\|/, $symbols))
		{
			addme("deleted ambiguous alias for symbol", $symbol);
			addme("deleted ambiguous alias for alias|symbol pair", "$alias|$symbol");
		}
	}
	else
	{
		addme("retained unambiguous alias for alias", $alias);
		foreach $symbol (split(/\|/, $symbols))
		{
			addme("retained unambiguous alias for symbol", $symbol);
			addme("retained unambiguous alias for alias|symbol pair", "$alias|$symbol");
		}
	}
	
	stepme(1000);
}
stopme();
stoptime();

state(commify($affected)." rows affected");

showmeall(1);

done();
