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

$mainquery = Query("SELECT marker_symbol, COUNT(id) AS c FROM mgi GROUP BY marker_symbol ORDER BY c DESC");
startme("Removing symbols that occur multiple times in table 'mgi'", 0, Numrows($mainquery));
starttime();
$affected = 0;
while (($symbol, $count) = Fetch($mainquery))
{
	addme("got $count occurrences for symbol", $symbol);

	if ($count > 1)
	{
		if (!switch('debug'))
		{
			$query = Query("DELETE FROM mgi WHERE marker_symbol='".esc($symbol)."'");
		}
		else
		{
			$query = Query("SELECT id FROM mgi WHERE marker_symbol='".esc($symbol)."'");
		}
		$affected += Numrows($query);
		
		addme("deleted repeated symbol", $symbol);
	}
	else
	{
		addme("retained unique symbol", $symbol);
	}
	
	stepme(1000);
}
stopme();
stoptime();

state(commify($affected)." rows affected");

showmeall(1);

done();
