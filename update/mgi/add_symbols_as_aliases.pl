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

$mainquery = Query("SELECT mgi_accession_id, marker_symbol FROM mgi GROUP BY marker_symbol ORDER BY marker_symbol");
startme("Getting all gene symbols from table 'mgi' and inserting them into table 'mgi_aliases' as pseudo-aliases (for easy lookup of e.g. MGI ID numbers via the 'alias' column)", 0, Numrows($mainquery));
starttime();
while (($mgiid, $symbol) = Fetch($mainquery))
{
	addme("total symbols", $symbol);
	

	# Make sure the symbol isn't also an alias for something else (and delete those incorrect aliases)
	$query = Query("SELECT symbol, alias FROM mgi_aliases WHERE alias='$symbol'");
	while (($othersymbol, $alias) = Fetch($query))
	{
		# die("Error: Official symbol '$symbol' is also recorded as an alias (for other symbol '$othersymbol')");
		if ($othersymbol eq $alias)
		{
			die("Error: Weird: Official symbol '$alias' is also recorded as an alias (for itself, '$othersymbol')") 
		}
		# elsif (lc($othersymbol) eq lc($alias))
		# {
		# 	# The only difference is the case, e.g. alias 'ADHFe1' for official symbol 'ADHFE1'. In these cases I am not inserting an alias (since we're not in the insert branch, below)
		# 	addme("got exactly one alias where the only difference to the symbol is the case (not adding another alias) for symbol", $symbol);
		# 	addme("got exactly one alias where the only difference to the symbol is the case (not adding another alias) for symbol|alias", "$symbol|$alias");
		# 	addme("got exactly one alias where the only difference to the symbol is the case (not adding another alias) for alias", $alias);
		# }
		else
		{
			addme("deleting alias that is also an official symbol for alias", $symbol);
			addme("deleting alias that is also an official symbol for other symbol", $othersymbol);

			if (!switch('debug'))
			{
				$query2 = Query("DELETE FROM mgi_aliases WHERE symbol='$othersymbol' AND alias='$symbol'");
			}
			else
			{
				$query2 = Query("SELECT id FROM mgi_aliases WHERE symbol='$othersymbol' AND alias='$symbol'");
			}
			$deleted += Numrows($query2);
		}
	}


	addme("adding gene symbol as pseudo-alias for gene symbol", $symbol);
	
	if (!switch('debug'))
	{
		$q = "INSERT INTO mgi_aliases SET mgiid='$mgiid', symbol='".esc($symbol)."', alias='".esc($symbol)."'";
		$q =~ s/=''/=NULL/g;
		$query2 = Query($q);
		$added += Numrows($query2);
	}
	else
	{
		$added++;
	}
	
	stepme(1000);
}
stopme();
stoptime();

state(commify($deleted)." rows deleted");
state(commify($added)." rows added");

showmeall(1);

done();
