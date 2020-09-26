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

Clear('mgi_aliases');

$mainquery = Query("SELECT DISTINCT mgi_accession_id, marker_symbol, marker_synonyms_pipe_separated FROM mgi ORDER BY marker_symbol");
startme("Getting all aliases from table 'mgi' and inserting them into table 'mgi_aliases'", 0, Numrows($mainquery));
starttime();
while (($mgiid, $symbol, $aliases) = Fetch($mainquery))
{
	die("Error: Couldn't parse MGI ID '$mgiid'") if ($mgiid !~ /^MGI:\d+$/);
	die("Error: Couldn't parse symbol '$symbol'") if ($symbol !~ /^[A-Za-z0-9_\.\- \(\)\/,\[\]:\+<>'=#*@;%~&`"]+$/);
	
	# Loop through aliases
	if (defined($aliases))
	{
		# Remove double quotes
		$aliases =~ s/^"//;
		$aliases =~ s/"$//;
		foreach $alias (unique(split(/\|/, $aliases)))
		{
			# Check aliases
			die("Error: Couldn't parse alias '$alias'") if ($alias !~ /^[A-Za-z0-9_\.\- \(\)\/,\[\]:\+<>'=#*@;%~&`"]+$/);
			
			if ($symbol eq $alias)
			{
				# die("Error: Weird: Official symbol '$alias' is also recorded as an alias (for itself, '$othersymbol')")
				addme("skipped alias that was identical to the gene name for symbol", $symbol);
				next;
			}

			$accs = '' if (!defined($accs));
			$ensg = '' if (!defined($ensg));
			$q = "INSERT INTO mgi_aliases SET species='MOUSE', mgiid='$mgiid', symbol='".esc($symbol)."', alias='".esc($alias)."'";
			$q =~ s/=''/=NULL/g;
			Query($q);
			
			addme("total aliases", $alias);
			addme("total alias|symbol pairs", "$alias|$symbol");
			addme("total symbols", $symbol);
			addme("total mgi ids", $mgiid);
		}
	}
	
	stepme(1000);
}
stopme();
stoptime();

showmeall(1);

Optimize('mgi_aliases');

done();
