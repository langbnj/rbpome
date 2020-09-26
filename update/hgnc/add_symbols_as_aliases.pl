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

$mainquery = Query("SELECT hgnc_id, symbol, ensembl_gene_id, uniprot_ids FROM hgnc GROUP BY symbol ORDER BY symbol");
startme("Getting all gene symbols from table 'hgnc' and inserting them into table 'hgnc_aliases' as pseudo-aliases (for easy lookup of e.g. HGNC ID numbers or ENSGs via the 'alias' column)", 0, Numrows($mainquery));
starttime();
while (($hgncid, $symbol, $ensg, $accs) = Fetch($mainquery))
{
	addme("total symbols", $symbol);
	
	# Check ENSG
	if (!defined($ensg))
	{
		addme("no ensg for symbol (kept)", $symbol);
		$ensg = '';
	}
	else
	{
		die("Error: Couldn't parse ENSG '$ensg'") if ($ensg !~ /^ENSG\d{11}$/);
		addme("exactly one ensg for symbol (kept)", $symbol);
	}
	
	# Check UniProt accessions
	if (!defined($accs))
	{
		addme("no uniprot accs for symbol (kept)", $symbol);
		$accs = '';
	}
	else
	{
		# Remove double quotes
		$accs =~ s/^"//;
		$accs =~ s/"$//;

		if ($accs =~ /\|/)
		{
			addme("multiple uniprot accs for symbol (kept)", $symbol);
		}
		else
		{
			addme("exactly one uniprot acc for symbol (kept)", $symbol);
		}
		
		# Check individual accs
		foreach $acc (nsort(split(/\|/, $accs)))
		{
			die("Error: Couldn't parse acc '$acc'") if (($acc !~ /^[A-Z0-9]{6}$/) and ($acc !~ /^[A-Z0-9]{10}$/));

			$query = Query("SELECT name FROM uniacc WHERE acc='$acc' AND species='human'");
			if (Numrows($query) == 0)
			{
				addme("no name in uniacc for acc (kept)", $acc);
			}
			else
			{
				($name) = FetchOne($query);
				$query = Query("SELECT gene FROM uniprot WHERE name='$name' AND species='human'");
				($uniprot_hgnc) = FetchOne($query);
				if ($uniprot_hgnc ne $symbol)
				{
					# die("Error: UniProt gene symbol '$uniprot_hgnc' doesn't match symbol '$symbol' for accession '$acc' (total accs: '$accs')");
					addme("uniprot gene symbol mismatch for symbol (kept)", $symbol);
					addme("uniprot gene symbol mismatch for acc (kept)", $acc);
				}
				else
				{
					addme("uniprot gene symbol match for symbol", $symbol);
					addme("uniprot gene symbol match for acc", $acc);
				}
				# {
				# 	addme("uniprot accession not found in table uniacc for acc (kept)", $acc);
				# }
				# else
				# {
				# 	addme("uniprot accession successfully found in table uniacc for acc (kept)", $acc);
				# }
			}
		}
	}


	# Make sure the symbol isn't also an alias for something else (and delete those aliases)
	$query = Query("SELECT symbol, alias FROM hgnc_aliases WHERE alias='$symbol'");
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
				$query2 = Query("DELETE FROM hgnc_aliases WHERE symbol='$othersymbol' AND alias='$symbol'");
			}
			else
			{
				$query2 = Query("SELECT id FROM hgnc_aliases WHERE symbol='$othersymbol' AND alias='$symbol'");
			}
			$deleted += Numrows($query2);
		}
	}


	addme("adding gene symbol as pseudo-alias for gene symbol", $symbol);
	
	if (!switch('debug'))
	{
		$q = "INSERT INTO hgnc_aliases SET hgncid='$hgncid', symbol='$symbol', alias='$symbol', ensg='$ensg', accs='$accs'";
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
