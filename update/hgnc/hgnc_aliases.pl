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

Clear('hgnc_aliases');

$mainquery = Query("SELECT hgnc_id, symbol, alias_symbol, ensembl_gene_id, uniprot_ids FROM hgnc ORDER BY symbol");
startme("Getting all aliases from table 'hgnc' and inserting them into table 'hgnc_aliases'", 0, Numrows($mainquery));
starttime();
while (($hgncid, $symbol, $aliases, $ensg, $accs) = Fetch($mainquery))
{
	die("Error: Couldn't parse HGNC ID '$hgncid'") if ($hgncid !~ /^HGNC:\d+$/);
	die("Error: Couldn't parse symbol '$symbol'") if ($symbol !~ /^[A-Za-z0-9\-\._@#αβγ]+$/);
	
	# Check ENSG
	if (!defined($ensg))
	{
		addme("no ensg for symbol (kept)", $symbol);
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
			addme("exactly one uniprot accs for symbol (kept)", $symbol);
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
	
	# Loop through aliases
	if (defined($aliases))
	{
		# Remove double quotes
		$aliases =~ s/^"//;
		$aliases =~ s/"$//;
		foreach $alias (nsort(split(/\|/, $aliases)))
		{
			# Check aliases
			die("Error: Couldn't parse alias '$alias'") if ($alias !~ /^[A-Za-z0-9\-\(\)\.\:\/_'#\+\[\]@\*αβγ]+$/);
			
			$accs = '' if (!defined($accs));
			$ensg = '' if (!defined($ensg));
			$q = "INSERT INTO hgnc_aliases SET hgncid='$hgncid', symbol='$symbol', alias='".esc($alias)."', ensg='$ensg', accs='$accs'";
			$q =~ s/=''/=NULL/g;
			Query($q);
			
			addme("total aliases", $alias);
			addme("total alias|symbol pairs", "$alias|$symbol");
			addme("total symbols", $symbol);
			addme("total hgnc ids", $hgncid);
			addme("total ensgs", $ensg);
			addme("total acc groups", $accs);
		}
	}
	
	stepme(1000);
}
stopme();
stoptime();

showmeall(1);

Optimize('hgnc_aliases');

done();
