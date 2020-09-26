#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

our $usage = "$0";
# ($var) = args(1);
args(0);

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

$mainquery = Query("SELECT symbol, ensembl_gene_id FROM hgnc ORDER BY symbol");
startme("Checking if the ENSGs mentioned in 'hgnc' match what the 'ensembl' and 'gencode_gff3_gene' tables say", 0, Numrows($mainquery));
starttime();
while (($symbol, $ensg) = Fetch($mainquery))
{
	stepme(100);
	
	addme("total symbols", $symbol);
	
	# Check HGNC's ENSG
	if (!defined($ensg))
	{
		addme("hgnc has a blank ensg field for symbol (skipped)", $symbol);
		next;
	}

	addme("total symbols with ensgs", $symbol);
	addme("total ensgs", $ensg);

	# Check ensembl
	$query = Query("SELECT DISTINCT hgnc FROM ensembl WHERE ensg='$ensg' AND hgnc IS NOT NULL");
	if (Numrows($query) == 0)
	{
		addme("'ensembl' table doesn't have a row for ensg (skipped)", $ensg);
		next;
	}
	($ensembl_symbol) = FetchOne($query);
	
	if ($symbol ne $ensembl_symbol)
	{
		# die("Error: Symbol in 'ensembl' table ('$ensembl_symbol') doesn't match HGNC's symbol ('$symbol')");
		addme("mismatch: symbol in 'ensembl' table vs. hgnc symbol for hgnc symbol", $symbol);
		addme("mismatch: symbol in 'ensembl' table vs. hgnc symbol for hgnc symbol|ensembl symbol", "$symbol|$ensembl_symbol");
	}
	else
	{
		addme("match: symbol in 'ensembl' table vs. hgnc symbol for hgnc symbol", $symbol);
	}
	
	# Check ensembl_gff3_gene
	$query = Query("SELECT DISTINCT symbol FROM ensembl_gff3_gene WHERE species='human' AND ensg='$ensg' AND symbol IS NOT NULL");
	if (Numrows($query) > 0)
	{
		($ensembl_gff3_symbol) = FetchOne($query);
	
		if ($symbol ne $ensembl_gff3_symbol)
		{
			# die("Error: Symbol in 'ensembl_gff3_gene' table ('$ensembl_gff3_symbol') doesn't match HGNC's symbol ('$symbol')");
			addme("mismatch: symbol in 'ensembl_gff3_gene' table vs. hgnc symbol for hgnc symbol", $symbol);
			addme("mismatch: symbol in 'ensembl_gff3_gene' table vs. hgnc symbol for hgnc symbol|ensembl gff3 symbol", "$symbol|$ensembl_gff3_symbol");
		}
		else
		{
			addme("match: symbol in 'ensembl_gff3_gene' table vs. hgnc symbol for hgnc symbol", $symbol);
		}
	}
	else
	{
		addme("no entry in 'ensembl_gff3_gene' table for ensg", $ensg);
	}

	# Check gencode_gff3_gene
	$query = Query("SELECT DISTINCT symbol FROM gencode_gff3_gene WHERE species='human' AND ensg='$ensg' AND symbol IS NOT NULL");
	if (Numrows($query) > 0)
	{
		($gencode_symbol) = FetchOne($query);
	
		if ($symbol ne $gencode_symbol)
		{
			# die("Error: Symbol in 'gencode_gff3_gene' table ('$gencode_symbol') doesn't match HGNC's symbol ('$symbol')");
			addme("mismatch: symbol in 'gencode_gff3_gene' table vs. hgnc symbol for hgnc symbol", $symbol);
			addme("mismatch: symbol in 'gencode_gff3_gene' table vs. hgnc symbol for hgnc symbol|gencode symbol", "$symbol|$gencode_symbol");
		}
		else
		{
			addme("match: symbol in 'gencode_gff3_gene' table vs. hgnc symbol for hgnc symbol", $symbol);
		}
	}
	else
	{
		addme("no entry in 'gencode_gff3_gene' table for ensg", $ensg);
	}
}
stopme();
stoptime();

# showmeall(1);
showmesome(100);

done();
