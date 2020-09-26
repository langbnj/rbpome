#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $myspecies = "'HUMAN', 'MOUSE', 'DANRE'";

our $usage = "$0";
# ($var) = args(1);
args(0);

$infile = "../../update/rbpome/input-aliases.txt";
# $outfile = "output.txt";
#
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");




# Gene alias mapping (symbol1 and symbol2 will be replaced by these in the 'rbpome' table)
%alias = ();

startme("Reading gene symbol alias mapping from '$infile'");
while (<IN>)
{
	# # Alias definitions
	#
	# # I'm going to follow UniProt, not HGNC. HGNC still prefers "C1orf35" over MMTAG2, even though it lists MMTAG2 as a synonym, which is really strange.
	#
	#
	#
	#
	#
	#
	#
	# # HUMAN
	# # Used this type of query to get these: https://www.uniprot.org/uniprot/?query=DBC1+AND+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
	#
	#
	#
	# HUMAN	ATP5B	ATP5F1B
	# HUMAN	ATP5C1	ATP5F1C
	
	chomp;
	
	# Skip empty lines
	next if /^$/;
	
	# Skip comments
	next if /^#/;
	
	@a = split(/\t/);
	die("Error: Line doesn't contain 3 fields:\n\n$_\n\n") if (scalar(@a) != 3);
	
	($species, $original, $alias) = @a;
	# $original = $a[1];
	# $alias = $a[2];
	
	# Convert to upper case
	$species = uc($species);
	$original = uc($original);
	$alias = uc($alias);
	
	$alias{"$species|$original"} = $alias;
	# # Discard species
	# $alias{$original} = $alias;
	
	stepme(10);
}
stopme();




# Manual ortholog fixes for the MGI homologene table (assigned from tmp-problematic-symbols.txt)
# delete($alias{"HUMAN|ATP5B"});
# $alias{"MOUSE|ATP5A1"} = "ATP5A1";





# Inverse %alias
startme("Creating hash of symbols that are the result of aliasing (aliases)");
%originals = ();
# foreach $original (keys(%alias))
foreach $tmp (keys(%alias))
{
	($species, $original) = split(/\|/, $tmp);
	$alias = $alias{"$species|$original"};
	
	if (exists($originals{"$species|$alias"}))
	{
		$originals{"$species|$alias"} .= "|$original";
	}
	else
	{
		$originals{"$species|$alias"} = $original;
	}

	# print "ORIGINAL $species|$alias = $original\n";
	# $is_alias{"$species|$alias"} = 1;
	# $original{$alias} = $original;
	
	addme("total species", $species);

	addme("total original species|symbols", "$species|$original");
	addme("total symbol species|aliases", "$species|$alias");

	addme("total original symbols", $original);
	addme("total symbol aliases", $alias);
	
	stepme(10);
}
stopme();







# start

startme("Checking how my manual gene symbol alias assignment affected the matchup between human, mouse and DANRE gene symbols in table 'rbpome'");
starttime();
# $mainquery = Query("SELECT symbol1, symbol2 FROM rbpome GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
# while (($symbol1, $symbol2) = Fetch($mainquery))
# $query = Query("SELECT symbol1, symbol2 FROM rbpome");
$query = Query("SELECT symbol1, symbol2, species1, species2 FROM rbpome");
@symbols  = ();
# while (($symbol1, $symbol2) = Fetch($query))
while (($symbol1, $symbol2, $species1, $species2) = Fetch($query))
{
	push(@symbols, "$species1|$symbol1");
	push(@symbols, "$species2|$symbol2");
	# push(@symbols, $symbol1);
	# push(@symbols, $symbol2);
	
	# print " >> $species1|$symbol1\n";
	# print " >> $species2|$symbol2\n";

	stepme(1000);
}
stopme();
@symbols = unique(@symbols);
# show(@symbols);

state("Got ".commify(scalar(@symbols))." unique species|symbols");
# state("Got ".scalar(@symbols)." unique symbols");




# Cycle through symbols
startme("Cycle through symbols, check if they fall into the same homologene group across species, and check if the alias assignment caused a difference for them");
foreach $tmp (@symbols)
{
	($species, $symbol) = split(/\|/, $tmp);
	
	# print " >> $species|$symbol\n";
	# if (exists($is_alias{"$species|$symbol"}))
	
	addme("total rbpome symbols", $symbol);
	addme("total rbpome species|symbols", "$species|$symbol");
	
	
	$aliased = 0;
	if (exists($originals{$symbol}))
	{
		$aliased = 1;
		$alias = $symbol;

		addme("total rbpome symbols that are aliases", $alias);
		
		# If the symbol has been aliased from something else, get its original symbol(s) from the fullsymbol field
		@fullsymbols = ();
		$query = Query("SELECT DISTINCT fullsymbol1 FROM rbpome WHERE symbol1='$alias'");
		while (($fullsymbol) = Fetch($query))
		{
			$fullsymbol = uc($fullsymbol);
			push(@fullsymbols, $fullsymbol);
		}
		$query = Query("SELECT DISTINCT fullsymbol2 FROM rbpome WHERE symbol2='$alias'");
		while (($fullsymbol) = Fetch($query))
		{
			$fullsymbol = uc($fullsymbol);
			push(@fullsymbols, $fullsymbol);
		}
		@fullsymbols = unique(@fullsymbols);
		
		# As an alternative, get the original symbols from %originals
		@originalsymbols = $originals{"$species|$alias"};
		@originalsymbols = unique(@originalsymbols);

		# Compare @fullsymbols and @originalsymbols
		die("Error: fullsymbols and originalsymbols don't match:\n\nfullsymbols: ".join("|", @fullsymbols)."\noriginalsymbols: ".join("|", @originalsymbols)."\n\n") if (join("|", @fullsymbols) ne join("|", @originalsymbols));
		exit;
		
		# Get the alns the original symbols fall into
		# $query =
		
		
		$query = Query("SELECT symbol1, symbol2, fullsymbol1, fullsymbol2 FROM rbpome WHERE symbol1='$alias' OR symbol2='$alias'");
		while (($symbol1, $symbol2, $fullsymbol1, $fullsymbol2) = Fetch($query))
		{
			if ($symbol1 eq $alias)
			{
				$original = uc($fullsymbol1);
			}
			if ($symbol2 eq $alias)
			{
				
				$original = uc($fullsymbol2);
			}
		}
		$original = $original{"$species|$alias"};
		
		# See what homologene cluster the symbol was in originally
		$query = Query("SELECT aln FROM homologene WHERE species='$species' AND symbol='$original'");
		($orialn) = FetchOne($query);
		
		# See which cluster the symbol got put in after assigning it the alias
		$query = Query("SELECT aln FROM homologene WHERE species='$species' AND symbol='$alias'");
		($alialn) = FetchOne($query);
		
		print " >> $orialn >> $alialn >> ";
		print "DIFFERENT\n" if ($orialn ne $alialn);
		print "SAME\n" if ($orialn eq $alialn);
		print "\n";
	}
	else
	{
		
	}
	
	
	$query = Query("SELECT species1 FROM rbpome WHERE symbol1='$symbol'");
	@tmpspecies = ();
	while (($tmpspecies) = Fetch($query))
	{
		push(@tmpspecies, $tmpspecies);
	}
	$query = Query("SELECT species2 FROM rbpome WHERE symbol2='$symbol'");
	while (($tmpspecies) = Fetch($query))
	{
		push(@tmpspecies, $tmpspecies);
	}
	@tmpspecies = unique(@tmpspecies);
	$tmpspecies = join("', '", @tmpspecies);
	$tmpspecies = "'".$tmpspecies."'";
	
	
	# Check if all species have the same homologene group for this symbol
	$query = Query("SELECT COUNT(DISTINCT aln), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR '|') FROM homologene WHERE symbol='$symbol' AND species IN ($tmpspecies)");
	($alncount, $speciescount, $specieslist) = FetchOne($query);
	if (switch('debug'))
	{
		print " >> $species|$symbol\n" if (switch('debug'));
		print "   >> $alncount alns\n" if (switch('debug'));
		print "   >> $speciescount species: $specieslist\n" if (switch('debug'));
	}
	if ($alncount == 0)
	{
		# die("Error: No alns for symbol '$symbol' in table 'homologene'");
		
		# Check if they are present for multiple species in the screen
		if (scalar(@tmpspecies) == 1)
		{
			addme("OK: no alns in table homologene for symbol, but we only have it for one species in table 'rbpome'", $symbol);
		}
		else
		{
			addme("BAD: no alns in table homologene for symbol and we have it in multiple species in table 'rbpome'", $symbol);
		}
	}
	elsif ($alncount == 1)
	{
		addme("GOOD: exactly one aln for $speciescount species for symbol", $symbol);
	}
	else
	{
		addme("GOOD: multiple alns for $speciescount species for symbol", $symbol);
	}
	
	
	stepme(100);
}
stopme();
stoptime();

# showme("BAD: no alns in table homologene for symbol and we have it in multiple species in table 'rbpome'");
showmeall(1);

done();
