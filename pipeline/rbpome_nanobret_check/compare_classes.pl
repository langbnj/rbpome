#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $table = 'rbpome';

# eCLIP parameters
$type = 'eclip_encode_12';
$hc = 0;
$minlog2fold = 0;
$mindist_threshold = 50;
$resamples = 100;

our $usage = "$0 [table (e.g. rbpome)] [score type] [threshold]\n\nExample: $0 rbpome SUM_IS 7.1";
($table, $scoretype, $threshold) = args(3);

$aliasfile = "../../update/rbpome/input-aliases.txt";
$speciesfile = "../../update/rbpome/input/symbol_to_species.txt";
$infile = "input/merged_nanobret_results_07-Feb-2020221237_ben_all.txt";
$classfile = "input/merged_refs_Silvia_03022020.txt";
$outfile = "output-$table-compare-classes.txt";
open(SPECIES, $speciesfile) or die("\nError: Couldn't open '$speciesfile'\n\n");
open(CLASS, $classfile) or die("\nError: Couldn't open '$classfile'\n\n");
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

print OUT "class\tpositives\tnegatives\ttested\tvalidation_ratio\tTP\tFP\tTN\tFN\tF1\tMCC\tpositive_list\tnegative_list\teclip_TP\teclip_FP\teclip_TN\teclip_FN\teclip_F1\teclip_MCC\teclip_pairs\tavg_jaccard\tavg_oddsratio\tavg_probability_mindist\tmedian_jaccard\tmedian_oddsratio\tmedian_probability_mindist\n";
# print OUT "class\tpositives\tpositive list\tnegatives\tnegative list\ttested\tvalidation_ratio\tTP\tFP\tTN\tFN\tF1\tMCC\teclip_TP\teclip_FP\teclip_TN\teclip_FN\teclip_F1\teclip_MCC\teclip_pairs\tavg_jaccard\tavg_oddsratio\tavg_probability_mindist\tmedian_jaccard\tmedian_oddsratio\tmedian_probability_mindist\n";


# start






# Read species mapping for gene symbols
startme("Reading symbol-to-species mapping from '$speciesfile'");
%species = ();
%realsymbol = ();	# with capitalisation
while (<SPECIES>)
{
	chomp;
	
	# 2700060E02Rik	M.musculus
	# A1CF	H.sapiens
	# ABCF2	H.sapiens
	
	@a = split(/\t/);
	
	$symbol = $a[0];
	$realsymbol = $symbol;	# With capitalisation etc.
	
	if (scalar(@a) == 1)
	{
		# Sebastian says to just skip all these. They're positive and negative controls, p53 is a control, and Lgals3 is a lectin, i.e. an extracellular carbohydrate binder. I'll exclude them from the output table.
		
		# Assigning some species manually (based on what they interact with)
		# GFP		# Jellyfish! UniProt species mnemonic 'AEQVI'.
		# Kif17v1	# human?
		# lam		# human?
		# larget	# Not in the screen (but assuming it's human)
		# Lgals3	# human?
		# p53		# Not in the screen (but assuming it's human)
	
		# $symbol = $a[0];
		# if ($symbol eq 'GFP')		{ $species = 'jellyfish'; }	# Interacts with exosc7 & fxr2, both H.sapiens
		# elsif ($symbol eq 'Kif17v1'){ $species = 'H.sapiens'; }	# Interacts with exosc8 & snrpa, both H.sapiens
		# elsif ($symbol eq 'lam')	{ $species = 'H.sapiens'; }	# Interacts with exosc7 & fxr2, both H.sapiens
		# elsif ($symbol eq 'larget')	{ $species = 'H.sapiens'; }	# Assuming they've used human. Not in the screen, doesn't matter
		# elsif ($symbol eq 'p53')	{ $species = 'H.sapiens'; }	# Assuming they've used human. Not in the screen, doesn't matter
		# # Lgals3:
		# # BAG4	H.sapiens
		# # BYSL	H.sapiens
		# # Pef1	M.musculus
		# # YTHDF3	H.sapiens
		# # ...so not sure. Assigning human for now.
		# elsif ($symbol eq 'Lgals3')	{ $species = 'H.sapiens'; }
		# else
		# {
			# addme("species mapping: no species assignment for symbol (skipped)", $symbol);
			next;
		# }
	}
	# Actually it's just a space
	# elsif (scalar(@a) == 3)
	# {
	# 	# RBM10 316F	H.sapiens
	# 	# RBM10 343G	H.sapiens
	# 	#
	# 	# These turn into "rbm10316f" / "rbm10343g" in the screen.
	#
	# 	$symbol = $a[0].$a[1];
	# 	$species = $a[2];
	#
	# 	addme("combined two columns for symbol")
	# }
	elsif (scalar(@a) != 2)
	{
		die("Error: Weird line:\n$_")
		# addme("species mapping: no species assignment in line (skipped)", $_);
		# next;
	}
	else
	{
		$species = $a[1];
	}
	
	# # Make symbol all upper case
	# $symbol = uc($symbol);
	
	# Make species all upper case
	$species = uc($species);
	
	# # Replace spaces in "RBM10 316F" >> "rbm10_316f", "eMIC SRRM3 mut 25" >> "emic_srrm3_mut_25" etc.
	# $symbol =~ s/ /_/g;
	# # Remove spaces in "RBM10 316F" >> "rbm10316f", "eMIC SRRM3 mut 25" >> "emicsrrm3mut25" etc. (as in Sebastian's table)
	# $symbol =~ s/ //g;
	
	# # Somehow "2700060E02Rik" turns into >> "x2700060e02rik" in the screen, add the x
	# # This is the only RIKEN ID.
	# $symbol = 'x'.$symbol if ($symbol eq '2700060e02rik');
	# if ($symbol eq '2700060e02rik')	{ $species = 'm.musculus'; $symbol = "x2700060e02rik"; }
	# die("Error: RIKEN ID '$symbol' needs to be handled") if ($symbol =~ /^\d+.+rik$/i);
	
	# if ($species eq 'H.SAPIENS')
	# {
	# 	$species = 'HUMAN';
	# }
	# elsif ($species eq 'M.MUSCULUS')
	# {
	# 	$species = 'MOUSE';
	# }
	# elsif ($species eq 'ZEBRAFISH')
	# {
	# 	$species = 'DANRE';
	# }
	# elsif ($species eq 'JELLYFISH')
	# {
	# 	$species = 'AEQVI';
	# }
	# else
	# {
	# 	# die("Error: Unhandled species '$species'");
	# 	addme("species mapping: skipped unhandled species for species", $species);
	# }
	
	if (exists($species{$symbol}))
	{
		die("Error: Duplicate species assignment for symbol '$symbol'");
	}
	$species{$symbol} = $species;
	$realsymbol{$symbol} = $realsymbol;
	
	stepme(100);
}
stopme();










# Gene alias mapping (symbol1 and symbol2 will be replaced by these in the 'rbpome' table)
%alias = ();

if (!switch('noalias'))
{
	open(ALIAS, $aliasfile) or die("\nError: Couldn't open '$aliasfile'\n\n");


	startme("Reading gene symbol alias mapping from '$aliasfile'");
	while (<ALIAS>)
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
	close(ALIAS);
}






# Get classes
startme("Reading pair class assignments from '$classfile'");
<CLASS>;	# Skip header
@classes = ();
%class = ();
%pairs_per_class = ();
%pair = ();
%pairids = ();
while (<CLASS>)
{
	chomp;
	
	# pair_no	gene in Halo vector	gene in Luc vector	interaction	processed plate
	# J103_01	csrp1	kif7	rbp_map
	# J103_02	kif16b_cter	ncbp3	rbp_map
	# J103_03	pura	kif3b	rbp_map
	# J103_04	TSNAX	Bicd2	rbp_map
	# J103_05	TSNAX	Kif6	rbp_map
	# J103_06	TSNAX	Kif2a	rbp_map
	
	@a = split(/\t/);

	$pairid = $a[0];

	$symbol1 = $a[1];
	$symbol2 = $a[2];

	$class = $a[3];
	
	$originalsymbol1 = $symbol1;
	$originalsymbol2 = $symbol2;
	
	
	
	
	# # Get species from table 'rbpome'
	#
	# # Try for symbol1 via fullsymbol1
	# $query = Query("SELECT DISTINCT species1 FROM $table WHERE fullsymbol1='$symbol1'");
	# die if (Numrows($query) > 1);
	# if (Numrows($query) == 0)
	# {
	# 	# Try for symbol1 via fullsymbol2
	# 	$query = Query("SELECT DISTINCT species2 FROM $table WHERE fullsymbol2='$symbol1'");
	# 	die if (Numrows($query) > 1);
	# 	if (Numrows($query) == 0)
	# 	{
	# 		# die("Error: Couldn't find fullsymbol '$symbol1' in table '$table'")
	# 		addme("classes: couldn't find fullsymbol in table '$table' for fullsymbol", $symbol1);
	# 	}
	# }
	# ($species1) = FetchOne($query);
	#
	# # Try for symbol2 via fullsymbol2
	# $query = Query("SELECT DISTINCT species1 FROM $table WHERE fullsymbol2='$symbol2'");
	# die if (Numrows($query) > 1);
	# if (Numrows($query) == 0)
	# {
	# 	# Try for symbol2 via fullsymbol2
	# 	$query = Query("SELECT DISTINCT species2 FROM $table WHERE fullsymbol2='$symbol2'");
	# 	die if (Numrows($query) > 1);
	# 	if (Numrows($query) == 0)
	# 	{
	# 		# die("Error: Couldn't find fullsymbol '$symbol2' in table '$table'")
	# 		addme("classes: couldn't find fullsymbol in table '$table' while getting species for fullsymbol", $symbol2);
	# 	}
	# }
	# ($species1) = FetchOne($query);
	


	
	
	# Make symbols all upper case
	$symbol1 = uc($symbol1);
	$symbol2 = uc($symbol2);
	



	if (!exists($species{$symbol1}))
	{
		# die("Error: No species assignment for symbol '$symbol1'");
		# addme("no species assignment for symbol pair (skipped)", "$symbol1|$symbol2");
		# addme("no species assignment for symbol (skipped)", $symbol1);
		# next;
		addme("no species assignment for symbol pair (kept)", "$symbol1|$symbol2");
		addme("no species assignment for symbol (kept)", $symbol1);
		$species1 = 'N/A';
	}
	else
	{
		$species1 = $species{$symbol1};
	}
	if (!exists($species{$symbol2}))
	{
		# die("Error: No species assignment for symbol '$symbol2'");
		# addme("no species assignment for symbol pair (skipped)", "$symbol1|$symbol2");
		# addme("no species assignment for symbol (skipped)", $symbol2);
		# next;
		addme("no species assignment for symbol pair (kept)", "$symbol1|$symbol2");
		addme("no species assignment for symbol (kept)", $symbol2);
		$species2 = 'N/A';
	}
	else
	{
		$species2 = $species{$symbol2};
	}
	
	addme("total species", $species1);
	addme("total species", $species2);
	
	# $same_species = 0;
	# $same_species = 1 if ($species1 eq $species2);




	
	# Retrieve original symbol (with capitalisation etc.) from the species mapping table
	# symbol1
	if (!exists($realsymbol{$symbol1}))
	{
		# If the realsymbol isn't set, the species shouldn't be set either
		die if (exists($species{$symbol1}));
		$fullsymbol1 = $originalsymbol1;
	}
	else
	{
		$fullsymbol1 = $realsymbol{$symbol1};

		# If the realsymbol is set, the species should be set too
		die if (!exists($species{$symbol1}));

		# The fullsymbol should always be the same as the originalsymbol
		# die("Error: fullsymbol '$fullsymbol1' from the species mapping file doesn't match the originalsymbol '$originalsymbol1'") if ($fullsymbol1 ne $originalsymbol1);
		if (uc($originalsymbol1) eq uc($fullsymbol1))
		{
			addme("case differs in the class assignment file '$classfile' for class fullsymbol | $table table fullsymbol (kept)", "$fullsymbol1|$originalsymbol1");
		}
		else
		{
			die("Error: Got '$originalsymbol1' instead of '$fullsymbol1' for fullsymbol1 for '$fullsymbol1'");
			# warn("Warning: Got '$originalsymbol1' instead of '$fullsymbol1' for fullsymbol1 for '$fullsymbol1'");
		}
	}
	# symbol2
	if (!exists($realsymbol{$symbol2}))
	{
		# If the realsymbol isn't set, the species shouldn't be set either
		die if (exists($species{$symbol2}));
		$fullsymbol2 = $originalsymbol2;
	}
	else
	{
		$fullsymbol2 = $realsymbol{$symbol2};

		# If the realsymbol is set, the species should be set too
		die if (!exists($species{$symbol2}));

		# The fullsymbol should always be the same as the originalsymbol
		# die("Error: fullsymbol '$fullsymbol2' from the species mapping file doesn't match the originalsymbol '$originalsymbol2'") if ($fullsymbol2 ne $originalsymbol2);
		if (uc($originalsymbol2) eq uc($fullsymbol2))
		{
			addme("case differs in the class assignment file '$classfile' for class fullsymbol | $table table fullsymbol (kept)", "$fullsymbol2|$originalsymbol2");
		}
		else
		{
			die("Error: Got '$originalsymbol2' instead of '$fullsymbol2' for fullsymbol2 for '$fullsymbol2'");
			# warn("Warning: Got '$originalsymbol2' instead of '$fullsymbol2' for fullsymbol2 for '$fullsymbol2'");
		}
	}
	
	
	
	
	# Make symbols all upper case
	$symbol1 = uc($symbol1);
	$symbol2 = uc($symbol2);




	# Apply alias mapping (e.g. ATP5B protein to ATP5F1B, which is its gene)
	# if (exists($alias{$symbol1}))
	if (exists($alias{"$species1|$symbol1"}))
	{
		addme("aliased from species|symbol", "$species1|$symbol1");
		# $symbol1 = $alias{$symbol1};
		$symbol1 = $alias{"$species1|$symbol1"};
		addme("aliased to species|symbol", "$species1|$symbol1");
	}
	# if (exists($alias{$symbol2}))
	if (exists($alias{"$species2|$symbol2"}))
	{
		addme("aliased from species|symbol", "$species2|$symbol2");
		# $symbol2 = $alias{$symbol2};
		$symbol2 = $alias{"$species2|$symbol2"};
		addme("aliased to species|symbol", "$species2|$symbol2");
	}
	
	# Check for pairs that are present more than once
	if (contains("$symbol1|$symbol2", returnme("classes: total symbol pairs before flipping")))
	{
		addme("classes: duplicate symbol pair before flipping for pair id", $pairid);
		addme("classes: duplicate symbol pair before flipping for symbol pair", "$symbol1|$symbol2");
		addme("classes: duplicate symbol pair before flipping for fullsymbol pair", "$fullsymbol1|$fullsymbol2");
	}

	addme("classes: total symbol pairs before flipping", "$symbol1|$symbol2");

	# # Flip symbols if required
	# die("NOTE: THIS IS NOT YET ADAPTED TO UNFLIPPED PAIRS");
	# #DEBUG
	# if ($symbol1 gt $symbol2)
	# {
	# 	addme("classes: pair flipped for pairid", $pairid);
	# 	addme("classes: pair flipped for original symbol1|symbol2", "$symbol1|$symbol2");
	# 	addme("classes: pair flipped for original fullsymbol1|fullsymbol2", "$fullsymbol1|$fullsymbol2");
	# 	($symbol1, $symbol2) = ($symbol2, $symbol1);
	# 	($fullsymbol1, $fullsymbol2) = ($fullsymbol2, $fullsymbol1);
	# 	($originalsymbol1, $originalsymbol2) = ($originalsymbol2, $originalsymbol1);
	# 	($species1, $species2) = ($species2, $species1);
	# }



	addme("classes: total classes", $class);
	addme("classes: total pair ids", $pairid);
	
	# Check for pairs that are present more than once
	if (contains("$symbol1|$symbol2", returnme("classes: total symbol pairs")))
	{
		addme("classes: duplicate symbol pair after flipping for pair id", $pairid);
		addme("classes: duplicate symbol pair after flipping for symbol pair", "$symbol1|$symbol2");
		addme("classes: duplicate symbol pair after flipping for fullsymbol pair", "$fullsymbol1|$fullsymbol2");
	}

	addme("classes: total symbol pairs", "$symbol1|$symbol2");
	addme("classes: total fullsymbol pairs", "$fullsymbol1|$fullsymbol2");
	addme("classes: total symbols", $symbol1);
	addme("classes: total symbols", $symbol2);
	addme("classes: total fullsymbols", $fullsymbol1);
	addme("classes: total fullsymbols", $fullsymbol2);
	
	
	# Now that I've gone to all this huge effort to get:
	# symbol1
	# symbol2
	# fullsymbol1
	# fullsymbol2
	# species1
	# species2
	# ...let's check if this matches up with the 'rbpome' table, for symbol1:
	# our $superloudmysql=1;
	# $query = Query("SELECT DISTINCT fullsymbol1, symbol1, species1 FROM $table WHERE fullsymbol1='$originalsymbol1' COLLATE latin1_general_cs");	# Case-sensitive query
	$query = Query("SELECT DISTINCT fullsymbol1, symbol1, species1 FROM $table WHERE fullsymbol1='$originalsymbol1'");
	if (Numrows($query) == 0)
	{
		# Try other column
		# $query = Query("SELECT DISTINCT fullsymbol2, symbol2, species2 FROM $table WHERE fullsymbol2='$originalsymbol1' COLLATE latin1_general_cs");	# Case-sensitive query
		$query = Query("SELECT DISTINCT fullsymbol2, symbol2, species2 FROM $table WHERE fullsymbol2='$originalsymbol1'");
		if (Numrows($query) == 0)
		{
			warn("Warning: fullsymbol not found in table '$table' for fullsymbol '$originalsymbol1' (skipped)");
			addme("fullsymbol not found in table '$table' for fullsymbol (skipped)", $originalsymbol1);
			next;
		}
	}
	($tmpfullsymbol1, $tmpsymbol1, $tmpspecies1) = FetchOne($query);
	if ($tmpfullsymbol1 ne $fullsymbol1)
	{
		if (uc($tmpfullsymbol1) eq uc($fullsymbol1))
		{
			addme("case differs in the class assignment file '$classfile' for class fullsymbol | $table table fullsymbol (kept)", "$fullsymbol1|$tmpfullsymbol1");
		}
		else
		{
			die("Error: Got '$tmpfullsymbol1' instead of '$fullsymbol1' for fullsymbol1 for '$fullsymbol1'");
			# warn("Warning: Got '$tmpfullsymbol1' instead of '$fullsymbol1' for fullsymbol1 for '$fullsymbol1'");
		}
	}
	# die("Error: Got '$tmpsymbol1' instead of '$symbol1' for symbol1 for '$fullsymbol1'") if ($tmpsymbol1 ne $symbol1);
	warn("Warning: Found '$tmpsymbol1' in table '$table' instead of '$symbol1' for symbol1 for '$fullsymbol1' (using '$tmpsymbol1')") if ($tmpsymbol1 ne $symbol1);
	$symbol1 = $tmpsymbol1;
	if ($tmpspecies1 ne $species1)
	{
		if ($species1 eq 'N/A')
		{
			# Use rbpome table species
			$species1 = $tmpspecies1;
		}
		else
		{
			die("Error: Got '$tmpspecies1' instead of '$species1' for species1 for '$fullsymbol1'") ;
		}
	}
	# ...let's check if this matches up with the 'rbpome' table, for symbol2:
	# $query = Query("SELECT DISTINCT fullsymbol1, symbol1, species1 FROM $table WHERE fullsymbol1='$originalsymbol2' COLLATE latin1_general_cs");	# Case-sensitive query
	$query = Query("SELECT DISTINCT fullsymbol1, symbol1, species1 FROM $table WHERE fullsymbol1='$originalsymbol2'");
	if (Numrows($query) == 0)
	{
		# Try other column
		# $query = Query("SELECT DISTINCT fullsymbol2, symbol2, species2 FROM $table WHERE fullsymbol2='$originalsymbol2' COLLATE latin1_general_cs");	# Case-sensitive query
		$query = Query("SELECT DISTINCT fullsymbol2, symbol2, species2 FROM $table WHERE fullsymbol2='$originalsymbol2'");
		if (Numrows($query) == 0)
		{
			warn("Warning: fullsymbol not found in table '$table' for fullsymbol '$originalsymbol2' (skipped)");
			addme("fullsymbol not found in table '$table' for fullsymbol (skipped)", $originalsymbol2);
			next;
		}
	}
	($tmpfullsymbol2, $tmpsymbol2, $tmpspecies2) = FetchOne($query);
	if ($tmpfullsymbol2 ne $fullsymbol2)
	{
		if (uc($tmpfullsymbol2) eq uc($fullsymbol2))
		{
			addme("case differs in the class assignment file '$classfile' for class fullsymbol | $table table fullsymbol (kept)", "$fullsymbol2|$tmpfullsymbol2");
		}
		else
		{
			die("Error: Got '$tmpfullsymbol2' instead of '$fullsymbol2' for fullsymbol2 for '$fullsymbol2'");
		}
	}
	# die("Error: Got '$tmpsymbol2' instead of '$symbol2' for symbol2 for '$fullsymbol2'") if ($tmpsymbol2 ne $symbol2);
	warn("Warning: Found '$tmpsymbol2' in table '$table' instead of '$symbol2' for symbol2 for '$fullsymbol2' (using '$tmpsymbol2')") if ($tmpsymbol2 ne $symbol2);
	$symbol2 = $tmpsymbol2;
	if ($tmpspecies2 ne $species2)
	{
		if ($species2 eq 'N/A')
		{
			# Use rbpome table species
			$species2 = $tmpspecies2;
		}
		else
		{
			die("Error: Got '$tmpspecies2' instead of '$species2' for species2 for '$fullsymbol2'") ;
		}
	}
	

	# Main assignment
	push(@classes, $class);
	$class{$pairid} = $class;
	$pair{$pairid} = "$symbol1|$symbol2";
	
	if (exists($pairids{"$symbol1|$symbol2"}))
	{
		$pairids{"$symbol1|$symbol2"} .= "|$pairid";
	}
	else
	{
		$pairids{"$symbol1|$symbol2"} = $pairid;
	}
	
	if (exists($pairs_per_class{$class}))
	{
		$pairs_per_class{$class} .= ";$symbol1|$symbol2";
	}
	else
	{
		$pairs_per_class{$class} = "$symbol1|$symbol2";
	}
	
	stepme(1);
}
stopme();

@classes = unique(@classes);

# showmesome(50);
# show(%class);
# show(%pair);
# show(%pairids);
# exit;
# >> OK, lots of pairs get tested more than once, but usually not in the same orientation.




# Get fullsymbol-symbol mapping
starttime();
$query = Query("SELECT DISTINCT fullsymbol1, symbol1, eclip1 FROM $table");
startme("Getting fullsymbol-symbol mapping from table '$table' (step 1/2)", 0, Numrows($query));
%symbol = ();
%eclip = ();
while (($fullsymbol, $symbol, $eclip) = Fetch($query))
{
	$fullsymbol = uc($fullsymbol);
	$symbol = uc($symbol);
	
	if (exists($symbol{$fullsymbol}))
	{
		if ($symbol{$fullsymbol} ne $symbol)
		{
			die("Error: Multiple symbols for fullsymbol '$fullsymbol': '$symbol{$fullsymbol}' and '$symbol'");
		}
	}

	$symbol{$fullsymbol} = $symbol;
	$eclip{$symbol} = 1 if ($eclip == 1);

	stepme(100);
}
stopme();

$query = Query("SELECT DISTINCT fullsymbol2, symbol2, eclip2 FROM $table");
startme("Getting fullsymbol-symbol mapping from table '$table' (step 2/2)", 0, Numrows($query));
while (($fullsymbol, $symbol, $eclip) = Fetch($query))
{
	$fullsymbol = uc($fullsymbol);
	$symbol = uc($symbol);

	if (exists($symbol{$fullsymbol}))
	{
		if ($symbol{$fullsymbol} ne $symbol)
		{
			die("Error: Multiple symbols for fullsymbol '$fullsymbol': '$symbol{$fullsymbol}' and '$symbol'");
		}
	}

	$symbol{$fullsymbol} = $symbol;
	$eclip{$symbol} = 1 if ($eclip == 1);

	stepme(100);
}
stopme();

state("Got mapping for ".commify(scalar(keys(%symbol)))." fullsymbols");
state("Got ".scalar(keys(%eclip))." eCLIP symbols");




# Getting NanoBRET data
%tested = ();
%positive = ();
%negative = ();
startme("Reading '$infile'");
<IN>;	# Skip header
while (<IN>)
{
	# pair_ID	p1	p2	nb_score	nb_std	merged_pairnames	classification_zscore095	Index
	# R_59	cpeb2	cpeb2	47.78505355	1.537508115	cpeb2_cpeb2	positive	1
	# J103_02	kif16b_cter	ncbp3	0.951010652	0.327433707	kif16b_cter_ncbp3	negative	2
	# R_48	srsf9	rbmx	39.77018046	1.315552926	srsf9_rbmx	positive	3
	# J103_04	tsnax	bicd2	1.247319423	0.282866362	tsnax_bicd2	negative	4
	
	@a = split(/\t/);
	$pairid = $a[0];
	$fullsymbol1 = $a[1];
	$fullsymbol2 = $a[2];
	$verdict = $a[6];

	# Skip homomers
	if ($fullsymbol1 eq $fullsymbol2)
	{
		addme("homomer skipped for fullsymbols", "$fullsymbol1|$fullsymbol2");
		next;
	}

	$fullsymbol1 = uc($fullsymbol1);
	$fullsymbol2 = uc($fullsymbol2);
	
	if (!exists($symbol{$fullsymbol1}))
	{
		# die("Error: No symbol for fullsymbol '$fullsymbol1'");
		addme("nanobret pair skipped because no symbol for fullsymbol (this also indicates the protein wasn't in the screen) (skipped)", $fullsymbol1);
		addme("nanobret pair skipped because no symbol for fullsymbol for pair (this also indicates the protein wasn't in the screen) (skipped)", "$fullsymbol1|$fullsymbol2");
		next;
	}
	if (!exists($symbol{$fullsymbol2}))
	{
		# die("Error: No symbol for fullsymbol '$fullsymbol2'");
		addme("nanobret pair skipped because no symbol for fullsymbol (this also indicates the protein wasn't in the screen) (skipped)", $fullsymbol2);
		addme("nanobret pair skipped because no symbol for fullsymbol for pair (this also indicates the protein wasn't in the screen) (skipped)", "$fullsymbol1|$fullsymbol2");
		next;
	}
	
	# Translate fullsymbol into symbol
	$symbol1 = $symbol{$fullsymbol1};
	$symbol2 = $symbol{$fullsymbol2};
	
	# # Flip if necessary
	# if ($symbol1 gt $symbol2)
	# {
	# 	($symbol1, $symbol2) = ($symbol2, $symbol1);
	# 	($fullsymbol1, $fullsymbol2) = ($fullsymbol2, $fullsymbol1);
	# }
	
	# Skip homomers
	if ($symbol1 eq $symbol2)
	{
		die("Error: Got a homomer for $symbol1|$symbol2 (only post-aliasing)");
	}
	
	# Skip mutants/fragments/etc. (containing an underscore, _)
	if ($symbol1 =~ /_/)
	{
		addme("skipped mutant/fragment for symbol", $symbol1);
		next;
	}
	if ($symbol2 =~ /_/)
	{
		addme("skipped mutant/fragment for symbol", $symbol2);
		next;
	}
	

	# Check if this pairid is defined in %pairs from above (otherwise it wasn't tested by NanoBRET, at least not successfully)
	if (!exists($pair{$pairid}))
	{
		# die("Error: No pairid for pair '$symbol1|$symbol2'") ;
		addme("nanobret: pairid not annotated in '$classfile' for pairid (skipped)", $pairid);
		next;
	}
	else
	{
		addme("nanobret: pairid was annotated in '$classfile' for pairid (kept)", $pairid);
	}
	
	# # Tested
	# if ($symbol1 lt $symbol2)
	# {
	# 	# "Correct pair"
	# 	$tested{"$symbol1|$symbol2"} = 1;
	# }
	# else
	# {
	# 	# "Wrong pair"
	# 	$tested{"$symbol2|$symbol1"} = 1;
	# }
	#
	# # Positives
	# if ($verdict eq 'positive')
	# {
	# 	if ($symbol1 lt $symbol2)
	# 	{
	# 		# "Correct pair"
	# 		$positive{"$symbol1|$symbol2"} = 1;
	# 	}
	# 	else
	# 	{
	# 		# "Wrong pair"
	# 		$positive{"$symbol2|$symbol1"} = 1;
	# 	}
	# }
	# elsif ($verdict eq 'negative')
	# {
	# 	if ($symbol1 lt $symbol2)
	# 	{
	# 		# "Correct pair"
	# 		$negative{"$symbol1|$symbol2"} = 1;
	# 	}
	# 	else
	# 	{
	# 		# "Wrong pair"
	# 		$negative{"$symbol2|$symbol1"} = 1;
	# 	}
	# }
	# else
	# {
	# 	die;
	# }

	# Tested
	$tested{$pairid} = 1;

	# Positives
	if ($verdict eq 'positive')
	{
		$positive{$pairid} = 1;
	}
	elsif ($verdict eq 'negative')
	{
		$negative{$pairid} = 1;
	}
	else
	{
		die;
	}

	
	
	
	stepme(100);
}
close(IN);
stopme();





# Cycle through classes and get NanoBRET validation ratios
state("Getting NanoBRET validation ratios:");
# $typequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='$scoretype'");
%screen_hit = ();
# while (($scoretype, $threshold) = Fetch($typequery))
print " >> $scoretype\n";
foreach $class (@classes)
{
	print "   >> $class\n";
	
	# $outfile = "output-$scoretype.txt";
	# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	
	# $mainquery = Query("SELECT DISTINCT symbol1, symbol2, eclip1, eclip2 FROM $table WHERE scoretype='$scoretype' AND hc=1 AND homodimer=0");
	
	@pairs = split(/;/, $pairs_per_class{$class});
	@pairs = unique(@pairs);
	
	print "     >> ".scalar(@pairs)." pairs tested\n";
	
	$total_tested = 0;
	$total_positive = 0;
	$total_negative = 0;
	$TP = 0;
	$FP = 0;
	$TN = 0;
	$FN = 0;
	$eclip_TP = 0;
	$eclip_FP = 0;
	$eclip_TN = 0;
	$eclip_FN = 0;
	@nanobret_positives = ();
	@nanobret_negatives = ();
	# while (($symbol1, $symbol2, $eclip1, $eclip2) = Fetch($mainquery))
	foreach $pair (@pairs)
	{
		($symbol1, $symbol2) = split(/\|/, $pair);
		
		# Get eCLIP data
		$eclip1 = 0;
		$eclip2 = 0;
		if (exists($eclip{$symbol1}))
		{
			$eclip1 = 1;
		}
		if (exists($eclip{$symbol2}))
		{
			$eclip2 = 1;
		}
		
		# Screen hit
		$screen_hit{"$symbol1|$symbol2"} = 1;
		
		# Init
		$tested = 0;
		$negative = 0;
		$positive = 0;
		
		if (exists($pairids{"$symbol1|$symbol2"}))
		{
			foreach $pairid (split(/\|/, $pairids{"$symbol1|$symbol2"}))
			{
				# Tested
				if (exists($tested{$pairid}))
				{
					$tested = 1;
				}
		
				# A single positive result overrides all negatives
				if (exists($positive{$pairid}))
				{
					$positive = 1;
					push(@nanobret_positives, $pair);
				}
			}
		}
		
		if (($tested == 1) and ($positive == 0))
		{
			$negative = 1;
			push(@nanobret_negatives, $pair);
		}
		
		
		if ($positive == 1)
		{
			$TP++;
			$eclip_TP++ if (($eclip1 == 1) and ($eclip2 == 1));
		}
		elsif ($negative == 1)
		{
			$FP++;
			$eclip_FP++ if (($eclip1 == 1) and ($eclip2 == 1));
		}
		
		$total_tested += $tested;
		$total_positive += $positive;
		$total_negative += $negative;
	}
	


	# Get TN/FN by cycling through the NanoBRET-tested pairs that weren't screen hits

	# Get all the tested pairs (via pairids)
	@testedpairs = ();
	foreach $pairid (keys(%tested))
	{
		push(@testedpairs, $pair{$pairid});
	}
	@testedpairs = unique(@testedpairs);
	
	# Cycle through all tested pairs
	foreach $pair (@testedpairs)
	{
		($symbol1, $symbol2) = split(/\|/, $pair);
		
		if (!exists($screen_hit{$pair}))
		{
			# Negative
			$negative = 0;
			if (exists($negative{$pair}))
			{
				$negative = 1;
			}
			
			# A positive result overrides negatives
			$positive = 0;
			if (exists($positive{$pair}))
			{
				$positive = 1;
				$negative = 0;
			}
			
			if ($positive == 1)
			{
				$FN++;
				$eclip_FN++ if (exists($eclip{$symbol1}) and exists($eclip{$symbol2}));
			}
			elsif ($negative == 1)
			{
				$TN++;
				$eclip_TN++ if (exists($eclip{$symbol1}) and exists($eclip{$symbol2}));
			}
		}
	}
	
	
	
	$validation = '0';
	if ($total_tested != 0)
	{
		$validation = $total_positive / $total_tested;
	}
	
	
	print "     >> $total_positive / $total_tested positive\n";
	
	@nanobret_positives = unique(@nanobret_positives);
	@nanobret_negatives = unique(@nanobret_negatives);
	
	foreach $pair (@nanobret_positives)
	{
		print "       >> $pair\n";
	}
	
	print "     >> ".sprintf("%.1f", $validation * 100)."%\n\n";


	# Calculate F1 score and Matthew's Correlation Coefficient
	$f1 = '0';
	if (((2 * $TP) + $FP + $FN) != 0)
	{
		$f1 = (2 * $TP) / ((2 * $TP) + $FP + $FN);
	}
	$mcc = '0';
	if (sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN)) != 0)
	{
		$mcc = (($TP * $TN) - ($FP * $FN)) / sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN));
	}

	# Calculate F1 score and Matthew's Correlation Coefficient for eCLIP pairs only
	$eclip_f1 = '0';
	if (((2 * $eclip_TP) + $eclip_FP + $eclip_FN) != 0)
	{
		$eclip_f1 = (2 * $eclip_TP) / ((2 * $eclip_TP) + $eclip_FP + $eclip_FN);
	}
	$eclip_mcc = '0';
	if (sqrt(($eclip_TP + $eclip_FP) * ($eclip_TP + $eclip_FN) * ($eclip_TN + $eclip_FP) * ($eclip_TN + $eclip_FN)) != 0)
	{
		$eclip_mcc = (($eclip_TP * $eclip_TN) - ($eclip_FP * $eclip_FN)) / sqrt(($eclip_TP + $eclip_FP) * ($eclip_TP + $eclip_FN) * ($eclip_TN + $eclip_FP) * ($eclip_TN + $eclip_FN));
	}

	# Round values
	$validation = round($validation, 3);
	$f1 = round($f1, 3);
	$mcc = round($mcc, 3);
	$eclip_f1 = round($eclip_f1, 3);
	$eclip_mcc = round($eclip_mcc, 3);


	# Get eclip_pairs number (how many screen hits are eCLIP-eCLIP interactions?)
	# $query = Query("SELECT COUNT(DISTINCT symbol1, symbol2) FROM $table WHERE scoretype='$scoretype' AND homodimer=0 AND hc=1 AND eclip1=1 AND eclip2=1");
	# $query = Query("SELECT COUNT(DISTINCT symbol1, symbol2) FROM $table WHERE scoretype='$scoretype' AND homodimer=0 AND hc=1 AND eclip1=1 AND eclip2=1 AND CONCAT_WS('|', symbol1, symbol2) IN ('".join("', '", @pairs)."')");
	# ($eclip_pairs) = FetchOne($query);
	$eclip_pairs = 0;
	foreach $pair (@pairs)
	{
		($symbol1, $symbol2) = split(/\|/, $pair);
		if (exists($eclip{$symbol1}) and exists($eclip{$symbol2}))
		{
			$eclip_pairs++;
		}
	}


	# Get Jaccard / oddsratio / probability_mindist from ~/pipeline/rbpome_analysis/output/ table files
	@list = @pairs;
	# @list = @nanobret_positives;
	# Average
	$mean_jaccard = get_cobind(\@list, 'jaccard', 'mean');
	$mean_oddsratio = get_cobind(\@list, 'oddsratio', 'mean');
	$mean_probability_mindist = get_cobind(\@list, 'probability_mindist', 'mean');
	# Median
	$median_jaccard = get_cobind(\@list, 'jaccard', 'median');
	$median_oddsratio = get_cobind(\@list, 'oddsratio', 'median');
	$median_probability_mindist = get_cobind(\@list, 'probability_mindist', 'median');

	# Round values
	# Average
	$mean_jaccard = round($mean_jaccard, 3);
	$mean_oddsratio = round($mean_oddsratio, 3);
	$mean_probability_mindist = round($mean_probability_mindist, 3);
	# Median
	$median_jaccard = round($median_jaccard, 3);
	$median_oddsratio = round($median_oddsratio, 3);
	$median_probability_mindist = round($median_probability_mindist, 3);
	
	$nanobret_positives = join(", ", @nanobret_positives);
	$nanobret_negatives = join(", ", @nanobret_negatives);
	
	print OUT "$class\t$total_positive\t$total_negative\t$total_tested\t$validation\t$TP\t$FP\t$TN\t$FN\t$f1\t$mcc\t$nanobret_positives\t$nanobret_negatives\t$eclip_TP\t$eclip_FP\t$eclip_TN\t$eclip_FN\t$eclip_f1\t$eclip_mcc\t$eclip_pairs\t$mean_jaccard\t$mean_oddsratio\t$mean_probability_mindist\t$median_jaccard\t$median_oddsratio\t$median_probability_mindist\n";
}
nl();
stoptime();



# Skipped:
showme("fullsymbol not found in table 'rbpome' for fullsymbol (skipped)");
showme("homomer skipped for fullsymbols");
showme("nanobret pair skipped because no symbol for fullsymbol for pair (this also indicates the protein wasn't in the screen) (skipped)");
showme("nanobret pair skipped because no symbol for fullsymbol (this also indicates the protein wasn't in the screen) (skipped)");
showme("skipped mutant/fragment for symbol");
showme("homomer");

showmeallsorted(1);
showme("classes: duplicate fullsymbol pairs");

done();




sub get_cobind
{
	my ($pairs, $mode, $median) = @_;
	
	my (@pairs) = @$pairs;
	# show(@pairs);
	# exit;
	
	$tmpfile = "../rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
	open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
	<TMP>;	# Skip header
	my %cobind = ();
	while (<TMP>)
	{
		chomp;
		
		# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	PCBP1	0.16151090849886
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	SF3B4	0.135072908672295
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FMR1	FXR2	0.420327027893556
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FXR1	FXR2	0.15188679245283
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	HNRNPK	U2AF2	0.234877126654064
		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	NONO	SFPQ	0.203641590800192
	
		@a = split(/\t/);
		
		$symbol1 = $a[10];
		$symbol2 = $a[11];
		
		my $pair = "$symbol1|$symbol2";
		
		# Keep only the desired pairs
		# show(@pairs);
		# exit;
		next if (!contains($pair, @pairs));
		
		# $set = $a[9];
		$cobind = $a[12];
		
		# Only use "screen_hits" (not positive controls, not background_resamples)
		# next if ($set ne 'screen_hits');
		
		# push(@cobind, $cobind);
		# Use a hash in order to unique the values (one per pair)
		$cobind{$pair} = $cobind;
	}
	close(TMP);
	
	my @cobind = ();
	# show(keys(%cobind));
	# exit;
	foreach my $pair (keys(%cobind))
	{
		push(@cobind, $cobind{$pair});
	}
	
	if (scalar(@cobind) == 0)
	{
		$res = 0;
	}
	else
	{
		if ($median eq 'median')
		{
			$res = median(@cobind);
		}
		else
		{
			$res = mean(@cobind);
		}
	}

	return($res);
}





