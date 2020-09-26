#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'rbpome';
# $hc_threshold = 3.25;	# Threshold above which 'hc' is set to '1' in the table (high confidence) (AVG_IS)

# Gene alias mapping (symbol1 and symbol2 will be replaced by these in the 'rbpome' table)
# Used this type of query to get these: https://www.uniprot.org/uniprot/?query=DBC1+AND+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
%alias = ();
$alias{uc("ATP5B")} = uc("ATP5F1B");
$alias{uc("ATP5C1")} = uc("ATP5F1C");
$alias{uc("C16orf80")} = uc("CFAP20");
$alias{uc("DBC1")} = uc("CCAR2");
$alias{uc("DGCR14")} = uc("ESS2");
$alias{uc("DIEXF")} = uc("UTP25");	# It's still DIEXF on UniProt, good to update already though. It'll throw a "not an official HGNC symbol" warning below.
$alias{uc("EIF2C1")} = uc("AGO1");
$alias{uc("EIF2C2")} = uc("AGO2");
# Skipping the eMIC mutants because we don't want them to pollute the screen, they're a specific experiment with Manu/Juan or something
# What's EXOSC3_5437?
$alias{uc("FAM190B")} = uc("CCSER2");
$alias{uc("HNRPDL")} = uc("HNRNPDL");
$alias{uc("HNRPLL")} = uc("HNRNPLL");
$alias{uc("Qk")} = uc("QKI");
$alias{uc("SF3B14")} = uc("SF3B6");
$alias{uc("SGOL2")} = uc("SGO2");
$alias{uc("SKIV2L2")} = uc("MTREX");
$alias{uc("UNR")} = uc("CSDE1");
$alias{uc("ZCCHC11")} = uc("TUT4");
$alias{uc("ZRSR1")} = uc("ZRSR2P1");

our $usage = "$0 [score type name] [input file] [hc threshold]\n\n 'hc threshold': Threshold above which 'hc' is set to '1' in the table (high confidence) (AVG_IS).\n\nExample: $0 AVG_IS input/RBPome_final_CF_0.txt 3.25\nExample: $0 F1 input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_F1_CF_0.96_14-Feb-2020134952_noHDs.txt 0.96\nExample: $0 MCC input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_MCC_CF_0.67_14-Feb-2020153225_noHDs.txt 0.67";
($scoretype, $infile, $hc_threshold) = args(3);

# $infile = "input/RBPome_final_noAAs.txt";
# $infile = "input/RBPome_final_noAAs_MCC_3_25.txt";
# $infile = "input/RBPome_final_CF_1.txt";
# $infile = "input/RBPome_final_CF_0.txt";
$speciesfile = "input/symbol_to_species.txt";
# $outfile = "output-table.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(SPECIES, $speciesfile) or die("\nError: Couldn't open '$speciesfile'\n\n");

# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# Clear($table);



# start

# Get eCLIP RBP symbols
startme("Getting eCLIP RBPs");
$query = Query("SELECT DISTINCT symbol FROM clip_gene WHERE type='eclip_encode'");
%eclip = ();
while (($symbol) = Fetch($query))
{
	$eclip{$symbol} = 1;
	
	stepme(100);
}
stopme();


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
	
	# Make symbol all lower case (like in the table)
	$symbol = lc($symbol);
	
	# Make species all lower case (like in the screen)
	$species = lc($species);
	
	# # Replace spaces in "RBM10 316F" >> "rbm10_316f", "eMIC SRRM3 mut 25" >> "emic_srrm3_mut_25" etc.
	# $symbol =~ s/ /_/g;
	# Remove spaces in "RBM10 316F" >> "rbm10316f", "eMIC SRRM3 mut 25" >> "emicsrrm3mut25" etc. (as in Sebastian's table)
	$symbol =~ s/ //g;
	
	# # Somehow "2700060E02Rik" turns into >> "x2700060e02rik" in the screen, add the x
	# # This is the only RIKEN ID.
	# $symbol = 'x'.$symbol if ($symbol eq '2700060e02rik');
	# if ($symbol eq '2700060e02rik')	{ $species = 'm.musculus'; $symbol = "x2700060e02rik"; }
	# die("Error: RIKEN ID '$symbol' needs to be handled") if ($symbol =~ /^\d+.+rik$/i);
	
	# if ($species eq 'h.sapiens')
	# {
	# 	$species = 'human';
	# }
	# elsif ($species eq 'm.musculus')
	# {
	# 	$species = 'mouse';
	# }
	# elsif ($species eq 'zebrafish')
	# {
	# 	$species = 'danre';
	# }
	# elsif ($species eq 'jellyfish')
	# {
	# 	$species = 'aeqvi';
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


startme("Reading edges from '$infile'", 0, chompme(`cat $infile | wc -l`));
starttime();
<IN>;	# Skip header
while (<IN>)
{
	chomp;
	
	# Protein_1	Protein_2	Times_detected	AVG_IS	AVG_S	AVG_NS	BG_alldirect	BG_all	Homodimer	Found_ADBD_BDAD
	# dusp14	pspc1	6	1341.496532	6.21979E-05	7.02258E-07	0	0	0	0
	# pspc1	rbm14	12	1218.777237	0.000270544	1.6835E-06	0	0	0	1
	# tia1	rbm38	3	1107.23183	1.5262E-06	5.9729E-07	0	0	0	0
	
	@a = split(/\t/);
	
	$symbol1 = $a[0];
	$symbol2 = $a[1];
	$times_detected = $a[2];
	$avg_is = $a[3];
	$avg_s = $a[4];
	$avg_ns = $a[5];
	$bg_alldirect = $a[6];
	$bg_all = $a[7];
	$homodimer = $a[8];
	$found_inverse = $a[9];
	
	if (!exists($species{$symbol1}))
	{
		# die("Error: No species assignment for symbol '$symbol1'");
		addme("no species assignment for symbol (skipped)", $symbol1);
		next;
	}
	if (!exists($species{$symbol2}))
	{
		# die("Error: No species assignment for symbol '$symbol2'");
		addme("no species assignment for symbol (skipped)", $symbol2);
		next;
	}
	
	$species1 = $species{$symbol1};
	$species2 = $species{$symbol2};
	
	$same_species = 0;
	$same_species = 1 if ($species1 eq $species2);
	
	# Retrieve original symbol (with capitalisation etc.) from the species mapping table
	$fullsymbol1 = $realsymbol{$symbol1};
	$fullsymbol2 = $realsymbol{$symbol2};
	
	# Make symbols all upper case
	$symbol1 = uc($symbol1);
	$symbol2 = uc($symbol2);
	
	# Apply alias mapping (e.g. ATP5B protein to ATP5F1B, which is its gene)
	if (exists($alias{$symbol1}))
	{
		$symbol1 = $alias{$symbol1};
	}
	if (exists($alias{$symbol2}))
	{
		$symbol2 = $alias{$symbol2};
	}
	
	addme("total symbol pairs", "$symbol1|$symbol2");
	addme("total symbols", $symbol1);
	addme("total symbols", $symbol2);
	addme("total 'full symbol' pairs", "$fullsymbol1|$fullsymbol2");
	addme("total 'full symbols'", $fullsymbol1);
	addme("total 'full symbols'", $fullsymbol2);
	
	# Check if the final symbols are real HGNC symbols
	if ($species1 eq 'human')
	{
		$query = Query("SELECT id FROM hgnc WHERE symbol='$symbol1'");
		if (Numrows($query) > 0)
		{
			addme("human symbol is a real hgnc symbol for symbol", $symbol1);
		}
		else
		{
			addme("human symbol is NOT a real hgnc symbol for symbol", $symbol1);

			# Check if it's an alias
			$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$symbol1'");
			if (Numrows($query) > 0)
			{
				addme("human symbol is NOT a real hgnc symbol for symbol, but it's an alias", $symbol1);
			}
			else
			{
				addme("human symbol is NOT a real hgnc symbol for symbol, nor is it an alias", $symbol1);
			}
		}
	}
	if ($species2 eq 'human')
	{
		$query = Query("SELECT id FROM hgnc WHERE symbol='$symbol2'");
		if (Numrows($query) > 0)
		{
			addme("human symbol is a real hgnc symbol for symbol", $symbol2);
		}
		else
		{
			addme("human symbol is NOT a real hgnc symbol for symbol", $symbol2);

			# Check if it's an alias
			$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$symbol2'");
			if (Numrows($query) > 0)
			{
				addme("human symbol is NOT a real hgnc symbol for symbol, but it's an alias", $symbol2);
			}
			else
			{
				addme("human symbol is NOT a real hgnc symbol for symbol, nor is it an alias", $symbol2);
			}
		}
	}
	
	# # Check if the 'full symbols' are real HGNC symbols (less relevant)
	# if ($species1 eq 'human')
	# {
	# 	$query = Query("SELECT id FROM hgnc WHERE symbol='$fullsymbol1'");
	# 	if (Numrows($query) > 0)
	# 	{
	# 		addme("human 'full symbol' is a real hgnc symbol for symbol", $fullsymbol1);
	# 	}
	# 	else
	# 	{
	# 		addme("human 'full symbol' is NOT a real hgnc symbol for symbol", $fullsymbol1);
	#
	# 		# Check if it's an alias
	# 		$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$fullsymbol1'");
	# 		if (Numrows($query) > 0)
	# 		{
	# 			addme("human 'full symbol' is NOT a real hgnc symbol for symbol, but it's an alias", $fullsymbol1);
	# 		}
	# 		else
	# 		{
	# 			addme("human 'full symbol' is NOT a real hgnc symbol for symbol, nor is it an alias", $fullsymbol1);
	# 		}
	# 	}
	# }
	# if ($species2 eq 'human')
	# {
	# 	$query = Query("SELECT id FROM hgnc WHERE symbol='$fullsymbol2'");
	# 	if (Numrows($query) > 0)
	# 	{
	# 		addme("human 'full symbol' is a real hgnc symbol for symbol", $fullsymbol2);
	# 	}
	# 	else
	# 	{
	# 		addme("human 'full symbol' is NOT a real hgnc symbol for symbol", $fullsymbol2);
	#
	# 		# Check if it's an alias
	# 		$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$fullsymbol2'");
	# 		if (Numrows($query) > 0)
	# 		{
	# 			addme("human 'full symbol' is NOT a real hgnc symbol for symbol, but it's an alias", $fullsymbol2);
	# 		}
	# 		else
	# 		{
	# 			addme("human 'full symbol' is NOT a real hgnc symbol for symbol, nor is it an alias", $fullsymbol2);
	# 		}
	# 	}
	# }
	
	# Get eCLIP statuses
	$eclip1 = 0;
	$eclip2 = 0;
	if (exists($eclip{$symbol1})) { $eclip1 = 1; }
	if (exists($eclip{$symbol2})) { $eclip2 = 1; }
	
	# Get hc status
	$hc = 0;
	# if ($avg_is > $hc_threshold)	# << wrong! Should be ≥:
	if ($avg_is >= $hc_threshold)
	{
		$hc = 1;
	}
	
	# eclip pairs above threshold (i.e. hc)
	if (($eclip1 == 1) and ($eclip2 == 1) and ($hc == 1))
	{
		if ($homodimer == 1)
		{
			addme("total eclip pairs (homodimers)", "$symbol1|$symbol2");
		}
		else
		{
			addme("total eclip pairs (heteromers)", "$symbol1|$symbol2");
		}
	}

	# Insert into MySQL	
	Query("INSERT INTO `$table` SET scoretype='$scoretype', source='".basename($infile)."', threshold='$hc_threshold', symbol1='$symbol1', symbol2='$symbol2', fullsymbol1='$fullsymbol1', fullsymbol2='$fullsymbol2', species1='$species1', species2='$species2', same_species='$same_species', times_detected='$times_detected', avg_is='$avg_is', avg_s='$avg_s', avg_ns='$avg_ns', bg_alldirect='$bg_alldirect', bg_all='$bg_all', homodimer='$homodimer', found_inverse='$found_inverse', eclip1='$eclip1', eclip2='$eclip2', hc='$hc'");
	
	stepme(100);
}
stopme();
stoptime();

# showmesome(50);
showmeall(1);

Optimize($table);

done();
