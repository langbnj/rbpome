#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'rbpome';
# $hc_threshold = 3.25;	# Threshold above which 'hc' is set to '1' in the table (high confidence) (AVG_IS)
$aliasfile = "input-aliases.txt";


our $usage = "$0 [score type name] [input file] [hc threshold] [-noalias]\n\n -noalias: Don't assign the aliases in $aliasfile (retain original names).\n\n 'hc threshold': Threshold above which 'hc' is set to '1' in the table (high confidence) (AVG_IS).\n\nExample: $0 AVG_IS input/RBPome_final_CF_0.txt 3.25\nExample: $0 F1 input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_F1_CF_0.96_14-Feb-2020134952_noHDs.txt 0.96\nExample: $0 MCC input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_MCC_CF_0.67_14-Feb-2020153225_noHDs.txt 0.67";
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
	
	# Make symbol all upper case
	$symbol = uc($symbol);
	
	# Make species all upper case
	$species = uc($species);
	
	# # Replace spaces in "RBM10 316F" >> "rbm10_316f", "eMIC SRRM3 mut 25" >> "emic_srrm3_mut_25" etc.
	# $symbol =~ s/ /_/g;
	# Remove spaces in "RBM10 316F" >> "rbm10316f", "eMIC SRRM3 mut 25" >> "emicsrrm3mut25" etc. (as in Sebastian's table)
	$symbol =~ s/ //g;
	
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


startme("Reading edges from '$infile'", 0, chompme(`cat $infile | wc -l`));
starttime();
<IN>;	# Skip header
while (<IN>)
{
	chomp;
	
	# # Protein_1	Protein_2	Times_detected	AVG_IS	AVG_S	AVG_NS	BG_alldirect	BG_all	Homodimer	Found_ADBD_BDAD
	# # dusp14	pspc1	6	1341.496532	6.21979E-05	7.02258E-07	0	0	0	0
	# # pspc1	rbm14	12	1218.777237	0.000270544	1.6835E-06	0	0	0	1
	# # tia1	rbm38	3	1107.23183	1.5262E-06	5.9729E-07	0	0	0	0
	#
	# Protein A (unique pairs!!)	Protein B (unique pairs!!)	Times detected (RIS >0) (Is still like to be able to filter for this. After all we have multi vector transformants that can give stochastic false positives)	uniqueABavgIS (max)	uniqueABsumIS (max)	Avg_RS (raw reads / # screens)	Avg_RIS (raw reads / # screens) => to know if interaction was sampled at all	Found in both orientations (1/0)	H47_avgIS	H47_sumIS	H47_found in both orientations	H47_times detected	Avg_RS (raw reads / # screens) in H47	Avg_RIS (raw reads / # screens) in H47	"Biogrid_all direct evidence (y2h, reconstituted complex, structure)"	Biogrid_all	Hippie	HuRI ?	Protein atlas (nuclear 1/0)	Youn-Gingras_data	Nanobret (pos/neg) unique	Nanobret MBU (mean)	Nanobret MBU std dev (from above)	Jaccard score	Cobinding prob	Cobinding odds ratio	Mean log min distance (bimodal)	Mean log min distance std dev	Meta RNA binding similarity A vs AB KS statistics	Meta RNA binding similarity B vs AB KS statistics
	# rbm12	rbm10776q	4	1.9	20.35	0.18	56.18	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10598f	3	1.74	18.82	0.09	19.18	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10580f	1	1.77	18.75	0	12.09	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10605f	2	0	17.82	0	7.27	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# hnrnpf	ddx3x	1	1.73	17.7	0	6	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10408v	2	0	17.58	0	6.27	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	
	@a = split(/\t/);
	
	# Get rid of N/A values
	foreach (@a)
	{
		$_ = '' if ($_ eq 'N/A');
	}
	
	$symbol1 = $a[0];
	$symbol2 = $a[1];
	$times_detected = $a[2];
	$avg_is = $a[4];
	# $avg_is = $a[3];
	$bg_alldirect = $a[14];
	$bg_all = $a[15];
	# $homodimer = $a[8];
	$found_inverse = $a[7];
	
	$homodimer = 0;
	$homodimer = 1 if ($symbol1 eq $symbol2);
	
	
	# Additional fields to retain
	$unique_ab_avg_is = $a[3];
	$avg_rs = $a[5];
	$avg_ris = $a[6];

	$h47_avg_is = $a[8];
	$h47_sum_is = $a[9];
	$h47_found_in_both_orientations = $a[10];
	$h47_times_detected = $a[11];	
	$h47_avg_rs = $a[12];
	$h47_avg_ris = $a[13];
	
	$hippie = $a[16];
	$huri = $a[17];
	
	# $hpa_nuclear = $a[18];	# blank
	# $bioid = $a[19];			# blank
	
	$nanobret = $a[20];
	$nanobret_mbu_avg = $a[21];
	$nanobret_mbu_sd = $a[22];
	
	
	
	
	# Make symbols upper case
	$symbol1 = uc($symbol1);
	$symbol2 = uc($symbol2);
	
	
	if (!exists($species{$symbol1}))
	{
		# die("Error: No species assignment for symbol '$symbol1'");
		addme("no species assignment for symbol pair (skipped)", "$symbol1|$symbol2");
		addme("no species assignment for symbol (skipped)", $symbol1);
		next;
	}
	if (!exists($species{$symbol2}))
	{
		# die("Error: No species assignment for symbol '$symbol2'");
		addme("no species assignment for symbol pair (skipped)", "$symbol1|$symbol2");
		addme("no species assignment for symbol (skipped)", $symbol2);
		next;
	}
	
	$species1 = $species{$symbol1};
	$species2 = $species{$symbol2};
	
	addme("total species", $species1);
	addme("total species", $species2);
	
	$same_species = 0;
	$same_species = 1 if ($species1 eq $species2);
	
	# Retrieve original symbol (with capitalisation etc.) from the species mapping table
	$fullsymbol1 = $realsymbol{$symbol1};
	$fullsymbol2 = $realsymbol{$symbol2};
	
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
	
	addme("total symbol pairs", "$symbol1|$symbol2");
	addme("total symbols", $symbol1);
	addme("total symbols", $symbol2);
	addme("total 'full symbol' pairs", "$fullsymbol1|$fullsymbol2");
	addme("total 'full symbols'", $fullsymbol1);
	addme("total 'full symbols'", $fullsymbol2);
	
	# Check if the final symbols are real HGNC symbols (human)
	# if ($species1 eq 'HUMAN')
	# {
		$query = Query("SELECT id FROM hgnc WHERE symbol='$symbol1'");
		if (Numrows($query) > 0)
		{
			addme("$species1 symbol is a real hgnc symbol for symbol", $symbol1);
		}
		else
		{
			addme("$species1 symbol is NOT a real hgnc symbol for symbol", $symbol1);

			# Check if it's an alias
			$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$symbol1'");
			if (Numrows($query) > 0)
			{
				addme("$species1 symbol is NOT a real hgnc symbol for symbol, but it's an alias", $symbol1);
			}
			else
			{
				addme("$species1 symbol is NOT a real hgnc symbol for symbol, nor is it an alias", $symbol1);
			}
		}
	# }
	# if ($species2 eq 'HUMAN')
	# {
		$query = Query("SELECT id FROM hgnc WHERE symbol='$symbol2'");
		if (Numrows($query) > 0)
		{
			addme("$species2 symbol is a real hgnc symbol for symbol", $symbol2);
		}
		else
		{
			addme("$species2 symbol is NOT a real hgnc symbol for symbol", $symbol2);

			# Check if it's an alias
			$query = Query("SELECT id FROM hgnc_aliases WHERE alias='$symbol2'");
			if (Numrows($query) > 0)
			{
				addme("$species2 symbol is NOT a real hgnc symbol for symbol, but it's an alias", $symbol2);
			}
			else
			{
				addme("$species2 symbol is NOT a real hgnc symbol for symbol, nor is it an alias", $symbol2);
			}
		}
	# }


	
	# Check if the final symbols are real MGI symbols (mouse)
	# I was going to check mouse using VGNC, but it turns out VGNC is missing mouse. Only MGI has mouse gene names (and synonyms). Got them now in table 'mgi'!
	# Best thing to use though: UniProt. Going to use that for all species (i.e. human & mouse).
	if ($species1 eq 'MOUSE')
	{
		$query = Query("SELECT id FROM mgi WHERE marker_symbol='$symbol1'");
		if (Numrows($query) > 0)
		{
			addme("mouse symbol is a real mgi symbol for symbol", $symbol1);
		}
		else
		{
			addme("mouse symbol is NOT a real mgi symbol for symbol", $symbol1);

			# Check if it's an alias
			$query = Query("SELECT id FROM mgi_aliases WHERE alias='$symbol1'");
			if (Numrows($query) > 0)
			{
				addme("mouse symbol is NOT a real mgi symbol for symbol, but it's an alias", $symbol1);
			}
			else
			{
				addme("mouse symbol is NOT a real mgi symbol for symbol, nor is it an alias", $symbol1);
			}
		}
	}
	if ($species2 eq 'MOUSE')
	{
		$query = Query("SELECT id FROM mgi WHERE marker_symbol='$symbol2'");
		if (Numrows($query) > 0)
		{
			addme("mouse symbol is a real mgi symbol for symbol", $symbol2);
		}
		else
		{
			addme("mouse symbol is NOT a real mgi symbol for symbol", $symbol2);

			# Check if it's an alias
			$query = Query("SELECT id FROM mgi_aliases WHERE alias='$symbol2'");
			if (Numrows($query) > 0)
			{
				addme("mouse symbol is NOT a real mgi symbol for symbol, but it's an alias", $symbol2);
			}
			else
			{
				addme("mouse symbol is NOT a real mgi symbol for symbol, nor is it an alias", $symbol2);
			}
		}
	}
	
	
	
	# Check if the final symbols are real UniProt symbols (human & mouse)
	$query = Query("SELECT id FROM uniprot WHERE species='$species1' AND gene='$symbol1'");
	if (Numrows($query) > 0)
	{
		addme("$species1 symbol is a real UniProt symbol for symbol", $symbol1);
	}
	else
	{
		addme("$species1 symbol is NOT a real UniProt symbol for symbol", $symbol1);
	}
	$query = Query("SELECT id FROM uniprot WHERE species='$species2' AND gene='$symbol2'");
	if (Numrows($query) > 0)
	{
		addme("$species2 symbol is a real UniProt symbol for symbol", $symbol2);
	}
	else
	{
		addme("$species2 symbol is NOT a real UniProt symbol for symbol", $symbol2);
	}
	
	


	
	# Check if the final symbols are real HUMAN UniProt symbols)
	$query = Query("SELECT id FROM uniprot WHERE species='human' AND gene='$symbol1'");
	if (Numrows($query) > 0)
	{
		addme("$species1 symbol is a real HUMAN UniProt symbol for symbol", $symbol1);
	}
	else
	{
		addme("$species1 symbol is NOT a real HUMAN UniProt symbol for symbol", $symbol1);
	}
	$query = Query("SELECT id FROM uniprot WHERE species='human' AND gene='$symbol2'");
	if (Numrows($query) > 0)
	{
		addme("$species2 symbol is a real HUMAN UniProt symbol for symbol", $symbol2);
	}
	else
	{
		addme("$species2 symbol is NOT a real HUMAN UniProt symbol for symbol", $symbol2);
	}
	
	


	
	# Check if the final symbols are real homologene symbols (human & mouse)
	$query = Query("SELECT id FROM homologene WHERE species='$species1' AND symbol='$symbol1'");
	if (Numrows($query) > 0)
	{
		addme("$species1 symbol is a real homologene symbol for symbol", $symbol1);
	}
	else
	{
		addme("$species1 symbol is NOT a real homologene symbol for symbol", $symbol1);
	}
	$query = Query("SELECT id FROM homologene WHERE species='$species2' AND symbol='$symbol2'");
	if (Numrows($query) > 0)
	{
		addme("$species2 symbol is a real homologene symbol for symbol", $symbol2);
	}
	else
	{
		addme("$species2 symbol is NOT a real homologene symbol for symbol", $symbol2);
	}
	
	
	
	# # Check if the 'full symbols' are real HGNC symbols (less relevant)
	# if ($species1 eq 'HUMAN')
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
	# if ($species2 eq 'HUMAN')
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
	# Query("INSERT INTO `$table` SET scoretype='$scoretype', source='".basename($infile)."', threshold='$hc_threshold', symbol1='$symbol1', symbol2='$symbol2', fullsymbol1='$fullsymbol1', fullsymbol2='$fullsymbol2', species1='$species1', species2='$species2', same_species='$same_species', times_detected='$times_detected', avg_is='$avg_is', bg_alldirect='$bg_alldirect', bg_all='$bg_all', homodimer='$homodimer', found_inverse='$found_inverse', eclip1='$eclip1', eclip2='$eclip2', hc='$hc'");
	$q = "INSERT INTO `$table` SET scoretype='$scoretype', source='".basename($infile)."', threshold='$hc_threshold', symbol1='$symbol1', symbol2='$symbol2', fullsymbol1='$fullsymbol1', fullsymbol2='$fullsymbol2', species1='$species1', species2='$species2', same_species='$same_species', times_detected='$times_detected', avg_is='$avg_is', bg_alldirect='$bg_alldirect', bg_all='$bg_all', homodimer='$homodimer', found_inverse='$found_inverse', eclip1='$eclip1', eclip2='$eclip2', hc='$hc', unique_ab_avg_is='$unique_ab_avg_is', avg_rs='$avg_rs', avg_ris='$avg_ris', h47_avg_is='$h47_avg_is', h47_sum_is='$h47_sum_is', h47_found_in_both_orientations='$h47_found_in_both_orientations', h47_times_detected='$h47_times_detected', h47_avg_rs='$h47_avg_rs', h47_avg_ris='$h47_avg_ris', hippie='$hippie', huri='$huri', nanobret='$nanobret', nanobret_mbu_avg='$nanobret_mbu_avg', nanobret_mbu_sd='$nanobret_mbu_sd'";
	$q =~ s/=''/=NULL/g;
	Query($q);
		
	# ALTER TABLE rbpome
	# ADD COLUMN `unique_ab_avg_is` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `avg_rs` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `avg_ris` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `h47_avg_is` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `h47_sum_is` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `h47_found_in_both_orientations` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `h47_times_detected` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `avg_rs` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `avg_ris` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `hippie` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `huri` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `nanobret` VARCHAR(50) NULL DEFAULT NULL,
	# ADD COLUMN `nanobret_mbu_avg` DOUBLE NULL DEFAULT NULL,
	# ADD COLUMN `nanobret_mbu_sd` DOUBLE NULL DEFAULT NULL;
		
	stepme(100);
}
stopme();
stoptime();

showmesome(100);
# showmeall(1);

Optimize($table);

done();
