#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $table = 'rbpome';

# eCLIP parameters (for getting Jaccard etc. from rbpome_analysis output files)
$scoretype = 'SUM_IS';
$threshold = 7.1;
$type = 'eclip_encode_12';
$hc = 0;
$minlog2fold = 0;
$mindist_threshold = 54;
$resamples = 10000;

# p-value alpha for Cytoscape pairs
$alpha = 0.05;



our $usage = "$0 [table (e.g. 'rbpome')] [-update]\n\n -update: Update the 'rbpome' table for some fields (HuRI, BioGRID, HIPPIE)\n\nExample: $0 rbpome";
($table) = args(1);

# $infile = input/_Final_Summary_Table.20200229.v3.txt";
$infile = "output-$table.txt";
$outfile = "output-$table-master.txt";
$mutfile = "output-$table-mutants-rbm10-srrms.txt";
$fragfile = "output-$table-fragments-etc.txt";
$edgefile = "output-$table-network-edges.txt";
$strictedgefile = "output-$table-network-edges-strict.txt";
$nodefile = "output-$table-network-nodes.txt";
$strictnodefile = "output-$table-network-nodes-strict.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
open(MUT, ">$mutfile") or die("\nError: Couldn't open '$mutfile'\n\n");
open(FRAG, ">$fragfile") or die("\nError: Couldn't open '$fragfile'\n\n");
open(EDGE, ">$edgefile") or die("\nError: Couldn't open '$edgefile'\n\n");
open(STRICTEDGE, ">$strictedgefile") or die("\nError: Couldn't open '$strictedgefile'\n\n");
open(NODE, ">$nodefile") or die("\nError: Couldn't open '$nodefile'\n\n");
open(STRICTNODE, ">$strictnodefile") or die("\nError: Couldn't open '$strictnodefile'\n\n");


print NODE "symbol\tnode_known_protein\tnode_hpa_nuclear\tnode_hpa_cytoplasmic\tnode_hpa_shuttling_loc\n";
print STRICTNODE "symbol\tnode_known_protein\tnode_hpa_nuclear\tnode_hpa_cytoplasmic\tnode_hpa_shuttling_loc\n";



# start



# These are impromptu remaps here, in order to get UniProt accessions. Most will be gotten via their gene symbol (because I made them match UniProt's):
# Query to help find these: SELECT a.primary_acc, p.* FROM uniprot p, uniacc a WHERE p.species='MOUSE' AND p.gene='KIF7' AND p.name=a.name GROUP BY a.primary_acc;
%acc = ();
# $acc{'DANRE|DDX39AB_ZEBRAFISH'} 		= 'Q803W0';
$acc{'DANRE|DDX39A'}			 		= 'Q803W0';
$acc{'DANRE|ESRP2_ZEBRAFISH'} 			= 'G1K2P5';
$acc{'HUMAN|AGO2_CTER'}					= 'Q9UKV8';
$acc{'HUMAN|DDX39A_NTER'}				= 'O00148';
$acc{'HUMAN|HNRNPDL_CTER'}				= 'O14979';
$acc{'HUMAN|RBM10_316F'} 				= 'P98175';
$acc{'HUMAN|RBM10_343G'} 				= 'P98175';
$acc{'HUMAN|RBM10_354E'} 				= 'P98175';
$acc{'HUMAN|RBM10_408V'} 				= 'P98175';
$acc{'HUMAN|RBM10_580F'} 				= 'P98175';
$acc{'HUMAN|RBM10_598F'} 				= 'P98175';
$acc{'HUMAN|RBM10_605F'} 				= 'P98175';
$acc{'HUMAN|RBM10_776Q'} 				= 'P98175';
$acc{'HUMAN|RBM10_784L'} 				= 'P98175';
$acc{'HUMAN|SRRM3_EMIC'} 				= 'A6NNA2';
$acc{'HUMAN|SRRM3_EMIC_MUT2'} 			= 'A6NNA2';
$acc{'HUMAN|SRRM3_EMIC_MUT5'} 			= 'A6NNA2';
$acc{'HUMAN|SRRM3_EMIC_MUT25'} 			= 'A6NNA2';
$acc{'HUMAN|SRRM4_EMIC'} 				= 'A7MD48';
$acc{'HUMAN|SRRM4_EMIC_MUT2'} 			= 'A7MD48';
$acc{'HUMAN|SRRM4_EMIC_MUT5'} 			= 'A7MD48';
$acc{'HUMAN|SRRM4_EMIC_MUT25'} 			= 'A7MD48';
$acc{'HUMAN|SRSF11_CONSERVED_PEPTIDE'} 	= 'Q05519';
$acc{'MOUSE|APC_CTER'} 					= 'Q61315';
$acc{'MOUSE|APC_NTER'} 					= 'Q61315';
$acc{'MOUSE|FAU_UBIM_RS30'} 			= 'P35545;P62862';
$acc{'MOUSE|FUBP3'} 					= 'Q3TIX6';
$acc{'MOUSE|KIF1B_CTER'} 				= 'Q60575';
$acc{'MOUSE|KIF1B_NTER'} 				= 'Q60575';
$acc{'MOUSE|KIF1C_CTER'} 				= 'O35071';
$acc{'MOUSE|KIF7_DEL'} 					= 'B7ZNG0';
$acc{'MOUSE|KIF11_CTER'} 				= 'Q6P9P6';
$acc{'MOUSE|KIF16B_CTER'} 				= 'B1AVY7';
$acc{'MOUSE|KIF17V2'} 					= 'A2AM72';
$acc{'MOUSE|KIF19A_CTER'} 				= 'Q99PT9';
$acc{'MOUSE|KIF20A_CTER'} 				= 'P97329';
$acc{'MOUSE|KIF23_CTER'} 				= 'E9Q5G3';
$acc{'MOUSE|KIFC2_NTER'} 				= 'O08672';
$acc{'MOUSE|KIFC3_NTER'} 				= 'O35231';
$acc{'MOUSE|MYO5A_CTER'} 				= 'Q99104';
$acc{'MOUSE|MYO5A_NTER'} 				= 'Q99104';
$acc{'MOUSE|SF3B2'} 					= 'Q3UJB0';
$acc{'MOUSE|SYT4_38_425'} 				= 'P40749';
$acc{'MOUSE|SYT4_38_STOP'} 				= 'P40749';
$acc{'MOUSE|TAF15'} 					= 'Q8BQ46';
$acc{'MOUSE|TUBB'} 						= 'P99024';
$acc{'MOUSE|ZC3H7B'} 					= 'F8VPP8';



# These are impromptu remaps for the fragments and mutants here, in order to get gene symbols to map by.
%symbol = ();
$symbol{'DDX39AB_ZEBRAFISH'} 			= 'DDX39A';
$symbol{'ESRP2_ZEBRAFISH'} 				= 'ESRP2';
$symbol{'AGO2_CTER'}					= 'AGO2';
$symbol{'DDX39A_NTER'}					= 'DDX39A';
$symbol{'HNRNPDL_CTER'}					= 'HNRNPDL';
$symbol{'RBM10_316F'} 					= 'RBM10';
$symbol{'RBM10_343G'} 					= 'RBM10';
$symbol{'RBM10_354E'} 					= 'RBM10';
$symbol{'RBM10_408V'} 					= 'RBM10';
$symbol{'RBM10_580F'} 					= 'RBM10';
$symbol{'RBM10_598F'} 					= 'RBM10';
$symbol{'RBM10_605F'} 					= 'RBM10';
$symbol{'RBM10_776Q'} 					= 'RBM10';
$symbol{'RBM10_784L'} 					= 'RBM10';
$symbol{'SRRM3_EMIC'} 					= 'SRRM3';
$symbol{'SRRM3_EMIC_MUT2'} 				= 'SRRM3';
$symbol{'SRRM3_EMIC_MUT5'} 				= 'SRRM3';
$symbol{'SRRM3_EMIC_MUT25'} 			= 'SRRM3';
$symbol{'SRRM4_EMIC'} 					= 'SRRM4';
$symbol{'SRRM4_EMIC_MUT2'} 				= 'SRRM4';
$symbol{'SRRM4_EMIC_MUT5'} 				= 'SRRM4';
$symbol{'SRRM4_EMIC_MUT25'} 			= 'SRRM4';
$symbol{'SRSF11_CONSERVED_PEPTIDE'} 	= 'SRSF11';
$symbol{'APC_CTER'} 					= 'APC';
$symbol{'APC_NTER'} 					= 'APC';
$symbol{'FAU_UBIM_RS30'} 				= 'FAU';
# $symbol{'FUBP3'} 						= 'FUBP3';
$symbol{'KIF1B_CTER'} 					= 'KIF1B';
$symbol{'KIF1B_NTER'} 					= 'KIF1B';
$symbol{'KIF1C_CTER'} 					= 'KIF1C';
$symbol{'KIF7_DEL'} 					= 'KIF7';
$symbol{'KIF11_CTER'} 					= 'KIF11';
$symbol{'KIF16B_CTER'} 					= 'KIF16B';
$symbol{'KIF17V2'} 						= 'KIF17';
$symbol{'KIF19A_CTER'} 					= 'KIF19';
$symbol{'KIF20A_CTER'} 					= 'KIF20A';
$symbol{'KIF23_CTER'} 					= 'KIF23';
$symbol{'KIFC2_NTER'} 					= 'KIFC2';
$symbol{'KIFC3_NTER'} 					= 'KIFC3';
$symbol{'MYO5A_CTER'} 					= 'MYO5A';
$symbol{'MYO5A_NTER'} 					= 'MYO5A';
# $symbol{'SF3B2'} 						= 'SF3B2';
$symbol{'SYT4_38_425'} 					= 'SYT4';
$symbol{'SYT4_38_STOP'} 				= 'SYT4';
# $symbol{'TAF15'} 						= 'TAF15';
# $symbol{'TUBB'} 						= 'TUBB';
# $symbol{'ZC3H7B'} 						= 'ZC3H7B';



# For mouse FAU:
# We've got, in H1001xH1057_current_with_uniprot_accs.xlsx:
# MQLFVRAQELHTLEVTGQETVAQIKDHVASLEGIAPEDQVVLLAGSPLEDEATLGQCGVEALTTLEVAGRMLGGKVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRMQYNRRFVNVVPTFGKKKGPNANS
#
# N-term:
# >sp|P35545|UBIM_MOUSE Ubiquitin-like protein FUBI OS=Mus musculus OX=10090 GN=Fau PE=1 SV=1
# MQLFVRAQELHTLEVTGQETVAQIKDHVASLEGIAPEDQVVLLAGSPLEDEATLGQCGVEALTTLEVAGRMLGG
#
# C-term:
# >sp|P62862|RS30_MOUSE 40S ribosomal protein S30 OS=Mus musculus OX=10090 GN=Fau PE=1 SV=1
# KVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRMQYNRRFVNVVPTFGKKKGPNANS
# 
# So...we have both concatenated, but they probably don't get cleaved in yeast.
# Should map it to both accessions.





# #DEBUG
# warn("DEBUG: Currently getting rbpome_analysis results from ../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/");
# # instead of ../../pipeline/rbpome_analysis/
# #END DEBUG



# Get jaccard values
$mode = 'jaccard';
%jaccard = ();
# $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$tmpfile = "../../pipeline/rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
<TMP>;	# Skip header
$res = '';
while (<TMP>)
{
	chomp;
	
	# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	PCBP1	0.16151090849886
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	SF3B4	0.135072908672295
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	FMR1	FXR2	0.420327027893556
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	FXR1	FXR2	0.15188679245283
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	HNRNPK	U2AF2	0.234877126654064
	# rbpome	eclip_encode_12	AVG_IS	jaccard	0	0	3.25	400	100	positives	NONO	SFPQ	0.203641590800192

	@a = split(/\t/);
	$set = $a[9];
	$rbp1 = $a[10];
	$rbp2 = $a[11];
	$cobind = $a[12];
	
	# Only use "screen_hits" (not positive controls, not background_resamples)
	next if ($set ne 'screen_hits');
	
	$jaccard{"$rbp1|$rbp2"} = $cobind;
	$jaccard{"$rbp2|$rbp1"} = $cobind;
}
close(TMP);

# Get probability values
$mode = 'probability';
%probability = ();
# $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$tmpfile = "../../pipeline/rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
<TMP>;	# Skip header
$res = '';
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
	$set = $a[9];
	$rbp1 = $a[10];
	$rbp2 = $a[11];
	$cobind = $a[12];
	
	# Only use "screen_hits" (not positive controls, not background_resamples)
	next if ($set ne 'screen_hits');
	
	$probability{"$rbp1|$rbp2"} = $cobind;
	# $probability{"$rbp2|$rbp1"} = $cobind;
}
close(TMP);

# Get probability_mindist values
$mode = 'probability_mindist';
%probability_mindist = ();
# $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
$tmpfile = "../../pipeline/rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
<TMP>;	# Skip header
$res = '';
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
	$set = $a[9];
	$rbp1 = $a[10];
	$rbp2 = $a[11];
	$cobind = $a[12];
	
	# Only use "screen_hits" (not positive controls, not background_resamples)
	next if ($set ne 'screen_hits');
	
	$probability_mindist{"$rbp1|$rbp2"} = $cobind;
	# $probability_mindist{"$rbp2|$rbp1"} = $cobind;
}
close(TMP);

# # Get oddsratio values
# $mode = 'oddsratio';
# %oddsratio = ();
# # $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
# $tmpfile = "../../pipeline/rbpome_analysis/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
# open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
# <TMP>;	# Skip header
# $res = '';
# while (<TMP>)
# {
# 	chomp;
#
# 	# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	PCBP1	0.16151090849886
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	SF3B4	0.135072908672295
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FMR1	FXR2	0.420327027893556
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FXR1	FXR2	0.15188679245283
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	HNRNPK	U2AF2	0.234877126654064
# 	# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	NONO	SFPQ	0.203641590800192
#
# 	@a = split(/\t/);
# 	$set = $a[9];
# 	$rbp1 = $a[10];
# 	$rbp2 = $a[11];
# 	$cobind = $a[12];
#
# 	# Only use "screen_hits" (not positive controls, not background_resamples)
# 	next if ($set ne 'screen_hits');
#
# 	$oddsratio{"$rbp1|$rbp2"} = $cobind;
# 	$oddsratio{"$rbp2|$rbp1"} = $cobind;
# }
# close(TMP);




# Get information from output-fit_vs_random-50.txt (resampling of the RBP pair binding distances)
%fit_vs_random = ();
# $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output-fit_vs_random-$table-$type-$scoretype-$mindist_threshold.txt";
$tmpfile = "../../pipeline/rbpome_analysis/output-fit_vs_random-$table-$type-$scoretype-$mindist_threshold.txt";
if (-s $tmpfile)
{
	open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
	<TMP>;	# Skip header
	while (<TMP>)
	{
		chomp;

		# # symbol1	symbol2	close_real	fraction_close_real	fraction_close_random	fraction_ratio	resampling_p	resampling_wilcox
		# # APOBEC3C	HNRNPK	24	0.04040404	0.03871021	1.04375667298111	0.45	1
		# # HNRNPK	APOBEC3C	24	0.01505646	0.01517115	0.99244025667138	0.57	1
		# # APOBEC3C	PCBP1	221	0.2036866	0.05850176	3.48171747311534	0	0
		# # PCBP1	APOBEC3C	277	0.1198615	0.02868432	4.17864185032101	0	0
		
		# symbol1	symbol2	median_dist_real	median_dist_random	close_real	fraction_close_real	fraction_close_random	fraction_close_delta	fraction_close_ratio	resampling_p	resampling_wilcox
		# APOBEC3C	HNRNPK	4358	3057	7	0.02071006	0.02091778	-0.000207720000000002	0.990069691908032	0.55	1
		# APOBEC3C	PCBP1	521	2131	70	0.1886792	0.02695708	0.16172212	6.99924472531891	0	0
		# APOBEC3C	PCBP2	1825	2518	65	0.1547619	0.02531133	0.12945057	6.1143329884285	0	0
		# CPEB4	CSTF2T	2286	2304	20	0.1724138	0.0891533	0.0832605	1.93390261493405	0	1
		
	
		@a = split(/\t/);
		$symbol1 = $a[0];
		$symbol2 = $a[1];

		$fit_vs_random{"$symbol1|$symbol2"} = $_;
	}
	close(TMP);
}
else
{
	warn("WARNING: '$tmpfile' didn't exist yet - Some columns will remain blank (incomplete)");
}





















# Copy table header to outfile unchanged
$_ = <IN>;
@head = split(/\t/, chompme($_));
$spliced_in = 0;

$head[0] = "Protein A";
$head[1] = "Protein B";
$head[8] = "Times detected (RIS >0)";
$head[9] = "avgIS";
$head[10] = "sumIS";
$head[20] = '"Biogrid_all direct evidence (Y2H, reconstituted complex, structure)"';
$head[23] = "HuRI";
$head[24] = "Known interaction";

# Splice in columns in the middle
# Known Complex Member
splice(@head, 25, 0, "Protein A has known interactions"); $spliced_in++;
splice(@head, 26, 0, "Protein B has known interactions"); $spliced_in++;
# HPA
splice(@head, 27, 0, "HPA nuclear A"); $spliced_in++;
splice(@head, 28, 0, "HPA nuclear B"); $spliced_in++;
splice(@head, 29, 0, "HPA cytoplasmic A"); $spliced_in++;
splice(@head, 30, 0, "HPA cytoplasmic B"); $spliced_in++;
splice(@head, 31, 0, "HPA main locations A"); $spliced_in++;
splice(@head, 32, 0, "HPA additional locations A"); $spliced_in++;
splice(@head, 33, 0, "HPA main locations B"); $spliced_in++;
splice(@head, 34, 0, "HPA additional locations B"); $spliced_in++;
splice(@head, 35, 0, "HPA shared locations"); $spliced_in++;
# eCLIP yes/no
splice(@head, 40, 0, "ENCODE eCLIP data A"); $spliced_in++;
splice(@head, 41, 0, "ENCODE eCLIP data B"); $spliced_in++;

$head[42] = "ENCODE eCLIP binding sites A";
$head[43] = "ENCODE eCLIP binding sites B";

splice(@head, 44, 0, "Jaccard index"); $spliced_in++;
splice(@head, 45, 0, "Cobinding prob p(A|B)"); $spliced_in++;
splice(@head, 46, 0, "Cobinding prob p(B|A)"); $spliced_in++;

splice(@head, 47, 0, "Cobinding prob p(A|B) (≤$mindist_threshold nt)"); $spliced_in++;
splice(@head, 48, 0, "Cobinding prob p(B|A) (≤$mindist_threshold nt)"); $spliced_in++;

$head[49] = "Close binding events (≤$mindist_threshold nt) A vs. B";
$head[50] = "Fraction of close binding events A vs. B";
$head[51] = "Fraction of close binding events in random data A vs. B";
$head[52] = "Ratio of fractions A vs. B";
$head[53] = "Resampling p-value A vs. B";
push(@head, "Resampling Wilcoxon p-value A vs. B");	# 54

push(@head, "Close binding events (≤$mindist_threshold nt) B vs. A");	# 55
push(@head, "Fraction of close binding events B vs. A");	# 56
push(@head, "Fraction of close binding events in random data B vs. A");	# 57
push(@head, "Ratio of fractions B vs. A");	# 58
push(@head, "Resampling p-value B vs. A");	# 59
push(@head, "Resampling Wilcoxon p-value B vs. A");	# 60

print OUT join("\t", @head)."\n";
print MUT join("\t", @head)."\n";
print FRAG join("\t", @head)."\n";
print EDGE join("\t", @head)."\tstrict\n";
print STRICTEDGE join("\t", @head)."\tstrict\n";


nl();
state("HEAD:");
show(@head);
nl();
nl();
# exit;


$clean = 0;
$mut = 0;
$frag = 0;
$edge = 0;
$strictedge = 0;

startme("Reading '$infile', filling in additional columns and splitting output into '$outfile', '$mutfile', '$fragfile', '$edgefile', and '$strictedgefile'", 0, chompme(`cat $infile | wc -l`) - 1);
@nodesymbols = ();
@strictnodesymbols = ();
%known_protein = ();
%node_known_protein = ();
%strict_node_known_protein = ();
%node_hpa_nuclear = ();
%node_hpa_cytoplasmic = ();
starttime();
# Table body
while (<IN>)
{
	chomp;
	
	# Protein A (unique pairs!!)	Protein B (unique pairs!!)	Times detected (RIS >0) (Is still like to be able to filter for this. After all we have multi vector transformants that can give stochastic false positives)	uniqueABavgIS (max)	uniqueABsumIS (max)	Avg_RS (raw reads / # screens)	Avg_RIS (raw reads / # screens) => to know if interaction was sampled at all	Found in both orientations (1/0)	H47_avgIS	H47_sumIS	H47_found in both orientations	H47_times detected	Avg_RS (raw reads / # screens) in H47	Avg_RIS (raw reads / # screens) in H47	"Biogrid_all direct evidence (y2h, reconstituted complex, structure)"	Biogrid_all	Hippie	HuRI ?	Protein atlas (nuclear 1/0)	Youn-Gingras_data	Nanobret (pos/neg) unique	Nanobret MBU (mean)	Nanobret MBU std dev (from above)	Jaccard score	Cobinding prob	Cobinding odds ratio	Mean log min distance (bimodal)	Mean log min distance std dev	Meta RNA binding similarity A vs AB KS statistics	Meta RNA binding similarity B vs AB KS statistics
	# rbm12	rbm10776q	4	1.9	20.35	0.18	56.18	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10598f	3	1.74	18.82	0.09	19.18	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10580f	1	1.77	18.75	0	12.09	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# rbm12	rbm10605f	2	0	17.82	0	7.27	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	# hnrnpf	ddx3x	1	1.73	17.7	0	6	0	N/A	N/A	N/A	N/A	N/A	N/A	0	0	0	0			N/A	N/A	N/A
	
	@a = split(/\t/);
	
	# #DEBUG
	# warn("DEBUG: Skipping to first NanoBRET positive");
	# next if ($_ !~ /\tpos\t/);
	# show(@a);
	# # exit;
	# #END DEBUG
	
	# Fill up the array, minus any fields spliced in below
	foreach $i (0 .. scalar(@head) - ($spliced_in - 1))
	{
		# if (undef($a[$i]))
		if (!exists($a[$i]))
		{
			$a[$i] = 'N/A';
		}
	}
	# show(@a);
	# $a[21] = 'asdfasdfasd';
	# show(@a);
	# exit;
	
	
	# Get UniProt accessions and species 1
	die if ($head[0] ne 'Protein A');
	$symbol1 = $a[0];
	$symbol1 = uc($symbol1);
	
	# Apply symbol alias if there is one (so I can look up additional information, like from HPA, on e.g. the RBM10 mutants)
	$originalsymbol1 = $symbol1;
	$symbol1 = $symbol{$symbol1} if (exists($symbol{$symbol1}));

	# $accs1 = '';

	print " >> SYMBOL $symbol1\n" if (switch('debug'));
	print "   >> ORIGINAL $symbol1\n" if (switch('debug'));

	# Get its original species
	@species = ();
	@names = ();
	@accs = ();
	$query = Query("SELECT species1 FROM $table WHERE symbol1='$originalsymbol1' GROUP BY species1 ORDER BY species1");
	while (($species) = Fetch($query))
	{
		push(@species, $species);

		print "     >> SPECIES $species\n" if (switch('debug'));
	}
	$query = Query("SELECT species2 FROM $table WHERE symbol2='$originalsymbol1' GROUP BY species2 ORDER BY species2");
	while (($species) = Fetch($query))
	{
		push(@species, $species);

		print "     >> SPECIES $species\n" if (switch('debug'));
	}
	# die if (scalar(@species) != unique(@species));
	@species = unique(@species);
	# # Don't unique @species, just sort it. There are cases like RBFOX2 that have two HUMAN proteins, for instance. Uniquing breaks these.
	# @species = nsort(@species);

	foreach $species (@species)
	{
		print "       >> SPECIES $species\n" if (switch('debug'));
		
		# Get UniProt name for this species
		$query = Query("SELECT name FROM uniprot WHERE species='$species' AND gene='$symbol1'");
		if (Numrows($query) == 0)
		{
			if (!exists($acc{"$species|$symbol1"}))
			{
				warn("Warning: No uniprot name in table 'uniprot' for species|symbol '$species|$symbol1' (needs to be mapped manually to an accession using the \%acc hash) (kept)");
				addme("no uniprot name in table 'uniprot' for species|symbol (needs to be mapped manually to an accession using the \%acc hash) (kept)", "$species|$symbol1");
			}
			else
			{
				# Use manual %acc hash for mapping (defined above)
				foreach $acc (split(/;/, $acc{"$species|$symbol1"}))
				{
					push(@accs, $acc);

					# Get name, too
					$query = Query("SELECT name FROM uniacc WHERE acc='$acc' AND species='$species'");
					if (Numrows($query) == 0)
					{
						# Not in Swiss-Prot
						$name = $acc.'_'.$species;
					}
					else
					{
						($name) = FetchOne($query);
					}
					push(@names, $name);
				}
			}
		}
		else
		{
			while (($name) = Fetch($query))
			{
				push(@names, $name);
				print "         >> NAME $name\n" if (switch('debug'));

				# Get UniProt accession
				$query = Query("SELECT DISTINCT primary_acc FROM uniacc WHERE name='$name' AND species='$species'");
				($acc) = FetchOne($query);
				push(@accs, $acc);
				print "         >> ACC $acc\n" if (switch('debug'));
			}
		}
	}

	$species1 = join(";", @species);
	$accs1 = join(";", @accs);
	$names1 = join(";", @names);
	print "           >> SPECIES1 $species1\n" if (switch('debug'));
	print "           >> ACCS1 $accs1\n" if (switch('debug'));
	print "           >> NAMES1 $names1\n" if (switch('debug'));
	$a[2] = $species1;	# 3
	$a[4] = $accs1;	# 5
	$a[6] = $names1;	# 7



	# Get UniProt accessions and species 2
	die if ($head[1] ne 'Protein B');
	$symbol2 = $a[1];
	$symbol2 = uc($symbol2);
	# $accs2 = '';

	# Apply symbol alias if there is one (so I can look up additional information, like from HPA, on e.g. the RBM10 mutants)
	$originalsymbol2 = $symbol2;
	$symbol2 = $symbol{$symbol2} if (exists($symbol{$symbol2}));
	
	print " >> SYMBOL $symbol2\n" if (switch('debug'));
	print "   >> ORIGINAL $symbol2\n" if (switch('debug'));

	# ...get its original species
	@species = ();
	@names = ();
	@accs = ();
	$query = Query("SELECT species1 FROM $table WHERE symbol1='$originalsymbol2' GROUP BY species1 ORDER BY species1");
	while (($species) = Fetch($query))
	{
		push(@species, $species);

		print "     >> SPECIES $species\n" if (switch('debug'));
	}
	$query = Query("SELECT species2 FROM $table WHERE symbol2='$originalsymbol2' GROUP BY species2 ORDER BY species2");
	while (($species) = Fetch($query))
	{
		push(@species, $species);

		print "     >> SPECIES $species\n" if (switch('debug'));
	}
	# die if (scalar(@species) != unique(@species));
	@species = unique(@species);
	# # Don't unique @species, just sort it. There are cases like RBFOX2 that have two HUMAN proteins, for instance. Uniquing breaks these.
	# @species = nsort(@species);

	foreach $species (@species)
	{
		print "       >> SPECIES $species\n" if (switch('debug'));

		# Get UniProt name for this species
		$query = Query("SELECT name FROM uniprot WHERE species='$species' AND gene='$symbol2'");
		if (Numrows($query) == 0)
		{
			if (!exists($acc{"$species|$symbol2"}))
			{
				warn("Warning: No uniprot name in table 'uniprot' for species|symbol '$species|$symbol2' (needs to be mapped manually to an accession using the \%acc hash) (kept)");
				addme("no uniprot name in table 'uniprot' for species|symbol (needs to be mapped manually to an accession using the \%acc hash) (kept)", "$species|$symbol2");
			}
			else
			{
				# Use manual %acc hash for mapping (defined above)
				foreach $acc (split(/;/, $acc{"$species|$symbol2"}))
				{
					push(@accs, $acc);

					# Get name, too
					$query = Query("SELECT name FROM uniacc WHERE acc='$acc' AND species='$species'");
					if (Numrows($query) == 0)
					{
						# Not in Swiss-Prot
						$name = $acc.'_'.$species;
					}
					else
					{
						($name) = FetchOne($query);
					}
					push(@names, $name);
				}
			}
		}
		else
		{
			while (($name) = Fetch($query))
			{
				push(@names, $name);
				print "         >> NAME $name\n" if (switch('debug'));

				# Get UniProt accession
				$query = Query("SELECT DISTINCT primary_acc FROM uniacc WHERE name='$name' AND species='$species'");
				($acc) = FetchOne($query);
				push(@accs, $acc);
				print "         >> ACC $acc\n" if (switch('debug'));
			}
		}
	}

	$species2 = join(";", @species);
	$accs2 = join(";", @accs);
	$names2 = join(";", @names);
	print "           >> SPECIES2 $species2\n" if (switch('debug'));
	print "           >> ACCS2 $accs2\n" if (switch('debug'));
	print "           >> NAMES2 $names2\n" if (switch('debug'));
	$a[3] = $species2;	# 4
	$a[5] = $accs2;		# 6
	$a[7] = $names2;	# 8






	# Update Biogrid_all direct evidence
	$i = 20;
	die if ($head[$i] ne "\"Biogrid_all direct evidence (Y2H, reconstituted complex, structure)\"");
	$query = Query("SELECT id FROM biogrid WHERE (gene1='$symbol1' AND gene2='$symbol2') OR (gene1='$symbol2' AND gene2='$symbol1') AND direct=1");
	$biogrid_direct_before = $a[$i];
	$biogrid_direct = '0';
	if (Numrows($query) > 0)
	{
		$biogrid_direct = '1';
	}
	if ($biogrid_direct ne $biogrid_direct_before)
	{
		addme("updated biogrid_direct assignment from $biogrid_direct_before to $biogrid_direct for pair", "$symbol1|$symbol2");

		# Update table
		if (switch('update'))
		{
			$query = Query("UPDATE $table SET bg_alldirect='$biogrid_direct' WHERE symbol1='$originalsymbol1' AND symbol2='$originalsymbol2'");
			die("Error: Couldn't update 'bg_alldirect' field in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) == 0);
		}
	}
	else
	{
		addme("biogrid_direct assignment remains $biogrid_direct_before for pair", "$symbol1|$symbol2");
	}
	$biogrid_direct = "biogrid_direct" if (switch('debug2'));
	$a[$i] = $biogrid_direct;





	
	# Update Biogrid
	$i = 21;
	die if ($head[$i] ne "Biogrid_all");
	$query = Query("SELECT id FROM biogrid WHERE (gene1='$symbol1' AND gene2='$symbol2') OR (gene1='$symbol2' AND gene2='$symbol1')");
	$biogrid_before = $a[$i];
	$biogrid = '0';
	if (Numrows($query) > 0)
	{
		$biogrid = '1';
	}
	if ($biogrid ne $biogrid_before)
	{
		addme("updated biogrid assignment from $biogrid_before to $biogrid for pair", "$symbol1|$symbol2");

		# Update table
		if (switch('update'))
		{
			$query = Query("UPDATE $table SET bg_all='$biogrid' WHERE symbol1='$originalsymbol1' AND symbol2='$originalsymbol2'");
			die("Error: Couldn't update 'bg_all' field in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) == 0);
		}
	}
	else
	{
		addme("biogrid assignment remains $biogrid_before for pair", "$symbol1|$symbol2");
	}
	$biogrid = "biogrid" if (switch('debug2'));
	$a[$i] = $biogrid;





	
	# Update HIPPIE
	$i = 22;
	die if ($head[$i] ne "Hippie");
	$query = Query("SELECT id FROM hippie WHERE (symbol1='$symbol1' AND symbol2='$symbol2') OR (symbol1='$symbol2' AND symbol2='$symbol1')");
	$hippie_before = $a[$i];
	$hippie = '0';
	if (Numrows($query) > 0)
	{
		$hippie = '1';
	}
	if ($hippie ne $hippie_before)
	{
		addme("updated hippie assignment from $hippie_before to $hippie for pair", "$symbol1|$symbol2");
		
		# Update table
		if (switch('update'))
		{
			$query = Query("UPDATE $table SET hippie='$hippie' WHERE symbol1='$originalsymbol1' AND symbol2='$originalsymbol2'");
			die("Error: Couldn't update 'hippie' field in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) == 0);
		}
	}
	else
	{
		addme("hippie assignment remains $hippie_before for pair", "$symbol1|$symbol2");
	}
	$hippie = "hippie" if (switch('debug2'));
	$a[$i] = $hippie;





	
	# Update HuRI (it's mostly correct, but it's missing a few pairs. Apparently older than Nov 2017, or not mapped successfully for all pairs)
	$i = 23;
	die if ($head[$i] ne "HuRI");
	# Using only type='huri', i.e. not the HI-Union dataset (which includes other databases)
	$query = Query("SELECT h.id FROM huri h, gencode_gff3_gene g1, gencode_gff3_gene g2 WHERE h.type='huri' AND g1.species='human' AND g2.species='human' AND ((g1.symbol='$symbol1' AND g2.symbol='$symbol2') OR (g1.symbol='$symbol2' AND g2.symbol='$symbol1')) AND h.ensg1=g1.ensg AND h.ensg2=g2.ensg");
	$huri_before = $a[$i];
	$huri = '0';
	if (Numrows($query) > 0)
	{
		$huri = '1';
	}
	if ($huri ne $huri_before)
	{
		addme("updated huri assignment from $huri_before to $huri for pair", "$symbol1|$symbol2");
		
		# Update table
		if (switch('update'))
		{
			$query = Query("UPDATE $table SET huri='$huri' WHERE symbol1='$originalsymbol1' AND symbol2='$originalsymbol2'");
			die("Error: Couldn't update 'huri' field in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) == 0);
		}
	}
	else
	{
		addme("huri assignment remains $huri_before for pair", "$symbol1|$symbol2");
	}
	$huri = "huri" if (switch('debug2'));
	$a[$i] = $huri;





	
	# Known interaction yes/no
	$i = 24;
	die if ($head[$i] ne "Known interaction");
	die if ($head[$i-4] ne "\"Biogrid_all direct evidence (Y2H, reconstituted complex, structure)\"");
	die if ($head[$i-3] ne "Biogrid_all");
	die if ($head[$i-2] ne "Hippie");
	die if ($head[$i-1] ne "HuRI");
	$known = 0;
	if (($a[$i-4] == 1) or ($a[$i-3] == 1) or ($a[$i-2] == 1) or ($a[$i-1] == 1))
	{
		$known = 1;
	}
	# Update table
	if (switch('update'))
	{
		$query = Query("UPDATE $table SET known='$known' WHERE symbol1='$originalsymbol1' AND symbol2='$originalsymbol2'");
		die("Error: Couldn't update 'known' field in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) == 0);

		$query = Query("SELECT id FROM $table WHERE symbol2='$originalsymbol1' AND symbol1='$originalsymbol2'");
		die("Error: Inverse pair exists in table '$table' for pair '$originalsymbol1|$originalsymbol2'") if (Numrows($query) != 0);
	}
	$known = "known" if (switch('debug2'));
	$a[$i] = $known;
	
	
	

	# Known complex member A (protein has at least one known interaction) yes/no
	$i = 25;
	die if ($head[$i] ne "Protein A has known interactions");
	$known_protein_1 = 0;
	if (($known == 1) or exists($known_protein{$originalsymbol1}))
	{
		$known_protein_1 = 1;
		$known_protein{$originalsymbol1} = 1;
	}
	# Update table
	if (switch('update'))
	{
		# Update where this protein is symbol1
		$query = Query("UPDATE $table SET knownprotein1='$known_protein_1' WHERE symbol1='$originalsymbol1'");
		die("Error: Couldn't update 'knownprotein1' field in table '$table' for symbol '$originalsymbol1'") if (Numrows($query) == 0);
		# Also update where this protein is symbol2
		$query = Query("UPDATE $table SET knownprotein2='$known_protein_1' WHERE symbol2='$originalsymbol1'");
	}
	$known_protein_1 = "known_protein_1" if (switch('debug2'));
	splice(@a, $i, 0, $known_protein_1);	# Add a column	
	
	

	# Known complex member B (protein has at least one known interaction) yes/no
	$i = 26;
	die if ($head[$i] ne "Protein B has known interactions");
	$known_protein_2 = 0;
	if (($known == 1) or exists($known_protein{$originalsymbol2}))
	{
		$known_protein_2 = 1;
		$known_protein{$originalsymbol2} = 1;
	}
	# Update table
	if (switch('update'))
	{
		# Update where this protein is symbol2
		$query = Query("UPDATE $table SET knownprotein2='$known_protein_2' WHERE symbol2='$originalsymbol2'");
		die("Error: Couldn't update 'knownprotein2' field in table '$table' for symbol '$originalsymbol2'") if (Numrows($query) == 0);
		# Also update where this protein is symbol1
		$query = Query("UPDATE $table SET knownprotein1='$known_protein_2' WHERE symbol1='$originalsymbol2'");
	}
	$known_protein_2 = "known_protein_2" if (switch('debug2'));
	splice(@a, $i, 0, $known_protein_2);	# Add a column	
	
	
	

	# HPA nuclear A
	$i = 27;
	die if ($head[$i] ne "HPA nuclear A");
	$query = Query("SELECT MAX(nuclear) FROM hpa_subcell_loc WHERE gene_name='$symbol1'");
	$hpa_nuclear_1 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_nuclear_1) = FetchOne($query);
		$hpa_nuclear_1 = 'N/A' if (!defined($hpa_nuclear_1));
	}
	$hpa_nuclear_1 = "hpa_nuclear_1" if (switch('debug2'));
	$a[$i] = $hpa_nuclear_1;

	# HPA nuclear B
	$i = 28;
	die if ($head[$i] ne "HPA nuclear B");
	$query = Query("SELECT MAX(nuclear) FROM hpa_subcell_loc WHERE gene_name='$symbol2'");
	$hpa_nuclear_2 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_nuclear_2) = FetchOne($query);
		$hpa_nuclear_2 = 'N/A' if (!defined($hpa_nuclear_2));
	}
	# show(@a);
	# state("SPLICE");
	$hpa_nuclear_2 = "hpa_nuclear_2" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_nuclear_2);	# Add a column
	# show(@a);
	# exit;

	# HPA cytoplasmic A
	$i = 29;
	die if ($head[$i] ne "HPA cytoplasmic A");
	# MAX() because ABCF2 has two ENSGs associated with it, and one is "Enhanced Cytosol" while the other is "Uncertain Cytosol", which doesn't count.
	$query = Query("SELECT MAX(cytoplasmic) FROM hpa_subcell_loc WHERE gene_name='$symbol1'");
	$hpa_cytoplasmic_1 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_cytoplasmic_1) = FetchOne($query);
		$hpa_cytoplasmic_1 = 'N/A' if (!defined($hpa_cytoplasmic_1));
	}
	$hpa_cytoplasmic_1 = "hpa_cytoplasmic_1" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_cytoplasmic_1);	# Add a column

	# HPA cytoplasmic B
	$i = 30;
	die if ($head[$i] ne "HPA cytoplasmic B");
	$query = Query("SELECT MAX(cytoplasmic) FROM hpa_subcell_loc WHERE gene_name='$symbol2'");
	$hpa_cytoplasmic_2 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_cytoplasmic_2) = FetchOne($query);
		$hpa_cytoplasmic_2 = 'N/A' if (!defined($hpa_cytoplasmic_2));
	}
	$hpa_cytoplasmic_2 = "hpa_cytoplasmic_2" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_cytoplasmic_2);	# Add a column
	# show(@a);
	
	# HPA main_location A
	$i = 31;
	die if ($head[$i] ne "HPA main locations A");
	$query = Query("SELECT DISTINCT main_location FROM hpa_subcell_loc WHERE gene_name='$symbol1' AND main_location IS NOT NULL");
	$hpa_main_locations_1 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_main_locations_1) = FetchOne($query);
	}
	$hpa_main_locations_1 = "hpa_main_locations_1" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_main_locations_1);	# Add a column

	# HPA additional_location A
	$i = 32;
	die if ($head[$i] ne "HPA additional locations A");
	$query = Query("SELECT DISTINCT additional_location FROM hpa_subcell_loc WHERE gene_name='$symbol1' AND additional_location IS NOT NULL");
	$hpa_additional_locations_1 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_additional_locations_1) = FetchOne($query);
	}
	$hpa_additional_locations_1 = "hpa_additional_locations_1" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_additional_locations_1);	# Add a column

	# HPA main_location B
	$i = 33;
	die if ($head[$i] ne "HPA main locations B");
	$query = Query("SELECT DISTINCT main_location FROM hpa_subcell_loc WHERE gene_name='$symbol2' AND main_location IS NOT NULL");
	$hpa_main_locations_2 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_main_locations_2) = FetchOne($query);
	}
	$hpa_main_locations_2 = "hpa_main_locations_2" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_main_locations_2);	# Add a column

	# HPA additional_location B
	$i = 34;
	die if ($head[$i] ne "HPA additional locations B");
	$query = Query("SELECT DISTINCT additional_location FROM hpa_subcell_loc WHERE gene_name='$symbol2' AND additional_location IS NOT NULL");
	$hpa_additional_locations_2 = 'N/A';
	if (Numrows($query) > 0)
	{
		($hpa_additional_locations_2) = FetchOne($query);
	}
	$hpa_additional_locations_2 = "hpa_additional_locations_2" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_additional_locations_2);	# Add a column

	# HPA shared locations
	$i = 35;
	die if ($head[$i] ne "HPA shared locations");
	$hpa_shared_locations = 'N/A';
	@hpa_locs_main_1 = split(/;/, $hpa_main_locations_1);
	@hpa_locs_main_2 = split(/;/, $hpa_main_locations_2);
	@hpa_locs_additional_1 = split(/;/, $hpa_additional_locations_1);
	@hpa_locs_additional_2 = split(/;/, $hpa_additional_locations_2);
	@hpa_locs_1 = (@hpa_locs_main_1, @hpa_locs_additional_1);
	@hpa_locs_2 = (@hpa_locs_main_2, @hpa_locs_additional_2);
	@hpa_locs_shared_raw = intersection(\@hpa_locs_1, \@hpa_locs_2);
	@hpa_locs_shared = ();
	foreach $tmploc (@hpa_locs_shared_raw)
	{
		push(@hpa_locs_shared, $tmploc) if ($tmploc ne 'N/A');
	}
	if (scalar(@hpa_locs_shared) == 0)
	{
		push(@hpa_locs_shared, "N/A");
	}
	@hpa_locs_shared = unique(@hpa_locs_shared);
	$hpa_shared_locations = join(';', @hpa_locs_shared);
	$hpa_shared_locations = "hpa_shared_locations" if (switch('debug2'));
	splice(@a, $i, 0, $hpa_shared_locations);	# Add a column

	# splice(@head, 28, 0, "HPA main locations A");
	# splice(@head, 29, 0, "HPA main locations B");
	# splice(@head, 30, 0, "HPA additional locations A");
	# splice(@head, 31, 0, "HPA additional locations B");
	# splice(@head, 32, 0, "HPA shared locations");
	



	# BioID
	$i = 36;
	die if ($head[$i] ne "Youn-Gingras_data");
	$query = Query("SELECT id FROM bioid WHERE ((symbol1='$symbol1' AND symbol2='$symbol2') OR (symbol1='$symbol2' AND symbol2='$symbol1')) AND fdr<0.05");
	$bioid = 0;
	if (Numrows($query) > 0)
	{
		$bioid = 1;
	}
	$bioid = 'bioid' if (switch('debug2'));
	splice(@a, $i, 0, $bioid);	# Add a column
	
	
	
	# # Store the NanoBRET pos/neg value
	# $nanobret = $a[37];
	
	
	
	# eCLIP data yes/no A
	$i = 40;
	die if ($head[$i] ne "ENCODE eCLIP data A");
	# $query = Query("SELECT id FROM $table WHERE (symbol1='$symbol1' AND eclip1=1) OR (symbol2='$symbol1' AND eclip2=1) LIMIT 1");
	$query = Query("SELECT COUNT(DISTINCT chr, start, stop, strand) FROM clip_raw_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$symbol1'");
	$eclip1 = 0;
	($eclip_sites_1) = FetchOne($query);
	if ($eclip_sites_1 > 0)
	{
		$eclip1 = 1;
	}
	$eclip1 = "eclip1" if (switch('debug2'));
	$a[$i] = $eclip1;
	
	# eCLIP data yes/no B
	$i = 41;
	die if ($head[$i] ne "ENCODE eCLIP data B");
	# $query = Query("SELECT id FROM $table WHERE (symbol1='$symbol2' AND eclip1=1) OR (symbol2='$symbol2' AND eclip2=1) LIMIT 1");
	$query = Query("SELECT COUNT(DISTINCT chr, start, stop, strand) FROM clip_raw_gene WHERE species='human' AND type='$type' AND map='gene' AND symbol='$symbol2'");
	$eclip2 = 0;
	($eclip_sites_2) = FetchOne($query);
	if ($eclip_sites_2 > 0)
	{
		$eclip2 = 1;
	}
	$eclip2 = "eclip2" if (switch('debug2'));
	$a[$i] = $eclip2;
	
	
	# eCLIP binding sites A
	$i = 42;
	die if ($head[$i] ne "ENCODE eCLIP binding sites A");
	$eclip_sites_1 = "eclip_sites_1" if (switch('debug2'));
	$a[$i] = $eclip_sites_1;
	
	# eCLIP binding sites B
	$i = 43;
	die if ($head[$i] ne "ENCODE eCLIP binding sites B");
	$eclip_sites_2 = "eclip_sites_2" if (switch('debug2'));
	$a[$i] = $eclip_sites_2;
	
	
	
	

	# Jaccard
	$i = 44;
	die if ($head[$i] ne "Jaccard index");
	# $a[30] = get_cobind("jaccard", $symbol1, $symbol2);
	$jaccard = 'N/A';
	$jaccard = $jaccard{"$symbol1|$symbol2"} if (exists($jaccard{"$symbol1|$symbol2"}));
	$a[$i] = $jaccard;
	
	# probability
	$i = 45;
	die if ($head[$i] ne "Cobinding prob p(A|B)");
	# $a[31] = get_cobind("probability", $symbol1, $symbol2);
	$probability = 'N/A';
	$probability = $probability{"$symbol1|$symbol2"} if (exists($probability{"$symbol1|$symbol2"}));
	$a[$i] = $probability;
	
	# probability
	$i = 46;
	die if ($head[$i] ne "Cobinding prob p(B|A)");
	# $a[32] = get_cobind("probability", $symbol1, $symbol2);
	$probability = 'N/A';
	$probability = $probability{"$symbol2|$symbol1"} if (exists($probability{"$symbol2|$symbol1"}));
	$a[$i] = $probability;
	
	# probability_mindist
	$i = 47;
	die if ($head[$i] ne "Cobinding prob p(A|B) (≤$mindist_threshold nt)");
	# $a[31] = get_cobind("probability_mindist", $symbol1, $symbol2);
	$probability_mindist = 'N/A';
	$probability_mindist = $probability_mindist{"$symbol1|$symbol2"} if (exists($probability_mindist{"$symbol1|$symbol2"}));
	$a[$i] = $probability_mindist;
	
	# probability_mindist
	$i = 48;
	die if ($head[$i] ne "Cobinding prob p(B|A) (≤$mindist_threshold nt)");
	# $a[32] = get_cobind("probability_mindist", $symbol1, $symbol2);
	$probability_mindist = 'N/A';
	$probability_mindist = $probability_mindist{"$symbol2|$symbol1"} if (exists($probability_mindist{"$symbol2|$symbol1"}));
	$a[$i] = $probability_mindist;
	
	# # Odds ratio
	# $i = 47;
	# die if ($head[$i] ne "Cobinding odds ratio");
	# # $a[33] = get_cobind("oddsratio", $symbol1, $symbol2);
	# $oddsratio = 'N/A';
	# $oddsratio = $oddsratio{"$symbol1|$symbol2"} if (exists($oddsratio{"$symbol1|$symbol2"}));
	# $a[$i] = $oddsratio;
	
	
	
	
	# Get information from output-fit_vs_random-$mindist_threshold.txt (resampling of the RBP pair binding distances)
	# symbol1	symbol2	close_real	fraction_close_real	fraction_close_random	fraction_ratio	resampling_p	resampling_wilcox
	# APOBEC3C	HNRNPK	24	0.04040404	0.03871021	1.04375667298111	0.45	1
	# HNRNPK	APOBEC3C	24	0.01505646	0.01517115	0.99244025667138	0.57	1
	# APOBEC3C	PCBP1	221	0.2036866	0.05850176	3.48171747311534	0	0
	# PCBP1	APOBEC3C	277	0.1198615	0.02868432	4.17864185032101	0	0
	
	# A vs. B
	
	$close_real_1 = 'N/A';
	$fraction_close_real_1 = 'N/A';
	$fraction_close_random_1 = 'N/A';
	$fraction_ratio_1 = 'N/A';
	$resampling_p_1 = 'N/A';
	$resampling_wilcox_1 = 'N/A';

	if (exists($fit_vs_random{"$symbol1|$symbol2"}))
	{
		$tmp = $fit_vs_random{"$symbol1|$symbol2"};
		@tmp_a = split(/\t/, $tmp);
		
		$close_real_1 = $tmp_a[30];
		$fraction_close_real_1 = $tmp_a[31];
		$fraction_close_random_1 = $tmp_a[32];
		$fraction_ratio_1 = $tmp_a[34];
		$resampling_p_1 = $tmp_a[35];
		$resampling_wilcox_1 = $tmp_a[36];
	}
	
	$a[49] = $close_real_1;
	$a[50] = $fraction_close_real_1;
	$a[51] = $fraction_close_random_1;
	$a[52] = $fraction_ratio_1;
	$a[53] = $resampling_p_1;
	$a[54] = $resampling_wilcox_1;
	

	# B vs. A
	
	$close_real_2 = 'N/A';
	$fraction_close_real_2 = 'N/A';
	$fraction_close_random_2 = 'N/A';
	$fraction_ratio_2 = 'N/A';
	$resampling_p_2 = 'N/A';
	$resampling_wilcox_2 = 'N/A';

	if (exists($fit_vs_random{"$symbol2|$symbol1"}))
	{
		$tmp = $fit_vs_random{"$symbol2|$symbol1"};
		@tmp_a = split(/\t/, $tmp);
	
		$close_real_2 = $tmp_a[30];
		$fraction_close_real_2 = $tmp_a[31];
		$fraction_close_random_2 = $tmp_a[32];
		$fraction_ratio_2 = $tmp_a[34];
		$resampling_p_2 = $tmp_a[35];
		$resampling_wilcox_2 = $tmp_a[36];
	}
	
	$a[55] = $close_real_2;
	$a[56] = $fraction_close_real_2;
	$a[57] = $fraction_close_random_2;
	$a[58] = $fraction_ratio_2;
	$a[59] = $resampling_p_2;
	$a[60] = $resampling_wilcox_2;
	
	
	
	# #DEBUG
	# show(@a);
	# exit;
	# #END DEBUG
	
	
	
	
	
	# Print
	# if (($symbol1 !~ /_/) and ($symbol2 !~ /_/))
	if (($originalsymbol1 !~ /_/) and ($originalsymbol2 !~ /_/))
	{
		# Normal proteins (main, "clean" master table)
		print OUT join("\t", @a)."\n";
		$clean++;
		


		# Additionally: Print pairs that have significant resampling_p and resampling_wilcox values to the Cytoscape edge file for Figure 5a
		if (($eclip1 == 1) and ($eclip2 == 1))
		{
			if (($resampling_p_1 ne 'N/A') and ($resampling_wilcox_1 ne 'N/A') and ($resampling_p_2 ne 'N/A') and ($resampling_wilcox_2 ne 'N/A'))
			{
				# # Stricter method: are all of the resampling and Wilcoxon p-values <0.05? (leads to 21 edges)
				# if (($resampling_p_1 < $alpha) and ($resampling_wilcox_1 < $alpha) and ($resampling_p_2 < $alpha) and ($resampling_wilcox_2 < $alpha))
				# Stricter method: are the resampling and Wilcoxon p-values <0.05 in at least one direction? (leads to ? edges)
				if ((($resampling_p_1 < $alpha) and ($resampling_wilcox_1 < $alpha)) or (($resampling_p_2 < $alpha) and ($resampling_wilcox_2 < $alpha)))
				{
					print STRICTEDGE join("\t", @a)."\t1\n";
					$strictedge++;

					print EDGE join("\t", @a)."\t1\n";
					$edge++;
					
					# Node annotation
					push(@nodesymbols, $symbol1);
					push(@nodesymbols, $symbol2);
					push(@strictnodesymbols, $symbol1);
					push(@strictnodesymbols, $symbol2);
					$node_hpa_nuclear{$symbol1} = $hpa_nuclear_1;
					$node_hpa_nuclear{$symbol2} = $hpa_nuclear_2;
					$node_hpa_cytoplasmic{$symbol1} = $hpa_cytoplasmic_1;
					$node_hpa_cytoplasmic{$symbol2} = $hpa_cytoplasmic_2;
					# node_known_protein
					if ($known == 1)
					{
						$strict_node_known_protein{$symbol1} = 1;
						$strict_node_known_protein{$symbol2} = 1;
						$node_known_protein{$symbol1} = 1;
						$node_known_protein{$symbol2} = 1;
					}
				}
				# More lenient method: is one of the resampling p-values <0.05? (leads to  edges)
				elsif (($resampling_p_1 < $alpha) or ($resampling_p_2 < $alpha))
				{
					print EDGE join("\t", @a)."\t0\n";
					$edge++;

					# Node annotation
					push(@nodesymbols, $symbol1);
					push(@nodesymbols, $symbol2);
					$node_hpa_nuclear{$symbol1} = $hpa_nuclear_1;
					$node_hpa_nuclear{$symbol2} = $hpa_nuclear_2;
					$node_hpa_cytoplasmic{$symbol1} = $hpa_cytoplasmic_1;
					$node_hpa_cytoplasmic{$symbol2} = $hpa_cytoplasmic_2;
					# node_known_protein
					if ($known == 1)
					{
						$node_known_protein{$symbol1} = 1;
						$node_known_protein{$symbol2} = 1;
					}
				}
			}
		}
	}
	elsif (($symbol1 =~ /^RBM10_/) or ($symbol1 =~ /^SRRM3_/) or ($symbol1 =~ /^SRRM4_/) or
		   ($symbol2 =~ /^RBM10_/) or ($symbol2 =~ /^SRRM3_/) or ($symbol2 =~ /^SRRM4_/))
	{
		# Mutants (supplementary table)
		print MUT join("\t", @a)."\n";
		$mut++;
	}
	else
	{
		# Fragments table (not used anywhere) (should only be the fragments)
		print FRAG join("\t", @a)."\n";
		$frag++;
	}
	
	
	
	# Also add the wild-type RBM10, SRRM3, and SRRM4 to the mutant table (they've already been written to the "clean" master table above)
	if (($symbol1 eq 'RBM10') or ($symbol1 eq 'SRRM3') or ($symbol1 eq 'SRRM4') or
	    ($symbol2 eq 'RBM10') or ($symbol2 eq 'SRRM3') or ($symbol2 eq 'SRRM4'))
	{
		# Mutants (supplementary table)
		print MUT join("\t", @a)."\n";
		$mut++;
	}
	
	
	
	# show(@a);
	# exit;
	
	
	stepme(100);
}
stopme();
stoptime();



# Write network node annotations
startme("Writing node annotation to '$nodefile'");
foreach $symbol (unique(@nodesymbols))
{
	# Set location
	$node_hpa_shuttling_loc = 'N/A';
	if ($node_hpa_nuclear{$symbol} == 1)
	{
		if ($node_hpa_cytoplasmic{$symbol} == 1)
		{
			$node_hpa_shuttling_loc = 'shuttling';
		}
		else
		{
			$node_hpa_shuttling_loc = 'nuclear';
		}
	}
	elsif ($node_hpa_cytoplasmic{$symbol} == 1)
	{
		$node_hpa_shuttling_loc = 'cytoplasmic';
	}
	
	# Set node_known_protein
	$node_known_protein = 0;
	if (exists($node_known_protein{$symbol}))
	{
		$node_known_protein = $node_known_protein{$symbol};
	}
	
	print NODE "$symbol\t$node_known_protein\t".$node_hpa_nuclear{$symbol}."\t".$node_hpa_cytoplasmic{$symbol}."\t".$node_hpa_shuttling_loc."\n";
	$node++;
	stepme(10);
}
stopme();



# Write network node annotations
startme("Writing strict node annotation to '$strictnodefile'");
foreach $symbol (unique(@strictnodesymbols))
{
	# Set location
	$node_hpa_shuttling_loc = 'N/A';
	if ($node_hpa_nuclear{$symbol} == 1)
	{
		if ($node_hpa_cytoplasmic{$symbol} == 1)
		{
			$node_hpa_shuttling_loc = 'shuttling';
		}
		else
		{
			$node_hpa_shuttling_loc = 'nuclear';
		}
	}
	elsif ($node_hpa_cytoplasmic{$symbol} == 1)
	{
		$node_hpa_shuttling_loc = 'cytoplasmic';
	}
	
	# Set node_known_protein
	$node_known_protein = 0;
	if (exists($strict_node_known_protein{$symbol}))
	{
		$node_known_protein = $strict_node_known_protein{$symbol};
	}
	
	print STRICTNODE "$symbol\t$node_known_protein\t".$node_hpa_nuclear{$symbol}."\t".$node_hpa_cytoplasmic{$symbol}."\t".$node_hpa_shuttling_loc."\n";
	$strictnode++;
	stepme(10);
}
stopme();






showmeallsorted(1);
# showme("no hgnc symbol for symbol (kept)");

nl();
state("Wrote ".commify($clean)." lines with clean proteins to '$outfile'", 1);
state("Wrote ".commify($mut)." lines with mutants (and the corresponding wild-type proteins) to '$mutfile'", 1);
state("Wrote ".commify($frag)." lines with fragments etc. to '$fragfile'", 1);
state("Wrote ".commify($edge)." lines with network edges for Cytoscape (Figure 5a) to '$edgefile'", 1);
state("Wrote ".commify($strictedge)." lines with strictly filtered network edges for Cytoscape (Figure 5a) to '$strictedgefile'", 1);
state("Wrote network annotation for ".commify($node)." genes for Cytoscape (Figure 5a) to '$nodefile'", 1);
state("Wrote network annotation for ".commify($strictnode)." strictly filtered genes for Cytoscape (Figure 5a) to '$strictnodefile'", 1);
nl();

done();





# sub get_cobind
# {
# 	my ($mode, $rbp1, $rbp2) = @_;
#
# 	my $tmpfile = "../../pipeline/rbpome_analysis/backup_2020_05_11_rbpome_results/output/output-txt-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$threshold-mindist_threshold$mindist_threshold-resamples$resamples.txt";
# 	open(TMP, $tmpfile) or die("Error: Couldn't open '$tmpfile'");
# 	<TMP>;	# Skip header
# 	my $res = '';
# 	while (<TMP>)
# 	{
# 		chomp;
#
# 		# table	type	scoretype	mode	hc	minlog2fold	minscore	mindist_threshold	resamples	set	rbp1	rbp2	cobind
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	PCBP1	0.16151090849886
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	EWSR1	SF3B4	0.135072908672295
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FMR1	FXR2	0.420327027893556
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	FXR1	FXR2	0.15188679245283
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	HNRNPK	U2AF2	0.234877126654064
# 		# rbpome	eclip_encode	AVG_IS	jaccard	0	0	3.25	400	100	positives	NONO	SFPQ	0.203641590800192
#
# 		my @a = split(/\t/);
# 		my $set = $a[9];
# 		my $this_rbp1 = $a[10];
# 		my $this_rbp2 = $a[11];
# 		my $cobind = $a[12];
#
# 		# Only use "screen_hits" (not positive controls, not background_resamples)
# 		next if ($set ne 'screen_hits');
#
# 		if (($this_rbp1 eq $rbp1) and ($this_rbp2 eq $rbp2))
# 		{
# 			$res = $cobind;
# 		}
# 	}
# 	close(TMP);
#
# 	return($res);
# }
