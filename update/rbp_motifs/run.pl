#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

#DEBUG
# waitforalljobs();
#END DEBUG

# args(0);
# $extend5 = 0;
# if (switch('extend5'))
# {
# 	$extend5 = 1;
# }

# Load grep motifs
Clear('rbp_motif_list');
run("Load Dominguez et al. (Chris Burge RNA Bind-N-Seq) grep motifs into table rbp_motif_list", "load_motifs_dominguez.pl");
run("Load ATTRACT grep motifs into table rbp_motif_list", "load_motifs_attract.pl");
run("Load my custom grep motifs (e.g. from UniProt for YBX3) into table rbp_motif_list", "load_motifs_custom.pl");

# run("Remove motif import files", "rm -f tmp-rbp_motifs*.txt");



# Get individual motif files
@fimosources = ('attract', 'cisbp', 'dominguez', 'rbpdb', 'rbpmap', 'rnacompete');
foreach $source (@fimosources)
{
	run("Get individual motif files (per RBP)", "get_motifs.pl $source");
}




# # Side analysis: Gene sequences
# # The question being: If I run FIMO on whole gene sequences, do I still get any hits at all?
#
# run("Get individual motif files for each RBP", "get_motifs.pl all");
# run("Get geneseqs", "get_geneseqs.pl");
# run("Grep (geneseqs)", "grep_geneseqs.pl");
# # run("Run FIMO (geneseqs)", "fimo_geneseqs.pl 0.05 -qvalue");
# run("Run FIMO (geneseqs)", "fimo_geneseqs.pl 0.001");
# # run("Parse FIMO (geneseqs)", "parse_fimo_geneseqs.pl fimo 0.05 -qvalue");
# run("Parse FIMO (geneseqs)", "parse_fimo_geneseqs.pl fimo 0.001");
#
# # output-fimo-geneseqs-AKAP1-0_05-qvalue.txt
# # parse_fimo.pl $method $type $source 0.001$tmpextend5
#
# # Conclusion:
# # "(Noch ein paar Zahlen um zu zeigen dass die Motif-Suche in den gesamten Genen keinen Sinn macht: FXR1, ein recht moderater Fall, hat ~400 ENCODE-Peaks mit Motifs. In den gesamten Genen findet FIMO mit p-value<0.001 ~318,000. Mit q-value<0.05 bleiben davon 0."
# # "Ja 😊 und mit „gesamte Gene“ meinte ich die komplette Sequenz nur von den Genen mit peaks!"
# # >> With a p-value threshold there are way too many hits, and none remain with a q-value threshold.





# Main analysis: Peak sequences


# foreach $extend5 (0, 1)
foreach $extend5 (1)
{
	if ($extend5 == 0)
	{
		$tmpextend5 = '';
		$tmpextend5_2 = '';
	}
	else
	{
		$tmpextend5 = ' -extend5';
		$tmpextend5_2 = '_extend5';
	}

	# @types = ('eclip_encode', 'eclip_tom', 'eclip_tom_local', 'eclip_tom_waldtest', 'eclip_tom_waldtest_local', 'eclip_encode_12', 'eclip_par_ohler');
	# @types = ('eclip_tom_waldtest', 'eclip_tom_waldtest_local');
	# @types = ('eclip_encode_12', 'eclip_par_ohler');
	# @types = ('eclip_encode_12');
	# @types = ('eclip_tom_width10_step10');
	# @types = ('eclip_tom_width10_step10', 'eclip_tom_width10_step5', 'eclip_tom_width20_step10', 'eclip_tom_width20_step5', 'eclip_tom_width30_step10', 'eclip_tom_width30_step5', 'eclip_tom_width40_step10', 'eclip_tom_width40_step5', 'eclip_tom_width50_step10', 'eclip_tom_width50_step5');
	# @types = ('eclip_tom_nocorr_width10_step10', 'eclip_tom_nocorr_width10_step5', 'eclip_tom_nocorr_width20_step10', 'eclip_tom_nocorr_width20_step5', 'eclip_tom_nocorr_width30_step10', 'eclip_tom_nocorr_width30_step5', 'eclip_tom_nocorr_width40_step10', 'eclip_tom_nocorr_width40_step5', 'eclip_tom_nocorr_width50_step10', 'eclip_tom_nocorr_width50_step5');
	# @types = (
	# 	'eclip_tom_width50_step10', 'eclip_tom_width50_step5', 'eclip_tom_width40_step10', 'eclip_tom_width40_step5', 'eclip_tom_width30_step10', 'eclip_tom_width30_step5', 'eclip_tom_width20_step10', 'eclip_tom_width20_step5', 'eclip_tom_width10_step10', 'eclip_tom_width10_step5',
	# 	'eclip_tom_nocorr_width50_step10', 'eclip_tom_nocorr_width50_step5', 'eclip_tom_nocorr_width40_step10', 'eclip_tom_nocorr_width40_step5', 'eclip_tom_nocorr_width30_step10', 'eclip_tom_nocorr_width30_step5', 'eclip_tom_nocorr_width20_step10', 'eclip_tom_nocorr_width20_step5', 'eclip_tom_nocorr_width10_step10', 'eclip_tom_nocorr_width10_step5',
	# 	'eclip_tom_relaxed_width50_step10', 'eclip_tom_relaxed_width50_step5', 'eclip_tom_relaxed_width40_step10', 'eclip_tom_relaxed_width40_step5', 'eclip_tom_relaxed_width30_step10', 'eclip_tom_relaxed_width30_step5', 'eclip_tom_relaxed_width20_step10', 'eclip_tom_relaxed_width20_step5', 'eclip_tom_relaxed_width10_step10', 'eclip_tom_relaxed_width10_step5',
	# 	'eclip_tom_lrt_nocorr_width50_step10', 'eclip_tom_lrt_nocorr_width50_step5', 'eclip_tom_lrt_nocorr_width40_step10', 'eclip_tom_lrt_nocorr_width40_step5', 'eclip_tom_lrt_nocorr_width30_step10', 'eclip_tom_lrt_nocorr_width30_step5', 'eclip_tom_lrt_nocorr_width20_step10', 'eclip_tom_lrt_nocorr_width20_step5', 'eclip_tom_lrt_nocorr_width10_step10', 'eclip_tom_lrt_nocorr_width10_step5',
	# 	'eclip_tom_lrt_relaxed_width50_step10', 'eclip_tom_lrt_relaxed_width50_step5', 'eclip_tom_lrt_relaxed_width40_step10', 'eclip_tom_lrt_relaxed_width40_step5', 'eclip_tom_lrt_relaxed_width30_step10', 'eclip_tom_lrt_relaxed_width30_step5', 'eclip_tom_lrt_relaxed_width20_step10', 'eclip_tom_lrt_relaxed_width20_step5', 'eclip_tom_lrt_relaxed_width10_step10', 'eclip_tom_lrt_relaxed_width10_step5'
	# );
	# @types = ('eclip_tom_width20_step20', 'eclip_tom_width30_step20', 'eclip_tom_width40_step20', 'eclip_tom_width50_step20', 'eclip_tom_nocorr_width20_step20', 'eclip_tom_nocorr_width30_step20', 'eclip_tom_nocorr_width40_step20', 'eclip_tom_nocorr_width50_step20', 'eclip_tom_relaxed_width20_step20', 'eclip_tom_relaxed_width30_step20', 'eclip_tom_relaxed_width40_step20', 'eclip_tom_relaxed_width50_step20', 'eclip_tom_lrt_nocorr_width20_step20', 'eclip_tom_lrt_nocorr_width30_step20', 'eclip_tom_lrt_nocorr_width40_step20', 'eclip_tom_lrt_nocorr_width50_step20', 'eclip_tom_lrt_relaxed_width20_step20', 'eclip_tom_lrt_relaxed_width30_step20', 'eclip_tom_lrt_relaxed_width40_step20', 'eclip_tom_lrt_relaxed_width50_step20');
	
	
	# @types = (
	# 'eclip_encode_12',
	# 't_lrt_ihw_nocorr_auto_w100_s5', # Best median peaks_with_motifs (FINAL)
	# 't_wald_ihw_corr_para_w100_s5', # Best compromise between high motif_fraction and high peak count
	# 't_lrt_bh_corr_para_w100_s5' # Best median motif_fraction
	# );

	# Varying window size 50/75/100:
	@types = (
	'eclip_encode_12',
	# Best median peaks_with_motifs (FINAL)
	't_lrt_ihw_nocorr_auto_w100_s5',
	't_lrt_ihw_nocorr_auto_w75_s5',
	't_lrt_ihw_nocorr_auto_w50_s5'
	);

	
	# # newparams
	# @types = ();
	# @bases =
	# (
	# 'n_lrt_bh_corr_auto',
	# 'n_lrt_bh_corr_para',
	# 'n_lrt_bh_nocorr_auto',
	# 'n_lrt_bh_nocorr_para',
	# 'n_lrt_ihw_corr_auto',
	# 'n_lrt_ihw_corr_para',
	# 'n_lrt_ihw_nocorr_auto',
	# 'n_lrt_ihw_nocorr_para',
	# 'n_wald_bh_corr_auto',
	# 'n_wald_bh_corr_para',
	# 'n_wald_bh_nocorr_auto',
	# 'n_wald_bh_nocorr_para',
	# 'n_wald_ihw_corr_auto',
	# 'n_wald_ihw_corr_para',
	# 'n_wald_ihw_nocorr_auto',
	# 'n_wald_ihw_nocorr_para'
	# );
	# # @widths = (10, 20, 30, 40, 50);
	# # @steps = (5, 10);
	# @widths = (50, 75, 100);
	# @steps = (5, 20);
	# foreach $base (@bases)
	# {
	# 	foreach $width (@widths)
	# 	{
	# 		foreach $step (@steps)
	# 		{
	# 			push(@types, "$base\_w$width\_s$step");
	# 			print "$base\_w$width\_s$step\n";
	# 		}
	# 	}
	# }
	# exit;
	# # @types = $types[0];
	
	# # All new types
	# @types = (
	# 't_lrt_bh_corr_auto_w50_s5',
	# 't_lrt_bh_corr_auto_w50_s20',
	# 't_lrt_bh_corr_auto_w75_s5',
	# 't_lrt_bh_corr_auto_w75_s20',
	# 't_lrt_bh_corr_auto_w100_s5',
	# 't_lrt_bh_corr_auto_w100_s20',
	# 't_lrt_bh_corr_para_w50_s5',
	# 't_lrt_bh_corr_para_w50_s20',
	# 't_lrt_bh_corr_para_w75_s5',
	# 't_lrt_bh_corr_para_w75_s20',
	# 't_lrt_bh_corr_para_w100_s5',
	# 't_lrt_bh_corr_para_w100_s20',
	# 't_lrt_bh_nocorr_auto_w50_s5',
	# 't_lrt_bh_nocorr_auto_w50_s20',
	# 't_lrt_bh_nocorr_auto_w75_s5',
	# 't_lrt_bh_nocorr_auto_w75_s20',
	# 't_lrt_bh_nocorr_auto_w100_s5',
	# 't_lrt_bh_nocorr_auto_w100_s20',
	# 't_lrt_bh_nocorr_para_w50_s5',
	# 't_lrt_bh_nocorr_para_w50_s20',
	# 't_lrt_bh_nocorr_para_w75_s5',
	# 't_lrt_bh_nocorr_para_w75_s20',
	# 't_lrt_bh_nocorr_para_w100_s5',
	# 't_lrt_bh_nocorr_para_w100_s20',
	# 't_lrt_ihw_corr_auto_w50_s5',
	# 't_lrt_ihw_corr_auto_w50_s20',
	# 't_lrt_ihw_corr_auto_w75_s5',
	# 't_lrt_ihw_corr_auto_w75_s20',
	# 't_lrt_ihw_corr_auto_w100_s5',
	# 't_lrt_ihw_corr_auto_w100_s20',
	# 't_lrt_ihw_corr_para_w50_s5',
	# 't_lrt_ihw_corr_para_w50_s20',
	# 't_lrt_ihw_corr_para_w75_s5',
	# 't_lrt_ihw_corr_para_w75_s20',
	# 't_lrt_ihw_corr_para_w100_s5',
	# 't_lrt_ihw_corr_para_w100_s20',
	# 't_lrt_ihw_nocorr_auto_w50_s5',
	# 't_lrt_ihw_nocorr_auto_w50_s20',
	# 't_lrt_ihw_nocorr_auto_w75_s5',
	# 't_lrt_ihw_nocorr_auto_w75_s20',
	# 't_lrt_ihw_nocorr_auto_w100_s5',
	# 't_lrt_ihw_nocorr_auto_w100_s20',
	# 't_lrt_ihw_nocorr_para_w50_s5',
	# 't_lrt_ihw_nocorr_para_w50_s20',
	# 't_lrt_ihw_nocorr_para_w75_s5',
	# 't_lrt_ihw_nocorr_para_w75_s20',
	# 't_lrt_ihw_nocorr_para_w100_s5',
	# 't_lrt_ihw_nocorr_para_w100_s20',
	# 't_wald_bh_corr_auto_w50_s5',
	# 't_wald_bh_corr_auto_w50_s20',
	# 't_wald_bh_corr_auto_w75_s5',
	# 't_wald_bh_corr_auto_w75_s20',
	# 't_wald_bh_corr_auto_w100_s5',
	# 't_wald_bh_corr_auto_w100_s20',
	# 't_wald_bh_corr_para_w50_s5',
	# 't_wald_bh_corr_para_w50_s20',
	# 't_wald_bh_corr_para_w75_s5',
	# 't_wald_bh_corr_para_w75_s20',
	# 't_wald_bh_corr_para_w100_s5',
	# 't_wald_bh_corr_para_w100_s20',
	# 't_wald_bh_nocorr_auto_w50_s5',
	# 't_wald_bh_nocorr_auto_w50_s20',
	# 't_wald_bh_nocorr_auto_w75_s5',
	# 't_wald_bh_nocorr_auto_w75_s20',
	# 't_wald_bh_nocorr_auto_w100_s5',
	# 't_wald_bh_nocorr_auto_w100_s20',
	# 't_wald_bh_nocorr_para_w50_s5',
	# 't_wald_bh_nocorr_para_w50_s20',
	# 't_wald_bh_nocorr_para_w75_s5',
	# 't_wald_bh_nocorr_para_w75_s20',
	# 't_wald_bh_nocorr_para_w100_s5',
	# 't_wald_bh_nocorr_para_w100_s20',
	# 't_wald_ihw_corr_auto_w50_s5',
	# 't_wald_ihw_corr_auto_w50_s20',
	# 't_wald_ihw_corr_auto_w75_s5',
	# 't_wald_ihw_corr_auto_w75_s20',
	# 't_wald_ihw_corr_auto_w100_s5',
	# 't_wald_ihw_corr_auto_w100_s20',
	# 't_wald_ihw_corr_para_w50_s5',
	# 't_wald_ihw_corr_para_w50_s20',
	# 't_wald_ihw_corr_para_w75_s5',
	# 't_wald_ihw_corr_para_w75_s20',
	# 't_wald_ihw_corr_para_w100_s5',
	# 't_wald_ihw_corr_para_w100_s20',
	# 't_wald_ihw_nocorr_auto_w50_s5',
	# 't_wald_ihw_nocorr_auto_w50_s20',
	# 't_wald_ihw_nocorr_auto_w75_s5',
	# 't_wald_ihw_nocorr_auto_w75_s20',
	# 't_wald_ihw_nocorr_auto_w100_s5',
	# 't_wald_ihw_nocorr_auto_w100_s20',
	# 't_wald_ihw_nocorr_para_w50_s5',
	# 't_wald_ihw_nocorr_para_w50_s20',
	# 't_wald_ihw_nocorr_para_w75_s5',
	# 't_wald_ihw_nocorr_para_w75_s20',
	# 't_wald_ihw_nocorr_para_w100_s5',
	# 't_wald_ihw_nocorr_para_w100_s20'
	# );


	# Search for motif PWMs
	foreach $type (@types)
	{
		run("Job", "~/scripts/qsub.sh job.pl $type")
		# run("Job", "~/scripts/qsub.sh job_grep_only.pl $type")
		# run("Job", "~/scripts/qsub.sh job_no_grep.pl $type")
	}
}

waitforjobs();

run("Plot logos of all the input PWMs", "logos.pl");

done();
