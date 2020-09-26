#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run

if (!switch('finish'))
{
	Clear('rbpome');


	# Final scores
	# run("Main", "main.pl SUM_IS input/_Final_Summary_Table.20200229.v3.txt 7.1 -noalias");
	run("Main", "main.pl SUM_IS input/_Final_Summary_Table.20200229.v3.txt 7.1");


	# # Main scores (Sebastian)
	# run("Main", "main_old.pl AVG_IS input/RBPome_final_CF_0.txt 3.25");
	#
	# # Alternative scores (Jae-Seong)
	# run("Main", "alternative.pl 11_IS 1.7");
	# run("Main", "alternative.pl SUM_IS 7.1");
	#
	# # Alternative scores (Sebastian)
	# run("Main", "main_old.pl 10_F1 input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_F1_CF_0.96_noHDs_noAAs.txt 0.96");
	# run("Main", "main_old.pl 10_MCC input/sebastian_alternative_scores/recY2H_RBPome_10x_wavgIS2_log10_MCC_CF_0.67_noHDs_noAAs.txt 0.67");
	# run("Main", "main_old.pl 10_4_1_F1 input/sebastian_alternative_scores/recY2H_RBPome_10x_4x_1x_wavgIS_log10_F1_CF_0.97_noHDs_noAAs.txt 0.97");
	# run("Main", "main_old.pl 10_4_1_MCC input/sebastian_alternative_scores/recY2H_RBPome_10x_4x_1x_wavgIS_log10_MCC_CF_0.53_noHDs_noAAs.txt 0.53");
	# run("Main", "main_old.pl noAA_AVG_IS input/sebastian_alternative_scores/RBPome_final_CF_0_noAAs_noHDs.txt 3.25");

	# Add UniProt accessions
	run("Add UniProt accessions", "add_accs.pl -addmeth");



	# Filter out autoactivators as baits (where Protein A is an autoactivator)
	run("Remove auto-activators", "remove_autoactivators.pl");



	# Deduplicate the table (collapsing inverse pair duplicates by keeping the higher score, and the alphanumerically ordered configuration (symbol1 lt symbol2))
	# This causes these two to be identical:
	# SELECT scoretype, COUNT(*) FROM rbpome GROUP BY scoretype;
	# SELECT scoretype, COUNT(DISTINCT symbol1, symbol2) FROM rbpome GROUP BY scoretype;
	# run("Deduplicate", "deduplicate_backup_before_one_query_merging.pl");
	# run("Deduplicate", "deduplicate.pl -debug");
	run("Deduplicate", "deduplicate.pl");


	# Dump the table
	run("Dump it", qq(~/scripts/query.pl "SELECT * FROM rbpome" -h > output-rbpome-complete.txt));
	# run("Get header", "head -n1 input/_Final_Summary_Table.20200229.v3.txt > output-rbpome.txt");
	run("Get header", "cat input/header.txt > output-rbpome.txt");
	# run("Dump it", q(~/scripts/query.pl "SELECT symbol1, symbol2, times_detected, unique_ab_avg_is, avg_is, avg_rs, avg_ris, found_inverse, h47_avg_is, h47_sum_is, h47_found_in_both_orientations, h47_times_detected, h47_avg_rs, h47_avg_ris, bg_alldirect, bg_all, hippie, huri, '' AS hpa_nuclear, '' AS bioid, nanobret, nanobret_mbu_avg, nanobret_mbu_sd, species1, species2 FROM rbpome" | perl -ne 's/\[NULL\]/N\/A/g; s/$/\t\t\t\t\t\t\t/; print' >> output-rbpome.txt));
	# run("Dump it", q(~/scripts/query.pl "SELECT symbol1, symbol2, times_detected, unique_ab_avg_is, avg_is, avg_rs, avg_ris, found_inverse, h47_avg_is, h47_sum_is, h47_found_in_both_orientations, h47_times_detected, h47_avg_rs, h47_avg_ris, bg_alldirect, bg_all, hippie, huri, '' AS hpa_nuclear, '' AS bioid, nanobret, nanobret_mbu_avg, nanobret_mbu_sd, NULL, NULL, NULL, NULL, NULL, NULL, NULL, species1, species2 FROM rbpome" | perl -ne 's/\[NULL\]/N\/A/g; print' >> output-rbpome.txt));
	run("Dump it", q(~/scripts/query.pl "SELECT symbol1, symbol2, species1, species2, NULL, NULL, NULL, NULL, times_detected, unique_ab_avg_is, avg_is, avg_rs, avg_ris, found_inverse, h47_avg_is, h47_sum_is, h47_found_in_both_orientations, h47_times_detected, h47_avg_rs, h47_avg_ris, bg_alldirect, bg_all, hippie, huri, '' AS hpa_nuclear, '' AS bioid, nanobret, nanobret_mbu_avg, nanobret_mbu_sd FROM rbpome" | perl -ne 's/\[NULL\]/N\/A/g; s/$/\t\t\t\t\t\t\t/; print' >> output-rbpome.txt));



	# Got eclip_pairs.txt by hand from main.pl showme output (total eclip pairs (heteromers) (71))
	# cat eclip_pairs.txt | perl -ne 'require("a"); chomp; ($symbol1, $symbol2) = split(/\t/); $query = Query("SELECT id FROM rbp_motifs_extend5_fimo_eclip_encode WHERE symbol=\"$symbol1\""); $motif1=0; $motif1=1 if (Numrows($query)>0); $query = Query("SELECT id FROM rbp_motifs_extend5_fimo_eclip_encode WHERE symbol=\"$symbol2\""); $motif2=0; $motif2=1 if (Numrows($query)>0); print "$symbol1\t$symbol2\t$motif1\t$motif2\n";' > eclip_pairs_with_motif_information.txt
	# Number of eCLIP pairs where both have known motifs:
	# cat eclip_pairs.txt | perl -ne 'chomp; ($s1, $s2) = split(/\t/); $m1=0; $m1=1 if ($s1 =~ /^(AKAP1|CPEB4|CSTF2|EIF4G2|EWSR1|FMR1|FUBP3|FUS|FXR1|FXR2|G3BP1|HNRNPA1|HNRNPC|HNRNPK|HNRNPL|HNRNPM|HNRNPU|IGF2BP1|IGF2BP2|IGF2BP3|KHDRBS1|KHSRP|MATR3|NONO|PABPC4|PABPN1|PCBP1|PCBP2|PTBP1|PUM1|PUM2|QKI|RBFOX2|RBM22|RBM5|SFPQ|SRSF1|SRSF7|SRSF9|SSB|SUPV3L1|TAF15|TARDBP|TIA1|TIAL1|TRA2A|U2AF2|YBX3|ZRANB2)$/); $m2=0; $m2=1 if ($s2 =~ /^(AKAP1|CPEB4|CSTF2|EIF4G2|EWSR1|FMR1|FUBP3|FUS|FXR1|FXR2|G3BP1|HNRNPA1|HNRNPC|HNRNPK|HNRNPL|HNRNPM|HNRNPU|IGF2BP1|IGF2BP2|IGF2BP3|KHDRBS1|KHSRP|MATR3|NONO|PABPC4|PABPN1|PCBP1|PCBP2|PTBP1|PUM1|PUM2|QKI|RBFOX2|RBM22|RBM5|SFPQ|SRSF1|SRSF7|SRSF9|SSB|SUPV3L1|TAF15|TARDBP|TIA1|TIAL1|TRA2A|U2AF2|YBX3|ZRANB2)$/); print "$s1\t$s2\t$m1\t$m2\n"'|g "\t1\t1"|wcl
}

if (switch('finish'))
{
	# -finish

	# Make master table (previously in ~/pipeline/rbpome_master_table)
	# NOTE: This also updates the HuRI, HIPPIE, BioGRID, known, and known_protein_1 and ..._2 columns
	run("Main (without UPDATE-ing)", "master.pl rbpome");
	# run("Main (with UPDATE-ing of the HuRI, HIPPIE, BioGRID, known, and known_protein_1 and ..._2 columns in table 'rbpome')", "master.pl rbpome -update");

	run("Import as rbpome_final", "~/scripts/import.pl output-rbpome-master.txt rbpome_final -allindices -overwrite");

	run("Make 'alternative' rbpome table", "rbpome_huri_hippie_biogrid.pl huri");
	run("Make 'alternative' rbpome table", "rbpome_huri_hippie_biogrid.pl hippie");
	run("Make 'alternative' rbpome table", "rbpome_huri_hippie_biogrid.pl biogrid");
	run("Make 'alternative' rbpome table", "rbpome_huri_hippie_biogrid.pl biogrid_direct");
}
else
{
	die("At this point: Please run ~/pipeline/rbpome_analysis to be able to fill the master table's eCLIP columns (Jaccard index etc.)\n\nThen, rerun:\n\nrun.pl -finish");
}


done();
