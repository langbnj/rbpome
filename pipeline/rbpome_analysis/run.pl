#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
args(0);

$max_connections = 128;		# Maximum number of jobs to submit simultaneously (the MySQL server is limited to ~150)

# $mindist_threshold = 50;
$mindist_threshold_strict = 54;	# "strict": This is the exact peak of the <=1000 positives curve in Figure 4e
$mindist_threshold_loose = 117;	# "looser": This is 95% of the peak of the <=1000 positives curve in Figure 4e

$num_resamples = 100;


$motif_method = 'fimo';
$motif_extend5 = 1;


# run

if (!switch('local'))
{
	$mine = mynodes() + queuednodes();
}

# @types = ('eclip_encode', 'eclip_encode_12');
# @types = ('eclip_encode');
# @types = ('t_lrt_ihw_nocorr_auto_w100_s5');
# @types = ('eclip_encode', 'eclip_tom', 'eclip_encode_12', 'eclip_par_ohler', 'eclip_tom_local', 'eclip_tom_waldtest', 'eclip_tom_waldtest_local');
# @types = ('eclip_encode', 'eclip_tom', 'eclip_encode_12', 'eclip_par_ohler');
# @types = ('clip_postar');
@types = ('eclip_encode_12');

@tables = ('rbpome');
# @tables = ('rbpome_huri');
# @tables = ('rbpome_hippie', 'rbpome_biogrid', 'rbpome_biogrid_direct');
# @tables = ('rbpome', 'rbpome_huri', 'rbpome_hippie', 'rbpome_biogrid', 'rbpome_biogrid_direct');






foreach $table (@tables)
{
	# Test if table exists
	# $query = Query("SELECT id FROM $table LIMIT 1");
	# ($tmpid) = FetchOne($query);
	die("Error: Table '$table' doesn't exist") if (!Exists($table));





	# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table");
	#DEBUG: get SUM_IS only (we've decided on this one)
	$scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
	#END DEBUG
	@scoretypes = ();
	while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
	{
		push(@scoretypes, "$scoretype|$minscore_adapted");
	}
	@scoretypes = unique(@scoretypes);




	cd("tmp");
	foreach $type (@types)
	{
		# foreach $scoretype_tmp ('avg_is|3.25', 'sum_is|7.1', '11_is|1.7', 'f1|0.96', 'mcc|0.67')
		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table ORDER BY scoretype");
		#DEBUG: get SUM_IS only (we've decided on this one)
		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
		#END DEBUG
		# while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
		# {
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			# foreach $mode ('probability_mindist', 'probability', 'fraction', 'jaccard', 'mindist', 'mediandist', 'mindist_thresh', 'mindist_thresh_frac')
			# foreach $mode ('jaccard', 'probability', 'mindist_thresh_frac', 'probability_mindist')
			# foreach $mode ('jaccard', 'probability', 'mindist_thresh_frac', 'probability_mindist', 'freq', 'probability_freqscaled', 'probability_mindist_freqscaled')
			# foreach $mode ('oddsratio', 'freq', 'probability_mindist', 'probability_freqscaled', 'probability_mindist_freqscaled')
			# foreach $mode ('jaccard', 'oddsratio')
			# foreach $mode ('oddsratio')
			# foreach $mode ('jaccard')
			# foreach $mode ('fraction', 'freq', 'jaccard', 'mediandist', 'mindist', 'mindist_thresh', 'mindist_thresh_frac', 'oddsratio', 'probability', 'probability_freqscaled', 'probability_mindist', 'probability_mindist_freqscaled')
			# foreach $mode ('fraction', 'freq', 'intersection', 'jaccard', 'mediandist', 'mindist', 'mindist_thresh', 'mindist_thresh_frac', 'oddsratio', 'probability', 'probability_freqscaled', 'probability_mindist', 'probability_mindist_freqscaled')
			foreach $mode ('fraction', 'freq', 'intersection', 'jaccard', 'mediandist', 'mindist', 'mindist_thresh', 'mindist_thresh_frac', 'oddsratio', 'probability', 'probability_mindist')
			{
				# foreach $hc (0, 1)
				foreach $hc (0)
				{
					# foreach $minlog2fold (0, 3, 5)
					foreach $minlog2fold (0)
					{
						# foreach $minscore (0, 1, 2, 3.25, 9, 20)
						# foreach $minscore (0, $minscore_adapted)
						foreach $minscore ($minscore_adapted)
						{
							# foreach $mindist_threshold (10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000)
							# foreach $mindist_threshold (mindist_threshold_loose, 400)
							# foreach $mindist_threshold (mindist_threshold_loose)
							foreach $mindist_threshold ($mindist_threshold_strict)
							{
								foreach $resamples (100)
								{
									# $pdffile = "../output/output-pdf-$table-$type-$scoretype-$mode-hc$hc-minlog2fold$minlog2fold-minscore$minscore-mindist_threshold$mindist_threshold-resamples$resamples.pdf";
									# if (-e $pdffile)
									# {
									# 	print " >> Output file '$pdffile' already exists, skipping!\n";
									# 	next;
									# }

									if (switch('local'))
									{
										# Run locally
										run("Main", "../main.pl $table $type $scoretype $mode $hc $minlog2fold $minscore $mindist_threshold $resamples");
									}
									else
									{
										# Run on the cluster
										# run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype $mode $hc $minlog2fold $minscore $mindist_threshold $resamples");
										run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype $mode $hc $minlog2fold $minscore $mindist_threshold $resamples -overwrite");

										$mine++;

										while ($mine >= $max_connections)
										{
											$mine = mynodes() + queuednodes();

											# $mine = 0 if (queuednodes() > 0);

											sleep(5);
										}
									}
								}
							}
						}
					}
				}
			}
		}


		waitforjobs() if (!switch('local'));

		# cd("..");

		# # Get top 20 pairs for Sebastian to test (across jaccard and mindist_thresh)
		# run("Best pairs (screen_hits)", "bestpairs.pl $table $type jaccard mindist_thresh positives 20");
		# run("Best pairs (screen_hits)", "bestpairs.pl $table $type jaccard mindist_thresh screen_hits 20 -skip_positives");
		# run("Best pairs (background_resamples)", "bestpairs.pl $table $type jaccard mindist_thresh background_resamples 20 -skip_screen_hits -skip_positives");
		#
		# # Get worst 20 pairs for Sebastian to test (across jaccard and mindist_thresh)
		# run("Worst pairs (screen_hits)", "worstpairs.pl $table $type jaccard mindist_thresh screen_hits 20");
		# run("Worst pairs (screen_hits)", "worstpairs.pl $table $type jaccard mindist_thresh background_resamples 20");


		# cd("tmp");
	}
	cd("..");


	waitforjobs() if (!switch('local'));








	# Figure 2 (the old figure2.pl)
	foreach $type (@types)
	{
		# foreach $scoretype_tmp ('avg_is|3.25', 'sum_is|7.1', '11_is|1.7', 'f1|0.96', 'mcc|0.67')
		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table ORDER BY scoretype");
		#DEBUG: get SUM_IS only (we've decided on this one)
		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
		#END DEBUG
		# 10_4_1_F1	0.97
		# 10_4_1_MCC	0.53
		# 10_F1	0.96
		# 10_MCC	0.67
		# 11_IS	1.7
		# AVG_IS	3.25
		# noAA_AVG_IS	3.25
		# SUM_IS	7.1
		cd("tmp");
		# while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
		# {
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			# Figures 2A-D
			run("Main 2A", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2B", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# run("Main 2C", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			foreach $mindist_threshold (10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000)
			{
				run("Main 2C Preparation", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold $num_resamples");
				run("Main 2D Preparation", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold $num_resamples");
				# run("Main 2D Preparation", "~/scripts/qsub.sh ~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 0 $mindist_threshold 100");
			}
			# Figure 2E I
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E II scaled
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist_freqscaled 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E III scaled
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype freq 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E IV odds ratio / jaccard
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E V oddsratio / probability_mindist
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E VI probability_mindist / oddsratio
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E VII oddsratio / fraction
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype fraction 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E VIII oddsratio / intersection
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype intersection 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figure 2E IX oddsratio / mindist_thresh
			# run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main 2E", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			# Figures 2F-G
			run("Main 2F", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			foreach $minscore (0, 1.07, 2.39, 3.072, 4.38, 5.4, 7.1, 10)
			{
				run("Main 2G Preparation", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore $mindist_threshold_strict $num_resamples");
				run("Main 2G Preparation", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore $mindist_threshold_strict $num_resamples");
				run("Main 2G Preparation", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh 0 0 $minscore $mindist_threshold_strict $num_resamples");
				# run("Main 2G Preparation", "~/scripts/qsub.sh ~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore 400 100");
			}


			# Distance plot

			# Run locally since the random pairs make things complicated (there could be overwrites, and the random pairs get written to a reproducible list at e.g. output/output-rbpome-randompairs.txt by the first -randompairs run)
			# run("Distances", "../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot");
			run("Distances (without bimodal fit, faster)", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -nofit");
			waitforjobs();
			run("Distances (without bimodal fit, faster) - Random pairs", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -nofit -randompairs");
			waitforjobs();

			# run("Distances (random data)", "../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -random");
			run("Distances (random data) (without bimodal fit, faster)", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -random -nofit");
			waitforjobs();
			run("Distances (random data) (without bimodal fit, faster) - Random pairs", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -random -nofit -randompairs");
			waitforjobs();


			# Distance plot (signed, retaining +/- for closest_per_peak mindists)
			
			# Run locally since the random pairs make things complicated (there could be overwrites, and the random pairs get written to a reproducible list at e.g. output/output-rbpome-randompairs.txt by the first -randompairs run)
			# run("Distances", "../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot");
			run("Distances (without bimodal fit, faster)", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -distance_signed -nofit");
			waitforjobs();
			run("Distances (without bimodal fit, faster) - Random pairs", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -distance_signed -nofit -randompairs");
			waitforjobs();

			# run("Distances (random data)", "../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -random");
			run("Distances (random data) (without bimodal fit, faster)", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -distance_signed -random -nofit");
			waitforjobs();
			run("Distances (random data) (without bimodal fit, faster) - Random pairs", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -distance_plot -distance_signed -random -nofit -randompairs");
			waitforjobs();
		}
		cd("..");
		waitforjobs();




		# Draw e.g. "Figure 2 - jaccard - eclip_encode_12 - SUM_IS.pdf"
		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
		cd("tmp");
		# while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
		# {
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			# Plot all the "modes"
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype fraction 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype freq 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype intersection 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mediandist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype oddsratio 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_freqscaled 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
			run("Main", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist_freqscaled 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples");
		}
		cd("..");

		waitforjobs();





		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
		# while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
		# {
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			# Draw prepared figures
			# Figures 2A-D
			cd("tmp");

			run("Main 2A", "~/scripts/qsub.sh ../main.pl $table $type $scoretype jaccard 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2A");
			run("Main 2B", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2B");
			run("Main 2C_old", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2C_old");
			run("Main 2C", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2C");
			run("Main 2D_old", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2D_old");
			run("Main 2D", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2D");

			run("Main 2F", "~/scripts/qsub.sh ../main.pl $table $type $scoretype mindist_thresh_frac 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2F");
			run("Main 2G", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2G");

			cd("..");

			# Compare

			foreach $mindist_threshold (10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000)
			# foreach $mindist_threshold (400, 1000)
			# foreach $mindist_threshold (400)
			{
				run("Compare", "~/scripts/qsub.sh compare.pl $table $type $scoretype mindist_thresh_frac 0.75 jaccard 0.35 0 0 $minscore_adapted $mindist_threshold $num_resamples");
				run("Compare", "~/scripts/qsub.sh compare.pl $table $type $scoretype probability_mindist 0.5 jaccard 0.35 0 0 $minscore_adapted $mindist_threshold $num_resamples");
				run("Compare", "~/scripts/qsub.sh compare.pl $table $type $scoretype probability_mindist 0.75 mindist_thresh_frac 0.75 0 0 $minscore_adapted $mindist_threshold $num_resamples");
				run("Compare", "~/scripts/qsub.sh compare.pl $table $type $scoretype mindist_thresh_frac 0.75 probability_mindist 0.75 0 0 $minscore_adapted $mindist_threshold $num_resamples");

				# Make motif plots
				cd("tmp");
				run("Make motif plots", "~/scripts/qsub.sh ../main.pl $table $type $scoretype probability_mindist 0 0 $minscore_adapted $mindist_threshold 100 -distance_plot -nofit -motif_plot");
				cd("..");
			}


			# Compare real vs. random
			foreach $mindist_threshold (10, 20, 40, $mindist_threshold_strict, $mindist_threshold_loose, 400, 1000, 10000)
			{
				# for I in 10 20 40 $mindist_threshold_strict $num_resamples 400 1000 10000; do q fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I; q fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I -randompairs; done
				# for I in 10 20 40 $mindist_threshold_strict $num_resamples 400 1000 10000; do le fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I; le fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I -randompairs; done
				# for I in 10 20 40 $mindist_threshold_strict $num_resamples 400 1000 10000; do lo fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I | tail; lo fit_vs_random.pl rbpome eclip_encode_12 SUM_IS $I -randompairs | tail; done
				run("Compare real distance data vs. random for real pairs", "~/scripts/qsub.sh fit_vs_random.pl $table $type $scoretype $mindist_threshold");
				run("Compare real distance data vs. random for random pairs", "~/scripts/qsub.sh fit_vs_random.pl $table $type $scoretype $mindist_threshold -randompairs");

				# Retaining sign (+/-) of the distances
				run("Compare real distance data vs. random for real pairs (retaining distance sign, +/-)", "~/scripts/qsub.sh fit_vs_random.pl $table $type $scoretype $mindist_threshold -distance_signed");
				run("Compare real distance data vs. random for random pairs (retaining distance sign, +/-)", "~/scripts/qsub.sh fit_vs_random.pl $table $type $scoretype $mindist_threshold -distance_signed -randompairs");
			}



			# Figure 2E oddsratio vs probability_mindist
			run("Compare 2E", "~/scripts/qsub.sh compare.pl $table $type $scoretype oddsratio 20 probability_mindist 0.4 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2E");

			# Figure 2E jaccard vs probability
			run("Compare 2E", "~/scripts/qsub.sh compare.pl $table $type $scoretype jaccard 0.3 probability 0.875 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2E");

			# Figure 2E jaccard vs probability_mindist (final)
			run("Compare 2E", "~/scripts/qsub.sh compare.pl $table $type $scoretype jaccard 0.3 probability_mindist 0.5 0 0 $minscore_adapted $mindist_threshold_strict $num_resamples -2E");
		}

		waitforjobs();



		# $scoretypequery = Query("SELECT DISTINCT scoretype, threshold FROM $table WHERE scoretype='SUM_IS'");
		# while (($scoretype, $minscore_adapted) = Fetch($scoretypequery))
		# {
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			# Make an overview plot combining the fit_vs_random plots
			run("Compare real distance data vs. random", "~/scripts/qsub.sh fit_vs_random_overview.pl $table $type $scoretype");
			run("Compare real distance data vs. random - Random pairs", "~/scripts/qsub.sh fit_vs_random_overview.pl $table $type $scoretype -randompairs");
			run("Consolidate real and random distance data into one table file", "~/scripts/qsub.sh fit_vs_random_consolidate.pl $table $type $scoretype");
		}

		waitforjobs();
		
		
		
		# Motif analysis (close binding)
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			run("Motif analysis (close binding)", "~/scripts/qsub.sh motif_analysis.pl $table $type $scoretype $mindist_threshold_strict $motif_method $motif_extend5");
		}
		waitforjobs();

		# Motif analysis (add non-close control, randomly sampled)
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			run("Motif analysis (distant binding, randomly sampled)", "~/scripts/qsub.sh motif_analysis_control.pl $table $type $scoretype $mindist_threshold_strict $motif_method $motif_extend5");
		}
		waitforjobs();

		# Motif analysis (add non-close control, randomly sampled, "solo"-bound genes only (where only one of the RBPs binds a gene))
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			run("Motif analysis (distant binding, randomly sampled)", "~/scripts/qsub.sh motif_analysis_control_solo.pl $table $type $scoretype $mindist_threshold_strict $motif_method $motif_extend5");
		}
		waitforjobs();

		# Motif ratio analysis (submits jobs, hence run locally)
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			run("Motif ratio analysis (motif pairs in close vs. distant binding)", "motif_ratios_submit.pl $table $type $scoretype $mindist_threshold_strict $motif_method $motif_extend5 $num_resamples");
		}
		waitforjobs();

		# Motif ratio analysis (combines & plots results)
		foreach $scoretype_tmp (@scoretypes)
		{
			($scoretype, $minscore_adapted) = split(/\|/, $scoretype_tmp);

			run("Motif ratio analysis (combine & plot results)", "~/scripts/qsub.sh motif_ratios_combine.pl $table $type $scoretype $mindist_threshold_strict $motif_method $motif_extend5");
		}
		waitforjobs();
	}












	run("Concatenate output files to make a combined p-value table & output-txt-$table.txt", "~/scripts/qsub.sh concat.pl $table");
	waitforjobs();
	run("Import output table", "~/scripts/qsub.sh ~/scripts/import.pl output-txt-$table.txt rbpome_analysis_$table -overwrite");
	waitforjobs();
	run("Import p-value table", "~/scripts/qsub.sh ~/scripts/import.pl output-pvalues-$table.txt rbpome_pvalues_$table -overwrite");
	waitforjobs();
	run("Import hits table", "~/scripts/qsub.sh ~/scripts/import.pl output-hits-$table.txt rbpome_hits_$table -overwrite");
	waitforjobs();


	# This part is a little broken if there's more than one type, the output files would need to be separated (or just run a single type above and uncomment)
	# run("Create output file without background_resamples", "cat output-txt-all.txt | grep -iP -v background_resample > output-txt-nobackground.txt");
	# # cat output-txt-all.txt | perl -ne 'chomp; @a = split(/\t/); print "$_\n" if ($a[0] eq "jaccard" and $a[9]>0.4)'>output-top-jaccard.txt
	# # cat output-txt-all.txt | perl -ne 'chomp; @a = split(/\t/); print "$_\n" if ($a[0] eq "mindist_thresh" and $a[9]>1000)'>output-top-mindist_thresh.txt
	# run("Sort output: Create header", "head -n1 output-txt-all.txt > output-txt-all-sorted.txt");
	# run("Sort output", "tail -n +2 output-txt-all.txt | sort -grk1,13 >> output-txt-all-sorted.txt");
	# run("Get Jaccard top 20", q(cat output-txt-all-sorted.txt | grep -iP "\tjaccard\t" | head -n 20 > output-jaccard-top20.txt));
	# run("Get Jaccard bottom 20", q(cat output-txt-all-sorted.txt | grep -iP "\tjaccard\t" | tail -n 20 > output-jaccard-bottom20.txt));
	# run("Get mindist_thresh top 20", q(cat output-txt-all-sorted.txt | grep -iP "\tmindist_thresh\t" | head -n 20 > output-mindist_thresh-top20.txt));
	# run("Get mindist_thresh bottom 20", q(cat output-txt-all-sorted.txt | grep -iP "\tmindist_thresh\t" | tail -n 20 > output-mindist_thresh-bottom20.txt));


	# # -- jaccard
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='jaccard' AND myset='positives' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC\" -h > output-pairs-jaccard-positives.txt");
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='jaccard' AND myset='screen_hits' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC\" -h > output-pairs-jaccard-screen_hits.txt");
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='jaccard' AND myset='background_resample' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC LIMIT 10000000\" -h > output-pairs-jaccard-background_resamples.txt");
	#
	# # -- mindist_thresh
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, mindist_threshold, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='mindist_thresh' AND myset='positives' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC\" -h > output-pairs-mindist_thresh-positives.txt");
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, mindist_threshold, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='mindist_thresh' AND myset='screen_hits' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC\" -h > output-pairs-mindist_thresh-screen_hits.txt");
	# run("Query", "~/scripts/query.pl \"SELECT DISTINCT rbp1, rbp2, scoretype, hc, minlog2fold, minscore, myset, type, mode, mindist_threshold, cobind FROM rbpome_analysis_$table WHERE scoretype='avg_is' AND mode='mindist_thresh' AND myset='background_resample' AND hc=0 AND minlog2fold=0 AND minscore=3.25 AND mindist_threshold=50 ORDER BY cobind DESC LIMIT 10000000\" -h > output-pairs-mindist_thresh-background_resamples.txt");
	
}


done();
