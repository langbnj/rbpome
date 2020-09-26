#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$max_connections = 128;		# Maximum number of jobs to submit simultaneously (the MySQL server is limited to ~150)

our $usage = "$0 [table] [eclip type] [screen score type] [mindist threshold for cobinding] [motif search method] [motif 5' extension 0/1] [number of resamples] [control to compare to] [-qvalue]\n\n -qvalue: Use motif hits with significant q-values, rather than p-values.\n\nExample: $0 rbpome eclip_encode_12 SUM_IS 54 fimo 1 100\n";
($table, $type, $scoretype, $mindist_threshold, $motif_method, $motif_extend5, $resamples) = args(7);

$tmpextend5 = '';
$tmpextend5 = '_extend5' if ($motif_extend5 == 1);

$tmpqvalue = '';
$tmpqvalue = '-qvalue' if (switch('qvalue'));

$tmpsig = 'psig=1';
$tmpsig = 'qsig=1' if (switch('qvalue'));


$infile = "output-motif_analysis-$table-$type-$scoretype-mindist_threshold$mindist_threshold-$motif_method-extend5_$motif_extend5$tmpqvalue.txt";
open(IN, $infile) or die("Error: Couldn't open '$infile'");



# start

startme("Reading RBP pairs with motifs from '$infile'");
@rbppairs = ();
while (<IN>)
{
	chomp;
	
	# symbol1	symbol2	bindtype	mindist	motifcount1	motifcount2	motifs1	motifs2
	# EWSR1	HNRNPK	close	5	0	9		attract|HNRNPK_2,attract|HNRNPK_5,attract|HNRNPK_7,attract|HNRNPK_8
	# EWSR1	HNRNPK	close	10	5	7	dominguez|EWSR1_1	attract|HNRNPK_2,attract|HNRNPK_4,attract|HNRNPK_6,attract|HNRNPK_9,cisbp|HNRNPK_1,rbpmap|HNRNPK_1,rnacompete|HNRNPK_1
	# EWSR1	HNRNPK	close	17	5	6	dominguez|EWSR1_1	attract|HNRNPK_1,attract|HNRNPK_13,attract|HNRNPK_2,attract|HNRNPK_21,attract|HNRNPK_5,attract|HNRNPK_8
	# EWSR1	HNRNPK	close	44	0	8		attract|HNRNPK_19,attract|HNRNPK_24,attract|HNRNPK_7
	# EWSR1	HNRNPK	close	7	0	0
	# EWSR1	HNRNPK	close	2	0	6		attract|HNRNPK_2,attract|HNRNPK_5,attract|HNRNPK_7,attract|HNRNPK_8
	# EWSR1	HNRNPK	close	19	0	1		dominguez|HNRNPK_1
	# EWSR1	HNRNPK	close	12	0	3		attract|HNRNPK_19,attract|HNRNPK_2
	# EWSR1	HNRNPK	close	40	0	3		attract|HNRNPK_19,attract|HNRNPK_2

	@a = split(/\t/);
	$symbol1 = $a[0];
	$symbol2 = $a[1];
	
	push(@rbppairs, "$symbol1|$symbol2");
}
@rbppairs = unique(@rbppairs);


cd('tmp');
startme("Submitting jobs to calculate close/control ratios");
$mine = mynodes() + queuednodes();
foreach $pair (@rbppairs)
{
	($symbol1, $symbol2) = split(/\|/, $pair);
	
	run("qsub", "~/scripts/qsub.sh ../motif_ratios_job.pl $table $type $scoretype $mindist_threshold $motif_method $motif_extend5 $resamples $symbol1 $symbol2");
	$mine++;
	
	addme("job submitted for RBP pair", "$symbol1|$symbol2");
	
	while ($mine >= $max_connections)
	{
		$mine = mynodes() + queuednodes();

		# $mine = 0 if (queuednodes() > 0);

		sleep(5);
	}
}
cd('..');

waitforjobs();


showmeall(1);

done();
