#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$scoretype = 'SUM_IS';
$threshold = 7.1;

# eCLIP parameters
$type = 'eclip_encode_12';

our $usage = "$0 [interaction source: biogrid/biogrid_direct/hippie/huri]\n\nExample: $0 huri";
($source) = args(1);
# args(0);

$new_table = "rbpome_$source";

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

# Get eCLIP proteins
startme("Getting eCLIP proteins from table 'clip_raw_gene'");
%eclip = ();
$query = Query("SELECT DISTINCT symbol FROM clip_raw_gene WHERE type='$type'");
while (($symbol) = Fetch($query))
{
	$eclip{$symbol} = 1;
	stepme(100);
}
stopme();




# Get positive controls (same as for rbpome)
state("Getting positive controls from 'rbpome' (>=2 direct studies in BioGRID)");
# state("Setting positive controls (manually)...");x
# warn("Warning: This uses manually set positive controls (6 of them)");
%positive = ();
# $positive{"FMR1|FXR2"} = 1;
# $positive{"FXR1|FXR2"} = 1;
# $positive{"PTBP1|PCBP1"} = 1;
# $positive{"PTBP1|QKI"} = 1;
# $positive{"SFPQ|NONO"} = 1;
# $positive{"U2AF1|U2AF2"} = 1;
# $query = Query("SELECT DISTINCT symbol1, symbol2 FROM rbpome WHERE eclip1=1 AND eclip2=1 AND bg_alldirect=1 AND homodimer=0 AND scoretype='$scoretype' AND threshold=$threshold AND species1='human' AND species2='human'");
# Require at least two direct studies in BioGRID for positive controls
$query = Query("SELECT DISTINCT r.symbol1, r.symbol2 FROM rbpome r, (SELECT gene1, gene2, COUNT(DISTINCT pmid) AS pmids FROM biogrid WHERE direct=1 GROUP BY gene1, gene2) g WHERE r.scoretype='$scoretype' AND r.threshold=$threshold AND r.symbol1!=r.symbol2 AND r.bg_alldirect=1 AND r.homodimer=0 AND r.eclip1=1 AND r.eclip2=1 AND r.hc=1 AND g.gene1=r.symbol1 AND g.gene2=r.symbol2 AND g.pmids>=2 ORDER BY symbol1, symbol2");
while (($symbol1, $symbol2) = Fetch($query))
{
	state(" >> positive >> $symbol1 >> $symbol2", 1);
	$positive{"$symbol1|$symbol2"} = 1;
	# stepme(100);
}
# stopme();




state("Dropping & recreating table '$new_table' using '$source' interactions based on table 'rbpome'...");
Query("DROP TABLE IF EXISTS $new_table");
Query("CREATE TABLE $new_table LIKE rbpome");

# id
# symbol1
# symbol2
# eclip1
# eclip2
# scoretype
# threshold
# hc
# homodimer
# bg_alldirect
# avg_is

if ($source eq 'biogrid')
{
	# query "SELECT * FROM biogrid LIMIT 5" -h
	# id	species	gene1	gene2	ncbigene1	ncbigene2	orderedlocus1	orderedlocus2	system	type	direct	pmid
	# 1	human	MAP2K4	FLNC	2318	112315	[NULL]	[NULL]	Two-hybrid	physical	1	9006895
	# 2	human	MYPN	ACTN2	88	124185	[NULL]	[NULL]	Two-hybrid	physical	1	11309420
	# 3	human	ACVR1	FNTA	2339	106605	[NULL]	[NULL]	Two-hybrid	physical	1	8599089

	$mainquery = Query("SELECT gene1, gene2 FROM biogrid GROUP BY gene1, gene2 ORDER BY gene1, gene2");
}
elsif ($source eq 'biogrid_direct')
{
	# query "SELECT * FROM biogrid WHERE direct=1 LIMIT 5" -h
	# id	species	gene1	gene2	ncbigene1	ncbigene2	orderedlocus1	orderedlocus2	system	type	direct	pmid
	# 1	human	MAP2K4	FLNC	2318	112315	[NULL]	[NULL]	Two-hybrid	physical	1	9006895
	# 2	human	MYPN	ACTN2	88	124185	[NULL]	[NULL]	Two-hybrid	physical	1	11309420
	# 3	human	ACVR1	FNTA	2339	106605	[NULL]	[NULL]	Two-hybrid	physical	1	8599089

	$mainquery = Query("SELECT gene1, gene2 FROM biogrid WHERE direct=1 GROUP BY gene1, gene2 ORDER BY gene1, gene2");
}
elsif ($source eq 'hippie')
{
	# query "SELECT * FROM hippie LIMIT 5" -h
	# id	acc1	acc2	name1	name2	symbol1	symbol2	species	ncbigene1	ncbigene2	systems	sources	pmids
	# 1	P00352	P00352	AL1A1_HUMAN	AL1A1_HUMAN	ALDH1A1	ALDH1A1	HUMAN	216	216	in vivo|Two-hybrid	BioGRID|HPRD|I2D|IntAct|MINT|Rual05	12081471|16189514|25416956
	# 2	Q13683	P02708	ITA7_HUMAN	ACHA_HUMAN	ITGA7	CHRNA1	HUMAN	3679	1134	Affinity Capture-Western|affinity chromatography technology|in vivo	BioGRID|HPRD|I2D	10910772
	# 3	Q9ULJ8	P63261	NEB1_HUMAN	ACTG_HUMAN	PPP1R9A	ACTG1	HUMAN	55607	71	in vitro|in vivo	HPRD	9362513|12052877
	
	$mainquery = Query("SELECT symbol1, symbol2 FROM hippie WHERE species='human' GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
}
elsif ($source eq 'huri')
{
	# query "SELECT * FROM huri WHERE type='huri' LIMIT 5" -h
	# id	ensg1	ensg2	symbol1	symbol2	type	species
	# 1	ENSG00000000005	ENSG00000099968	TNMD	BCL2L13	huri	human
	# 2	ENSG00000000005	ENSG00000104765	TNMD	BNIP3L	huri	human
	# 3	ENSG00000000005	ENSG00000105383	TNMD	CD33	huri	human
		
	$mainquery = Query("SELECT symbol1, symbol2 FROM huri WHERE type='huri' AND species='human' GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2");
}
else
{
	die("Error: Unhandled source '$source'");
}
startme("Filling table '$new_table' with interactions from source '$source'", Numrows($mainquery));
starttime();
while (($symbol1, $symbol2) = Fetch($mainquery))
{
	if ($symbol1 eq $symbol2)
	{
		addme("homomer skipped for symbol1|symbol2", "$symbol1|$symbol2");
		next;
	}
	
	# id
	# symbol1
	# symbol2
	# eclip1
	# eclip2
	# scoretype
	# threshold
	# hc
	# homodimer
	# bg_alldirect
	# avg_is
	
	# Get eCLIP proteins
	$eclip1 = 0;
	$eclip2 = 0;
	if (exists($eclip{$symbol1}))
	{
		$eclip1 = 1;
		addme("eclip protein inserted", $symbol1);
	}
	if (exists($eclip{$symbol2}))
	{
		$eclip2 = 1;
		addme("eclip protein inserted", $symbol2);
	}
	
	if (($eclip1 == 1) and ($eclip2 == 1))
	{
		addme("eclip-eclip pair inserted", "$symbol1|$symbol2");
	}
	
	# For the positive control: use the same set of proteins as for 'rbpome' (bg_alldirect=1 proteins in 'rbpome')
	$bg_alldirect = 0;
	if (exists($positive{"$symbol1|$symbol2"}))
	{
		$bg_alldirect = 1;
		addme("positive control inserted", "$symbol1|$symbol2");
	}
	
	$avg_is = $threshold;
	$hc = 1;
	$homodimer = 0;
	
	Query("INSERT INTO $new_table SET symbol1='$symbol1', symbol2='$symbol2', eclip1='$eclip1', eclip2='$eclip2', scoretype='$scoretype', threshold='$threshold', hc='$hc', homodimer='$homodimer', bg_alldirect='$bg_alldirect', avg_is='$avg_is'");
	
	addme("pair inserted", "$symbol1|$symbol2");
	
	stepme(1000);
}
stopme();
stoptime();

showmeall(1);

done();
