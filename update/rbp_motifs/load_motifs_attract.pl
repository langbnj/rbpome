#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# CREATE TABLE `rbp_motif_list` (
#   `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
#   `symbol` varchar(25) DEFAULT NULL,
#   `acc` varchar(25) DEFAULT NULL,
#   `species` varchar(25) DEFAULT NULL,
#   `source` varchar(25) DEFAULT NULL,
#   `motif` varchar(25) DEFAULT NULL,
#   `first` tinyint(1) DEFAULT NULL,
#   `logonum` tinyint(1) DEFAULT NULL,
#   `score` float DEFAULT NULL,
#   PRIMARY KEY (`id`),
#   KEY `Symbol` (`symbol`),
#   KEY `Acc` (`acc`),
#   KEY `Species` (`species`),
#   KEY `Source` (`source`),
#   KEY `First` (`first`),
#   KEY `Logonum` (`logonum`)
# ) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='RBP binding motifs';
#
#
# CREATE TABLE `rbp_motif_instances` (
#   `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
#   `symbol` varchar(25) DEFAULT NULL,
#   `acc` varchar(25) DEFAULT NULL,
#   `species` varchar(25) DEFAULT NULL,
#   `celltype` varchar(50) DEFAULT NULL,
#   `rep` tinyint(1) unsigned DEFAULT NULL,
#   `type` varchar(25) DEFAULT NULL,
#   `map` varchar(25) DEFAULT NULL,
#   `chr` varchar(50) DEFAULT NULL,
#   `start` int(11) unsigned DEFAULT NULL,
#   `stop` int(11) unsigned DEFAULT NULL,
#   `strand` char(1) DEFAULT NULL,
#   `ensgv` varchar(25) DEFAULT NULL,
#   `ensg` varchar(25) DEFAULT NULL,
#   `source` varchar(25) DEFAULT NULL,
#   `motif` varchar(25) DEFAULT NULL,
#   `firstonly` tinyint(1) unsigned DEFAULT NULL,
#   `mode` varchar(25) DEFAULT NULL,
#   `hits` int(11) unsigned DEFAULT NULL,
#   `pvalue` double DEFAULT NULL,
#   `qvalue` double DEFAULT NULL,
#   PRIMARY KEY (`id`),
#   KEY `Symbol` (`symbol`),
#   KEY `Acc` (`acc`),
#   KEY `Species` (`species`),
#   KEY `Celltype` (`celltype`),
#   KEY `Rep` (`rep`),
#   KEY `Type` (`type`),
#   KEY `Map` (`map`),
#   KEY `Chr` (`chr`),
#   KEY `Strand` (`strand`),
#   KEY `Source` (`source`),
#   KEY `Motif` (`motif`),
#   KEY `Hits` (`hits`),
#   KEY `Ensgv` (`ensgv`),
#   KEY `Ensg` (`ensg`),
#   KEY `Firstonly` (`firstonly`),
#   KEY `Mode` (`mode`)
# ) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='RBP motif instances within eCLIP peaks';



$table = 'rbp_motif_list';
$source = 'attract';
$species = 'human';

our $usage = "$0 [-firstonly]\n\n -firstonly: Only use the first (top) motif per RBP (skip all lower-ranked significant motifs)\n\nExample: $0 -firstonly";
# ($var) = args(1);
args(0);

$infile = "input/attract/attract_db.txt";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

Query("DELETE FROM `$table` WHERE source='$source'");
state("Cleared '$source' rows from table '$table'");



# start

startme("Reading '$infile' and inserting motifs into table '$table'");
starttime();
<IN>;	# Skip header
$inserted = 0;
while (<IN>)
{
	chomp;

	# Gene_name	Gene_id	Mutated	Organism	Motif	Len	Experiment_description	Database	Pubmed	Experiment_description	Family	Matrix_id	Score
	# 3IVK	3IVK	no	Mus_musculus	GAAACA	6	X-RAY DIFFRACTION	PDB	19965478	X-RAY DIFFRACTION	N/A	519	1.000000**
	# 3IVK	3IVK	no	Mus_musculus	UGGG	4	X-RAY DIFFRACTION	PDB	19965478	X-RAY DIFFRACTION	N/A	574	1.000000**
	# 4KZD	4KZD	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	464	1.000000**
	# 4KZE	4KZE	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	437	1.000000**
	# 4Q9Q	4Q9Q	no	Mus_musculus	GAAAC	5	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	423	1.000000**
	# 4Q9R	4Q9R	no	Mus_musculus	CGAAAC	6	X-RAY DIFFRACTION	PDB	24952597	X-RAY DIFFRACTION	N/A	433	1.000000**
	# A1CF	ENSG00000148584	no	Homo_sapiens	UGAUCAGUAUA	11	UV cross-linking	R	10669759	UV cross-linking	RRM	110	1.000000**
	# A1CF	ENSG00000148584	no	Homo_sapiens	AUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.126952
	# A1CF	ENSG00000148584	no	Homo_sapiens	UUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.126411
	# A1CF	ENSG00000148584	no	Homo_sapiens	AUAAUUG	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.189114**
	# A1CF	ENSG00000148584	no	Homo_sapiens	UUAAUUG	7	RNAcompete	C	23846655	RNAcompete	RRM	M001_0.6	0.188308
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	AUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.029649
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	UUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.04527
	# A1CF	ENSGALG00000003765	no	Gallus_gallus	GUAAUUA	7	RNAcompete	C	23846655	RNAcompete	RRM	M244_0.6	0.075703
	
	# Gene_name	( no need to explain :-) right?)
	# Gene_id	( no need to explain :-) right?)
	# Mutated	(if the target gene is mutated)
	# Organism	( no need to explain :-) right?)
	# Motif	( no need to explain :-) right?)
	# Len	(lenght of the motif)
	# Experiment_description(when available)
	# Database (Database from where the motifs were extracted PDB: Protein data bank, C: Cisbp-RNA, R:RBPDB, S: Spliceaid-F, AEDB:ASD)
	# Pubmed (pubmed ID)
	# Experiment (type of experiment; short description)
	# Family (domain)
	# Matrix_id (linked to the file PWM.txt)
	# Score (Qscore refer to the paper)
	#
	# The field Matrix_id refers to the pwm id that you can find in the pwm.txt file.
	# The position weight matrices are annotated in fasta format.
	
	@a = split(/\t/);
	
	$symbol = $a[0];
	# $ensg = $a[1];
	$mutated = $a[2];
	$fullspecies = $a[3];
	$motif = $a[4];
	# $motiflen = $a[5];
	# $exptype = $a[6];
	# $motifid = $a[11];
	$score = $a[12];
	
	# Remove whitespace from symbol (necessary for NOVA1)
	$symbol =~ s/^\s+//;
	$symbol =~ s/\s+$//;
	
	# From the paper: "Score" is a quality score for SELEX experiments etc. (i.e. where there is data on multiple motifs). Where there is only one motif, the score is 1. Also, the top score seems to be highlighted with two asterisks: **.
	# I'll keep only lines with the asterisks (since I only care about the mapping).
	
	# Skip mutated motifs
	next if ($mutated ne 'no');
	
	# Process score
	$first = 1;
	if ($score !~ /\*\*$/)
	{
		# # Skip non-"top" motifs (they'd be represented as one PWM anyway)
		# next;
		$first = 0;
	}
	# Remove asterisks for the "top" motifs
	$score =~ s/\*\*$//;
	
	# Skip non-human motifs
	if ($fullspecies ne 'Homo_sapiens')
	{
		addme("skipped non-human motif for fullspecies", $fullspecies);
		next;
	}
	
	addme("total symbols", $symbol);
	
	# Get UniProt accession
	# $query = Query("SELECT DISTINCT name FROM uniprot WHERE gene='$symbol' AND species='$species'");
	# ($name) = FetchOne($query);
	# $query = Query("SELECT DISTINCT acc FROM uniacc WHERE name='$name' AND species='$species'");
	# ($acc) = FetchOne($query);
	$query = Query("SELECT DISTINCT acc FROM clip_gene WHERE symbol='$symbol' AND species='$species'");
	$acc = '';
	if (Numrows($query) > 0)
	{
		($acc) = FetchOne($query);
	}
	
	addme("total symbol|motifs", "$symbol|$motif");
	addme("total species", $fullspecies);
	
	$s = "INSERT INTO `$table` SET symbol='$symbol', acc='$acc', species='$species', source='$source', motif='$motif', first='$first', logonum='', score='$score'";
	$s =~ s/=''/=NULL/g;
	Query($s);
	$inserted++;

	stepme(100);
}
stopme();
stoptime();

showmeall(1);

state("Rows inserted: ".commify($inserted));

Optimize($table);

done();
