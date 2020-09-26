#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


our $usage = "$0 [eCLIP type]\n\nExample: $0 s_lrt_bh_corr_auto_w50_s5";
($type) = args(1);


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

	@grepsources = ('attract', 'dominguez', 'custom');
	@fimosources = ('attract', 'cisbp', 'dominguez', 'rbpdb', 'rbpmap', 'rnacompete');
	
	# @fimomethods = ('fimo', 'fimobg', 'fimobgi');
	@fimomethods = ('fimo');


	# Grep for non-PWM motifs (not available for all sources)
	# Includes parsing into rbp_motifs
	$table = "rbp_motifs$tmpextend5_2\_grep_$type";
	# if (Exists($table))
	# {
	# 	# Not clearing existing tables, keeping them & not rerunning
	# 	state(" >> Table $table already exists, skipping");
	# 	next;
	# }
	# Clear($table);
	Query("DROP TABLE IF EXISTS `$table`");
	Query("CREATE TABLE IF NOT EXISTS `$table` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `symbol` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `celltype` varchar(50) DEFAULT NULL,
  `rep` tinyint(1) unsigned DEFAULT NULL,
  `method` varchar(25) DEFAULT NULL,
  `extend5` tinyint(1) unsigned DEFAULT NULL,
  `type` varchar(50) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int(11) unsigned DEFAULT NULL,
  `stop` int(11) unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `motif` varchar(25) DEFAULT NULL,
  `motifstart` int(11) unsigned DEFAULT NULL,
  `motifstop` int(11) unsigned DEFAULT NULL,
  `score` double DEFAULT NULL,
  `pvalue` double DEFAULT NULL,
  `qvalue` double DEFAULT NULL,
  `psig` tinyint(1) DEFAULT NULL,
  `qsig` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Symbol` (`symbol`),
  KEY `Acc` (`acc`),
  KEY `Species` (`species`),
  KEY `Celltype` (`celltype`),
  KEY `Chr` (`chr`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Source` (`source`),
  KEY `Type` (`type`),
  KEY `Motif` (`motif`),
  KEY `Method` (`method`),
  KEY `Rep` (`rep`),
  KEY `Extend5` (`extend5`),
  KEY `Psig` (`psig`),
  KEY `Qsig` (`qsig`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='RBP motif instances within eCLIP peaks';
");
	foreach $source (@grepsources)
	{
		# Find motifs in eCLIP peak regions using MEME FIMO (only for RBPs with $source motifs)
		run("Get peak sequences", "get_peakseqs.pl grep $type $source$tmpextend5");

		# All motifs (stored with a significant p-value)
		run("Grep for non-PWM motifs", "grep.pl $type $source$tmpextend5");

		# Top motifs only (stored with a significant q-value)
		run("Grep for non-PWM motifs", "grep.pl $type $source -firstonly$tmpextend5");
	}


	# Search for motif PWMs
	$table = "rbp_motifs$tmpextend5_2\_fimo_$type";
	# if (Exists($table))
	# {
	# 	# Not clearing existing tables, keeping them & not rerunning
	# 	state(" >> Table $table already exists, skipping");
	# 	next;
	# }



	# Clear tables & DISABLE KEYS for LOADing
	foreach $method (@fimomethods)
	{
		$table = "rbp_motifs$tmpextend5_2\_$method\_$type";
		# if (Exists($table))
		# {
		# 	# Not clearing existing tables, keeping them & not rerunning
		# 	state(" >> Table $table already exists, skipping");
		# 	next;
		# }
		# Clear($table);
		Query("DROP TABLE IF EXISTS `$table`");
		Query("CREATE TABLE IF NOT EXISTS `$table` (
			  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
			  `symbol` varchar(25) DEFAULT NULL,
			  `acc` varchar(25) DEFAULT NULL,
			  `species` varchar(25) DEFAULT NULL,
			  `celltype` varchar(50) DEFAULT NULL,
			  `rep` tinyint(1) unsigned DEFAULT NULL,
			  `method` varchar(25) DEFAULT NULL,
			  `extend5` tinyint(1) unsigned DEFAULT NULL,
			  `type` varchar(50) DEFAULT NULL,
			  `chr` varchar(50) DEFAULT NULL,
			  `start` int(11) unsigned DEFAULT NULL,
			  `stop` int(11) unsigned DEFAULT NULL,
			  `strand` char(1) DEFAULT NULL,
			  `source` varchar(25) DEFAULT NULL,
			  `motif` varchar(25) DEFAULT NULL,
			  `motifstart` int(11) unsigned DEFAULT NULL,
			  `motifstop` int(11) unsigned DEFAULT NULL,
			  `score` double DEFAULT NULL,
			  `pvalue` double DEFAULT NULL,
			  `qvalue` double DEFAULT NULL,
			  `psig` tinyint(1) DEFAULT NULL,
			  `qsig` tinyint(1) DEFAULT NULL,
			  PRIMARY KEY (`id`),
			  KEY `Symbol` (`symbol`),
			  KEY `Acc` (`acc`),
			  KEY `Species` (`species`),
			  KEY `Celltype` (`celltype`),
			  KEY `Chr` (`chr`),
			  KEY `Start` (`start`),
			  KEY `Stop` (`stop`),
			  KEY `Strand` (`strand`),
			  KEY `Source` (`source`),
			  KEY `Type` (`type`),
			  KEY `Motif` (`motif`),
			  KEY `Method` (`method`),
			  KEY `Rep` (`rep`),
			  KEY `Extend5` (`extend5`),
			  KEY `Psig` (`psig`),
			  KEY `Qsig` (`qsig`)
			) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='RBP motif instances within eCLIP peaks';
			");
		Query("ALTER TABLE $table DISABLE KEYS");
	}


	# Run FIMO & Parse
	foreach $source (@fimosources)
	{
		# Find motifs in eCLIP peak regions using MEME FIMO (only for RBPs with $source motifs)
		run("Get peak sequences", "get_peakseqs.pl fimo $type $source$tmpextend5");

		# Run FIMO (p-value threshold)
		run("Run MEME FIMO", "fimo.pl $type $source 0.001$tmpextend5");

		# Run FIMO (q-value threshold)
		run("Run MEME FIMO", "fimo.pl $type $source 0.05 -qvalue$tmpextend5");

		# Parse
		foreach $method (@fimomethods)
		{
			$table = "rbp_motifs$tmpextend5_2\_$method\_$type";

			run("Parse MEME FIMO output", "parse_fimo.pl $method $type $source 0.001$tmpextend5");
			run("Parse MEME FIMO output", "parse_fimo.pl $method $type $source 0.05 -qvalue$tmpextend5");
		}
	}
	
	# Re-enable keys & optimize
	foreach $method ('grep', @fimomethods)
	{
		$table = "rbp_motifs$tmpextend5_2\_$method\_$type";

		state("Re-enabling keys on table '$table'...", 1);
		starttime();
		Query("ALTER TABLE $table ENABLE KEYS");
		stoptime();

		starttime();
		Optimize($table);
		stoptime();
	}
}


done();
