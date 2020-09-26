#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'hippie';

our $usage = "$0";
# ($var) = args(1);
args(0);

$infile = "input/hippie_current.txt";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

Clear($table);


# start

startme("Reading '$infile' and inserting into table '$table'");
starttime();
while (<IN>)
{
	chomp;
	
	@a = split(/\t/);
	
	# AL1A1_HUMAN	216	AL1A1_HUMAN	216	0.76	experiments:in vivo,Two-hybrid;pmids:12081471,16189514,25416956;sources:HPRD,BioGRID,IntAct,MINT,I2D,Rual05
	# ITA7_HUMAN	3679	ACHA_HUMAN	1134	0.73	experiments:in vivo,Affinity Capture-Western,affinity chromatography technology;pmids:10910772;sources:HPRD,BioGRID,I2D
	# NEB1_HUMAN	55607	ACTG_HUMAN	71	0.65	experiments:in vitro,in vivo;pmids:9362513,12052877;sources:HPRD
	# SRGN_HUMAN	5552	CD44_HUMAN	960	0.63	experiments:in vivo;pmids:9334256,16189514,16713569;sources:HPRD,I2D,Rual05,Lim06
	# GRB7_HUMAN	2886	ERBB2_HUMAN	2064	0.9	experiments:in vitro,in vivo,Reconstituted Complex,protein array,affinity chromatography technology,nuclear magnetic resonance,x-ray crystallography,competition binding,pull down;pmids:9079677,16273093,12975581,12061724,17875712;sources:HPRD,BioGRID,MINT,I2D,IntAct
	# PAK1_HUMAN	5058	ERBB2_HUMAN	2064	0.73	experiments:in vivo,Affinity Capture-Western,affinity chromatography technology;pmids:9774445;sources:HPRD,BioGRID,I2D,STRING

	$names1 = $a[0];
	$ncbigene1 = $a[1];
	$names2 = $a[2];
	$ncbigene2 = $a[3];
	
	$names1 = join('|', unique(split(/,/, $names1)));
	$names2 = join('|', unique(split(/,/, $names2)));
	
	$score = $a[4];
	
	# Score filter
	# A threshold on the HIPPIE confidence score can be chosen. The user can either specify a custom value between 0 and 1 or choose a predefined confidence level:
	# medium confidence (0.63 - second quartile of the HIPPIE score distribution)
	 # or high confidence (0.73 - third quartile).
	$hc = 0;
	$hc = 1 if ($score >= 0.63);
	$hc = 2 if ($score >= 0.73);
	
	$anns = $a[5];
	$systems = '';
	$pmids = '';
	$sources = '';
	$species = 'HUMAN';
	if (defined($anns))
	{
		foreach $ann (split(/;/, $anns))
		{
			if ($ann =~ /^experiments:(.+)/)
			{
				$systems = $1;
			}
			elsif ($ann =~ /^pmids:(.+)/)
			{
				$pmids = $1;
			}
			elsif ($ann =~ /^sources:(.+)/)
			{
				$sources = $1;
			}
			elsif ($ann =~ /^species:(.+)/)
			{
				if ($1 eq 'Mus musculus (Mouse)')
				{
					$species = 'MOUSE';
				}
				elsif ($1 eq 'Bos taurus (Bovine)')
				{
					$species = 'MOUSE';
				}
				else
				{
					# die("Error: Unhandled species '$1'")
					addme("unhandled species", $1);
				}
			}
			else
			{
				die("Error: Unhandled annotation type: '$ann'");
			}
		}
	}
	
	# # Replace separators
	# $systems =~ s/,/|/g;
	# $pmids =~ s/,/|/g;
	# $sources =~ s/,/|/g;
	
	# Replace separators & natsort
	$systems = join('|', unique(split(/,/, $systems)));
	$pmids = join('|', unique(split(/,/, $pmids)));
	$sources = join('|', unique(split(/,/, $sources)));
	
	
	
	# Look up UniProt accessions
	# our $superloudmysql = 1;
	$query1 = Query("SELECT DISTINCT primary_acc, name FROM uniacc WHERE name IN ('".join("', '", split(/\|/, $names1))."') AND species='human'");
	if (Numrows($query1) == 0)
	{
		addme("no acc in uniacc for names (skipped)", $names1);
		$acc1 = '';
		$name1 = $names1;
	}
	elsif (Numrows($query1) > 1)
	{
		# die("Error: Multiple primary_accs in uniacc for names '$names1'");
		addme("multiple primary_accs in uniacc for names (skipped)", $names1);
	}
	# else
	# {
	# 	($acc1, $name1) = FetchOne($query1);
	# }
	while (($acc1, $name1) = Fetch($query1))
	{
		$query2 = Query("SELECT DISTINCT primary_acc, name FROM uniacc WHERE name IN ('".join("', '", split(/\|/, $names2))."') AND species='human'");
		if (Numrows($query2) == 0)
		{
			addme("no acc in uniacc for names (skipped)", $names2);
			$acc2 = '';
			$name2 = $names2;
		}
		elsif (Numrows($query2) > 1)
		{
			# die("Error: Multiple primary_accs in uniacc for names '$names2'");
			addme("multiple primary_accs in uniacc for names (skipped)", $names2);
		}
		# else
		# {
		# 	($acc2, $name2) = FetchOne($query2);
		# }
		while (($acc2, $name2) = Fetch($query2))
		{
			# Get gene symbols
			$query = Query("SELECT DISTINCT gene FROM uniprot WHERE name='$name1'");
			($symbol1) = FetchOne($query);
			$query = Query("SELECT DISTINCT gene FROM uniprot WHERE name='$name2'");
			($symbol2) = FetchOne($query);
			
			
			# Insert into table
			$q = "INSERT INTO $table SET acc1='$acc1', acc2='$acc2', name1='$name1', name2='$name2', symbol1='$symbol1', symbol2='$symbol2', species='$species', ncbigene1='$ncbigene1', ncbigene2='$ncbigene2', systems='$systems', sources='$sources', pmids='$pmids'";
			$q =~ s/=''/=NULL/g;

			Query($q);
			
			addme("acc pair inserted", "$acc1|$acc2");
			addme("name pair inserted", "$name1|$name2");
			addme("names pair inserted", "$names1|$names2");
			addme("symbol pair inserted", "$symbol1|$symbol2");
			addme("ncbigene pair inserted", "$ncbigene1|$ncbigene2");

			stepme(1000);
		}
	}
}
stopme();
stoptime();

showmeall(1);

done();
