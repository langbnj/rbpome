#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

$table = "biogrid";

# initialize

# our $superloudmysql = 1;

our $usage = "[input file]";
($infile) = args(1);

#$outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
#open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start
Clear($table);

startme("Inserting interactions from '$infile' into table '$table'");
starttime();

while (<IN>)
{
	chomp;
	
	# #BioGRID Interaction ID	Entrez Gene Interactor A	Entrez Gene Interactor B	BioGRID ID Interactor A	BioGRID ID Interactor B	Systematic Name Interactor A	Systematic Name Interactor B	Official Symbol Interactor A	Official Symbol Interactor B	Synonyms Interactor A	Synonyms Interactor B	Experimental System	Experimental System Type	Author	Pubmed ID	Organism Interactor A	Organism Interactor B	Throughput	Score	Modification	Phenotypes	Qualifications	Tags	Source Database
	# 103	6416	2318	112315	108607	-	-	MAP2K4	FLNC	JNKK|JNKK1|MAPKK4|MEK4|MKK4|PRKMK4|SAPKK-1|SAPKK1|SEK1|SERK1|SKK1	ABP-280|ABP280A|ABPA|ABPL|FLN2|MFM5|MPD4	Two-hybrid	physical	Marti A (1997)	9006895	9606	9606	Low Throughput	-	-	-	-	-	BIOGRID
	# 117	84665	88	124185	106603	-	-	MYPN	ACTN2	CMD1DD|CMH22|MYOP|RCM4	CMD1AA	Two-hybrid	physical	Bang ML (2001)	11309420	9606	9606	Low Throughput	-	BIOGRID

	@a = split(/\t/);
	
	# next if ($a[5] eq '-');
	# next if ($a[6] eq '-');
	die("Error: Gene 1 is '' for line:\n\n$_\n\n") if ($a[7] eq '-');
	die("Error: Gene 2 is '' for line:\n\n$_\n\n") if ($a[8] eq '-');
	
	$ncbigene1 = $a[2];
	$ncbigene2 = $a[3];
	
	$a[5] = '' if ($a[5] eq '-');
	$a[6] = '' if ($a[6] eq '-');
	
	$orderedlocus1 = $a[5];
	$orderedlocus2 = $a[6];
	
	$gene1 = $a[7];
	$gene2 = $a[8];
	
	$system = $a[11];
	$type = $a[12];
	
	$pmid = $a[14];
	
	$tax1 = $a[15];
	$tax2 = $a[16];
	
	$species = '';
	
	# NCBI Taxon IDs
    # if (($tax1 eq $tax2) and ($tax1 eq '3702')) { $species = 'arath'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '224308')) { $species = 'bsub'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '9031')) { $species = 'chick'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '9913')) { $species = 'cow'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '9615')) { $species = 'dog'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '511145')) { $species = 'ecoli'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '7227')) { $species = 'fly'; }
    if (($tax1 eq $tax2) and ($tax1 eq '9606')) { $species = 'human'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '9544')) { $species = 'macaque'; }
    if (($tax1 eq $tax2) and ($tax1 eq '10090')) { $species = 'mouse'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '10116')) { $species = 'rat'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '4896')) { $species = 'schpo'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '6239')) { $species = 'worm'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '8355')) { $species = 'xenopus'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '559292')) { $species = 'yeast'; }
    # if (($tax1 eq $tax2) and ($tax1 eq '7955')) { $species = 'zebrafish'; }
	
	# Direct interaction detected? (these are Sebastian Maurer's choices for "all_direct" interactions)
	$direct = 0;
	$direct = 1 if ($system eq 'Two-hybrid');
	$direct = 1 if ($system eq 'Reconstituted Complex');
	$direct = 1 if ($system eq 'Co-crystal Structure');
	
	
	# Insert into table
	if ($species ne '')
	{
		$q = "INSERT INTO `$table` SET species='$species', gene1='$gene1', gene2='$gene2', ncbigene1='$ncbigene1', ncbigene2='$ncbigene2', orderedlocus1='$orderedlocus1', orderedlocus2='$orderedlocus2', system='$system', type='$type', direct='$direct', pmid='$pmid'";
		$q =~ s/=''/=NULL/g;
        Query($q);
	}
	
	stepme(100000);
}
stopme();

Optimize($table);

stoptime();
done();


