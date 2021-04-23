#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_cds';
$ensgtable = 'gencode_gff3_transcript';

our $usage = "$0";
args(0);

Clear($table);





# start
state("Inserting CDSs from GENCODE GFF3 files into table '$table'...");
starttime();
open(DIR, "ls input/*.gff3 -1|");
while (<DIR>)
{
    chomp;
    
    # DEBUG
    # $_ = 'Homo_sapiens.GRCh38.84.chr_patch_hapl_scaff.gff3';
    # $_ = 'Drosophila_melanogaster.BDGP6.84.gff3';
    
    /^input\/gencode\.(\w+)\.[\w\.\-]+\.gff3$/ or die("Error: Couldn't match file name '$_'");
    $species = uc($1);
    startme(" >> $species >> $_", 1);

    # Read GFF3
    open(IN, $_) or die("Error: Couldn't open '$_'\n");
    while (<IN>)
    {
        chomp;
        
        # skip comments
		next if /^#/;
      
		@a = split(/\t/);
        # state($_);
        # show(@a);
        
        # ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = @a;
        $chr = $a[0];
        $source = $a[1];
        $feature = $a[2];
        $start = $a[3];
        $stop = $a[4];
        $strand = $a[6];
        $phase = $a[7];
        $att = $a[8];
        
        # next if (($feature ne 'gene') and ($feature ne 'transcript') and ($feature ne 'CDS') and ($feature ne 'exon'));
        next if ($feature ne 'CDS');

        # Parse attributes
        # ID=CDS:ENSP00000334393;Parent=transcript:ENST00000335137;protein_id=ENSP00000334393
        %att = ();
        foreach $s (split(/;/, $att))
        {
            ($a, $v) = split(/=/, $s);
            die("Error: Duplicate attribute '$a' in '$att'") if (exists($att{$a}));
            $att{$a} = $v;
        }
        
        $att{'ID'} =~ /^CDS:(.+)/ or die("Error: Couldn't parse ID in '$_'");
        $id = $1;
        $enstv = $att{'Parent'};
		$enst = $enstv;
		if ($enst =~ /_PAR_Y$/) { $enst =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $enst =~ s/\.\d+$//; }
        $enspv = $att{'protein_id'};
		$ensp = $enspv;
		if ($ensp =~ /_PAR_Y$/) { $ensp =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $ensp =~ s/\.\d+$//; }
        if ($id ne $enstv)
        {
            # die("Error: enstv mismatch between ID and protein_ID in '$_'");
            addme("mismatch between id and enstv for enstv (kept enstv)", $enstv);
            addme("mismatch between id and enstv for id|enstv (kept enstv)", "$id|$enstv");
            addme("mismatch between id and enstv for species (kept enstv)", $species);
        }
		
		# Get ENSG and ENSGV
		$query = Query("SELECT DISTINCT ensg, ensgv FROM gencode_gff3 WHERE species='$species' AND enst='$enst' AND enstv='$enstv'");
		($ensg, $ensgv) = FetchOne($query);
        
        Query("INSERT INTO `$table` SET ensg='$ensg', ensgv='$ensgv', enst='$enst', enstv='$enstv', ensp='$ensp', enspv='$enspv', species='$species', source='$source', chr='$chr', start='$start', stop='$stop', strand='$strand', phase='$phase'");

        addme("total ensgs inserted", $ensg);
        addme("total ensgvs inserted", $ensgv);
        addme("total ensts inserted", $enst);
        addme("total enstvs inserted", $enstv);
        addme("total ensps inserted", $ensp);
        addme("total enspvs inserted", $enspv);

        stepme(10000, 1);
    }
    close(IN);
    stopme(1);
}
nl();
stoptime();

showmesome(20);
nl();

state("Note - There are slightly fewer ENSPs than ENSTs exclusively because of _PAR_Y transcripts, i.e. ENST A from chromosome X and ENST B from chromosome Y produce the same ENSP.");
# SELECT ensp, GROUP_CONCAT(DISTINCT enst), COUNT(DISTINCT enst) AS c FROM blang.gencode_gff3_cds GROUP BY ensp ORDER BY c DESC;
# SELECT enspv, GROUP_CONCAT(DISTINCT enstv), COUNT(DISTINCT enstv) AS c FROM blang.gencode_gff3_cds GROUP BY enspv ORDER BY c DESC;

Optimize($table);

done();
