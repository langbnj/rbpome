#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_exon';

our $usage = "$0";
args(0);

Clear($table);





# start
state("Inserting exons from GENCODE GFF3 files into table '$table'...");
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

        # ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = @a;   # GTF
        # ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = @a;        # GFF3
        $chr = $a[0];
        $source = $a[1];
        $feature = $a[2];
        $start = $a[3];
        $stop = $a[4];
        $strand = $a[6];
        # $phase = $a[7];
        $att = $a[8];

        die("Error: Seqid is . in '$_'") if ($chr eq '.');
        die("Error: Start is . in '$_'") if ($start eq '.');
        if ($start eq '###')
        {
            # die("Error: start is blank in '$_'");
            addme("start is ### for line (kept)", $_);
            addme("start is ### for species (kept)", $species);

            $start = '';
        }
        die("Error: Stop is . in '$_'") if ($stop eq '.');
        if ($stop eq '###')
        {
            # die("Error: stop is blank in '$_'");
            addme("stop is ### for line (kept)", $_);
            addme("stop is ### for species (kept)", $species);

            $stop = '';
        }
        if ($strand eq '.')
        {
            # die("Error: Strand is blank in '$_'");
            addme("strand is . for line (kept)", $_);
            addme("strand is . for species (kept)", $species);

            $strand = '';
        }
        if ($strand eq '###')
        {
            # die("Error: Strand is blank in '$_'");
            addme("strand is ### for line (kept)", $_);
            addme("strand is ### for species (kept)", $species);

            $strand = '';
        }

        # next if (($feature ne 'gene') and ($feature ne 'transcript') and ($feature ne 'CDS') and ($feature ne 'exon'));
        next if ($feature ne 'exon');

        # Parse attributes
        # Parent=transcript:ENSAMET00000017883;Name=ENSAMEE00000172154;constitutive=1;gencode_end_phase=0;gencode_phase=0;exon_id=ENSAMEE00000172154;rank=1;version=1
        # cat *.gff3 | g "\texon\t" | cut -f 9 | perl -ne 's/;/;\n/g; print' | perl -ne 's/=.+$/$1/g; print' | suq
        # 17279126 version
        # 17716605 Name
        # 17716605 Parent
        # 17716605 constitutive
        # 17716605 gencode_end_phase
        # 17716605 gencode_phase
        # 17716605 exon_id
        # 17716605 rank
        %att = ();
        foreach $s (split(/;/, $att))
        {
            ($a, $v) = split(/=/, $s);
            if (exists($att{$a}))
            {
                # die("Error: Duplicate attribute '$a' in '$att'") if ($a ne 'tag');
                $att{$a} .= "|$v";
            }
            else
            {
                $att{$a} = $v;
            }
        }

        $att{'ID'} =~ /^exon:(.+)/ or die("Error: Couldn't parse ID in '$_'");
        $id = $1;
        $ensev = $att{'exon_id'};
        # if ($id ne $enstv)
        # {
        #     # die("Error: ENSE mismatch between ID and exon_id in '$_'");
        #     addme("mismatch between id and ense for ense (kept ense)", $ense);
        #     addme("mismatch between id and ense for id|ense (kept ense)", "$id|$ense");
        #     addme("mismatch between id and ense for species (kept ense)", $species);
        # }
        # if (exists($att{'version'})) { $ensev = $ense.".".$att{'version'}; } else { addme("unversioned ense used as ensev for species", $species); $ensev = $ense; }
        # Workaround for GRCh37 remap (see gene.pl):
    		if ($ensev =~ /_\d$/)
    		{
          addme("_2, _3 etc. enses inserted before replace (GRCh37 remap workaround)", $ensev);
    			$ensev =~ s/_\d$//;
    			addme("_2, _3 etc. enses inserted after replace (GRCh37 remap workaround)", $ensev);
    		}
		$ense = $ensev;
		if ($ense =~ /_PAR_Y$/) { $ense =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $ense =~ s/\.\d+$//; }
        $enstv = $att{'Parent'};
		$enst = $enstv;
		if ($enst =~ /_PAR_Y$/) { $enst =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $enst =~ s/\.\d+$//; }
        # # $constitutive = $att{'constitutive'}; # Field doesn't exist in GENCODE
        # $constitutive = '';
        # # $gencode_end_phase = $att{'gencode_end_phase'}; # Field doesn't exist in GENCODE
        # $gencode_end_phase = '';
        # # $gencode_phase = $att{'gencode_phase'}; # Field doesn't exist in GENCODE
        # $gencode_phase = '';
        $exon_number = $att{'exon_number'};
		$biotype = $att{'transcript_type'};
        # Get CDS coordinates (like for Ensembl)
        $query = Query("SELECT start, stop FROM gencode_gff3_cds WHERE species='$species' AND enstv='$enstv' AND start>=$start AND stop<=$stop");
        if (Numrows($query) == 0)
        {
            # die("Error: No overlapping CDS found in table 'ensembl_gff3_cds' for exon '$ensev')");
            # die("Error: Biotype is 'protein_coding', but no CDS found in exon '$ensev'") if ($biotype eq 'protein_coding');
            addme("no cds in exon for ensev (kept)", $ensev);
            if ($biotype eq 'protein_coding')
            {
                addme("no cds in exon of protein_coding transcript for ensev (kept)", $ensev);
            }
            $coding = 0;
            $cds_start = '';
            $cds_stop = '';
        }
        else
        {
            addme("cds successfully found in exon for ensev", $ensev);
            addme("cds successfully found in exon for biotype", $biotype);
            ($cds_start, $cds_stop) = FetchOne($query);
            $coding = 1;

            die if ($cds_start !~ /^\d+$/);
            die if ($cds_stop !~ /^\d+$/);
        }


		# Get ENSG and ENSGV
		$query = Query("SELECT DISTINCT ensg, ensgv FROM gencode_gff3 WHERE species='$species' AND enst='$enst' AND enstv='$enstv'");
		($ensg, $ensgv) = FetchOne($query);
        
        $q = "INSERT INTO `$table` SET ensg='$ensg', ensgv='$ensgv', enst='$enst', enstv='$enstv', biotype='$biotype', ense='$ense', ensev='$ensev', exonid='$id', species='$species', source='$source', exon_number='$exon_number', chr='$chr', start='$start', stop='$stop', strand='$strand', coding='$coding', cds_start='$cds_start', cds_stop='$cds_stop'";
        $q =~ s/=''/=NULL/g;
        Query($q);

        addme("total ensgs inserted", $ensg);
        addme("total ensts inserted", $enst);
        addme("total enses inserted", $ense);
        addme("total ensevs inserted", $ensev);

        stepme(10000, 1);
    }
    close(IN);
    stopme(1);
}
nl();
stoptime();

showmesome(20);
nl();

Optimize($table);

done();
