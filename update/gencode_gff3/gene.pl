#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_gene';

our $usage = "$0";
args(0);

Clear($table);





# start
state("Inserting genes from GENCODE GFF3 files into table '$table'...");
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
        # $phase = $a[7];
        $att = $a[8];

        # next if (($feature ne 'gene') and ($feature ne 'transcript') and ($feature ne 'CDS') and ($feature ne 'exon'));
        next if ($feature ne 'gene');

        # Parse attributes
        # ID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2
        # cat *.gff3 | g "\tgene\t" | cut -f 9 | perl -ne 's/;/;\n/g; print' | perl -ne 's/=.+$/$1/g; print' | suq

        # Ensembl
        # 101931 havana_gene
        # 101931 havana_version
        # 1060786 Name
        # 1154564 description
        # 1266601 version
        # 1323360 ID
        # 1323360 biotype
        # 1323360 gene_id
        # 1323360 logic_name

        # GENCODE
        # 20928 tag
        # 91486 havana_gene
        # 114229 ID
        # 114229 gene_id
        # 114229 gene_name
        # 114229 gene_status
        # 114229 gene_type
        # 114229 level

        %att = ();
        foreach $s (split(/;/, $att))
        {
            ($a, $v) = split(/=/, $s);
            if (exists($att{$a}))
            {
                die("Error: Duplicate attribute '$a' in '$att'") if ($a ne 'tag');
                $att{$a} .= "|$v";
            }
            else
            {
                $att{$a} = $v;
            }
        }

        $id = $att{'ID'};
        if (exists($att{'gene_name'})) { $symbol = $att{'gene_name'}; } else { $symbol = ''; }
        if (exists($att{'description'})) { $desc = $att{'description'}; } else { $desc = ''; }
        $biotype = $att{'gene_type'};
		# Field doesn't exist in GENCODE
        # $status = $att{'gene_status'};
        # $status = '';
        $level = $att{'level'};
        $ensgv = $att{'gene_id'};
        if ($ensgv =~ /_\d$/)
    		{
          addme("_2, _3 etc. ensgs inserted before replace (GRCh37 remap workaround)", $ensgv);
    			$ensgv =~ s/_\d$//;
    			addme("_2, _3 etc. ensgs inserted after replace (GRCh37 remap workaround)", $ensgv);
    		}
		$ensg = $ensgv;
		if ($ensg =~ /_PAR_Y$/) { $ensg =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $ensg =~ s/\.\d+$//; }
        if (exists($att{'tag'})) { $tags = $att{'tag'}; } else { $tags = ''; }
        if (contains('basic', split(/\|/, $tags))) { $basic = 1; $tags =~ s/^basic\|//; $tags =~ s/\|basic\|/\|/; $tags =~ s/\|basic$//; } else { $basic = ''; }

		# Workaround for GRCh37 remap:
		# It contains e.g. ENSG00000223972.5_2, which has ID=ENSG00000223972.5 (so they don't match because of the _2).
		# Removing all _2, _3 etc. here for that reason.
		# Funnily, a _3 case (ENSG00000243485.5) doesn't have a corresponding _2.
        if ($id ne $ensgv)
        {
            # addme("mismatch between id and ensgv for ensgv (kept ensgv)", $ensgv);
            # addme("mismatch between id and ensgv for id|ensgv (kept ensgv)", "$id|$ensgv");
            # addme("mismatch between id and ensgv for species (kept ensgv)", $species);

            # addme("mismatch between id and ensgv due to _PAR_Y for ensgv (kept id, for _PAR_Y)", $ensgv);
            # addme("mismatch between id and ensgv due to _PAR_Y for id|ensgv (kept id, for _PAR_Y)", "$id|$ensgv");
            # addme("mismatch between id and ensgv due to _PAR_Y for species (kept id, for _PAR_Y)", $species);

			die("Error: Mismatch between id and ensgv not due to _PAR_Y for id '$id', ensgv '$ensgv', species '$species'") if ($id ne $ensgv.'_PAR_Y');
			$ensgv = $id;
			$ensg = $ensgv;
			if ($ensg =~ /_PAR_Y$/) { $ensg =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $ensg =~ s/\.\d+$//; }

            addme("_PAR_Y ensgs inserted", $ensg);
            # addme("_PAR_Y ensgvs inserted", $ensgv);
        }

        $q = "INSERT INTO `$table` SET ensg='$ensg', ensgv='$ensgv', species='$species', source='$source', symbol='".esc($symbol)."', biotype='$biotype', level='$level', chr='$chr', start='$start', stop='$stop', strand='$strand'";
        $q =~ s/=''/=NULL/g;
        Query($q);

        addme("total ensgs inserted", $ensg);
        addme("total ensgvs inserted", $ensgv);

        stepme(10000, 1);
    }
    close(IN);
    stopme(1);
}
nl();
stoptime();

showmesome(20);
# showmesome(50);
nl();

Optimize($table);

done();
