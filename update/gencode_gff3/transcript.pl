#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_transcript';

our $usage = "$0";
args(0);

Clear($table);





# start
state("Inserting transcripts from GENCODE GFF3 files into table '$table'...");
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
        next if ($feature ne 'transcript');

        # Parse attributes
        # ID=
        # cat *.gff3 | g "\ttranscript\t" | cut -f 9 | perl -ne 's/;/;\n/g; print' | perl -ne 's/=.+$/$1/g; print' | suq

        # Ensembl
        # 69283 ccdsid
        # 120571 tag
        # 161007 transcript_support_level
        # 175333 havana_transcript
        # 175333 havana_version
        # 1381998 Name
        # 1577444 version
        # 1670462 ID
        # 1670462 Parent
        # 1670462 biotype
        # 1670462 transcript_id

        # GENCODE
        # 24841 ont
        # 69916 ccdsid
        # 160026 protein_id
        # 236563 tag
        # 295337 havana_transcript
        # 309281 havana_gene
        # 326579 transcript_support_level
        # 333109 ID
        # 333109 Parent
        # 333109 gene_id
        # 333109 gene_name
        # 333109 gene_status
        # 333109 gene_type
        # 333109 level
        # 333109 transcript_id
        # 333109 transcript_name
        # 333109 transcript_status
        # 333109 transcript_type

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
        $ensgv = $att{'Parent'};
        $ensg = $ensgv;
		if ($ensg =~ /_PAR_Y$/) { $ensg =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $ensg =~ s/\.\d+$//; }
		# In GENCODE, every transcript has a gene_name and a transcript_name (which has the 3-digit number)
        # if (exists($att{'Name'}))
        # {
        #     $symbolv = $att{'Name'};
        #     $symbol = $symbolv;
        #     # Check if transcript symbol always has a 3 digit number after the gene symbol
        #     if ($symbol !~ /-\d{3}$/)
        #     {
        #         # die("Error: Couldn't parse transcript symbol '$symbolv'");
        #         addme("non-numbered transcript symbol for symbolv", $symbolv);
        #         addme("non-numbered transcript symbol for species", $species);
        #     }
        #     $symbol =~ s/-\d{3}$//;
        # }
        # else
        # {
        #     $symbolv = '';
        #     $symbol = '';
        # }
		$symbol = $att{'gene_name'};
		$symbolv = $att{'transcript_name'};
        $biotype = $att{'transcript_type'};
        $enst = $att{'transcript_id'};
        # if (exists($att{'version'})) { $enstv = $enst.".".$att{'version'}; } else { addme("unversioned enst used as enstv for species", $species); $enstv = $enst; }		# GENCODE always has versioned ENSTs (mouse/human)
        # Workaround for GRCh37 remap (see gene.pl):
    		if ($enst =~ /_\d$/)
    		{
          addme("_2, _3 etc. ensts inserted before replace (GRCh37 remap workaround)", $enst);
    			$enst =~ s/_\d$//;
    			addme("_2, _3 etc. ensts inserted after replace (GRCh37 remap workaround)", $enst);
    		}
		$enstv = $enst;
		if ($enst =~ /_PAR_Y$/) { $enst =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $enst =~ s/\.\d+$//; }
        if (exists($att{'transcript_support_level'})) { $tsl = $att{'transcript_support_level'}; } else { $tsl = ''; }
        $tsl =~ s/ \(assigned to previous version \d+\)//;
        $tsl = '' if ($tsl eq 'NA');
        if (exists($att{'ccdsid'}))
        {
            $ccdsv = $att{'ccdsid'};
            $ccds = $ccdsv;
            die("Error: Couldn't match versioned CCDS ID '$ccdsv") if ($ccds !~ /\.\d+$/);
            $ccds = s/\.\d+$//;
        }
        else
        {
            $ccdsv = '';
            $ccds = '';
        }
        if (exists($att{'tag'})) { $tags = $att{'tag'}; } else { $tags = ''; }


		# GENCODE BASIC
        # if (contains('basic', split(/\|/, $tags))) { $basic = 1; $tags =~ s/^basic\|//; $tags =~ s/\|basic\|/\|/; $tags =~ s/\|basic$//; $tags =~ s/^basic$//; } else { $basic = ''; }
        if (contains('basic', split(/,/, $tags))) { $basic = 1; $tags =~ s/^basic,//; $tags =~ s/,basic,/,/; $tags =~ s/,basic$//; $tags =~ s/^basic$//; } else { $basic = ''; }
		
		# APPRIS levels
		# 
		# http://appris.bioinfo.cnio.es/#/downloads:
		# Principal Isoform flags
		#
		# APPRIS selects a single CDS variant for each gene as the 'PRINCIPAL' isoform based on the range of protein features. Principal isoforms are tagged with the numbers 1 to 5, with 1 being the most reliable. The definition of the flags are as follows:
		#
		# PRINCIPAL:1
		# Transcript(s) expected to code for the main functional isoform based solely on the core modules in the APPRIS database. The APPRIS core modules map protein structural and functional information and cross-species conservation to the annotated variants.
		#
		# PRINCIPAL:2
		# Where the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes), the database chooses two or more of the CDS variants as "candidates" to be the principal variant.
		#
		# If one of these candidates has a distinct CCDS identifier it is selected as the principal variant for that gene. A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.
		#
		# PRINCIPAL:3
		# Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants have distinct CCDS identifiers, APPRIS selects the variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated.
		#
		# Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS identifiers are not included in this flag, since they will have been annotated in the same release of CCDS. These are distinguished with the next flag.
		#
		# In addition, there is more than one variant with a distinct (but consecutive) CCDS identifiers, APPRIS choose the variant when all splice junctions are supported by at least one non-suspect mRNA. This information is reported by the method Transcript Support Level (TSL), which is a method to highlight the well-supported and poorly-supported transcript models for users. The method relies on the primary data that can support full-length transcript structure: mRNA and EST alignments supplied by UCSC and Ensembl.
		#
		# PRINCIPAL:4
		# Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one variant with a distinct (but consecutive) CCDS identifiers and all the splice junctions are not well-supported, APPRIS selects the longest CCDS isoform as the principal variant.
		#
		# PRINCIPAL:5
		# Where the APPRIS core modules are unable to choose a clear principal variant and none of the candidate variants are annotated by CCDS, APPRIS selects the longest of the candidate isoforms as the principal variant.
		#
		# For genes in which the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes) the "candidate" variants not chosen as principal are labeled in the following way:
		#
		# ALTERNATIVE:1
		# Candidate transcript(s) models that are conserved in at least three tested non-primate species.
		#
		# ALTERNATIVE:2
		# Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.
		#
		# Non-candidate transcripts are not flagged and are considered as "MINOR" transcripts.
		# 
		# 
		# [ant-login7] ~/update/gencode_gff3/input >> g -o "appris_\w+" gencode.human.v27.chr_patch_hapl_scaff.annotation.gff3 | suq
		#  478546 appris_principal_1
		#   55497 appris_principal_2
		#   93136 appris_principal_3
		#   26286 appris_principal_4
		#    7785 appris_principal_5
		#   50636 appris_alternative_1
		#  224984 appris_alternative_2
		$appris = '';
        if (contains('appris_principal_1', split(/,/, $tags))) { $appris = 'P1'; $tags =~ s/^appris_principal_1,//; $tags =~ s/,appris_principal_1,/,/; $tags =~ s/,appris_principal_1$//; $tags =~ s/^appris_principal_1$//; }
        if (contains('appris_principal_2', split(/,/, $tags))) { $appris = 'P2'; $tags =~ s/^appris_principal_2,//; $tags =~ s/,appris_principal_2,/,/; $tags =~ s/,appris_principal_2$//; $tags =~ s/^appris_principal_2$//; }
        if (contains('appris_principal_3', split(/,/, $tags))) { $appris = 'P3'; $tags =~ s/^appris_principal_3,//; $tags =~ s/,appris_principal_3,/,/; $tags =~ s/,appris_principal_3$//; $tags =~ s/^appris_principal_3$//; }
        if (contains('appris_principal_4', split(/,/, $tags))) { $appris = 'P4'; $tags =~ s/^appris_principal_4,//; $tags =~ s/,appris_principal_4,/,/; $tags =~ s/,appris_principal_4$//; $tags =~ s/^appris_principal_4$//; }
        if (contains('appris_principal_5', split(/,/, $tags))) { $appris = 'P5'; $tags =~ s/^appris_principal_5,//; $tags =~ s/,appris_principal_5,/,/; $tags =~ s/,appris_principal_5$//; $tags =~ s/^appris_principal_5$//; }
        if (contains('appris_alternative_1', split(/,/, $tags))) { $appris = 'ALT1'; $tags =~ s/^appris_alternative_1,//; $tags =~ s/,appris_alternative_1,/,/; $tags =~ s/,appris_alternative_1$//; $tags =~ s/^appris_alternative_1$//; }
        if (contains('appris_alternative_2', split(/,/, $tags))) { $appris = 'ALT2'; $tags =~ s/^appris_alternative_2,//; $tags =~ s/,appris_alternative_2,/,/; $tags =~ s/,appris_alternative_2$//; $tags =~ s/^appris_alternative_2$//; }


        if ($id ne $enstv)
        {
            # die("Error: ENST mismatch between ID and transcript_id in '$_'");
            # addme("mismatch between id and enst for enst (kept enst)", $enst);
            # addme("mismatch between id and enst for id|enst (kept enst)", "$id|$enst");
            # addme("mismatch between id and enst for species (kept enst)", $species);

            # addme("mismatch between id and enstv due to _PAR_Y for enstv (kept id, for _PAR_Y)", $enstv);
            # addme("mismatch between id and enstv due to _PAR_Y for id|enstv (kept id, for _PAR_Y)", "$id|$enstv");
            # addme("mismatch between id and enstv due to _PAR_Y for species (kept id, for _PAR_Y)", $species);

			die("Error: Mismatch between id and enstv not due to _PAR_Y for id '$id', enstv '$enstv', species '$species'") if ($id ne $enstv.'_PAR_Y');
			$enstv = $id;
			$enst = $enstv;
			if ($enst =~ /_PAR_Y$/) { $enst =~ s/\.\d+_PAR_Y$/_PAR_Y/; } else { $enst =~ s/\.\d+$//; }

            # addme("_PAR_Y ensts inserted", $enst);
            addme("_PAR_Y enstvs inserted", $enstv);
        }

        $q = "INSERT INTO `$table` SET ensg='$ensg', ensgv='$ensgv', enst='$enst', enstv='$enstv', species='$species', source='$source', symbol='".esc($symbol)."', symbolv='".esc($symbolv)."', biotype='$biotype', tsl='$tsl', basic='$basic', appris='$appris', ccds='$ccds', ccdsv='$ccdsv', tags='$tags', chr='$chr', start='$start', stop='$stop', strand='$strand'";
        $q =~ s/=''/=NULL/g;
        Query($q);

        addme("total ensgs inserted", $ensg);
        addme("total ensts inserted", $enst);
        addme("total enstvs inserted", $enstv);

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
