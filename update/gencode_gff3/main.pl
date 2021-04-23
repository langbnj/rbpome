#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3';
$ensembltable = 'ensembl';

our $usage = "$0";
args(0);




Clear($table);


# start

startme("Getting UniProt species codes for Ensembl species from table '$ensembltable'");
$query = Query("SELECT DISTINCT fullspecies, species FROM `$ensembltable`");
%unispec = ();
while (($species, $unispec) = Fetch($query))
{
    $unispec{$species} = $unispec;
    # state(" >> $species >> $unispec", 1);
    stepme(10);
}
stopme();



state("Reading transcript annotation from Ensembl GFF3 file (Gene Symbol, Transcript Biotype, Transcript Support Level, GENCODE Basic, CCDS ID and other Ensembl tags)...");
open(DIR, "ls input -1|");
while (<DIR>)
{
  chomp;

    # DEBUG
    # $_ = 'Homo_sapiens.GRCh38.84.chr_patch_hapl_scaff.gff3';

    /^(\w+)\.[\w\.\-]+\.gff3$/ or die("Error: Couldn't match file name '$_'");
    # $species = lc($1);
  $species = $unispec{$1};
  
  startme(" >> $species >> $_", 1);
    # Read GFF3
  open(IN, "input/$_") or die("Error: Couldn't open '$_'\n");
  while (<IN>)
  {
      chomp;
      
      # skip comments
      next if /^#/;
      
      @a = split(/\t/);
        # state($_);
        # show(@a);
        
        # ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = @a;
        $source = $a[1];
        $feature = $a[2];
        $attributes = $a[8];
        
        next if (($feature ne 'gene') and ($feature ne 'transcript') and ($feature ne 'CDS') and ($feature ne 'exon'));
        
        # gene:
        # 

        # transcript:
        # 

        # CDS:
        # ID=CDS:ENSP00000334393;Parent=transcript:ENST00000335137;protein_id=ENSP00000334393

        # exon:
        # 

        
        $symbol = '';
        $ensg = '';
        $enst = '';
        $ensp = '';
        $biotype = '';
        $tsl = '';
        @tags = ();
        $basic = '';
        $ccds = '';
        foreach $attribute (split/; ?/, $attributes)
        {
            # state($attribute, 1);

            $attribute =~ /^(\w+) "([^"]+)"$/ or die("Error: Couldn't parse attribute '$attribute'");
            $tag = $1;
            $value = $2;
            
            # addme("tags", $tag);
            # addme("tag '$tag' values", $value);

            $symbol = $value if ($tag eq 'gene_name');
            $ensg = $value if ($tag eq 'gene_id');
            $enst = $value if ($tag eq 'transcript_id');
            $ensp = $value if ($tag eq 'protein_id');
            $biotype = $value if ($tag eq 'transcript_biotype');
            if ($tag eq 'transcript_support_level')
            {
                $tsl = $value;
                $tsl =~ /^(\d+|NA)( \(assigned to previous version \d+\))?$/ or die("Error: Couldn't parse transcript support level '$tsl' in line:\n\n$_\n\n");
                $tsl = $1;
                $tsl = '' if ($tsl eq 'NA');
            }
            if ($tag eq 'tag')
            {
                if ($value eq 'basic')
                {
                    $basic = 1;
                }
                elsif ($value eq 'CCDS')
                {
                    # CCDS ID gets parsed below
                }
                else
                {
                    push(@tags, $value) if ($tag eq 'tag');
                }
            }
            $ccds = $value if ($tag eq 'ccds_id');
        }
        die("Error: Duplicate tags in line '$_'") if (scalar(unique(@tags)) != scalar(@tags));
        $tags = join("|", @tags);
        
        # if (($ensg eq '') or ($enst eq '') or ($ensp eq ''))
        # {
            die("Error: No ensg in line:\n\n$_\n\n") if ($ensg eq '');
            die("Error: No enst in line:\n\n$_\n\n") if ($enst eq '');
            die("Error: No ensp in line:\n\n$_\n\n") if ($ensp eq '');
            # if ($ensp eq '')
            # {
            #     addme("line skipped because ensp is blank for enst", $enst);
            # }
        #     next;
        # }
        
        addme("successfully parsed symbols in ensembl gff3", $symbol);
        addme("successfully parsed ensgs", $ensg);
        addme("successfully parsed ensts", $enst);
        addme("successfully parsed ensps", $ensp);
        addme("successfully parsed biotypes", $biotype);
        addme("successfully parsed tsls", $tsl);
        addme("successfully parsed sources", $source);
        addme("successfully parsed ccds ids", $ccds);
        foreach $tag (split(/\|/, $tags)) { addme("successfully parsed tags", $tag); }
        
        $q = "INSERT INTO `$table` SET ensg='$ensg', enst='$enst', ensp='$ensp', species='$species', source='$source', symbol='".esc($symbol)."', biotype='$biotype', tsl='$tsl', basic='$basic', ccds='$ccds', gencode_tags='$tags'";
        $q =~ s/=''/=NULL/g;
        $query = Query($q);
        # state($q);

        stepme(10000, 1);
        
        # DEBUG
        # last if (getme() % 100000 == 0);
    }
    close(IN);
    stopme(1);

    # DEBUG
    # last;
}


# # Workaround: 9 Flybase transcripts are missing from the Ensembl GTF, but they are in the FASTA (along with their biotype and gene symbol) (and in the GFF3). Here, I'm getting their biotypes from the FASTA (biotype is a field that's always annotated in the GTF).
# $mainquery = Query("SELECT ensgv, enstv, enspv FROM `$table` WHERE biotype IS NULL");
# if (Numrows($mainquery) != 0)
# {
#     nl();
#     state("Using workaround to map biotypes for ".Numrows($mainquery)." rows from Ensembl FASTA (they're missing from the GTF):");
#     warn("Warning: Using workaround to map biotypes for ".Numrows($mainquery)." rows from Ensembl FASTA (they're missing from the GTF)");
#     # >FBpp0083480 pep:known chromosome:BDGP6:3R:21368380:21377399:1 gene:FBgn0002781 transcript:FBtr0084081 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083483 pep:known chromosome:BDGP6:3R:21366004:21377399:1 gene:FBgn0002781 transcript:FBtr0084084 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083484 pep:known chromosome:BDGP6:3R:21361398:21377399:1 gene:FBgn0002781 transcript:FBtr0084085 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083481 pep:known chromosome:BDGP6:3R:21367910:21377399:1 gene:FBgn0002781 transcript:FBtr0084082 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083482 pep:known chromosome:BDGP6:3R:21366450:21377399:1 gene:FBgn0002781 transcript:FBtr0084083 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0300270 pep:known chromosome:BDGP6:3R:21367363:21377399:1 gene:FBgn0002781 transcript:FBtr0307759 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0300271 pep:known chromosome:BDGP6:3R:21366974:21377399:1 gene:FBgn0002781 transcript:FBtr0307760 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083478 pep:known chromosome:BDGP6:3R:21360390:21377399:1 gene:FBgn0002781 transcript:FBtr0084079 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     # >FBpp0083479 pep:known chromosome:BDGP6:3R:21369462:21377401:1 gene:FBgn0002781 transcript:FBtr0084080 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:mod(mdg4) description:modifier of mdg4 [Source:FlyBase gene name;Acc:FBgn0002781]
#     while (($ensgv, $enstv, $enspv) = Fetch($mainquery))
#     {
#         $_ = chompme(`grep -iP 'transcript:$enstv' ~/update/ensembl/input/Drosophila_melanogaster.BDGP6.pep.all.fa`);
#         
#         /transcript_biotype:(\S+)/ or die("Error: Couldn't get transcript_biotype from line:\n\n$_\n\n");
#         $biotype = $1;
# 
#         /gene_symbol:(\S+)/ or die("Error: Couldn't get gene_symbol from line:\n\n$_\n\n");
#         $symbol = $1;
#         
#         /Source:(\S+)/ or die("Error: Couldn't get source from line:\n\n$_\n\n");
#         $source = $1;
#         
#         /gene:(\S+)/ or die("Error: Couldn't get source from line:\n\n$_\n\n");
#         die("Error: Gene mismatch (should be '$ensgv', but is '$1') in line:\n\n$_\n\n") if ($1 ne $ensgv);
#         
#         state(" >> $ensgv|$enstv|$enspv >> $source >> $symbol >> $biotype", 1);
#         $query = Query("SELECT id FROM `$table` WHERE ensgv='$ensgv' AND enstv='$enstv' AND enspv='$enspv'");
#         die("Error: Multiple or zero lines in table '$table' for '$ensgv|$enstv|$enspv'") if (Numrows($query) != 1);
#         $query = Query("UPDATE `$table` SET source='$source', symbol='".esc($symbol)."', biotype='$biotype' WHERE ensgv='$ensgv' AND enstv='$enstv' AND enspv='$enspv'");
#         addme("Workaround: Successfully added source ('ensembl' column), gene symbol and biotype information from Ensembl FASTA for enstv", $enstv);
#     }
# }



# showmesome(10000);
showmesome(50);
# showmeall();

Optimize($table);

done();
