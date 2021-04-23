#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download GENCODE human and mouse
run("Create input subdirectory for GENCODE", "mkdir -p input/sorted");

foreach $species ('human', 'mouse')
{
    $rel = 27 if ($species eq 'human');
    $rel = 'M16' if ($species eq 'mouse');

    if ($species eq 'human')
    {
		# GRCh37-mapped release (for Joehanes et al. eQTLs, which were still done on GRCh37)
		# Only available for human
		run("Download", "wget -O input/gencode.".$species."37.v$rel.chr_patch_hapl_scaff.annotation.gff3.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/GRCh37_mapping/gencode.v".$rel."lift37.annotation.gff3.gz'");
		run("Unpack", "gunzip -f input/gencode.".$species."37.v$rel.chr_patch_hapl_scaff.annotation.gff3.gz");

		# Sort (for bedtools intersect)
		foreach $type ('transcript', 'exon')
		{
			$infile = $outfile = "input/gencode.".$species."37.v$rel.chr_patch_hapl_scaff.annotation.gff3";
			$outfile =~ s/\.gff3$/\.$type\.sorted\.gff3/;
			$outfile =~ s/^input\//input\/sorted\//;
			run("Sort: Get header", qq(grep '^#' $infile > $outfile));
			run("Sort: Sort", qq(grep -v '^#' $infile | grep -P '\\t$type\\t' | sort -k1,1 -k4,4n -k5,5n >> $outfile));
		}
	
	    # run("Download", "wget -O input/GRCh38.p10.genome.fa.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz'");
	    run("Download", "wget -O input/GRCh38.p10.genome.fa.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/GRCh38.p10.genome.fa.gz'");
	    run("Unpack", "gunzip -f input/GRCh38.p10.genome.fa.gz");
    }

    run("Download", "wget -O input/gencode.$species.v$rel.chr_patch_hapl_scaff.annotation.gff3.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/gencode.v$rel.chr_patch_hapl_scaff.annotation.gff3.gz'");
    run("Unpack", "gunzip -f input/gencode.$species.v$rel.chr_patch_hapl_scaff.annotation.gff3.gz");

	# Sort (for bedtools intersect)
	foreach $type ('transcript', 'exon')
	{
		$infile = $outfile = "input/gencode.$species.v$rel.chr_patch_hapl_scaff.annotation.gff3";
		$outfile =~ s/\.gff3$/\.$type\.sorted\.gff3/;
		$outfile =~ s/^input\//input\/sorted\//;
		run("Sort: Get header", "grep '^#' $infile > $outfile");
		run("Sort: Sort", "grep -v '^#' $infile | grep -P '\\t$type\\t' | sort -k1,1 -k4,4n -k5,5n >> $outfile");
	}

    run("Download", "wget -O input/gencode.$species.v$rel.metadata.Annotation_remark.txt.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/gencode.v$rel.metadata.Annotation_remark.gz'");
    run("Unpack", "gunzip -f input/gencode.$species.v$rel.metadata.Annotation_remark.txt.gz");

    run("Download", "wget -O input/gencode.$species.v$rel.metadata.Transcript_supporting_feature.txt.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/gencode.v$rel.metadata.Transcript_supporting_feature.gz'");
    run("Unpack", "gunzip -f input/gencode.$species.v$rel.metadata.Transcript_supporting_feature.txt.gz");

    run("Download", "wget -O input/gencode.$species.v$rel.pc_translations.fa.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/gencode.v$rel.pc_translations.fa.gz'");
    run("Unpack", "gunzip -f input/gencode.$species.v$rel.pc_translations.fa.gz");

    run("Download", "wget -O input/gencode.$species.v$rel.transcripts.fa.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$species/release_$rel/gencode.v$rel.transcripts.fa.gz'");
    run("Unpack", "gunzip -f input/gencode.$species.v$rel.transcripts.fa.gz");
	

    state("Downloaded '$species' GENCODE release '$rel'");
}

# run("Add '.txt' to file names", "ls -1 input | xargs -I{} mv input/{} input/{}.txt");



# Note: GENCODE human 27 is equivalent (nearly identical) to Ensembl releases 90 and 91.
# Only GENCODE has 45 _PAR_Y ENSGVs.
# Meanwhile, the following 87 ENSGVs are only in Ensembl:
# ENSG00000274081.4
# ENSG00000277603.4
# ENSG00000280467.2
# ENSG00000280606.3
# ENSG00000280674.2
# ENSG00000280873.3
# ENSG00000281178.2
# ENSG00000281277.3
# ENSG00000281434.2
# ENSG00000281557.1
# ENSG00000281702.3
# ENSG00000281897.2
# ENSG00000282243.1
# ENSG00000282513.1
# ENSG00000282833.2
# ENSG00000282838.1
# ENSG00000282848.1
# ENSG00000282873.1
# ENSG00000282878.1
# ENSG00000282892.1
# ENSG00000282903.1
# ENSG00000282908.1
# ENSG00000282928.1
# ENSG00000282932.2
# ENSG00000282956.1
# ENSG00000282962.1
# ENSG00000282971.1
# ENSG00000282972.1
# ENSG00000282974.1
# ENSG00000283007.1
# ENSG00000283017.1
# ENSG00000283021.1
# ENSG00000283022.1
# ENSG00000283024.1
# ENSG00000283026.1
# ENSG00000283030.1
# ENSG00000283034.1
# ENSG00000283035.1
# ENSG00000283067.1
# ENSG00000283079.1
# ENSG00000283087.1
# ENSG00000283102.1
# ENSG00000283109.1
# ENSG00000283244.1
# ENSG00000283530.1
# ENSG00000283651.1
# ENSG00000283658.1
# ENSG00000283668.1
# ENSG00000283681.1
# ENSG00000283715.1
# ENSG00000283730.1
# ENSG00000283746.1
# ENSG00000283747.1
# ENSG00000283802.1
# ENSG00000283806.1
# ENSG00000283825.1
# ENSG00000283837.1
# ENSG00000283852.1
# ENSG00000283951.1
# ENSG00000283953.1
# ENSG00000283964.1
# ENSG00000283965.1
# ENSG00000283975.1
# ENSG00000283997.1
# ENSG00000284004.1
# ENSG00000284046.1
# ENSG00000284053.1
# ENSG00000284063.1
# ENSG00000284096.1
# ENSG00000284101.1
# ENSG00000284113.1
# ENSG00000284126.1
# ENSG00000284137.1
# ENSG00000284208.1
# ENSG00000284211.1
# ENSG00000284245.1
# ENSG00000284302.1
# ENSG00000284367.1
# ENSG00000284384.1
# ENSG00000284389.1
# ENSG00000284400.1
# ENSG00000284409.1
# ENSG00000284462.1
# ENSG00000284470.1
# ENSG00000284539.1
# ENSG00000284581.1
# ENSG00000284590.1




# show directory
run("ls", "ls -1 input");

done();
