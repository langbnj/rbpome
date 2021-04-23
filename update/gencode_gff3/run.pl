#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run
run("Parse genes", "gene.pl");
run("Parse transcripts", "transcript.pl");
run("Parse transcript sequences and translated (protein) sequences from GENCODE FASTA files", "transcript_add_sequences.pl");
run("Parse RefSeq_NT (NM_...) identifiers from GENCODE metadata files", "transcript_add_ncbiids.pl"); # Note: Some (very few, ~20) ENSTVs have multiple ncbiids (RefSeq_NT), hence the ncbiids field is concatenated (...|...) - also adding them to the 'ncbi' table for easier lookup (see ~/update/ncbi)
run("Flag longest transcripts per ENSGV", "transcript_flag_longest.pl ensgv");
run("Flag longest transcripts per gene symbol", "transcript_flag_longest.pl symbol");
run("Parse CDSs", "cds.pl");
run("Parse exons", "exon.pl");
run("Parse exon sequences and translated (protein) sequences from GENCODE FASTA files", "exon_add_sequences.pl");

# run("Combine into 'gencode_gtf'", "main.pl");
# -- query "SELECT ensg, enst, ensp, species, source, symbol, biotype, tsl, basic, ccds, gencode_tags FROM blang.gencode_gtf ORDER BY species, ensg, enst, ensp;" > gtf.txt
# -- query "SELECT g.ensg, t.enst, c.ensp, g.species, g.source, g.symbol, t.biotype, t.tsl, t.basic, t.ccds, t.tags AS gencode_tags FROM gencode_gff3_gene g, gencode_gff3_transcript t, gencode_gff3_cds c WHERE g.ensg=t.ensg AND t.enst=c.enst ORDER BY species, g.ensg, t.enst, c.ensp;" > gff3.txt
# -- diff gtf.txt gff3.txt

done();
