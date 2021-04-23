#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize




# Download

run("Download latest data", "download.pl");
run("Import metadata file (comes with the download)", "~/scripts/import.pl input/metadata.tsv clip_encode_metadata -overwrite");




# run

#DEBUG
# waitforalljobs();
#END DEBUG

# Main (clip_raw)

Clear('clip_raw_gene');
# Clear('clip_raw_transcript');
# Clear('clip_raw_exon');

run("Main Genes",	"raw.pl human GRCh38 gene 3 5");
# run("Main Transcripts",	"raw.pl human GRCh38 transcript 3 5");
# run("Main Exons", "raw.pl human GRCh38 exon 3 5");
# run("Add region type annotation to entries in clip_raw_gene (requires clip_raw_exon)", "add_regions.pl eclip_encode");

# Weaker p-value threshold: -log10p 3
run("Combined replicates only",	"raw.pl human GRCh38 gene 3 3 -12");
# run("Combined replicates only",	"raw.pl human GRCh38 exon 3 3 -12");
# run("Combined replicates only",	"raw.pl human GRCh38 transcript 3 3 -12");

done();
