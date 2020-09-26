#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download Table S3 from Dominguez et al. (2018) (Chris M. Burge) (https://www.sciencedirect.com/science/article/pii/S1097276518303514#app2)
run("Download", "wget -O input/dominguez-s3.xlsx 'https://ars.els-cdn.com/content/image/1-s2.0-S1097276518303514-mmc4.xlsx'");
# Use sheet 1 & ignore sheets 2-4. These have "proportions" with which the motifs are represented in the logos they drew, but they're not real "enrichment scores" (these are only in sheet 1, "stepwise R-1").
# Exported sheet 1 to TSV using Excel: dominguez-s3.txt
run("dos2unix", "dos2unix input/dominguez-s3.txt");
# Raw RNA-Bind-N-Seq data (not used here): https://www.encodeproject.org/search/?type=Experiment&assay_title=RNA+Bind-N-Seq&assay_title=RNA+Bind-n-Seq



# MEME
run("Download MEME Suite", "wget -O bin/meme-5.0.5.tar.gz 'http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz'");
cd("bin");
run("Unpack", "tar -xavpf meme-5.0.5.tar.gz");
cd("..");
# Get MEME Databases
run("Download MEME Databases", "wget -O input/motif_databases.12.18.tgz 'http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz'");
cd("input");
run("Unpack", "tar -xavpf motif_databases.12.18.tgz");
cd("..");
# MEME databases to use:
#
# 1) RNA/Ray2013_rbp_Homo_sapiens.meme 102 motifs, between 7 and 8 in width (average width 7.1).
# Homo sapiens RNA-binding motifs from Ray et al., Nature 11:172-177, 2013. They are from in vitro experiments using the RNAcompete method for the rapid and systematic analysis of the RNA sequence preferences of RBPs. The motifs are converted from the data in this archive file titled "Top10align PFMs learned from all data" that is supplied by the Ray et al. at their website.
# 
# Not using this one, it's identical to the one above (it only differs in the header, where it specifies ACGT instead of ACGU)
# # RNA/Ray2013_rbp_Homo_sapiens.dna_encoded.meme 102 motifs, between 7 and 8 in width (average width 7.1).
# # Homo sapiens RNA-binding motifs from Ray et al., Nature 11:172-177, 2013. They are from in vitro experiments using the RNAcompete method for the rapid and systematic analysis of the RNA sequence preferences of RBPs. The motifs are converted from the data in this archive file titled "Top10align PFMs learned from all data" that is supplied by the Ray et al. at their website.
# 
# 2) CISBP-RNA/Homo_sapiens.meme 98 motifs, between 4 and 18 in width (average width 7.2).
# Direct and inferred motifs for Homo sapiens from the CISBP-RNA (Catalog of Inferred Sequence Binding Preferences of RNA binding proteins) database. To reduce redundancy, for each RNA-binding protein (RBP) in this species that has a CISBP-RNA motif, we selected a single motif according to the following precedence rules. We chose the direct motif if there is one, otherwise we chose the inferred motif with the highest RNA binding domain (RBD) similarity (according to CISBP-RNA) to an RBP in another species that has a direct motif. If there is more than one direct motif or inferred motif with the highest RBD similarity, we chose among them according to their provenance (CISBP-RNA's "Motif_Type" attribute) in the following order: CLIP-seq, PAR-clip, RIP-chip, RNAcompete, SELEX, yeast three-hybrid system. We then linked each motif thus determined to a single RBP in the CISBP-RNA database, following the same precedence rules.




# RBPDB (last updated in 2012)
# Some information on the file names (e.g. 1171_19561594.pwm: internal ID followed by PMID)
# http://rbpdb.ccbr.utoronto.ca/docs/RBPDB_README_v1.2.txt
# Downloads:
# http://rbpdb.ccbr.utoronto.ca/download.php
# I got the "Human database content - TDT format" from there, but I'm not using it.
# Only using:
# All PWMs
# http://rbpdb.ccbr.utoronto.ca/downloads/PWMDir.zip
# I took the matrix_list.txt file from there and filtered it to only include the human PWMs:
# Human PWMs and PFMs
# http://rbpdb.ccbr.utoronto.ca/downloads/matrices_human.zip
# Then, in matrix_list.txt, I updated a few gene symbols (those that don't appear in the 'symbol' column in table 'hgnc').
# Then I renamed the .pwm files to something useful (gene symbols):
# [ant-login7] ~/update/rbp_motifs/input/rbpdb/matrices_human/PWMDir >> cat matrix_list.txt | cut -f1,2 | perl -ne 'chomp; s/(\w+)/$1.pwm/g; $i = 1; ($from, $to) = split(/\t/); $to =~ s/(_\d+)?\.pwm$/_$i.pwm/; while (-s $to) { $i++; $to =~ s/(_\d+)?\.pwm$/_$i.pwm/; } ; system("mv -v $from $to");'
# Got PWMs for 29 RBPs
# Most of these also have simple motif strings listed in RBPDB_v1.3.1_experiments_human_2012-11-21.tbt.
# Many others have motifs, but a lot of them are really short and degenerate, so I don't think I should use them.
# The PWMs range between 4 and 10+ nt.
# [ant-login7] ~/update/rbp_motifs/input/rbpdb >> ../../masta2meme.pl all_PWMs_human.masta
run("Convert MASTA to MEME", "download_masta2meme.pl");




# ATTRACT
# https://attract.cnic.es/attract/static/ATtRACT.zip
# Turn it into a MEME file:
run("Convert ATTRACT to MEME", "download_attract2meme.pl");



# RBPMAP
# http://rbpmap.technion.ac.il/download.html#download
# http://rbpmap.technion.ac.il/source/RBPmap_1.1.tar.gz
run("Convert RBPmap to MEME", "download_rbpmap2meme.pl");



# Dominguez et al. (Chris Burge RNA Bind-N-Seq)
# Got this by email from Daniel Dominguez
# Turn PWMs into a MEME file:
run("Convert PWMs to MEME", "download_dominguez2meme.pl");



# show directory
run("ls", "ls -lah input");

done();
