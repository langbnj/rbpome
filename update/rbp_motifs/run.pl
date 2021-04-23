#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

#DEBUG
# waitforalljobs();
#END DEBUG

# Load grep motifs
Clear('rbp_motif_list');
run("Load Dominguez et al. (Chris Burge RNA Bind-N-Seq) grep motifs into table rbp_motif_list", "load_motifs_dominguez.pl");
run("Load ATTRACT grep motifs into table rbp_motif_list", "load_motifs_attract.pl");
run("Load my custom grep motifs (e.g. from UniProt for YBX3) into table rbp_motif_list", "load_motifs_custom.pl");

# run("Remove motif import files", "rm -f tmp-rbp_motifs*.txt");



# Get individual motif files
@fimosources = ('attract', 'cisbp', 'dominguez', 'rbpdb', 'rbpmap', 'rnacompete');
foreach $source (@fimosources)
{
	run("Get individual motif files (per RBP)", "get_motifs.pl $source");
}


# Main analysis: Peak sequences
@types = ('eclip_encode_12');


# Search for motif PWMs
foreach $type (@types)
{
	run("Job", "~/scripts/qsub.sh job.pl $type")	# This currently is no_grep (and no background models)
	# run("Job", "~/scripts/qsub.sh job_grep_only.pl $type")
	# run("Job", "~/scripts/qsub.sh job_no_grep.pl $type")
}

waitforjobs();

run("Plot logos of all the input PWMs", "logos.pl");

done();
