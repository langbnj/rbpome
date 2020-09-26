#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download

run("Download", "wget -O input/hgnc_complete_set.txt 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt'");

# Mapping to mouse gene symbols
# (not actually used)
run("Download", "wget -O input/custom_hgnc_to_mgi.txt 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=family.name&col=gd_mgd_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'");

# Remove quotes (they're only around some fields, which is why the import script doesn't catch them)
run("cat", "cat input/hgnc_complete_set.txt | perl -ne 's/\\t\"/\\t/g; s/\"\\t/\\t/g; print' > input/hgnc_complete_set_noquotes.txt");





# show directory
run("ls", "ls -1 input");

done();
