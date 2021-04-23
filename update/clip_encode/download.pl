#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download

# Got files.txt from the ENCODE data matrix page (hit Download):
# https://www.encodeproject.org/matrix/?type=Experiment

# Currently using this one (eCLIP and GRCh38, and no "adrenal gland"):
# https://www.encodeproject.org/matrix/?type=Experiment&target.investigated_as=RNA+binding+protein&assay_slims=RNA+binding&assay_title=eCLIP&assembly=GRCh38&files.file_type=bed+narrowPeak&biosample_ontology.classification=cell+line
# >> hit Download >> files.txt

cd("input");
run("Download", "xargs -n 1 curl -O -L < ../files.txt");
run("Unpack", "gunzip -vf *.gz");
cd("..");





# show directory
run("ls", "ls -1 input");

done();
