#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download
# http://interactome.baderlab.org/download
# run("Download", "wget -O input/huri.tsv 'http://interactome.baderlab.org/data/HuRI.tsv'");
# run("Download", "wget -O input/hi-union.tsv 'http://interactome.baderlab.org/data/HI-union.tsv'");

# http://www.interactome-atlas.org/download
run("Download", "wget -O input/huri.tsv 'http://www.interactome-atlas.org/data/HuRI.tsv'");
run("Download", "wget -O input/hi-union.tsv 'http://www.interactome-atlas.org/data/HI-union.tsv'");

run("dos2unix", "dos2unix input/huri.tsv");
run("dos2unix", "dos2unix input/hi-union.tsv");





# show directory
run("ls", "ls -lah input");

done();
