#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');

# $rel = "3.4.154";	# Nov 2017
$rel = "3.5.185";	# May 2020



# Download
run("Download", "wget -O input/BIOGRID-ALL-$rel.tab2.zip 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-$rel/BIOGRID-ALL-$rel.tab2.zip'");
run("Unpack", "unzip input/BIOGRID-ALL-$rel.tab2.zip -d input");
run("Rename", "mv input/BIOGRID-ALL-3.5.185.tab2.txt input/BIOGRID-ALL.tab2.txt");
run("Clean up", "rm -f input/BIOGRID-ALL-$rel.tab2.zip");





# show directory
run("ls", "ls input");

done();
