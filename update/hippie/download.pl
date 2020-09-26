#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download
run("Download", "wget -O input/hippie_current.txt 'http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt'");





# show directory
run("ls", "ls -lah input");

done();
