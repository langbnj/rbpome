#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run
run("Get validation ratios", "main.pl 0.05");
run("cat", "cat output-table.txt");

done();
