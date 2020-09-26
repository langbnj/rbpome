#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$table = 'rbpome';

# run

run("Compare validation ratios for different score types (SUM_IS, AVG_IS, etc.)", "compare_scoretypes.pl $table");
run("Show output", "cat output-$table-compare-scoretypes.txt");

run("Compare validation ratios for different classes (control, high_interest, ben_screen_top, etc.) for SUM_IS", "compare_classes.pl $table SUM_IS 7.1");
# run("Show output", "cat output-$table-compare-classes.txt");

done();
