#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table (e.g. rbpome)]\n\nExample: $0 rbpome";
($table) = args(1);
#args(0);

# Test if table exists
# $query = Query("SELECT id FROM $table LIMIT 1");
# ($tmpid) = FetchOne($query);
die("Error: Table '$table' doesn't exist") if (!Exists($table));



# run

# p-values table
# run("Concatenate p-value tables into one", "cat output/output-pvalues-*.txt | head -n 1 > output-pvalues-table.txt; cat output/output-pvalues-*.txt | grep -vP 'mode\\thc\\tminlog2fold\\tminscore\\tmindist_threshold\\tresamples\\tvalues\\tp_value\\tresampling_p_value' | natsort >> output-pvalues-table.txt");
#DEBUG
# Some output-pvalues files currently have 11 columns, others have 12 (adding scalar(@fgv) as "values"). Here, I'm filtering to include only the 12-column tables in the final one.
# Later, rerun everything and remove the 12-filtering step to make sure I get an error if some files don't have 12 columns
run("Concatenate p-value tables into one", "cat output/output-pvalues-$table-*.txt | head -n 1 > output-pvalues-$table.txt; cat output/output-pvalues-$table-*.txt | perl -ne '\@a = split(/\t/); print if (scalar(\@a) == 12)' | grep -vP 'mode\\thc\\tminlog2fold\\tminscore\\tmindist_threshold\\tresamples\\tvalues\\tp_value\\tresampling_p_value' | natsort >> output-pvalues-$table.txt");
#END DEBUG

# txt output file
run("Concatenate txt output files into one", "cat output/output-txt-$table-*.txt | head -n 1 > output-txt-$table-tmp.txt; cat output/output-txt-$table-*.txt | grep -vP 'mode\\thc\\tminlog2fold\\tminscore\\tmindist_threshold\\tresamples\\tset\\trbp1\\trbp2\\tcobind' | natsort >> output-txt-$table-tmp.txt");
run("Collapse background_resample_1 etc. into one", "cat output-txt-$table-tmp.txt | perl -ne 's/background_resample_\\d+/background_resample/; print' > output-txt-$table.txt");
# run("Remove background_resample_1 etc. (and positives, too - keep only screen_hits) (much smaller table)", "cat output-txt-all-tmp.txt | grep '\tscreen_hits\t' > output-txt-all.txt");

# hits output file
run("Concatenate hits output files into one", "cat output/output-hits-$table-*.txt | head -n 1 > output-hits-$table-tmp.txt; cat output/output-hits-$table-*.txt | grep -vP 'table\\ttype\\tscoretype\\tmode\\thc\\tminlog2fold\\tminscore\\tmindist_threshold\\tresamples\\tset\\trbp1\\trbp2\\tensgv\\tmindist' | natsort >> output-hits-$table-tmp.txt");
run("Simply copy to final", "cat output-hits-$table-tmp.txt > output-hits-$table.txt");
# run("Collapse background_resample_1 etc. into one", "cat output-hits-table-tmp.txt | perl -ne 's/background_resample_\\d+/background_resample/; print' > output-hits-table.txt");
# run("Remove background_resample_1 etc. (and positives, too - keep only screen_hits) (much smaller table)", "cat output-hits-table-tmp.txt | head -n 1 > output-hits-table.txt; cat output-hits-table-tmp.txt | grep -iP '\\tscreen_hits\\t' >> output-hits-table.txt");



done();

