#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

Clear('huri');

# run
# huri & hi-union
# Jaccard index: 0.821
# huri (52,547) 100.00% in intersection
# hi-union (63,994) 82.11% in intersection
# >> huri is a subset of hi-union. hi-union adds in ~11500 more interactions from other databases.
run("Main", "main.pl huri");

# Skipping HI-Union now
# run("Main", "main.pl hi-union");

done();
