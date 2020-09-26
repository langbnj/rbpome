#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run

# Import
run("Import", "~/scripts/import.pl input/hgnc_complete_set_noquotes.txt hgnc -overwrite -quiet");


# Delete withdrawn entries
state("Deleting withdrawn entries...");
$query = Query("DELETE FROM hgnc WHERE status!='Approved'");
$affected = Numrows($query);
state(commify($affected)." rows affected", 1);
done();


# Make hgnc_aliases table for easy lookup (and only retaining unambiguous aliases)
run("Make hgnc_aliases table for easy lookup", "hgnc_aliases.pl");
run("Remove ambiguous aliases from table 'hgnc_aliases'", "remove_ambiguous_aliases.pl");
run("Add symbols as 'aliases' to table 'hgnc_aliases' (so I can simply look up ENSGs via the 'alias' column alone)", "add_symbols_as_aliases.pl");

run("Check ENSG-symbol mapping", "check_ensgs.pl");

done();


# Stats
$query = Query("SELECT COUNT(DISTINCT symbol), COUNT(DISTINCT alias), COUNT(DISTINCT CONCAT(symbol, alias)) FROM hgnc_aliases");
($symbols, $aliases, $pairs) = FetchOne($query);

state("Unique symbols: ".commify($symbols)."\nUnique aliases: ".commify($aliases)."\nUnique pairs:   ".commify($pairs));

# Note: Each symbol has exactly one ENSG, but 3 ENSGs have multiple symbols.
# SELECT ensg, GROUP_CONCAT(DISTINCT symbol ORDER BY symbol), COUNT(DISTINCT symbol) AS c FROM hgnc_aliases GROUP BY ensg ORDER BY c DESC;
# SELECT symbol, GROUP_CONCAT(DISTINCT ensg ORDER BY ensg), COUNT(DISTINCT ensg) AS c FROM hgnc_aliases GROUP BY symbol ORDER BY c DESC;
# >> Can look up ENSGs fine via the symbol, but looking up symbols via ENSGs will give multiple symbols in 3 cases.

done();
