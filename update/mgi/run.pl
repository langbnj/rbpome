#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run
run("Import", "~/scripts/import.pl input/MRK_List2.txt mgi -overwrite -quiet");

# Make mgi_aliases table for easy lookup (and only retaining unambiguous aliases)
run("Remove repeated symbols from table 'mgi'", "remove_repeated_symbols.pl");
run("Make mgi_aliases table for easy lookup", "mgi_aliases.pl");
run("Remove ambiguous aliases from table 'mgi_aliases'", "remove_ambiguous_aliases.pl");
run("Add symbols as 'aliases' to table 'mgi_aliases' (so I can simply look up symbols via the 'alias' column alone)", "add_symbols_as_aliases.pl");

done();


# Stats
$query = Query("SELECT COUNT(DISTINCT symbol), COUNT(DISTINCT alias), COUNT(DISTINCT CONCAT(symbol, alias)) FROM mgi_aliases");
($symbols, $aliases, $pairs) = FetchOne($query);

state("Unique symbols: ".commify($symbols)."\nUnique aliases: ".commify($aliases)."\nUnique pairs:   ".commify($pairs));

# Note: Each symbol has exactly one ENSG, but 3 ENSGs have multiple symbols.
# SELECT ensg, GROUP_CONCAT(DISTINCT symbol ORDER BY symbol), COUNT(DISTINCT symbol) AS c FROM mgi_aliases GROUP BY ensg ORDER BY c DESC;
# SELECT symbol, GROUP_CONCAT(DISTINCT ensg ORDER BY ensg), COUNT(DISTINCT ensg) AS c FROM mgi_aliases GROUP BY symbol ORDER BY c DESC;
# >> Can look up ENSGs fine via the symbol, but looking up symbols via ENSGs will give multiple symbols in 3 cases.


# Get a couple of problematic symbols by hand:
# cat ~/update/mgi/input/HOM_AllOrganism.tsv | grep -iP "\thuman\t9606\t(ATP5F1B|ATP5F1C|CFAP20|ESS2|H1-2|H1-4|H2AX|KARS1|MMTAG2|MTREX|QARS1|RARS1|SARS1|SGO2|TUT4|ZRSR2P1)\t"|cut -f4
# cat ~/update/mgi/input/HOM_AllOrganism.tsv | grep -iP "\tmouse, laboratory\t10090\t(ATP5F1A|Ddx39a|Dynlt1|NCBP3|NOLC1|Nop53|PNISR|PUM3|RACK1|Ramac|RRP1|RTRAF|Snu13|SPOUT1|Srp54|Txn|Yars1|ZRANB2)\t"|cut -f4
# cat ~/update/mgi/input/HOM_AllOrganism.tsv | grep -iP "\tzebrafish\t7955\thuman\t9606\tDDX39A\t"|cut -f4

done();
