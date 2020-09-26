#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run

run("Get problematic mouse-human symbol pairs from tables 'hgnc' and 'mgi' (where they're mapped to each other by HGNC, but the symbols don't match)", qq(~/scripts/query.pl "SELECT UPPER(m.marker_symbol) AS mouse_symbol, UPPER(h.symbol) AS human_symbol FROM hgnc h, mgi m WHERE h.mgd_id=m.mgi_accession_id AND UPPER(h.symbol)!=UPPER(m.marker_symbol) LIMIT 1000000;
" -h > tmp-problematic-symbols.txt));

run("Main", "main.pl");

# This was intended for checking if deduplicate.pl in update-rbpome changes ortholog group assignment (it doesn't).
# Outcome:
# Before ~/update/rbpome/deduplicate.pl: No clashes between different alns (different orthology groups). Already all good as far as homologene is concerned.
# After: All good, too!

done();



# Before log:
#
# --------------------------------------------------------------------------------
# > Get problematic mouse-human symbol pairs from tables 'hgnc' and 'mgi' (where they're mapped to each other by HGNC, but the symbols don't match) (~/scripts/query.pl "SELECT UPPER(m.marker_symbol) AS mouse_symbol, UPPER(h.symbol) AS human_symbol FROM hgnc h, mgi m WHERE h.mgd_id=m.mgi_accession_id AND UPPER(h.symbol)!=UPPER(m.marker_symbol) LIMIT 1000000;
# " -h > tmp-problematic-symbols.txt)
# --------------------------------------------------------------------------------
#
#
# --------------------------------------------------------------------------------
# > Main (main.pl)
# --------------------------------------------------------------------------------
#
#
# Reading gene symbol alias mapping from '../../update/rbpome/input-aliases.txt':
# 0..10..20..30..40..50..60..65
#
#
# Creating hash of symbols that are the result of aliasing (aliases):
# 0..10..20..30..40..50..60..65
#
#
# Checking how my manual gene symbol alias assignment affected the matchup between human, mouse and DANRE gene symbols in table 'rbpome':
# 0..1,000..2,000..3,000..4,000..5,000..6,000..7,000..8,000..9,000..10,000..11,000..12,000..13,000..14,000..15,000..16,000..17,000..18,000..19,000..20,000..21,000..22,000..23,000..24,000..25,000..26,000..27,000..28,000..29,000..30,000..31,000..32,000..33,000..34,000..34,342
#
#
# Got 1,024 unique species|symbols
#
#
# Cycle through symbols, check if they fall into the same homologene group across species, and check if the alias assignment caused a difference for them:
# 0..100..200..300..400..500..600..700..800..900..1,000..1,024
#
# Time elapsed:	3.270 sec
#
# total species:   ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 2 (e.g. HUMAN)
# total original species|symbols:   ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 65 (e.g. HUMAN|ATP5B)
# total symbol species|aliases:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 64 (e.g. HUMAN|AGO1)
# total original symbols:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 65 (e.g. ATP5A1)
# total symbol aliases:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 64 (e.g. AGO1)
# total rbpome symbols:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 1,024 (e.g. A1CF)
# total rbpome species|symbols:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 1,024 (e.g. DANRE|DDX39AB)
# OK: no alns in table homologene for symbol, but we only have it for one species in table 'rbpome':  · 60 (e.g. APC_CTER)
# GOOD: exactly one aln for 1 species for symbol:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 964 (e.g. A1CF)
#
# Done!
#
#
# Done!
#
#





# After log:


--------------------------------------------------------------------------------
> Get problematic mouse-human symbol pairs from tables 'hgnc' and 'mgi' (where they're mapped to each other by HGNC, but the symbols don't match) (~/scripts/query.pl "SELECT UPPER(m.marker_symbol) AS mouse_symbol, UPPER(h.symbol) AS human_symbol FROM hgnc h, mgi m WHERE h.mgd_id=m.mgi_accession_id AND UPPER(h.symbol)!=UPPER(m.marker_symbol) LIMIT 1000000;
" -h > tmp-problematic-symbols.txt)
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
> Main (main.pl)
--------------------------------------------------------------------------------


Reading gene symbol alias mapping from '../../update/rbpome/input-aliases.txt':
0..10..20..30..40..50..60..65


Creating hash of symbols that are the result of aliasing (aliases):
0..10..20..30..40..50..60..65


Checking how my manual gene symbol alias assignment affected the matchup between human, mouse and DANRE gene symbols in table 'rbpome':
0..1,000..2,000..3,000..4,000..5,000..6,000..7,000..8,000..9,000..10,000..11,000..12,000..13,000..14,000..15,000..16,000..17,000..18,000..19,000..20,000..21,000..22,000..23,000..24,000..25,000..26,000..27,000..28,000..29,000..30,000..31,000..32,000..33,000..34,000..34,079


Got 1,020 unique species|symbols


Cycle through symbols, check if they fall into the same homologene group across species, and check if the alias assignment caused a difference for them:
0..100..200..300..400..500..600..700..800..900..1,000..1,020

Time elapsed:	2.658 sec

total species:   ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 2 (e.g. HUMAN)
total original species|symbols:   ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 65 (e.g. HUMAN|ATP5B)
total symbol species|aliases:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 64 (e.g. HUMAN|AGO1)
total original symbols:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 65 (e.g. ATP5A1)
total symbol aliases:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 64 (e.g. AGO1)
total rbpome symbols:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 1,017 (e.g. A1CF)
total rbpome species|symbols:  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 1,020 (e.g. DANRE|DDX39AB)
OK: no alns in table homologene for symbol, but we only have it for one species in table 'rbpome':  · 71 (e.g. APC_CTER)
GOOD: exactly one aln for 1 species for symbol:    ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 944 (e.g. A1CF)
GOOD: exactly one aln for 2 species for symbol:   ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  · 2 (e.g. AGO2)

Done!


Done!



