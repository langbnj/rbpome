# rbpome
Scripts used for our RNA-binding protein interactome (RBPome) study

## Contents

**analyses**: R scripts used to perform analyses and produce figures (mentioned in the folder names)

**include**: Perl utility functions and R MySQL connection functions

**mysql_tables**: MySQL table structures (CREATE statements)

**pipeline**: Perl scripts used to produce various figures (see "Note" within folders)

**scripts**: Additional utility Perl scripts

**update**: Perl scripts used to fill MySQL tables. "run.pl" is the main script in each folder

For the Perl scripts, first run "download.pl" to download datasets, followed by "run.pl", the main script in each folder.

## Note

Please do not hesitate to contact us for assistance with these scripts — they have grown historically as this analysis grew far bigger than originally envisioned and aren't polished, though hopefully the comments make it possible to follow them. Thank you!