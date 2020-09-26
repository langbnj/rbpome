#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0";
($var) = args(1);
# args(0);

$infile = "input/MRK_List2.txt";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

startme("Starting");
starttime();
<IN>;	# Skip header
while (<IN>)
{
	chomp;
	
	# MGI Accession ID	Chr	cM Position	genome coordinate start	genome coordinate end	strand	Marker Symbol	Status	Marker Name	Marker Type	Feature Type	Marker Synonyms (pipe-separated)
	# MGI:1341858	5	  syntenic				03B03F	O	DNA segment, 03B03F (Research Genetics)	BAC/YAC end	BAC/YAC end
	# MGI:1341869	5	  syntenic				03B03R	O	DNA segment, 03B03R (Research Genetics)	BAC/YAC end	BAC/YAC end
	# MGI:1337005	11	  syntenic				03.MMHAP34FRA.seq	O	DNA segment, 03.MMHAP34FRA.seq	DNA Segment	DNA segment
	# MGI:1918911	7	  29.36	45567795	45575176	-	0610005C13Rik	O	RIKEN cDNA 0610005C13 gene	Gene	lncRNA gene
	# MGI:1923503	7	  syntenic	74818818	74853813	-	0610006L08Rik	O	RIKEN cDNA 0610006L08 gene	Gene	lncRNA gene
	
	
	stepme(100);
}
stopme();
stoptime();

showmeall(1);

done();
