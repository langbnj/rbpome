#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_transcript';

our $usage = "$0";
args(0);

Query("UPDATE `$table` SET seq=NULL, aaseq=NULL");
state("Cleared 'seq' and 'aaseq' fields in table '$table'");





# start
state("Inserting transcript sequences and translations from GENCODE FASTA files into table '$table'...");
starttime();
open(DIR, "ls input/*.transcripts.fa input/*.pc_translations.fa -1|");
while (<DIR>)
{
    chomp;

    # DEBUG
    # $_ = 'gencode.human.v27.transcripts.fa';
    # $_ = 'gencode.human.v27.pc_translations.fa';
    # $_ = 'gencode.mouse.vM16.transcripts.fa';
    # $_ = 'gencode.mouse.vM16.pc_translations.fa';

    /^input\/gencode\.(\w+)\.[\w\.\-]+\.(transcripts|pc_translations)\.fa$/ or die("Error: Couldn't match file name '$_'");
    $species = uc($1);
    $seqtype = $2;
    startme(" >> $species >> $_", 1);
	die("Error: seqtype is '$seqtype' in '$_'") if (($seqtype ne 'transcripts') and ($seqtype ne 'pc_translations'));

    # Read FASTA
	if ($seqtype eq 'transcripts')
	{
		# RNA sequences
	    open(IN, $_) or die("Error: Couldn't open '$_'\n");
		fastabreak();
	    while (<IN>)
	    {
	        ($title, $seq) = getfasta();
			
			# >ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|
			$title =~ /^([^\|]+)\|/ or die("Error: Couldn't parse title '$title'");
			$enstv = $1;
		
			$q = "UPDATE `$table` SET seq='$seq' WHERE enstv='$enstv' AND species='$species' AND seq IS NULL";
			$query = Query($q);
			die("Error: ".Numrows($query)." rows affected by query (should be 1):\n\n$q\n\n") if (Numrows($query) != 1);

	        addme("total '$species' enstvs updated with $seqtype sequence", $enstv);

	        stepme(10000, 1);
	    }
		normalbreak();
	    close(IN);
	}
	elsif ($seqtype eq 'pc_translations')
	{
		# Protein sequences
	    open(IN, $_) or die("Error: Couldn't open '$_'\n");
		fastabreak();
	    while (<IN>)
	    {
	        ($title, $seq) = getfasta();
			
			# >ENSP00000493376.1|ENST00000641515.1|ENSG00000186092.5|OTTHUMG00000001094.3|OTTHUMT00000003223.3|OR4F5-202|OR4F5|305
			$title =~ /^[^\|]+\|([^\|]+)\|/ or die("Error: Couldn't parse title '$title'");
			$enstv = $1;
		
			$q = "UPDATE `$table` SET aaseq='$seq' WHERE enstv='$enstv' AND species='$species' AND aaseq IS NULL";
			$query = Query($q);
			die("Error: ".Numrows($query)." rows affected by query (should be 1):\n\n$q\n\n") if (Numrows($query) != 1);

	        addme("total '$species' enstvs updated with $seqtype sequence", $enstv);

	        stepme(10000, 1);
	    }
		normalbreak();
	    close(IN);
	}
    stopme(1);
}
nl();
stoptime();

showmesome(20);
nl();

Optimize($table);

done();
