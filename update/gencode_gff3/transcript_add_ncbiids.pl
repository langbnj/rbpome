#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'gencode_gff3_transcript';

our $usage = "$0";
args(0);

Query("UPDATE `$table` SET ncbiids=NULL");
state("Cleared 'ncbiid' field in table '$table'");





# start
state("Inserting ncbiidds from GENCODE metadata files into table '$table'...");
starttime();
open(DIR, "ls input/*.metadata.Transcript_supporting_feature.txt -1|");
while (<DIR>)
{
    chomp;

    # DEBUG
    # $_ = 'gencode.human.v27.transcripts.fa';
    # $_ = 'gencode.human.v27.pc_translations.fa';
    # $_ = 'gencode.mouse.vM16.transcripts.fa';
    # $_ = 'gencode.mouse.vM16.pc_translations.fa';

    /^input\/gencode\.(\w+)\.[\w\.\-]+\.metadata\.Transcript_supporting_feature\.txt$/ or die("Error: Couldn't match file name '$_'");
    $species = uc($1);
    startme(" >> $species >> $_", 1);

	# RNA sequences
    open(IN, $_) or die("Error: Couldn't open '$_'\n");
    while (<IN>)
    {
		chomp;
        @a = split(/\t/);
		
		# ENST00000618181.4	NM_152486.2	RefSeq_dna
		# ENST00000618181.4	I7FC29.1	Uniprot/SPTREMBL
		# ENST00000622503.4	NM_152486.2	RefSeq_dna
		
		die("Error: Expected 3 columns in line:\n\n$_\n\n") if (scalar(@a) != 3);
		
		$enstv = $a[0];
		$ncbiid = $a[1];
		$type = $a[2];

		# Skip everything that isn't RefSeq_NT
		next if ($type ne 'RefSeq_dna');
		
		# Verify ENSTV
		$enstv =~ /^(ENST|ENSMUST)\d+\.\d+$/ or die("Error: Couldn't parse enstv '$enstv'");
		
		# Get previous NCBIIDs, if there are multiple ones for this ENSTV
		$query = Query("SELECT ncbiids FROM `$table` WHERE enstv='$enstv' AND species='$species'");
		($ncbiids) = FetchOne($query);
		# Add new ID
		if (defined($ncbiids))
		{
			die if ($ncbiids eq '');
			$ncbiids .= "|$ncbiid";
		}
		else
		{
			$ncbiids = "$ncbiid";
		}
		$ncbiids = join('|', unique(split(/\|/, $ncbiids)));
		
		# Update table
		$q = "UPDATE `$table` SET ncbiids='$ncbiids' WHERE enstv='$enstv' AND species='$species'";
		$query = Query($q);
		die("Error: ".Numrows($query)." rows affected by query (should be 1):\n\n$q\n\n") if (Numrows($query) != 1);

        addme("total '$species' enstvs updated with ncbiid(s)", $enstv);
        addme("total '$species' ncbiids", $ncbiid);

        stepme(1000, 1);
    }
    close(IN);
    stopme(1);
}
nl();
stoptime();

showmesome(20);
nl();

Optimize($table);

done();
