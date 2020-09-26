#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

$table = 'rbpome';

our $usage = "$0 [-addmeth]\n\n -addmeth: Add an initiator methionine at the start of the CDS if there isn't one";
# ($var) = args(1);
args(0);

$infile = "input/OrfSeqDNA.txt";
$tablefile = "output-table.txt";
$outfilep = "output-fasta-protein.txt";
$outfiled = "output-fasta-dna.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(TABLE, ">$tablefile") or die("\nError: Couldn't open '$tablefile'\n\n");
open(OUTP, ">$outfilep") or die("\nError: Couldn't open '$outfilep'\n\n");
open(OUTD, ">$outfiled") or die("\nError: Couldn't open '$outfiled'\n\n");

print TABLE "ProtSeq\tUniProtAccs\tUniProtNames\tUniProtSpecies\n";


# start

# Read species mapping for gene symbols
startme("Reading ORF DNA sequences from '$infile', translating them to amino acids, and checking if they match UniProt proteins");
$i = 0;
starttime();
while (<IN>)
{
	chomp;
	
	$i++;
	$seqname = "protein_$i";
	
	$accs = '';
	$names = '';
	$speciess = '';

	$dna = $_;
	
	if ($dna =~ /^[acgt]+$/)
	{
		addme("sequence processing: all-lower case sequences (converted)", $seqname);
		$dna = uc($dna);
	}
	elsif ($dna =~ /^[ACGT]+$/)
	{
		addme("sequence processing: all-upper case sequences (kept)", $seqname);
	}
	else
	{
		$dna =~ tr/acgt/ACGT/;
		
		if ($dna =~ /^[ACGT]+$/)
		{
			addme("sequence processing: there were mixed upper and lower case characters in sequence (fixed)", $seqname);
		}
		elsif ($dna =~ / /)
		{
			$dna =~ s/ //g;
			if ($dna =~ /^[ACGT]+$/)
			{
				addme("sequence processing: there were spaces in sequence (fixed)", $seqname);
			}
		}
		else
		{
			die("Error: Non-ACGT characters in sequence:\n\n$dna\n\n");
		}
	}
	
	# If there is no initiator methionine, add one
	if ($dna !~ /^ATG/)
	{
		if (switch('addmeth'))
		{
			# Add an initiator methionine artificially (apparently they left them out very often)
			$dna = 'ATG'.$dna;
		}
	}
	
	if ($dna !~ /^ATG/) 
	{
		addme("sequence processing: sequence doesn't start with ATG for sequence (skipped to the first ATG codon)", $seqname);
		$dna =~ s/^[ACGT]+?ATG/ATG/;
		
		
		# Do a test translation to check for internal stop codons
		$tmp_dna = $dna;
		if (length($tmp_dna) % 3 == 1)
		{
			$tmp_dna = substr($tmp_dna, 0, length($tmp_dna) - 1);
		}
		elsif (length($tmp_dna) % 3 == 2)
		{
			$tmp_dna = substr($tmp_dna, 0, length($tmp_dna) - 2);
		}
		die if (length($tmp_dna) % 3 != 0);
	
		$tmp_protein = translate($tmp_dna);
	
		$skip = 0;
		while ($tmp_protein =~ /[ACDEFGHIKLMNPQRSTVWY]\*[ACDEFGHIKLMNPQRSTVWY]/)
		{
			$skip++;
			addme("sequence processing: sequence doesn't start with ATG for sequence and the first codon produced a frameshift (skipped to ATG codon #$skip until there were no internal stop codons)", $seqname);
			$dna =~ s/^[ACGT]+?ATG/ATG/;


			# Do a test translation to check for internal stop codons
			$tmp_dna = $dna;
			if (length($tmp_dna) % 3 == 1)
			{
				$tmp_dna = substr($tmp_dna, 0, length($tmp_dna) - 1);
			}
			elsif (length($tmp_dna) % 3 == 2)
			{
				$tmp_dna = substr($tmp_dna, 0, length($tmp_dna) - 2);
			}
			die if (length($tmp_dna) % 3 != 0);
	
			$tmp_protein = translate($tmp_dna);
		}
	}
	
	if ($dna !~ /^ATG/)
	{
		die("Error: sequence still doesn't start with a start codon for sequence\n\n$dna\n\n");
	}
	
	if (length($dna) % 3 == 1)
	{
		addme("sequence processing: sequence length is not a multiple of 3 (removed one base)", $seqname);
		$dna = substr($dna, 0, length($dna) - 1);
	}
	elsif (length($dna) % 3 == 2)
	{
		addme("sequence processing: sequence length is not a multiple of 3 (removed two bases)", $seqname);
		$dna = substr($dna, 0, length($dna) - 2);
	}
	die if (length($dna) % 3 != 0);
	
	$seq = translate($dna);
	
	if ($seq =~ /\*$/)
	{
		addme("sequence processing: got a stop codon at the end (good, removed it)", $seqname);
		$seq =~ s/\*$//;
	}
	else
	{
		addme("sequence processing: no stop codon at the end (weird)", $seqname);
	}
	
	# foreach $species ('human', 'mouse')
	# {
		# $query = Query("SELECT name FROM uniseq WHERE type='UniProt' AND species='$species' AND seq='$seq'");
		$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq='$seq'");
		if (Numrows($query) == 1)
		{
			addme("uniseq matching: found a perfect match for seq", $seqname);
		}
		elsif (Numrows($query) > 1)
		{
			addme("uniseq matching: found multiple perfect matches for seq", $seqname);
		}
		elsif (Numrows($query) == 0)
		{
			# Try to add an initiator methione
			$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq='M$seq'");
			if (Numrows($query) > 0)
			{
				# die("Error: Found a match by dropping the initiator methionine for seq", $seqname);
				addme("uniseq matching: found perfect match(es) by adding an initiator methionine for seq", $seqname);
			}
			else
			{
				# Try to remove the initiator methione
				$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq='".substr($seq, 1)."'");
				if (Numrows($query) > 0)
				{
					# die("Error: Found a match by dropping the initiator methionine for seq", $seqname);
					addme("uniseq matching: found perfect match(es) by dropping the initiator methionine for seq", $seqname);
				}
				else
				{
					# Try LIKE query at the start (since they all start with a start codon)
					$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq LIKE '$seq\%'");
					if (Numrows($query) > 0)
					{
						# die("Error: Found a match by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
						addme("uniseq matching: found match(es) using LIKE ...% for seq", $seqname);
					}
					else
					{
						# Try LIKE query at the start (since they all start with a start codon), and drop the initiator methionine
						$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq LIKE '".substr($seq, 1)."\%'");
						if (Numrows($query) > 0)
						{
							# die("Error: Found a match by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
							addme("uniseq matching: found match(es) by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
						}
						else
						{
							# Try to find the sequence internally
							$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq LIKE '\%$seq\%'");
							if (Numrows($query) > 0)
							{
								# die("Error: Found an internal match for seq", $seqname);
								addme("uniseq matching: found internal match(es) using LIKE %...% for seq", $seqname);
							}
							else
							{
								# Try to find the sequence internally, dropping the initiator methionine
								$query = Query("SELECT name, species, '' AS acc FROM uniseq WHERE type='UniProt' AND species IN ('human', 'mouse') AND seq LIKE '\%".substr($seq, 1)."\%'");
								if (Numrows($query) > 0)
								{
									# die("Error: Found an internal match for seq", $seqname);
									addme("uniseq matching: found internal match(es) by dropping the initiator methionine using LIKE %...% for seq", $seqname);
								}
								else
								{


									# From here on: try uniiso instead of uniseq
									# Try to remove the initiator methione
									$query = Query("SELECT DISTINCT name, species, acc FROM uniiso WHERE species IN ('human', 'mouse') AND seq='".substr($seq, 1)."'");
									if (Numrows($query) > 0)
									{
										# die("Error: Found a match by dropping the initiator methionine for seq", $seqname);
										addme("uniiso matching: found perfect match(es) by dropping the initiator methionine for seq", $seqname);
									}
									else
									{
										# Try LIKE query at the start (since they all start with a start codon)
										$query = Query("SELECT DISTINCT name, species, acc FROM uniiso WHERE species IN ('human', 'mouse') AND seq LIKE '$seq\%'");
										if (Numrows($query) > 0)
										{
											# die("Error: Found a match by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
											addme("uniiso matching: found match(es) using LIKE ...% for seq", $seqname);
										}
										else
										{
											# Try LIKE query at the start (since they all start with a start codon), and drop the initiator methionine
											$query = Query("SELECT DISTINCT name, species, acc FROM uniiso WHERE species IN ('human', 'mouse') AND seq LIKE '".substr($seq, 1)."\%'");
											if (Numrows($query) > 0)
											{
												# die("Error: Found a match by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
												addme("uniiso matching: found match(es) by dropping the initiator methionine and using LIKE ...% for seq", $seqname);
											}
											else
											{
												# Try to find the sequence internally
												$query = Query("SELECT DISTINCT name, species, acc FROM uniiso WHERE species IN ('human', 'mouse') AND seq LIKE '\%$seq\%'");
												if (Numrows($query) > 0)
												{
													# die("Error: Found an internal match for seq", $seqname);
													addme("uniiso matching: found internal match(es) using LIKE %...% for seq", $seqname);
												}
												else
												{
													# Try to find the sequence internally, dropping the initiator methionine
													$query = Query("SELECT DISTINCT name, species, acc FROM uniiso WHERE species IN ('human', 'mouse') AND seq LIKE '\%".substr($seq, 1)."\%'");
													if (Numrows($query) > 0)
													{
														# die("Error: Found an internal match for seq", $seqname);
														addme("uniiso matching: found internal match(es) by dropping the initiator methionine using LIKE %...% for seq", $seqname);
													}
												}
											}
										}
									}								
								
								
								}
							}
						}
					}
				}
			}
		}
		if (Numrows($query) == 1)
		{
			addme("uniseq matching: eventually found one match for seq", $seqname);
		}
		elsif (Numrows($query) > 1)
		{
			addme("uniseq matching: eventually found multiple matches for seq", $seqname);
		}
		elsif (Numrows($query) == 0)
		{
			addme("uniseq matching: couldn't find a match in uniseq by any method for sequence", $seqname);
		}
		
		while (($name, $species, $acc) = Fetch($query))
		{
			if ($acc eq '')
			{
				# If acc isn't already set from uniiso (the match came from uniseq), get it here
				$accquery = Query("SELECT DISTINCT primary_acc FROM uniacc WHERE name='$name' AND species='$species'");
				# while (($acc) = Fetch($query))
				($acc) = FetchOne($accquery);
			}

			$names .= "|$name";
			$accs .= "|$acc";
			$speciess .= "|$species";
			
			print OUTP ">$acc|$name|$species\n".split60($seq)."\n";
			print OUTD ">$acc|$name|$species\n".split60($dna)."\n";
			
			addme("uniseq matching: successfully found an acc in uniseq/uniacc for sequence", $seqname);
		}
	# }
	
	# Strip initial spacers
	$names =~ s/^\|//;
	$accs =~ s/^\|//;
	$speciess =~ s/^\|//;

	# Reorder them
	$names = join("|", unique(split(/\|/, $names)));
	$accs = join("|", unique(split(/\|/, $accs)));
	$speciess = join("|", unique(split(/\|/, $speciess)));

	# Convert them to semicolons
	$names =~ s/\|/;/;
	$accs =~ s/\|/;/;
	$speciess =~ s/\|/;/;
	
	# Print table columns for copying into the large Excel table
	print TABLE "$seq\t$accs\t$names\t$speciess\n";

	stepme(1);
}
stopme();
stoptime();

showmeallsorted(1);

Optimize($table);

nl();
state("Wrote $i proteins to '$outfilep'", 1);
state("Wrote $i proteins to '$outfiled'", 1);
state("Wrote $i proteins to '$tablefile'", 1);
nl();

done();
