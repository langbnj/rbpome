#!/usr/bin/perl

# turn off buffering on stdout
$| = 1;

# prettier progress numbers, but gigantic log files. disable when logging.
our $awesomeprogress = 0;

# Use an EXTREME amount of nodes (2/3 of the free ones) (regardless of time of day). Extreme!
our $extrememode = 1;

# Use a POLITE amount of nodes (1/3 of the free ones) (regardless of time of day). Polite!
our $politemode = 0;

# Use .instakill file as workaround since bkill is so slow on the LSF cluster
our $instakill_enabled = 1;

# Features
# use feature qw(say);		# Would only work within this scope.

# Install modules
# sudo apt-get install libmysqlclient-dev # For DBD::mysql's 'mysql_config' test
# cpan
# install Bio::Phylo Statistics::Descriptive Statistics::R Text::CSV List::Compare Sort::Key::Natural Carp::Always Math::Round Statistics::Robust::Scale List::Vectorize Statistics::Multtest Statistics::RankCorrelation DBD::mysql

# Modules
use POSIX qw(ceil floor);
use Data::Dumper qw(Dumper);
use Digest::MD5;
use Time::HiRes;
use Statistics::Descriptive;
use Statistics::R;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Set::Scalar;
use Bio::Phylo::IO;
use Text::CSV;
use List::Compare;
# Natural sorting
# use Sort::Naturally qw(nsort);
# use Sort::Key::Natural qw( natkeysort );
# *nsort = \&natkeysort;
# sub nsort {return natkeysort(@_);}
use Sort::Key::Natural qw(natkeysort); sub nsort { return natkeysort { lc $_ } @_; }
# Crash traceback
use Carp::Always;
# use Carp::Always::Dump;
use autodie;
use Term::ReadKey;
use Statistics::Robust::Scale;
use Statistics::Multtest;
use Scalar::Util qw(looks_like_number);
use Statistics::RankCorrelation;
use List::Vectorize qw (sample);


#
# Exit if /users/gt/blang/.instakill exists (bkill alternative)
#
# Also in run() and state()
#

sub instakill
{
	our $instakill_enabled;
	# If this is enabled:
	if ($instakill_enabled == 1)
	{
		if (-e '/users/gt/blang/.instakill')
		{
			print("\n==================================================================================\n"
			."              ****  WARNING: EXITING BECAUSE .INSTAKILL FOUND  ****               "
			."\n==================================================================================\n\n");
			warn("\n==================================================================================\n"
			."              ****  WARNING: EXITING BECAUSE .INSTAKILL FOUND  ****               "
			."\n==================================================================================\n\n");
			die("\n");
		}
	}
}
# Reset if running from one of the mainframes
if (qtype() ne 'node')
# if (qtype() eq 'lsf')
{
  system("rm -f /users/gt/blang/.instakill");
}
else
{
  instakill();
}

#
# Format number with commas
#
# Example:
# my $number = 987654321;
# print commify($number);
#

sub commify
{
    my ($s) = @_;

    if (defined($s))
    {
        $s = reverse($s);
        $s =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
        return scalar reverse $s;
    }
    else
    {
        return '';
    }
}


# Free disk space functions

# Example:
#
# print "Got .getspace()." free space!";		# Returns free space in current working directory in GB (with up to 2 trailing decimals)
# minspace(2);									# will exit if there is less than 2 GB free space in the current working directory
#

sub minspace
{
	my ($minfree) = @_;

	if ($minfree !~ /^\d+$/)
	{
		die("Error: minspace() requires just a number (in gigabytes)");
	}

	$free = readpipe(q(df -H . | perl -ne '@a = split(/ +/); print $a[3] if ($a[3] ne "Avail");'));
	$free =~ /^([\d\.]+)(\w)$/ or die("Error: Couldn't parse free space '$free'");
	if (($1 <= $minfree) or ($2 ne 'G'))
	{
		$pwd = `pwd`;
		if ($pwd =~ /Volumes\/([^\/]+)\//)
		{
			$pwd = $1;
		}
		else
		{
			$pwd = "Macintosh HD";
		}
		die("\nError: Not enough free space (only $free available on Volume '$pwd', need at least ".$minfree."G)\n\n");
	}
}

sub getspace
{
	# in gigabytes

	my $free;

	$free = readpipe(q(df -k . | perl -ne '@a = split(/ +/); print $a[3] if ($a[3] ne "Available");'));

	$free =~ /^\d+$/ or die("Error: Couldn't parse free space '$free'");

	# $free /= 1024;
	# $free /= 1024;

	$free *= 1024;
	$free /= 1000;

	$free /= 1000;
	$free /= 1000;

	$free =~ /\d+\.(\d+)/ or die("Error: Couldn't parse free space '$free'");
	my $trailing = length($1);

	if ($trailing > 3)
	{
		$free = sprintf("%.3f", $free);
	}

	return $free;
}


# Break key functions

# Example:
#
# while (...)
# {
# 	last if break('x');
# }
# breakoff();					# << don't forget! This returns the terminal to normal (i.e. handy visible text)
#

sub break
{
	my ($break) = @_;

	my $key;
	ReadMode 3;
	if (defined($key = ReadKey(-1)))
	{
		if ($key eq $break)
		{
			print "\nExiting: Break key '$break' encountered!\n\n";

			ReadMode 0;    # Reset tty mode before exiting
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}

sub breakoff
{
	ReadMode 0;    # Reset tty mode before exiting
}


# Print out a complex variable
sub show
{
    print Dumper(@_);
}

# Print out a complex variable on STDERR
sub showe
{
    warn(Dumper(@_));
}

# # Print out a simple variable
# sub report
# {
# 	my ($var) = @_;
# 	$$var;
#
# 	state(uc($var)."=['$$var']");
# }

sub Log
{
	my ($value, @a) = @_;

	Query("INSERT INTO `output` SET name='".join("|", @a)."', pipeline='".locale()."', value='$value'");
}

sub LogRetrieve
{
	my (@a) = @_;

	$query = Query("SELECT value FROM `output` WHERE name='".join("|", @a)."'");

	return FetchOne($query);
}

sub locale
{
	return `/users/gt/blang/scripts/locale.sh`;
}

sub round
{
	# Standard decimal notation
	my ($num, $len) = @_;

	die("Error: Length not defined") if (!defined($len));

	my $res = sprintf("%.".$len."f", $num);

	# remove trailing zeroes
	# just use sciround instead?
	# $res =~ s/0+$//;
	# $res =~ s/\.$//;

	return $res;
}

sub sciround
{
	# Flexible scientific notation. Note that 'len' here specifies the number of digits before and after the period/exponent E.
	my ($num, $len) = @_;

	die("Error: Length not defined") if (!defined($len));

	my $res = sprintf("%.".$len."g", $num);

	return $res;
}

sub base
# return file basename
{
	my ($s) = @_;

	$s =~ /^([^\.]+)\./;

	return $1;
}

sub md5
{
	my ($s) = @_;

	$res = Digest::MD5::md5($s);

	return $res;
}

sub truefalse
# translate 1/0 to true/false
{
	my ($i) = @_;

	$res = 'n/a';

	if ($i eq '1') { $res = 'TRUE'; }
	if ($i eq '0') { $res = 'FALSE'; }

	return $res;
}


sub comp
# compare two numbers and give +, 0 or - as a result.
# results are trends. if two is bigger than one, >> +
{
	my ($one, $two) = @_;

	if ($two > $one) { return '+'; }
	elsif ($two < $one) { return '-'; }
	else { return '0'; }
}

sub comptext
# compare two numbers and give "up", "equal" or "down" as a result.
# results are trends. if two is bigger than one, >> up
{
	my ($one, $two) = @_;

	if ($two > $one) { return 'up'; }
	elsif ($two < $one) { return 'down'; }
	else { return 'equal'; }
}

sub strident
# returns the number of positions that are identical between two strings
{
	my ($one, $two) = @_;

	die("Error: Length of '$one' and '$two' is not equal") if (length($one) != length($two));

	my $id = 0;
	for ($i = 0; $i < length($one); $i++)
	{
		if (substr($one, $i, 1) eq substr($two, $i, 1))
		{
			$id++;
		}
	}

	return $id;
}

sub stridentnoblank
# returns the number of positions that are identical between two strings, but ignores mismatches involving whitespace for the count
# (i.e. if a position in one of the strings is a space, the position will be reported as an identical match)
{
	my ($one, $two) = @_;

	die("Error: Length of '$one' and '$two' is not equal") if (length($one) != length($two));

	my $id = 0;
	for ($i = 0; $i < length($one); $i++)
	{
		if ((substr($one, $i, 1) eq substr($two, $i, 1)) or (substr($one, $i, 1) =~ /^\s$/) or (substr($two, $i, 1) =~ /^\s$/))
		{
			$id++;
		}
	}

	return $id;
}

sub consensus
# takes an array of strings
# returns a consensus string
# mismatches will be spaces (' ')
{
	my (@a) = @_;

	my $len = length($a[0]);
	foreach my $s (@a)
	{
		die("Error: Length mismatch") if (length($s) != $len);
	}

	my $con = '';

	$i = 0;
	while ($i < $len)
	{
		my $c = substr($a[0], $i, 1);
		foreach $s (@a)
		{
			if (substr($s, $i, 1) ne $c)
			{
				$c = ' ';
			}
		}

		$con .= $c;

		$i++;
	}

	return $con;
}

# Checks if string consists only of the 20 standard amino acids (and at least one)
sub aa
{
	my ($s) = @_;

	if ($s =~ /^[ACDEFGHIKLMNPQRSTVWY]+$/) { return 1; }
	else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + U
sub aau
{
	my ($s) = @_;

	if ($s =~ /^[ACDEFGHIKLMNPQRSTUVWY]+$/) { return 1; }
	else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + X
sub aax
{
    my ($s) = @_;

    if ($s =~ /^[ACDEFGHIKLMNPQRSTVWYX]+$/) { return 1; }
    else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + X + U
sub aaxu
{
    my ($s) = @_;

    if ($s =~ /^[ACDEFGHIKLMNPQRSTVWYXU]+$/) { return 1; }
    else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + X + U + B + Z
sub aaxubz
{
    my ($s) = @_;

    if ($s =~ /^[ACDEFGHIKLMNPQRSTVWYXUBZ]+$/) { return 1; }
    else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + X + U
sub aaux
{
    my ($s) = @_;

    if ($s =~ /^[ACDEFGHIKLMNPQRSTVWYXU]+$/) { return 1; }
    else { return 0; }
}

# Checks if string consists only of the 20 standard amino acids (and at least one) + X + *
sub aaxstar
{
    my ($s) = @_;

    if ($s =~ /^[ACDEFGHIKLMNPQRSTVWYX\*]+$/) { return 1; }
    else { return 0; }
}

# Returns non-standard amino acids, if any
sub nonaa
{
	my ($s) = @_;

	while ($s =~ /([^ACDEFGHIKLMNPQRSTVWY])/g)
	{
		$res .= $1;
	}

	return $res;
}

# Returns non-standard amino acids - X, if any
sub nonaax
{
	my ($s) = @_;

	while ($s =~ /([^ACDEFGHIKLMNPQRSTVWYX])/g)
	{
		$res .= $1;
	}

	return $res;
}

# Checks if string consists only of the 4 DNA bases (and at least one)
sub dna
{
	my ($s) = @_;

	if ($s =~ /^[ACGT]+$/) { return 1; }
	else { return 0; }
}

# Checks if string consists only of the 4 DNA bases + N (and at least one)
sub dnan
{
	my ($s) = @_;

	if ($s =~ /^[ACGTN]+$/) { return 1; }
	else { return 0; }
}

# Checks if string consists only of the 4 RNA bases (and at least one)
sub rna
{
	my ($s) = @_;

	if ($s =~ /^[ACGU]+$/) { return 1; }
	else { return 0; }
}

# Checks if string consists only of the 4 RNA bases + N (and at least one)
sub rnan
{
	my ($s) = @_;

	if ($s =~ /^[ACGUN]+$/) { return 1; }
	else { return 0; }
}

# Checks if it's a yeast ordered locus name
sub orderedlocus
{
	my ($s) = @_;
	
	# This is a regular expression I built by hand for RNAct. It works for "YMR193W", "Q0010", "tT(UGU)G2", "15S_rRNA", etc. (all yeast ENSTs)
	# This returns 0, i.e. they all match: SELECT * FROM ensembl WHERE species='yeast' AND enst NOT REGEXP '^(T|SN)?[A-Z0-9\-\(\)_]{2,9}([A-Z]|_RRNA)?$';
	$s = uc($s);
	if ($s =~ /^(T|SN)?[A-Z0-9\-\(\)_]{2,9}([A-Z]|_RRNA)?$/) { return 1; }
	else { return 0; }
}

# Checks if it's a UniProt accession
sub acc
{
	my ($s) = @_;
	
	# This is a regular expression I got straight from the UniProt help section.
	# https://www.uniprot.org/help/accession_numbers
	$s = uc($s);
	if ($s =~ /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$/) { return 1; }
	else { return 0; }
}

# Split CSV string into an array of fields
sub splitcsv
{
	instakill();
	my ($s) = @_;

	# state("STRING '$s'");

	my $csv = Text::CSV->new ( { binary => 1 } ) or die("Error: Couldn't use module Text::CSV: ".Text::CSV->error_diag());

	$status = $csv->parse($s);
	if ($status != 1) { die("Error: Something went wrong with Text::CSV while parsing '$s'"); }
	# state("STATUS '$status'");

	my @a = $csv->fields();

	# Remove empty fields
	# @a = grep(!/^$/, @a);

	# Remove trailing empty fields at end of array
	# @a = split(/\0/, join("\0", @a));


	# state("\@A '".join("', '", @a)."'");

	return @a;
}

sub removetrailingfields
{
	my (@a) = @_;

	@a = split(/\0/, join("\0", @a));

	return @a;
}

sub numberduplicates
# Add "|1", "|2" etc. to duplicates in the given array of strings
# Only works with ADJACENT duplicates, i.e. they must be next to each other. Sort string first if in doubt.
# Useful for e.g. microarray expression data with sample replicates of the same name
{
	my (@a) = @_;

	$i = 0;
	while ($i < scalar(@a))
	{
		$i2 = 1;
		if (!exists($a[$i]) or !exists($a[$i + $i2]))
		{
			$i++;
			next;
		}

		if ($a[$i + $i2] eq $a[$i])
		{
			while ((($i + $i2) < scalar(@a)) and ($a[$i + $i2] eq $a[$i]))
			{
				$a[$i + $i2] .= "|".($i2 + 1);
				$i2++;
			}
			$a[$i] .= "|1";
		}

		$i++;
	}

	return @a;
}

# Phylogenetics functions
#
# Examples:
#
# $happytree = fixtree($inputtree)				# Just reads a tree, collapses unbranched internal nodes and returns it (fixes TreeBeST trees)
#
# $happytree = fixtree_compara($inputtree)		# Reads a tree, collapses unbranched internal nodes and removes _taxonids from node names (fixes TreeBeST trees)
#
# $comparatree = comparaprune($inputtree)		# Prunes off all terminal nodes that aren't Compara species & collapses unbranched internal nodes
# $maffttree = mafftprune($inputtree)			# Prunes off all terminal nodes that aren't numbers         & collapses unbranched internal nodes
#
# showbranchlengths($happytree)					# Lists all leaves with the branch lengths leading to them from the root
#
# @leaves = returnterminals($happytree)			# Return all leaves as an array
#

sub fixtree
{
	my ($treestring) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	my $old = '';
	while ($old ne $tree->to_newick)
	{
		$old = $tree->to_newick;

		# Collapse unbranched internal nodes
		$tree->remove_unbranched_internals;
	}

	return $tree->to_newick;
}

sub fixtree_compara
{
	my ($treestring) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	my $old = '';
	while ($old ne $tree->to_newick)
	{
		$old = $tree->to_newick;

		foreach my $node (@{$tree->get_root->get_terminals()})
		{
			my $name = $node->get_name;

			$name =~ s/_\d+$//;

			$node->set_name($name);
		}

		# Collapse unbranched internal nodes
		$tree->remove_unbranched_internals;
	}

	return $tree->to_newick;
}

sub comparaprune
{
	my ($treestring, $silent) = @_;

	if (!defined($silent))
	{
		my $silent = 0;
	}

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	startme("Pruning all non-Compara species from tree") if (!$silent);
	my $old = '';
	while ($old ne $tree->to_newick)
	{
		$old = $tree->to_newick;

		foreach my $node (@{$tree->get_root->get_terminals()})
		{
			# print " >> ".ref($node)." >> ".$node->get_name;

			my $display = $node->get_name;
			$display =~ s/_/ /;
			my $query = Query("SELECT species FROM comparaspecies WHERE display='$display'");
			die if (Numrows($query) > 1);
			if (Numrows($query) == 1)
			{
				($species) = FetchOne($query);
				$node->set_name($species);
			}
			elsif (Numrows($query) == 0)
			{
				addme("comparaprune: node id not a comparaspecies for species", $display) if (!$silent);
				$node->get_parent->prune_child($node);
			}
		}
		stepme(10) if (!$silent);

		# Collapse unbranched internal nodes
		$tree->remove_unbranched_internals;
	}
	stopme() if (!$silent);
	showme("comparaprune: node id not a comparaspecies for species") if (!$silent);
	clearme("comparaprune: node id not a comparaspecies for species") if (!$silent);

	return $tree->to_newick;
}

sub mafftprune
{
	my ($treestring) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	my $old = '';
	while ($old ne $tree->to_newick)
	{
		$old = $tree->to_newick;
		foreach my $node (@{$tree->get_root->get_terminals()})
		{
			if ($node->get_name !~ /^\d+$/)
			{
				if (defined($node->get_parent))
				{
					$node->get_parent->prune_child($node);
				}
			}
		}

		# Collapse unbranched internal nodes
		# $tree->remove_unbranched_internals;		# This function has problems with the modified tree here. They go away once it is converted back to Newick (see below).
		foreach my $node (@{$tree->get_root->get_internals()})
		{
			if (scalar(@{$node->get_children} == 1))
			{
				$node->collapse();
			}
		}
	}

	# fixtree: necessary to fix root
	$tree = $tree->to_newick;
	$tree = fixtree($tree);

	return $tree;
}

sub showbranchlengths
{
	my ($treestring, $silent) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	foreach my $node (@{$tree->get_root->get_terminals()})
	{
		print " >> ".$node->get_name." >> ".$node->calc_path_to_root()."\n";
	}
}

sub showterminaldistances
{
	my ($treestring, $silent) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	foreach my $node1 (@{$tree->get_root->get_terminals()})
	{
		foreach my $node2 (@{$tree->get_root->get_terminals()})
		{
			next if ($node2 == $node1);
			print " >> ".$node1->get_name." >> ".$node2->get_name." >> ".$node1->calc_patristic_distance($node2)."\n";
		}
	}
}

sub returnterminals
{
	my ($treestring) = @_;

	my $tree = Bio::Phylo::IO->parse(-string => $treestring, -format => 'newick')->first;

	my @a = ();
	foreach my $node (@{$tree->get_root->get_terminals()})
	{
		push(@a, $node->get_name);
	}

	return unique(@a);
}

# Statistics functions

# Example:
#
# $mean = mean(@a);
# $average = avg(@a);		# mean() / avg() = same function
#
# $variance = var(@a);
# $stdev = stdev(@a);
#
# $TissueSpecificityIndex = tissuespecificity(@a);
#
# characterise(@a);			# Prints a characterisation of an array of numbers, using most of the functions above
#

sub mean
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> mean();
}

sub avg
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> mean();
}

sub var
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> variance();
}

sub sd
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> standard_deviation();
}

sub std
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> standard_deviation();
}

sub stdev
{
	my $stat = Statistics::Descriptive::Full -> new(); $stat -> add_data(@_);
	return $stat -> standard_deviation();
}

sub mad
{
	# Median absolute deviation
	return Statistics::Robust::Scale::MAD(\@_);
}

sub median
{
	my @values = sort { $a <=> $b } (@_);
	if (scalar(@values) % 2)
	{
		return $values[int(scalar(@values) / 2)];
	}
	else
	{
		return ($values[scalar(@values) / 2] + $values[(scalar(@values) / 2) - 1]) / 2;
	}
}

sub characterise
{
	my @a = @_;

	if (scalar(@a) > 0)
	{
		print "NUM ".scalar(@a)."\t"."MIN ".sprintf("%.2f", min(@a))."\t"."AVG ".sprintf("%.2f", avg(@a))."\t"."MED ".sprintf("%.2f", median(@a))."\t"."MAX ".sprintf("%.2f", max(@a))."\t"."STD ".sprintf("%.2f", stdev(@a))."\n";
	}
	else
	{
		print "NUM N/A\tMIN N/A\tAVG N/A\tMED N/A\tMAX N/A\tSTD N/A\n";
	}
}

sub tissuespecificity
{
	# tissue specificity coefficient

	my @a = @_;

	my $sum = 0;

	foreach (@a)
	{
		$sum += (1 - ($_ / max(@a)));
	}

	return ($sum / (scalar(@a) - 1));
}

# sub yekutieli
sub fdr
{
    # Benjamini-Hochberg-Yekutieli correction (FDR with positive dependency)
    my $adj = Statistics::Multtest::BY(\@_);
	my @adj = @$adj;

	return @adj;
}

# FASTA functions

#
# $flippedseq = strandflip($dnaseq);	# Flips DNA strand
# $aaseq = translate($dnaseq);			# Translates bases into amino acids!
# $aaseq = translate_loosely($dnaseq);	# Translates bases into amino acids! Tolerant to gaps and purine/pyrimidine/wobble X symbols.
# $aaseq = translate_loosely_mitochondrial($dnaseq);	# Translates bases into amino acids! Tolerant to gaps and purine/pyrimidine/wobble X symbols. Uses mitochondrial genetic code.
#
# Example:
#
# fastabreak();
# while (<IN>)
# {
#	($title, $seq) = getfasta();
#
#	[...]
#
#	print OUT ">$title\n".split60($seq)."\n";
# }
# normalbreak();
#
#
#

sub strandflip
{
    my ($in) = @_;
    my $out = '';

    die("Error while strand-flipping: Sequence contains characters other than A, C, G, T: '$in'") if ($in !~ /^[ACGT]+$/);

    $out = reverse($in);
    $out =~ tr/ACGT/TGCA/;

    return $out;
}

sub strandflipu
{
	# For RNA
    my ($in) = @_;
    my $out = '';

    die("Error while strand-flipping: Sequence contains characters other than A, C, G, U: '$in'") if ($in !~ /^[ACGU]+$/);

    $out = reverse($in);
    $out =~ tr/ACGU/UGCA/;

    return $out;
}

sub strandflip_loosely
{
    my ($in) = @_;
    my $out = '';

    die("Error while strand-flipping: Sequence contains characters other than A, C, G, T, N, -: '$in'") if ($in !~ /^[ACGTN-]+$/);

    $out = reverse($in);
    $out =~ tr/ACGT/TGCA/;

    return $out;
}

sub strandflip_looselyu
{
    my ($in) = @_;
    my $out = '';

    die("Error while strand-flipping: Sequence contains characters other than A, C, G, U, N, -: '$in'") if ($in !~ /^[ACGUN-]+$/);

    $out = reverse($in);
    $out =~ tr/ACGU/UGCA/;

    return $out;
}

sub translate
{
	my ($in) = @_;
	my $out = '';

	# AAA=K|AAC=N|AAG=K|AAU=N|ACA=T|ACC=T|ACG=T|ACU=T|AGA=R|AGC=S|AGG=R|AGU=S|AUA=I|AUC=I|AUG=M|AUU=I|CAA=Q|CAC=H|CAG=Q|CAU=H|CCA=P|CCC=P|CCG=P|CCU=P|CGA=R|CGC=R|CGG=R|CGU=R|CUA=L|CUC=L|CUG=L|CUU=L|GAA=E|GAC=D|GAG=E|GAU=D|GCA=A|GCC=A|GCG=A|GCU=A|GGA=G|GGC=G|GGG=G|GGU=G|GUA=V|GUC=V|GUG=V|GUU=V|UAA=*|UAC=Y|UAG=*|UAU=Y|UCA=S|UCC=S|UCG=S|UCU=S|UGA=*|UGC=C|UGG=W|UGU=C|UUA=L|UUC=F|UUG=L|UUU=F|

	# die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'") if (length($in) % 3 != 0);
	die("Error while translating: Sequence contains characters other than A, C, G, T, N, -: '$in'") if ($in !~ /^[ACGTN\-]+$/);

	# DNA >> AA
	foreach $triplet ($in =~ /.../g)
	{
		   if ($triplet eq 'AAA') { $out .= 'K'; }
		elsif ($triplet eq 'AAC') { $out .= 'N'; }
		elsif ($triplet eq 'AAG') { $out .= 'K'; }
		elsif ($triplet eq 'AAT') { $out .= 'N'; }
		elsif ($triplet eq 'ACA') { $out .= 'T'; }
		elsif ($triplet eq 'ACC') { $out .= 'T'; }
		elsif ($triplet eq 'ACG') { $out .= 'T'; }
		elsif ($triplet eq 'ACT') { $out .= 'T'; }
		elsif ($triplet eq 'AGA') { $out .= 'R'; }
		elsif ($triplet eq 'AGC') { $out .= 'S'; }
		elsif ($triplet eq 'AGG') { $out .= 'R'; }
		elsif ($triplet eq 'AGT') { $out .= 'S'; }
		elsif ($triplet eq 'ATA') { $out .= 'I'; }
		elsif ($triplet eq 'ATC') { $out .= 'I'; }
		elsif ($triplet eq 'ATG') { $out .= 'M'; }
		elsif ($triplet eq 'ATT') { $out .= 'I'; }
		elsif ($triplet eq 'CAA') { $out .= 'Q'; }
		elsif ($triplet eq 'CAC') { $out .= 'H'; }
		elsif ($triplet eq 'CAG') { $out .= 'Q'; }
		elsif ($triplet eq 'CAT') { $out .= 'H'; }
		elsif ($triplet eq 'CCA') { $out .= 'P'; }
		elsif ($triplet eq 'CCC') { $out .= 'P'; }
		elsif ($triplet eq 'CCG') { $out .= 'P'; }
		elsif ($triplet eq 'CCT') { $out .= 'P'; }
		elsif ($triplet eq 'CGA') { $out .= 'R'; }
		elsif ($triplet eq 'CGC') { $out .= 'R'; }
		elsif ($triplet eq 'CGG') { $out .= 'R'; }
		elsif ($triplet eq 'CGT') { $out .= 'R'; }
		elsif ($triplet eq 'CTA') { $out .= 'L'; }
		elsif ($triplet eq 'CTC') { $out .= 'L'; }
		elsif ($triplet eq 'CTG') { $out .= 'L'; }
		elsif ($triplet eq 'CTT') { $out .= 'L'; }
		elsif ($triplet eq 'GAA') { $out .= 'E'; }
		elsif ($triplet eq 'GAC') { $out .= 'D'; }
		elsif ($triplet eq 'GAG') { $out .= 'E'; }
		elsif ($triplet eq 'GAT') { $out .= 'D'; }
		elsif ($triplet eq 'GCA') { $out .= 'A'; }
		elsif ($triplet eq 'GCC') { $out .= 'A'; }
		elsif ($triplet eq 'GCG') { $out .= 'A'; }
		elsif ($triplet eq 'GCT') { $out .= 'A'; }
		elsif ($triplet eq 'GGA') { $out .= 'G'; }
		elsif ($triplet eq 'GGC') { $out .= 'G'; }
		elsif ($triplet eq 'GGG') { $out .= 'G'; }
		elsif ($triplet eq 'GGT') { $out .= 'G'; }
		elsif ($triplet eq 'GTA') { $out .= 'V'; }
		elsif ($triplet eq 'GTC') { $out .= 'V'; }
		elsif ($triplet eq 'GTG') { $out .= 'V'; }
		elsif ($triplet eq 'GTT') { $out .= 'V'; }
		elsif ($triplet eq 'TAA') { $out .= '*'; }
		elsif ($triplet eq 'TAC') { $out .= 'Y'; }
		elsif ($triplet eq 'TAG') { $out .= '*'; }
		elsif ($triplet eq 'TAT') { $out .= 'Y'; }
		elsif ($triplet eq 'TCA') { $out .= 'S'; }
		elsif ($triplet eq 'TCC') { $out .= 'S'; }
		elsif ($triplet eq 'TCG') { $out .= 'S'; }
		elsif ($triplet eq 'TCT') { $out .= 'S'; }
		elsif ($triplet eq 'TGA') { $out .= '*'; }
		elsif ($triplet eq 'TGC') { $out .= 'C'; }
		elsif ($triplet eq 'TGG') { $out .= 'W'; }
		elsif ($triplet eq 'TGT') { $out .= 'C'; }
		elsif ($triplet eq 'TTA') { $out .= 'L'; }
		elsif ($triplet eq 'TTC') { $out .= 'F'; }
		elsif ($triplet eq 'TTG') { $out .= 'L'; }
		elsif ($triplet eq 'TTT') { $out .= 'F'; }
		else { die("Error while translating: Unrecognised triplet '$triplet': '$in'\n\nGot as far as: '$out'"); }
	}

	die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'\n\nGot as far as: '$out'") if (length($in) % 3 != 0);

	return $out;
}

sub translate_loosely
{
    my ($in) = @_;
    my $out = '';

    # AAA=K|AAC=N|AAG=K|AAU=N|ACA=T|ACC=T|ACG=T|ACU=T|AGA=R|AGC=S|AGG=R|AGU=S|AUA=I|AUC=I|AUG=M|AUU=I|CAA=Q|CAC=H|CAG=Q|CAU=H|CCA=P|CCC=P|CCG=P|CCU=P|CGA=R|CGC=R|CGG=R|CGU=R|CUA=L|CUC=L|CUG=L|CUU=L|GAA=E|GAC=D|GAG=E|GAU=D|GCA=A|GCC=A|GCG=A|GCU=A|GGA=G|GGC=G|GGG=G|GGU=G|GUA=V|GUC=V|GUG=V|GUU=V|UAA=*|UAC=Y|UAG=*|UAU=Y|UCA=S|UCC=S|UCG=S|UCU=S|UGA=*|UGC=C|UGG=W|UGU=C|UUA=L|UUC=F|UUG=L|UUU=F|

    # die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'") if (length($in) % 3 != 0);
    die("Error while translating: Sequence contains characters other than A, C, G, T, N, -: '$in'") if ($in !~ /^[ACGTN\-]+$/);

    # DNA >> AA
    foreach $triplet ($in =~ /.../g)
    {
           if ($triplet eq 'AAA') { $out .= 'K'; }
        elsif ($triplet eq 'AAC') { $out .= 'N'; }
        elsif ($triplet eq 'AAG') { $out .= 'K'; }
        elsif ($triplet eq 'AAT') { $out .= 'N'; }
        elsif ($triplet eq 'ACA') { $out .= 'T'; }
        elsif ($triplet eq 'ACC') { $out .= 'T'; }
        elsif ($triplet eq 'ACG') { $out .= 'T'; }
        elsif ($triplet eq 'ACT') { $out .= 'T'; }
        elsif ($triplet eq 'AGA') { $out .= 'R'; }
        elsif ($triplet eq 'AGC') { $out .= 'S'; }
        elsif ($triplet eq 'AGG') { $out .= 'R'; }
        elsif ($triplet eq 'AGT') { $out .= 'S'; }
        elsif ($triplet eq 'ATA') { $out .= 'I'; }
        elsif ($triplet eq 'ATC') { $out .= 'I'; }
        elsif ($triplet eq 'ATG') { $out .= 'M'; }
        elsif ($triplet eq 'ATT') { $out .= 'I'; }
        elsif ($triplet eq 'CAA') { $out .= 'Q'; }
        elsif ($triplet eq 'CAC') { $out .= 'H'; }
        elsif ($triplet eq 'CAG') { $out .= 'Q'; }
        elsif ($triplet eq 'CAT') { $out .= 'H'; }
        elsif ($triplet eq 'CCA') { $out .= 'P'; }
        elsif ($triplet eq 'CCC') { $out .= 'P'; }
        elsif ($triplet eq 'CCG') { $out .= 'P'; }
        elsif ($triplet eq 'CCT') { $out .= 'P'; }
        elsif ($triplet eq 'CGA') { $out .= 'R'; }
        elsif ($triplet eq 'CGC') { $out .= 'R'; }
        elsif ($triplet eq 'CGG') { $out .= 'R'; }
        elsif ($triplet eq 'CGT') { $out .= 'R'; }
        elsif ($triplet eq 'CTA') { $out .= 'L'; }
        elsif ($triplet eq 'CTC') { $out .= 'L'; }
        elsif ($triplet eq 'CTG') { $out .= 'L'; }
        elsif ($triplet eq 'CTT') { $out .= 'L'; }
        elsif ($triplet eq 'GAA') { $out .= 'E'; }
        elsif ($triplet eq 'GAC') { $out .= 'D'; }
        elsif ($triplet eq 'GAG') { $out .= 'E'; }
        elsif ($triplet eq 'GAT') { $out .= 'D'; }
        elsif ($triplet eq 'GCA') { $out .= 'A'; }
        elsif ($triplet eq 'GCC') { $out .= 'A'; }
        elsif ($triplet eq 'GCG') { $out .= 'A'; }
        elsif ($triplet eq 'GCT') { $out .= 'A'; }
        elsif ($triplet eq 'GGA') { $out .= 'G'; }
        elsif ($triplet eq 'GGC') { $out .= 'G'; }
        elsif ($triplet eq 'GGG') { $out .= 'G'; }
        elsif ($triplet eq 'GGT') { $out .= 'G'; }
        elsif ($triplet eq 'GTA') { $out .= 'V'; }
        elsif ($triplet eq 'GTC') { $out .= 'V'; }
        elsif ($triplet eq 'GTG') { $out .= 'V'; }
        elsif ($triplet eq 'GTT') { $out .= 'V'; }
        elsif ($triplet eq 'TAA') { $out .= '*'; }
        elsif ($triplet eq 'TAC') { $out .= 'Y'; }
        elsif ($triplet eq 'TAG') { $out .= '*'; }
        elsif ($triplet eq 'TAT') { $out .= 'Y'; }
        elsif ($triplet eq 'TCA') { $out .= 'S'; }
        elsif ($triplet eq 'TCC') { $out .= 'S'; }
        elsif ($triplet eq 'TCG') { $out .= 'S'; }
        elsif ($triplet eq 'TCT') { $out .= 'S'; }
        elsif ($triplet eq 'TGA') { $out .= '*'; }
        elsif ($triplet eq 'TGC') { $out .= 'C'; }
        elsif ($triplet eq 'TGG') { $out .= 'W'; }
        elsif ($triplet eq 'TGT') { $out .= 'C'; }
        elsif ($triplet eq 'TTA') { $out .= 'L'; }
        elsif ($triplet eq 'TTC') { $out .= 'F'; }
        elsif ($triplet eq 'TTG') { $out .= 'L'; }
        elsif ($triplet eq 'TTT') { $out .= 'F'; }
        # Gap
        elsif ($triplet eq '---') { $out .= '-'; }
        # Wobble
        elsif ($triplet =~ /^CT.$/) { $out .= 'L'; }
        elsif ($triplet =~ /^GT.$/) { $out .= 'V'; }
        elsif ($triplet =~ /^TC.$/) { $out .= 'S'; }
        elsif ($triplet =~ /^CC.$/) { $out .= 'P'; }
        elsif ($triplet =~ /^AC.$/) { $out .= 'T'; }
        elsif ($triplet =~ /^GC.$/) { $out .= 'A'; }
        elsif ($triplet =~ /^CG.$/) { $out .= 'R'; }
        elsif ($triplet =~ /^GG.$/) { $out .= 'G'; }
        # Ambiguity
        elsif ($triplet eq 'AAY') { $out .= 'N'; }
        elsif ($triplet eq 'AMA') { $out .= 'X'; }
        elsif ($triplet eq 'AMT') { $out .= 'X'; }
        elsif ($triplet eq 'ARA') { $out .= 'X'; }
        elsif ($triplet eq 'ATK') { $out .= 'X'; }
        elsif ($triplet eq 'ATY') { $out .= 'I'; }
        elsif ($triplet eq 'AYA') { $out .= 'X'; }
        elsif ($triplet eq 'CAK') { $out .= 'X'; }
        elsif ($triplet eq 'CMA') { $out .= 'X'; }
        elsif ($triplet eq 'CRT') { $out .= 'X'; }
        elsif ($triplet eq 'GAK') { $out .= 'X'; }
        elsif ($triplet eq 'GAR') { $out .= 'E'; }
        elsif ($triplet eq 'GAW') { $out .= 'X'; }
        elsif ($triplet eq 'GKT') { $out .= 'X'; }
        elsif ($triplet eq 'GMA') { $out .= 'X'; }
        elsif ($triplet eq 'GMT') { $out .= 'X'; }
        elsif ($triplet eq 'GYG') { $out .= 'X'; }
        elsif ($triplet eq 'KAG') { $out .= 'X'; }
        elsif ($triplet eq 'KAT') { $out .= 'X'; }
        elsif ($triplet eq 'KCT') { $out .= 'X'; }
        elsif ($triplet eq 'KTG') { $out .= 'X'; }
        elsif ($triplet eq 'MGT') { $out .= 'X'; }
        elsif ($triplet eq 'RAG') { $out .= 'X'; }
        elsif ($triplet eq 'RAT') { $out .= 'X'; }
        elsif ($triplet eq 'RGT') { $out .= 'X'; }
        elsif ($triplet eq 'RTG') { $out .= 'X'; }
        elsif ($triplet eq 'SCA') { $out .= 'X'; }
        elsif ($triplet eq 'SCT') { $out .= 'X'; }
        elsif ($triplet eq 'STA') { $out .= 'X'; }
        elsif ($triplet eq 'TKG') { $out .= 'X'; }
        elsif ($triplet eq 'TTK') { $out .= 'X'; }
        elsif ($triplet eq 'TYA') { $out .= 'X'; }
        elsif ($triplet eq 'WAA') { $out .= 'X'; }
        elsif ($triplet eq 'WAC') { $out .= 'X'; }
        elsif ($triplet eq 'YAT') { $out .= 'X'; }
        elsif ($triplet eq 'YCC') { $out .= 'X'; }
        elsif ($triplet eq 'YCT') { $out .= 'X'; }
        elsif ($triplet eq 'YTA') { $out .= 'L'; }
        # else { $out .= "X"; print "\n[unknown codon $triplet >> X]\n"; }
        else { $out .= "X"; }

        # # Stop on first stop codon
        # last if ($out =~ /\*$/);
    }

    die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'\n\nGot as far as: '$out'") if (length($in) % 3 != 0);
    # die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'\n\nGot as far as: '$out'") if ((length($in) % 3 != 0) and ($out !~ /\*$/));

    return $out;
}

sub translate_loosely_mitochondrial
{
    my ($in) = @_;
    my $out = '';

    # Mitochondrial changes:
    # AGA >> * instead of R
    # AGG >> * instead of R
    # ATA >> M instead of I
    # ATT >> M instead of I, but only for initiator methionines
    # TGA >> W instead of *
    # No changes to any wobble codons.

    # AAA=K|AAC=N|AAG=K|AAU=N|ACA=T|ACC=T|ACG=T|ACU=T|AGA=R|AGC=S|AGG=R|AGU=S|AUA=I|AUC=I|AUG=M|AUU=I|CAA=Q|CAC=H|CAG=Q|CAU=H|CCA=P|CCC=P|CCG=P|CCU=P|CGA=R|CGC=R|CGG=R|CGU=R|CUA=L|CUC=L|CUG=L|CUU=L|GAA=E|GAC=D|GAG=E|GAU=D|GCA=A|GCC=A|GCG=A|GCU=A|GGA=G|GGC=G|GGG=G|GGU=G|GUA=V|GUC=V|GUG=V|GUU=V|UAA=*|UAC=Y|UAG=*|UAU=Y|UCA=S|UCC=S|UCG=S|UCU=S|UGA=*|UGC=C|UGG=W|UGU=C|UUA=L|UUC=F|UUG=L|UUU=F|

    # die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'") if (length($in) % 3 != 0);
    die("Error while translating: Sequence contains characters other than A, C, G, T, N, -: '$in'") if ($in !~ /^[ACGTN\-]+$/);

    # DNA >> AA
    my $initiator = 1;
    foreach $triplet ($in =~ /.../g)
    {
           if ($triplet eq 'AAA') { $out .= 'K'; }
        elsif ($triplet eq 'AAC') { $out .= 'N'; }
        elsif ($triplet eq 'AAG') { $out .= 'K'; }
        elsif ($triplet eq 'AAT') { $out .= 'N'; }
        elsif ($triplet eq 'ACA') { $out .= 'T'; }
        elsif ($triplet eq 'ACC') { $out .= 'T'; }
        elsif ($triplet eq 'ACG') { $out .= 'T'; }
        elsif ($triplet eq 'ACT') { $out .= 'T'; }
        elsif ($triplet eq 'AGA') { $out .= '*'; }  # mitochondria-specific
        elsif ($triplet eq 'AGC') { $out .= 'S'; }
        elsif ($triplet eq 'AGG') { $out .= '*'; }  # mitochondria-specific
        elsif ($triplet eq 'AGT') { $out .= 'S'; }
        elsif ($triplet eq 'ATA') { $out .= 'M'; }  # mitochondria-specific
        elsif ($triplet eq 'ATC') { $out .= 'I'; }
        elsif ($triplet eq 'ATG') { $out .= 'M'; }
        elsif (($triplet eq 'ATT') and ($initiator == 1)) { $out .= 'M'; }  # mitochondria-specific
        elsif (($triplet eq 'ATT') and ($initiator == 0)) { $out .= 'I'; }
        elsif ($triplet eq 'CAA') { $out .= 'Q'; }
        elsif ($triplet eq 'CAC') { $out .= 'H'; }
        elsif ($triplet eq 'CAG') { $out .= 'Q'; }
        elsif ($triplet eq 'CAT') { $out .= 'H'; }
        elsif ($triplet eq 'CCA') { $out .= 'P'; }
        elsif ($triplet eq 'CCC') { $out .= 'P'; }
        elsif ($triplet eq 'CCG') { $out .= 'P'; }
        elsif ($triplet eq 'CCT') { $out .= 'P'; }
        elsif ($triplet eq 'CGA') { $out .= 'R'; }
        elsif ($triplet eq 'CGC') { $out .= 'R'; }
        elsif ($triplet eq 'CGG') { $out .= 'R'; }
        elsif ($triplet eq 'CGT') { $out .= 'R'; }
        elsif ($triplet eq 'CTA') { $out .= 'L'; }
        elsif ($triplet eq 'CTC') { $out .= 'L'; }
        elsif ($triplet eq 'CTG') { $out .= 'L'; }
        elsif ($triplet eq 'CTT') { $out .= 'L'; }
        elsif ($triplet eq 'GAA') { $out .= 'E'; }
        elsif ($triplet eq 'GAC') { $out .= 'D'; }
        elsif ($triplet eq 'GAG') { $out .= 'E'; }
        elsif ($triplet eq 'GAT') { $out .= 'D'; }
        elsif ($triplet eq 'GCA') { $out .= 'A'; }
        elsif ($triplet eq 'GCC') { $out .= 'A'; }
        elsif ($triplet eq 'GCG') { $out .= 'A'; }
        elsif ($triplet eq 'GCT') { $out .= 'A'; }
        elsif ($triplet eq 'GGA') { $out .= 'G'; }
        elsif ($triplet eq 'GGC') { $out .= 'G'; }
        elsif ($triplet eq 'GGG') { $out .= 'G'; }
        elsif ($triplet eq 'GGT') { $out .= 'G'; }
        elsif ($triplet eq 'GTA') { $out .= 'V'; }
        elsif ($triplet eq 'GTC') { $out .= 'V'; }
        elsif ($triplet eq 'GTG') { $out .= 'V'; }
        elsif ($triplet eq 'GTT') { $out .= 'V'; }
        elsif ($triplet eq 'TAA') { $out .= '*'; }
        elsif ($triplet eq 'TAC') { $out .= 'Y'; }
        elsif ($triplet eq 'TAG') { $out .= '*'; }
        elsif ($triplet eq 'TAT') { $out .= 'Y'; }
        elsif ($triplet eq 'TCA') { $out .= 'S'; }
        elsif ($triplet eq 'TCC') { $out .= 'S'; }
        elsif ($triplet eq 'TCG') { $out .= 'S'; }
        elsif ($triplet eq 'TCT') { $out .= 'S'; }
        elsif ($triplet eq 'TGA') { $out .= 'W'; }  # mitochondria-specific
        elsif ($triplet eq 'TGC') { $out .= 'C'; }
        elsif ($triplet eq 'TGG') { $out .= 'W'; }
        elsif ($triplet eq 'TGT') { $out .= 'C'; }
        elsif ($triplet eq 'TTA') { $out .= 'L'; }
        elsif ($triplet eq 'TTC') { $out .= 'F'; }
        elsif ($triplet eq 'TTG') { $out .= 'L'; }
        elsif ($triplet eq 'TTT') { $out .= 'F'; }
        # Gap
        elsif ($triplet eq '---') { $out .= '-'; }
        # Wobble
        elsif ($triplet =~ /^CT.$/) { $out .= 'L'; }
        elsif ($triplet =~ /^GT.$/) { $out .= 'V'; }
        elsif ($triplet =~ /^TC.$/) { $out .= 'S'; }
        elsif ($triplet =~ /^CC.$/) { $out .= 'P'; }
        elsif ($triplet =~ /^AC.$/) { $out .= 'T'; }
        elsif ($triplet =~ /^GC.$/) { $out .= 'A'; }
        elsif ($triplet =~ /^CG.$/) { $out .= 'R'; }
        elsif ($triplet =~ /^GG.$/) { $out .= 'G'; }
        # Ambiguity
        elsif ($triplet eq 'AAY') { $out .= 'N'; }
        elsif ($triplet eq 'AMA') { $out .= 'X'; }
        elsif ($triplet eq 'AMT') { $out .= 'X'; }
        elsif ($triplet eq 'ARA') { $out .= 'X'; }
        elsif ($triplet eq 'ATK') { $out .= 'X'; }
        elsif ($triplet eq 'ATY') { $out .= 'I'; }
        elsif ($triplet eq 'AYA') { $out .= 'X'; }
        elsif ($triplet eq 'CAK') { $out .= 'X'; }
        elsif ($triplet eq 'CMA') { $out .= 'X'; }
        elsif ($triplet eq 'CRT') { $out .= 'X'; }
        elsif ($triplet eq 'GAK') { $out .= 'X'; }
        elsif ($triplet eq 'GAR') { $out .= 'E'; }
        elsif ($triplet eq 'GAW') { $out .= 'X'; }
        elsif ($triplet eq 'GKT') { $out .= 'X'; }
        elsif ($triplet eq 'GMA') { $out .= 'X'; }
        elsif ($triplet eq 'GMT') { $out .= 'X'; }
        elsif ($triplet eq 'GYG') { $out .= 'X'; }
        elsif ($triplet eq 'KAG') { $out .= 'X'; }
        elsif ($triplet eq 'KAT') { $out .= 'X'; }
        elsif ($triplet eq 'KCT') { $out .= 'X'; }
        elsif ($triplet eq 'KTG') { $out .= 'X'; }
        elsif ($triplet eq 'MGT') { $out .= 'X'; }
        elsif ($triplet eq 'RAG') { $out .= 'X'; }
        elsif ($triplet eq 'RAT') { $out .= 'X'; }
        elsif ($triplet eq 'RGT') { $out .= 'X'; }
        elsif ($triplet eq 'RTG') { $out .= 'X'; }
        elsif ($triplet eq 'SCA') { $out .= 'X'; }
        elsif ($triplet eq 'SCT') { $out .= 'X'; }
        elsif ($triplet eq 'STA') { $out .= 'X'; }
        elsif ($triplet eq 'TKG') { $out .= 'X'; }
        elsif ($triplet eq 'TTK') { $out .= 'X'; }
        elsif ($triplet eq 'TYA') { $out .= 'X'; }
        elsif ($triplet eq 'WAA') { $out .= 'X'; }
        elsif ($triplet eq 'WAC') { $out .= 'X'; }
        elsif ($triplet eq 'YAT') { $out .= 'X'; }
        elsif ($triplet eq 'YCC') { $out .= 'X'; }
        elsif ($triplet eq 'YCT') { $out .= 'X'; }
        elsif ($triplet eq 'YTA') { $out .= 'L'; }
        # else { $out .= "X"; print "\n[unknown codon $triplet >> X]\n"; }
        else { $out .= "X"; }

        # # Stop on first stop codon
        # last if ($out =~ /\*$/);

        $initiator = 0;
    }

    die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'\n\nGot as far as: '$out'") if (length($in) % 3 != 0);
    # die("Error while translating: Sequence length is not a multiple of three (".length($in)." / 3 = ".sprintf("%.3f", (length($in) / 3))."): '$in'\n\nGot as far as: '$out'") if ((length($in) % 3 != 0) and ($out !~ /\*$/));

    return $out;
}

sub fastabreak
{
	our $superreservedlinebreakvariablethingy = $/;
	our $/ = "\n>";
}

sub getfasta
{
	/^>?(.+?)\n(.+)/s or die("Error: Couldn't read FASTA block\n\n['$_']\n\n");

	my $title = $1;
	my $seq = $2;
	#$seq =~ s/[\n> \*]//g;
	$seq =~ s/[\n> ]//g;

	return ($title, $seq);
}

sub split60
{
	my ($seq) = @_;
	$seq =~ s/(.{60})/$1\n/g;

	$seq =~ s/\n$//;

	return $seq;
}

sub normalbreak
{
	our $/ = our $superreservedlinebreakvariablethingy;
}

# Reporting functions

# Example:
#
# addme('things', 'thing A!');			# add to list
# addme('things', 'thing B!');			# add to list
# showme('things');						# show list								unique()s the list.
# showmesome(1000);						# show all lists <1000 (default) items	unique()s the lists.
# showmeall();							# show all lists, ever!					unique()s the lists.
# countme('things');					# return number of list elements		unique()s the list.
# returnme('things');					# return all list elements				unique()s the list.
# findme('thing A!', 'things');			# return 1 if list contains thing A		unique()s the list.
# foreach (returnme('things'))			# return the list						unique()s the list.
# {
# 	[...]
# }
# clearme('things');					# clear list
# clearall();							# clear all lists
#


sub addme
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($things, $thing) = @_;

	die("Error in addme: undefined value") if (!defined($things) or !defined($thing));

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# push(@{$superreservedreportingvariablethingy{$things}}, $thing);
		$superreservedreportingvariablethingy{$things}{$thing} = 1;

			#         # Intermediate uniquing to save memory
			#         if (scalar(@{$superreservedreportingvariablethingy{$things}}) % 100000)
			#         {
			# uniqueme($things);
			# # warn("Warning: addme('$things') is getting gigantic (>100,000)");
			#         }
	}
	else
	{
		# $superreservedreportingvariablethingy{$things} = [$thing];
		$superreservedreportingvariablethingy{$things}{$thing} = 1;
		if (@superreservedreportingvariablethingylist)
		{
			push(@superreservedreportingvariablethingylist, $things);
		}
		else
		{
			@superreservedreportingvariablethingylist = ($things);
		}
	}
}

sub showme
{
	our %superreservedreportingvariablethingy;

	my ($things, $tiny) = @_;

	if (!defined($tiny)) { $tiny = 0; }

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# @{$superreservedreportingvariablethingy{$things}} = unique(@{$superreservedreportingvariablethingy{$things}});

		if ($tiny)
		# if ($tiny == 1)
		{
			# tiny display
			# print "$things: ".scalar(@{$superreservedreportingvariablethingy{$things}})."\n";

			my $max1 = 0;
			foreach my $key (keys(%superreservedreportingvariablethingy))
			{
				if (length($key) > $max1) { $max1 = length($key); }
			}
			my $max2 = 0;
			foreach my $key (keys(%superreservedreportingvariablethingy))
			{
				# if (length(commify(scalar(@{$superreservedreportingvariablethingy{$key}}))) > $max2) { $max2 = length(commify(scalar(@{$superreservedreportingvariablethingy{$key}}))); }
				if (length(commify(scalar(keys(%{$superreservedreportingvariablethingy{$key}})))) > $max2) { $max2 = length(commify(scalar(keys(%{$superreservedreportingvariablethingy{$key}})))); }
			}
			my $space = '';
			my $i = 0;
			# while ($i < (($max1 + $max2) - length($things) - length(commify(scalar(@{$superreservedreportingvariablethingy{$things}})))))
			while ($i < (($max1 + $max2) - length($things) - length(commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}}))))))
			{
				$space .= ' ';
				$i++;
			}

			# Add visual guides
			$space =~ s/   / . /g;
            # $space =~ s/   $/  ./g;   # add them closer to the things string
			$space = reverse($space);
			# Use interpunct
			$space =~ s/\./·/g;

			# print "$things: $space".commify(scalar(@{$superreservedreportingvariablethingy{$things}}))." (e.g. ".${$superreservedreportingvariablethingy{$things}}[0].")\n";
			# print "$things: $space".commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}})))." (e.g. ".(keys(%{$superreservedreportingvariablethingy{$things}}))[0].")\n";
			if (defined($superfastsort))
			{
				# Sort alphabetically
				print "$things: $space".commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}})))." (e.g. ".(sort(keys(%{$superreservedreportingvariablethingy{$things}})))[0].")\n";
			}
			else
			{
				# Sort naturally
				print "$things: $space".commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}})))." (e.g. ".(nsort(keys(%{$superreservedreportingvariablethingy{$things}})))[0].")\n";
			}
		}
		# elsif ($tiny == 2)
		# {
		# 	# extended display, good for multi-line items
		# 	print "\n$things (".scalar(@{$superreservedreportingvariablethingy{$things}})."):\n\n >> ['".join("']\n\n\n >> ['", @{$superreservedreportingvariablethingy{$things}})."']\n\n";
		# }
		else
		{
			# default display
			# print "\n$things (".commify(scalar(@{$superreservedreportingvariablethingy{$things}}))."):\n\n".join("\n", @{$superreservedreportingvariablethingy{$things}})."\n\n";
			if (defined($superfastsort))
			{
				# Sort alphabetically
				print "\n$things (".commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}})))."):\n\n".join("\n", sort(keys(%{$superreservedreportingvariablethingy{$things}})))."\n\n";
			}
			else
			{
				# Sort naturally
				print "\n$things (".commify(scalar(keys(%{$superreservedreportingvariablethingy{$things}})))."):\n\n".join("\n", nsort(keys(%{$superreservedreportingvariablethingy{$things}})))."\n\n";
			}
		}
	}
	else
	{
		# if ($tiny)
		# {
			# tiny display
			print "$things: "."0"."\n";
		# }
		# else
		# {
		# 	# default display
		# 	print "\n$things ("."0"."):\n\n".""."\n\n";
		# }
	}
}

sub showmeall
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($tiny) = @_;

	if (!defined($tiny)) { $tiny = 0; }

	foreach my $key (@superreservedreportingvariablethingylist)					# not sorted
	# foreach my $key (unique(@superreservedreportingvariablethingylist))		# sorted
	{
		showme($key, $tiny);
	}
}

sub showmeallsorted
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($tiny) = @_;

	if (!defined($tiny)) { $tiny = 0; }

	# foreach my $key (@superreservedreportingvariablethingylist)					# not sorted
	foreach my $key (unique(@superreservedreportingvariablethingylist))		# sorted
	{
		showme($key, $tiny);
	}
}

sub showmesome
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($limit) = @_;

	if (!defined($limit)) { $limit = 1000; }

	foreach my $key (@superreservedreportingvariablethingylist)					# not sorted
	# foreach my $key (unique(@superreservedreportingvariablethingylist))		# sorted
	{
		# if there are less than $limit unique items in list...
		# if (scalar(unique(@{$superreservedreportingvariablethingy{$key}})) <= $limit)
		if (scalar(keys(%{$superreservedreportingvariablethingy{$key}})) <= $limit)
		{
			# normal
			showme($key);
		}
		else
		{
			# tiny
			showme($key, 1);
		}
	}
}

sub showmesomesorted
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($limit) = @_;

	if (!defined($limit)) { $limit = 1000; }

	# foreach my $key (@superreservedreportingvariablethingylist)					# not sorted
	foreach my $key (unique(@superreservedreportingvariablethingylist))		# sorted
	{
		# if there are less than $limit unique items in list...
		# if (scalar(unique(@{$superreservedreportingvariablethingy{$key}})) <= $limit)
		if (scalar(keys(%{$superreservedreportingvariablethingy{$key}})) <= $limit)
		{
			# normal
			showme($key);
		}
		else
		{
			# tiny
			showme($key, 1);
		}
	}
}

sub countme
{
	our %superreservedreportingvariablethingy;

	my ($things) = @_;

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# @{$superreservedreportingvariablethingy{$things}} = unique(@{$superreservedreportingvariablethingy{$things}});
		# return(scalar(@{$superreservedreportingvariablethingy{$things}}));

		return(scalar(keys(%{$superreservedreportingvariablethingy{$things}})));
	}
	else
	{
		# warn("Warning: countme(): '$things' is 0");

		return 0;
	}
}

sub findme
{
	our %superreservedreportingvariablethingy;

	my ($thing, $things) = @_;

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# @{$superreservedreportingvariablethingy{$things}} = unique(@{$superreservedreportingvariablethingy{$things}});
		# foreach (@{$superreservedreportingvariablethingy{$things}})

		foreach (keys(%{$superreservedreportingvariablethingy{$things}}))
		{
			if ($_ eq $thing)
			{
				return 1;
			}
		}
	}

	return 0;
}

sub returnme
{
	our %superreservedreportingvariablethingy;

	# my ($things, $nounique) = @_;
	# if (!defined($nounique)) { $nounique = 0; }

	my ($things) = @_;

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# if (!$nounique)
		# {
		# 	@{$superreservedreportingvariablethingy{$things}} = unique(@{$superreservedreportingvariablethingy{$things}});
		# }

		return(keys(%{$superreservedreportingvariablethingy{$things}}));
	}
	else
	{
		return 0;
	}
}

sub uniqueme
{
	our %superreservedreportingvariablethingy;

	my ($things) = @_;

	if (exists($superreservedreportingvariablethingy{$things}))
	{
		# @{$superreservedreportingvariablethingy{$things}} = unique(@{$superreservedreportingvariablethingy{$things}});
		# Already uniqued now.

		return 1;
	}
	else
	{
		return 0;
	}
}

sub clearme
{
	our %superreservedreportingvariablethingy;
	our @superreservedreportingvariablethingylist;

	my ($thing) = @_;

	if (exists($superreservedreportingvariablethingy{$thing}))
	{
		delete($superreservedreportingvariablethingy{$thing});

		@superreservedreportingvariablethingylist = grep { $_ ne $thing } @superreservedreportingvariablethingylist;

		return 1;
	}
	else
	{
		return 0;
	}
}

sub clearall
{
	our %superreservedreportingvariablethingy = ();
	our @superreservedreportingvariablethingylist = ();

	return 1;
}

# Progress functions

# Example:
#
# startme("Doing things!!!!!", 0, 1000);					# 1000 max things!
# stepme(100);
# stopme();
#


sub startme
{
	instakill();
	my ($s, $tiny, $max) = @_;

	our $awesomeprogress;
	our $superreservedprogressvariablethingy = 0;

	if (defined($tiny) and ($tiny == 1))
	{
		# tiny

		print "$s: ";
		# print "[$max] " if (defined($max));
		print "[".commify($max)."] " if (defined($max));
		print "\n" if ($awesomeprogress);
	}
	else
	{
		# big

		# print "\n$s:\n";

		print "\n$s";
		# print " [$max]" if (defined($max));
		print " [".commify($max)."]" if (defined($max));
		print ":\n";
	}

    $| = 0; $| = 1; print "0.."; $| = 0; $| = 1;
}

sub stepme
{
	instakill();
	my ($stepsize, $tiny) = @_;
	our $superreservedprogressvariablethingy;
	our $awesomeprogress;

	$superreservedprogressvariablethingy++;

	if ($awesomeprogress)
	{
        $| = 0; $| = 1; print "\r".commify($superreservedprogressvariablethingy).".."; $| = 0; $| = 1;
	}
	else
	{
        if (($superreservedprogressvariablethingy % $stepsize) == 0) { $| = 0; $| = 1; print commify($superreservedprogressvariablethingy).".."; $| = 0; $| = 1; }
	}
}

sub stopme
{
	my ($tiny) = @_;
	our $superreservedprogressvariablethingy;
	our $awesomeprogress;

	if (defined($tiny) and ($tiny == 1))
	{
		# tiny

		if ($awesomeprogress)
		{
			# print "\r$superreservedprogressvariablethingy  \n";
			print "\r".commify($superreservedprogressvariablethingy)."  \n";
		}
		else
		{
			# print "$superreservedprogressvariablethingy\n";
			print commify($superreservedprogressvariablethingy)."\n";
		}
	}
	else
	{
		# big

		if ($awesomeprogress)
		{
			# print "\r$superreservedprogressvariablethingy  \n\n";
			print "\r".commify($superreservedprogressvariablethingy)."  \n\n";
		}
		else
		{
			# print "$superreservedprogressvariablethingy\n\n";
			print commify($superreservedprogressvariablethingy)."\n\n";
		}
	}
}

sub getme
{
	our $superreservedprogressvariablethingy;
	return $superreservedprogressvariablethingy;
}

# Secondary progress functions: None of the printing, all of the stepping and getting.

# Example:
#
# startme2();
# stepme2();
# getme2();
#
# Same with startme3(), startme4() and startme5().
#

sub startme2 { our $superreservedprogressvariablethingy2 = 0; }
sub stepme2 { our $superreservedprogressvariablethingy2; $superreservedprogressvariablethingy2++; }
sub getme2 { our $superreservedprogressvariablethingy2; return $superreservedprogressvariablethingy2; }

sub startme3 { our $superreservedprogressvariablethingy3 = 0; }
sub stepme3 { our $superreservedprogressvariablethingy3; $superreservedprogressvariablethingy3++; }
sub getme3 { our $superreservedprogressvariablethingy3; return $superreservedprogressvariablethingy3; }

sub startme4 { our $superreservedprogressvariablethingy4 = 0; }
sub stepme4 { our $superreservedprogressvariablethingy4; $superreservedprogressvariablethingy4++; }
sub getme4 { our $superreservedprogressvariablethingy4; return $superreservedprogressvariablethingy4; }

sub startme5 { our $superreservedprogressvariablethingy5 = 0; }
sub stepme5 { our $superreservedprogressvariablethingy5; $superreservedprogressvariablethingy5++; }
sub getme5 { our $superreservedprogressvariablethingy5; return $superreservedprogressvariablethingy5; }

sub startme6 { our $superreservedprogressvariablethingy6 = 0; }
sub stepme6 { our $superreservedprogressvariablethingy6; $superreservedprogressvariablethingy6++; }
sub getme6 { our $superreservedprogressvariablethingy6; return $superreservedprogressvariablethingy6; }

# Easy Timer functions

# Example:
#
# starttime();
# stoptime();
#
# starttime2();
# stoptime2();
#

sub starttime
{
	instakill();
	our $superreservedtimervariablethingy = Time::HiRes::gettimeofday;
}

sub starttime2
{
	our $superreservedtimervariablethingy2 = Time::HiRes::gettimeofday;
}

sub stoptime
{
	our $superreservedtimervariablethingy;

	my ($tiny) = @_;

	if (!defined($tiny))
	{
		$tiny = 0;
	}

	my $s = (Time::HiRes::gettimeofday - $superreservedtimervariablethingy);

	my $min = 0; my $h = 0; my $d = 0;

	while ($s > 60) { $min++; $s -= 60; }
	while ($min > 60) { $h++; $min -= 60; }
	while ($h > 24) { $d++; $h -= 24; }

	if ($d != 0) { $d = $d.' days  '; } else { $d = ''; }
	if ($h != 0) { $h = $h.' h  '; } else { $h = ''; }
	if ($min != 0) { $min = $min.' min  '; } else { $min = ''; }
	if ($s != 0) { $s = sprintf("%.3f", $s).' sec'; } else { $s = ''; }

	if ($tiny == 1)
	{
		print "Time elapsed:\t$d$h$min$s\n";
	}
	elsif ($tiny == 2)
	{
		return "Time elapsed:\t$d$h$min$s\n";
	}
	else
	{
		print "Time elapsed:\t$d$h$min$s\n\n";
	}
}

sub stoptime2
{
	our $superreservedtimervariablethingy2;

	my ($tiny) = @_;

	if (!defined($tiny))
	{
		$tiny = 0;
	}

	my $s = (Time::HiRes::gettimeofday - $superreservedtimervariablethingy2);

	my $min = 0; my $h = 0; my $d = 0;

	while ($s > 60) { $min++; $s -= 60; }
	while ($min > 60) { $h++; $min -= 60; }
	while ($h > 24) { $d++; $h -= 24; }

	if ($d != 0) { $d = $d.' days  '; } else { $d = ''; }
	if ($h != 0) { $h = $h.' h  '; } else { $h = ''; }
	if ($min != 0) { $min = $min.' min  '; } else { $min = ''; }
	if ($s != 0) { $s = sprintf("%.3f", $s).' sec'; } else { $s = ''; }

	if ($tiny == 1)
	{
		print "Time elapsed:\t$d$h$min$s\n";
	}
	elsif ($tiny == 2)
	{
		return "Time elapsed:\t$d$h$min$s\n";
	}
	else
	{
		print "Time elapsed:\t$d$h$min$s\n\n";
	}
}

# Advanced Timer functions

# Example:
#
# $timer = startdis();
# showdis(stopdis($timer));
#

sub startdis
{
	$res = Time::HiRes::gettimeofday;

	return $res;
}

sub stopdis
{
	my ($timer) = @_;

	my $res = (Time::HiRes::gettimeofday - $timer);

	return $res;
}

sub showdis
{
	my ($s, $text) = @_;

	if (!defined($text)) { $text = ''; }

	my $min = 0; my $h = 0; my $d = 0;

	while ($s > 60) { $min++; $s -= 60; }
	while ($min > 60) { $h++; $min -= 60; }
	while ($h > 24) { $d++; $h -= 24; }

	if ($d != 0) { $d = $d.' days  '; } else { $d = ''; }
	if ($h != 0) { $h = $h.' h  '; } else { $h = ''; }
	if ($min != 0) { $min = $min.' min  '; } else { $min = ''; }
	if ($s != 0) { $s = sprintf("%.3f", $s).' sec'; } else { $s = ''; }

	print "Time elapsed:\t$d$h$min$s\t$text\n\n";
}

# Grid Engine
# Job Control Functions

# Example:
#
# $free = freenodes();
# while (($free > 0)
# {
# 	run("Submitting $this", "qsub -wd '$rwd/running/$this' -N 'bl$this' -e '$rwd/running/$this/_ERROR_.txt' -o '$rwd/running/$this/_OUTPUT_.txt' '$lwd/running/$this/_SCRIPT_.sh' '$this'", !switch('debug'));
# }
#
# 3307043 0.34775 j18196.sh  qwang        dr    11/20/2010 02:43:45 all.q@fmb23.lmb.internal           1
# 3317941 0.26470 nd2        blang        qw    11/25/2010 18:47:28                                    1
#

sub mynodes			# get number of BLADE farm CPUs currently running my jobs.
{
	my $n = 0;

	if (qtype() eq 'qsub')
	{
		# GridEngine
		
		my $qstat = `qstat -u \$USER 2>&1`;
		foreach (split(/\n/, $qstat))
		{
			# skip headers
			next if (/^job-ID/);
			next if (/^------/);
			
			# Workaround for qstat failures
			if (/^error:/)
			{
				$n += 9999;
			}

	        # if (/\s+(\d+)\s+$/) { $n += $1; }
	        if (/\s+r\s+\S+\s+\S+\s+\S+\s+(\d+)\s+$/) { $n += $1; }
			# else { die("Error: Couldn't parse '$_'"); }
		}
		
		# # open(QSTAT, "qstat -u \$USER -U \$USER|");
		# open(QSTAT, "qstat -u \$USER|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	# skip headers
		# 	next if (/^job-ID/);
		# 	next if (/^------/);
		#
		# 	        # if (/\s+(\d+)\s+$/) { $n += $1; }
		# 	        if (/\s+r\s+\S+\s+\S+\s+\S+\s+(\d+)\s+$/) { $n += $1; }
		# }
		# close(QSTAT);


		# open(QSTAT, "qsummary|grep \$USER|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		#   if (/^\S+\s+(\d+)\s+\d+$/)
		# 	{
		# 		$n += $1;
		# 	}
		# }
		# close(QSTAT);
	}
	elsif (qtype() eq 'lsf')
	{
		# LSF (Note: It's actually just the job number here, not number of CPUs)
		$tmp = `busers \$USER | awk '{print \$4}' | grep -v 'NJOBS'`;
		chomp($tmp);
		$n += $tmp;
	}
	else
	{
        # Local (iMac)
        open(BJOBS, "ps -lax|");
        $n = 0;
        while (<BJOBS>)
        {
        	chomp;

        	# Remove leading whitespace before split
        	s/^\s+//;

        	# Skip empty lines
        	next if /^$/;

            # Split input
        	@a = split(/ +/, $_);

            if (($a[14] eq '/users/gt/blang/bin/perl') and ($a[15] eq '-w') and ($a[16] =~ /\.pl$/))
        	{
                # Count jobs
                $n++;
        	}
        }
        close(BJOBS);
	}

	return $n;
}

sub thesenodes			# get number of BLADE farm CPUs currently running jobs from a particular pipeline
{
	my $n = 0;

	if (qtype() eq 'qsub')
	{
		# GridEngine
		# open(QSTAT, "qstat -u \$USER -U \$USER -r|");
		open(QSTAT, "qstat -u \$USER -r|");
		while (<QSTAT>)
		{
			chomp;

			# "       Full jobname:     pipeline_evolutionary_rate_analysis_main_pl_human_Capra3_DISOPRED2"
			# get e.g. 'pipeline_evolutionary_rate_analysis'
			$locale = locale();
			$locale =~ s/\s+$//;
			#state("LOCALE '$locale'");
	        if (/Full jobname:\s+$locale\_/) { $n += 1; }
		}
		close(QSTAT);
	}
    elsif (qtype() eq 'lsf')
    {
        # LSF (Note: It's actually just the job number here, not number of CPUs)
        open(BJOBS, "bjobs -w -u \$USER | awk '{print \$7}' | grep -v 'JOB_NAME'|");
        while (<BJOBS>)
        {
            chomp;

            # "update_systemage_run_pl"
            # get e.g. 'pipeline_evolutionary_rate_analysis'
            $locale = locale();
            $locale =~ s/\s+$//;
            # state("LOCALE '$locale'");
            if (/^$locale\_/) { $n += 1; }
        }
        close(BJOBS);
    }
    elsif (qtype() eq 'lsf')
    {
        # LSF (Note: It's actually just the job number here, not number of CPUs)
        open(BJOBS, "bjobs -w -u \$USER | awk '{print \$7}' | grep -v 'JOB_NAME'|");
        while (<BJOBS>)
        {
            chomp;

            # "update_systemage_run_pl"
            # get e.g. 'pipeline_evolutionary_rate_analysis'
            $locale = locale();
            $locale =~ s/\s+$//;
            # state("LOCALE '$locale'");
            if (/^$locale\_/) { $n += 1; }
        }
        close(BJOBS);
    }
    else
    {
        # Local (iMac)
        # open(BJOBS, "ps -lax|");
        $bjobs = `ps -lax`;
        %parents = ();
        # while (<BJOBS>)
        foreach (split(/\n/, $bjobs))
        {
        	chomp;

        	# Remove leading whitespace before split
        	s/^\s+//;

        	# Skip empty lines
        	next if /^$/;

            # Split input
        	@a = split(/ +/, $_);

            if (($a[14] eq '/users/gt/blang/bin/perl') and ($a[15] eq '-w') and ($a[16] =~ /\.pl$/))
            {
        	    # Process output
                @b = @a;
                splice(@b, 0, 16);
                $b[0] =~ s/^\/users\/gt\/blang\///;

                # "update_splice_elm"
                $locale = locale();
                $locale =~ s/\s+$//;

                $b[0] =~ s/\//_/g;

                if (($b[0] =~ /^$locale\_/) or (exists($parents{$a[1]})))
                {
                    # Count jobs
                    $n += 1;
                    # Count jobs that are children of these, too
                    $parents{$a[0]} = 1;
                }
            }
        }
        # close(BJOBS);
    }

	return $n;
}

sub busynodes			# get number of BLADE farm CPUs currently in use by anybody.
{
	my $n = 0;

	if (qtype() eq 'qsub')
	{
		# GridEngine

		# open(QSTAT, "qstat -g c -U \$USER|");
		# open(QSTAT, "qstat -g c -q long-sl7|");
		# open(QSTAT, "qstat -g c -q long-sl7|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	# skip headers
		# 	next if (/^CLUSTER QUEUE/);
		# 	next if (/^------/);
		#
		# 	@a = split(/\s+/);
		#
		# 	$n = $a[2];
		# }
		# close(QSTAT);
		# open(QSTAT, "qhost -q | grep -P '(long-sl7|mem_256|mem_512|mem_1tb)' | perl -ane 'next if ((\$F[3] ne \"\") and (\$F[3] ne \"a\")); print \$F[2].\"\n\"' |");		# 'a' for alarm is still okay, it means no new jobs but existing ones continue to run normally
		open(QSTAT, "qhost -q | grep -P '(long-sl7)' | perl -ane 'next if ((\$F[3] ne \"\") and (\$F[3] ne \"a\")); print \$F[2].\"\n\"' |");		# 'a' for alarm is still okay, it means no new jobs but existing ones continue to run normally
		   # long-sl7             BIP   0/0/16        d
		   # long-sl7             BIP   0/16/16
		   # long-sl7             BIP   0/0/16        adu");
		while (<QSTAT>)
		{
			# 0/16/16
			chomp;

			@a = split(/\//);

			$n += $a[1];
		}
		close(QSTAT);

		# open(QSTAT, "qsummary|grep 'CPUs in use'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^#CPUs in use=\s+-(\d+)$/)
		# 	{
		# 		$n = $1;
		# 	}
		# }
		# close(QSTAT);
	}
	else
	{
		# LSF (Here it's the CPU numbers, like for GridEngine)

		# open(BHOSTS, q(bhosts -w | perl -ane 'if (($F[1] eq 'ok')) {print join("\t", @F)."\n";}'|));
		# while (<BHOSTS>)
		# {
		# 	chomp;
		# 	@a = split(/\t/);
		# 	$n += $a[4];
		# }
		# close(BHOSTS);

        # $n = readpipe(q(bhosts -w | perl -ane 'if (($F[1] eq 'ok')) {print join("\t", @F)."\n";}' | datamash sum 5));
		$n = readpipe(q(bhosts -w | perl -ane 'if (($F[1] eq 'ok') or ($F[1] eq 'closed_Full')) {print join("\t", @F)."\n";}' | awk '{sum+=$5} END {print sum}'));
		chomp($n);
	}

	return $n;
}

sub queuednodes
{
	my $n = 0;

	if (qtype() eq 'qsub')
	{
		# GridEngine
		my $qstat = `qstat -u \$USER -s p 2>&1`;
		foreach (split(/\n/, $qstat))
		{
			# skip headers
			next if (/^job-ID/);
			next if (/^------/);
			
			# Workaround for qstat failures
			if (/^error:/)
			{
				$n += 9999;
			}
			
	        if (/\s+(\d+)\s+$/) { $n += $1; }
			# else { die("Error: Couldn't parse '$_'"); }
		}

		# open(QSTAT, "qstat -u \$USER -s p|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	# skip headers
		# 	next if (/^job-ID/);
		# 	next if (/^------/);
		#
		# 	        if (/\s+(\d+)\s+$/) { $n += $1; }
		# }
		# close(QSTAT);
	}
	else
	{
		# # LSF (Here it's the CPU numbers, like for GridEngine)
		# $n = readpipe(q(bhosts -w | perl -ane 'if (($F[1] eq 'ok') or ($F[1] eq 'closed_Full')) {print join("\t", @F)."\n";}' | awk '{sum+=$5} END {print sum}'));
		# chomp($n);
		die("Error: Not implemented");
	}

	return $n;
}

sub freenodes
{
	my ($talk) = @_;

	our $extrememode;

	my $avail = 0;
	my $availpercent = 0;
	my $cpus = 0;
	my $nonerror = 0;

	my $freeday = 0;
	my $freenight = 0;

	if (qtype() eq 'qsub')
	{
		# GridEngine

		# # open(QSTAT, "qstat -g c -U \$USER|");
		# open(QSTAT, "qstat -g c -q long-sl7|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	# skip headers
		# 	next if (/^CLUSTER QUEUE/);
		# 	next if (/^------/);
		#
		# 	@a = split(/\s+/);
		#
		# 	$avail = $a[4];
		# 	$cpus = $a[5];
		# 	# $nonerror = $cpus - $a[7];
		# 	$nonerror = $cpus - ($a[6] + $a[7]);
		# }
		# close(QSTAT);

		# open(QSTAT, "qhost -q | grep -P '(long-sl7|mem_256|mem_512|mem_1tb)' | perl -ane 'next if ((\$F[3] ne \"\") and (\$F[3] ne \"a\")); print \$F[2].\"/\$F[3]\n\"' |");
		open(QSTAT, "qhost -q | grep -P '(long-sl7)' | perl -ane 'next if ((\$F[3] ne \"\") and (\$F[3] ne \"a\")); print \$F[2].\"/\$F[3]\n\"' |");
		   # long-sl7             BIP   0/0/16        d
		   # long-sl7             BIP   0/16/16
		   # long-sl7             BIP   0/0/16        adu");
		while (<QSTAT>)
		{
			# 0/16/16
			chomp;

			@a = split(/\//);

			$a[3] = '' if (!exists($a[3]));

			if ($a[3] eq '')
			{
				# No error
				$nonerror += $a[2];
				$cpus += $a[2];
				$avail += $a[2] - $a[1];
			}
			elsif ($a[3] eq 'a')
			{
				 # Alarm state (won't add new jobs, but will run normally)
				$cpus += $a[2];
				$avail += $a[2] - $a[1];
			}
		}
		close(QSTAT);

		# open(QSTAT, "qsummary|grep 'CPUs on cluster'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^#CPUs on cluster=\s+(\d+)$/)
		# 	{
		# 		$cpus = $1;
		# 	}
		# }
		# close(QSTAT);
		#
		# open(QSTAT, "qsummary|grep 'CPUs unavailable'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^#CPUs unavailable=\s+-(\d+)$/)
		# 	{
		# 		$nonerror = $cpus - $1;
		# 	}
		# }
		# close(QSTAT);
		#
		# open(QSTAT, "qsummary|grep 'CPUs available'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^#CPUs available=\s+(\d+)\s+\((\d+)%\)$/)
		# 	{
		# 		$avail = $1;
		# 		# $availpercent = $2;
		# 	}
		# }
		# close(QSTAT);

		$mine = mynodes();
		$busy = busynodes() - $mine;

		$free = $avail;

		# open(QSTAT, "qsummary|grep 'From 8am to 8pm you may now claim'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^From 8am to 8pm you may now claim (\d+) CPUs$/)
		# 	{
		# 		$freeday = $1;
		# 	}
		# }
		# close(QSTAT);
		#
		# open(QSTAT, "qsummary|grep 'From 8pm to 8am you may now claim'|");
		# while (<QSTAT>)
		# {
		# 	chomp;
		#
		# 	if (/^From 8pm to 8am you may now claim (\d+) CPUs$/)
		# 	{
		# 		$freenight = $1;
		# 	}
		# }
		# close(QSTAT);

		$freeday = round(($free + $mine) / 3, 0) - $mine;
		$freenight = round(($free + $mine) / 2, 0) - $mine;
	}
	else
	{
		# LSF (Here it's the CPU numbers, like for GridEngine)

		$cpus = readpipe(q(bhosts -w | awk '{sum+=$4} END {print sum}'));
		chomp($cpus);

		$nonerror = readpipe(q(bhosts -w | perl -ane 'if (($F[1] eq 'ok') or ($F[1] eq 'closed_Full')) { print join("\t", @F)."\n"; }' | awk '{sum+=$4} END {print sum}'));
		chomp($nonerror);

		$avail = $nonerror - busynodes();

		$mine = mynodes();
		$busy = busynodes() - $mine;

		$free = $avail;

		$freeday = round(($free + $mine) / 3, 0) - $mine;
		$freenight = round(($free + $mine) / 2, 0) - $mine;
	}

	my $time = '';
	my $submit = 0;

	if ($extrememode == 1)
	{
		$time = 'extreme mode, 2/3';
		$submit = round(($free + $mine) * 2 / 3, 0) - $mine;		# extreme mode: use up to two thirds of available nodes
	}
	elsif ($politemode == 1)
	{
		$time = 'polite mode, 1/3';
		$submit = round(($free + $mine) / 3, 0) - $mine;			# polite mode: use up to one third of available nodes (also at nighttime)
	}
	else
	{
		# ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
		my @t = localtime();
		if (($t[2] >= 10) and ($t[2] < 19) and ($t[6] != 6) and ($t[6] != 0))
		{
			$time = 'daytime, 1/3';
			# $submit = round(($free + $mine) / 3, 0) - $mine;		# daytime: use up to one third of available nodes
            # $submit = $freeday - $mine;       # daytime: use up to one third of available nodes
			$submit = $freeday;		# daytime: use up to one third of available nodes
		}
		else
		{
			$time = 'nighttime or weekend, 1/2';
			# $submit = round(($free + $mine) / 2, 0) - $mine;		# nighttime: use up to half of available nodes

            # $submit = $freenight - $mine;     # nighttime: use up to half of available nodes
			$submit = $freenight;		# nighttime: use up to half of available nodes

			# $time = 'nighttime or weekend, 1/3';
			# $submit = round(($free + $mine) / 3, 0) - $mine;		# nighttime: be highly polite and still only use up to one third of available nodes
		}
	}

	$queued = queuednodes();

	if ($talk)
	{
		if (qtype() eq 'qsub')
		{
			print " >> Gridengine\n";
            print "   >> $nonerror of $cpus cores available, $busy in use by others, $mine are mine ($queued queued), $avail free\n";
            print "   >> $submit to submit ($time)\n\n";
		}
		else
		{
            print " >> LSF cluster\n";
            print "   >> $nonerror of $cpus cores available, $busy in use by others, $mine are mine ($queued queued), $avail free\n";
            print "   >> $submit to submit ($time)\n\n";
            # system("~/scripts/bsummary.sh");
		}
	}

	return $submit;
}

sub qtype
{
	$hostname = `hostname`;
	chomp($hostname);
	$hostname =~ s/\.lmb\.internal$//;
	$hostname =~ s/\.embl\.de$//;
	$hostname =~ s/\.linux\.crg\.es$//;

    # @mainframes_qsub = qw(hex hal alpha beta delta eta epsilon sigma);
	@mainframes_qsub = qw(hex hal alpha beta delta eta epsilon sigma ant-login1 ant-login2 ant-login3 ant-login4 ant-login5 ant-login6 ant-login7 ant-login8 ant-login9);
    # @mainframes_lsf = qw(submaster build000);
    # @mainframes_lsf = qw(submaster1 compute1-001);
	@mainframes_lsf = qw(submaster1);

    # "compute-0-10.local" is what a sigma node looks like and I'm returning qsub for them because they run qsub.
	if (contains($hostname, @mainframes_qsub) or ($hostname =~ /^compute-\d+-\d+\.local$/))
	{
		return "qsub";
	}
	elsif (contains($hostname, @mainframes_lsf))
	{
		return "lsf";
	}
	else
	{
		# Assume compute node if unknown hostname
		return "node";
	}
}

sub waitforjobs			# Wait for jobs from a particular script pipeline directory to finish
{
	state("Waiting for jobs to finish...");

	if (qtype() ne 'node')
	{
		# if script is running locally on a mainframe:
		$minnodes = 0;
	}
	else
	{
		# if script is running on its own node:
		$minnodes = 1;
	}

	while (thesenodes() > $minnodes)
	{
		instakill();
		sleep(3);
	}
}

sub waitforalljobs			# Wait for all of my jobs to finish
{
	state("Waiting for all jobs to finish...");

	if (qtype() ne 'node')
	{
		# if script is running locally on a mainframe:
		$minnodes = 0;
	}
	else
	{
		# if script is running on its own node:
		$minnodes = 1;
	}

	while (mynodes() > $minnodes)
	{
		sleep(3);
	}
}

sub getnumjobs			# LEGACY == get number of BLADE farm CPUs currently in use. slightly dodgy but works. == LEGACY
{
	my $j = 0;

	open(QSTAT, "qstat -u \"*\"|");
	while (<QSTAT>)
	{
		chomp;
		if (/\s+(\d+)\s+$/) { $j += $1; }
	}
	close(QSTAT);

	return $j;
}

# MySQL
# Query count / wait functions

# Example:
#
# state(getnumqueries()." queries running");
# waitforallqueries();
#

sub getnumqueries         # get number of MySQL queries currently running
{
    my $n = 0;

    $n = `~/scripts/p.sh | grep -vP "Sleep\\s+\\d+\\s+\\[NULL\\]" | wc -l`;
    $n -= 2;
    die("Error: Couldn't parse number of queries from ~/scripts/p.sh") if ($n !~ /^\d+$/);

    return $n;
}

sub waitforallqueries 		# Wait for all of my MySQL queries to finish
{
	state("Waiting for all queries to finish...");

    # Debug: Sleep a minimum of 5 seconds to ensure freshly started jobs are counted
	sleep(5);

	while (getnumqueries() > 0)
	{
        # state(getnumqueries());
		sleep(5);
	}
}

# Argument parsing (@ARGV)

# Example:
#
# script.pl human oct4 -aples
#
# ($needed1, $needed2) = args(2);
# if (switch('aples')) { print "APLE DAAY! :D" }
#

sub args
{
	my ($needed, $silent) = @_;
	our @ARGV;
	our $usage;
	our %superreservedswitchvariablethingy = ();

	if (!defined($silent)) { $silent = 0; }

	my @r = ();
	my $i;
	my $s = '';

	if (@ARGV < $needed)
	{
		# print "\nError: More parameters needed. ($needed required)\n";
		$s = "\nError: More parameters needed. ($needed required)\n";

		$i = 1;
		foreach $a (@ARGV)
		{
			# print " >> parameter $i: '$a'\n";
			$s .= " >> parameter $i: '$a'\n";
			$i++;
		}

		# print "\n";
		$s .= "\n";
		if (defined($usage))
		{
			# print "Usage:\n$usage";
			$s .= "Usage:\n$usage";
		}
		else
		{
			# print "Usage: $0";
			$s .= "Usage: $0";
			$i = 1;
			while ($i <= $needed)
			{
				# print " [parameter $i]";
				$s .= " [parameter $i]";
				$i++;
			}
		}

		# print "\n\n";
		$s .= "\n\n";
		# exit;
		die($s);
	}

	$i = 1;
	foreach my $a (@ARGV)
	{
		if ($i <= $needed)
		{
		    if ($a =~ /^-/)
		    {
        		# print "\nError: Parameter looks like a switch (starts with '-') ($needed required)\n\n";
        		# exit;
        		die("\nError: Parameter looks like a switch (starts with '-') ($needed required)\n\n");
		    }

			push(@r, $a);
		}
		else
		{
			$a =~ s/^-//;
			$superreservedswitchvariablethingy{$a} = 1;

			print "Switch active: -$a\n" unless ($silent == 1);
		}

		$i++;
	}

	return @r;
}

sub switch
{
	my ($switch) = @_;
	our %superreservedswitchvariablethingy;

	if (exists($superreservedswitchvariablethingy{$switch}) and ($superreservedswitchvariablethingy{$switch} == 1)) { return 1 }
	else { return 0 }
}

sub switchstring
{
	my ($switch) = @_;
	our %superreservedswitchvariablethingy;

	if (exists($superreservedswitchvariablethingy{$switch}) and ($superreservedswitchvariablethingy{$switch} == 1)) { return "-$switch" }
	else { return "" }
}

sub setswitch
{
	my ($switch) = @_;
	our %superreservedswitchvariablethingy;
	
	$superreservedswitchvariablethingy{$switch} = 1;
}

# New R functions (statistics) (using Statistics::R)

# Examples:
#
#
# startr();															# initialize
#
# setr('a', \@a);													# set a to @a!
# runr('w = wilcox.test(a, b)')										# run custom R code (like a Mann-Whitney U test!)
# $p = getr('w$p.value');											# get that p-value into $p!
#
# setr('x', 2);														# set x to 2!
# $x = getr('x');													# get that 2 into $x!
#
# stopr();															# shut down (not necessary)


sub startr
{
	instakill();

	our $superreservedrvariablethingy = Statistics::R -> new( r_bin => "/users/gt/blang/bin/R" );

	runr("options(max.print=1000000000)");
}

sub setr
{
	instakill();

    # Note: This only works for single text variables and arrays of numbers
	our $superreservedrvariablethingy;

	startr() if (!defined($superreservedrvariablethingy));

	my ($var, @a) = @_;

    # Simple version (sometimes very slow): just use the built-in set function (broken)
    # $superreservedrvariablethingy -> set($var, @{$a[0]});


    # More elaborate version (should be much faster): sometimes use temporary files to set values (much faster for large arrays)

    if ((scalar(@a) == 1) and (ref($a[0]) eq ''))
    {
            # If it's just one value and not a reference:
            # Use the built-in set function from Statistics::R
            # state("SIMPLE VALUE: SETTING '$var' TO '$a[0]'");
          $superreservedrvariablethingy -> set($var, $a[0]);
    }
    elsif ((scalar(@a) == 1) and (ref($a[0]) eq 'SCALAR'))
    {
            # If it's a scalar reference (not directly supported by Statistics::R):
            # Use the built-in set function from Statistics::R
            # state("SCALAR REFERENCE: SETTING '$var' TO '".${$a[0]}."'");
          $superreservedrvariablethingy -> set($var, ${$a[0]});
    }
    elsif ((scalar(@a) == 1) and (ref($a[0]) ne 'ARRAY'))
    {
            # If it's a non-scalar, non-array reference:
            # Use the built-in set function from Statistics::R
            # state("NON-SCALAR, NON-ARRAY REFERENCE: SETTING '$var' TO '$a[0]'");
          $superreservedrvariablethingy -> set($var, $a[0]);
    }
    else
    {
            # Otherwise, import it via a temporary file

            # Turn reference into array if it's an array reference
          if ((scalar(@a) == 1) and (ref($a[0]) eq 'ARRAY'))
            {
                # state("ARRAY REFERENCE:");
                @a = @{$a[0]};
            }
            # else
            # {
            #     state("ARRAY:");
            # }

            # Check if these are numbers
            my $numeric = 1;
            foreach my $s (@a)
            {
                if (!looks_like_number($s))
                {
                    $numeric = 0;
                }
            }

            # Generate temporary file ($$ means process id)
          $rand = int(rand(1000000000));
          my $rin = "tmp-R-set-$$-$rand.txt";

            # Write temporary file
            # c(1, 2, 3, 4, 5)
          open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");
          if ($numeric == 1)
          {
                # state("NUMERIC:");
              print R "c(".join(", ", @a).")\n";
          }
          else
          {
                # state("NON-NUMERIC:");
              print R 'c("'.join('", "', @a).'")'."\n";
          }
          close(R);

            # run("cat", "cat $rin");

            # Read temporary file
      $superreservedrvariablethingy -> run("$var <- dget(file=\"$rin\")");

            # Remove temporary file
          system("rm -f '$rin'");
    }
}

sub getr
{
	instakill();

	our $superreservedrvariablethingy;

	startr() if (!defined($superreservedrvariablethingy));

	return $superreservedrvariablethingy -> get(@_);
}

sub runr
{
	instakill();

	startr() if (!defined($superreservedrvariablethingy));

	my ($command, $loud) = @_;

	if (!defined($loud)) { $loud = 0; }

	our $superreservedrvariablethingy;

	my $result = '';

	if ($loud == 1)
	{
		foreach my $line (split(/\n/, $command))
		{
			# print command
			print "R >> $line\n";

			# strip whitespace at beginning and end
			$line =~ s/^\s+//g;
			$line =~ s/\s+$//g;

			# # strip comments
			$line =~ s/\s+#.+$//g;

			# skip empty lines
			next if ($line =~ /^\s*$/);
			next if ($line =~ /^\s*;$/);

			# skip comment lines
			next if ($line =~ /^\s*#/);

			# run command
			state("\n\nR COMMAND >> ['$line']\n\n");

			$res = $superreservedrvariablethingy -> run($line);

			print "R >> $res\n\n";

			$result .= $res;
		}
	}
	else
	{
		# strip comments
		$command =~ s/\s+#.+$//mg;
		$command =~ s/^\s*#.+$//mg;

		# strip trailing whitespace/newlines and multiple newlines
		$command =~ s/^\s+//mg;
		$command =~ s/\s+$//mg;
		$command =~ s/\n{2,}/\n/g;
		$command =~ s/^\s*$//mg;
		$command =~ s/^\s*;$//mg;
		$command =~ s/^[\s\n]+//sg;
		$command =~ s/[\s\n;]+$//sg;

		# strip comments
		$command =~ s/\s+#.+$//mg;
		$command =~ s/^\s*#.+$//mg;

		if ($command ne '')
		{
			#DEBUG
            # state("\n\nR COMMAND >> ['$command']\n\n");
			#END DEBUG

			$result = $superreservedrvariablethingy -> run($command);
		}
		
		if ($loud == 2)
		{
			if ($result ne '')
			{
				state("R >> $result");
			}
		}
	}

	return $result;
}

sub stopr
{
	our $superreservedrvariablethingy;

	$superreservedrvariablethingy -> stop();
}




# Old R functions (statistics) (old clunky method involving temporary files)

# Examples:
#
#
# $result = runr("a <- 1"."\n"."a");								# arbitrary R code (old clunky method)
#
# $p_value = fisher(1, 9, 11, 3);									# small numbers - 2x2 contingency table
# $p_value = fisher(1, 9, 11, 3, 1);								# small numbers - 2x2 contingency table with output (loud)
# $p_value = chisquare(18239, 10845, 1044, 10933);					# large numbers - 2x2 contingency table
# $p_value = chisquare(18239, 10845, 1044, 10933, 104, 10933);		# large numbers - 3x2 contingency table (one more row)
#
# $p_value = wilcoxon(\@a, \@b);									# Wilcoxon rank-sum test (= Mann-Whitney U) (two independent samples, not the paired Wilcoxon signed-rank test)
# $p_value = ks(\@a, \@b);											# KS test
# $p_value = ksboot(\@a, \@b);										# KS test with bootstrapping
# $p_value = ttest(\@a, \@b);										# two-sample t-test
# $p_value = ttest(\@a, 5);											# one-sample t-test, comparing vs. mu = 5 (theoretical mean of the population, given null hypothesis)
#

# Contingency tables - 2x2, small values (total N < 400)

sub runR
{
	my($s, $loud) = @_;

	if (!defined($loud)) { $loud = 0; }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";
	my $res = '';

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");
	print R "$s\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		$res .= $_;
	}
	print "\n" if ($loud == 1);

	if ($res eq '') { die("Error: Empty R output in '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $res;
}

sub fisher
{
	if (scalar(@_) % 2 == 1)
	{
		$loud = pop(@_);
	}
	else
	{
		$loud = 0;
	}

	# # Warning if numbers are too large for Fisher's Exact
	# # OK actually - it's an exact test after all. Can use it for anything.
	# if (sum(@_) > 400)
	# {
	# 	warn("Warning: Numbers seem too large for Fisher's Exact Test") if ($loud == 1);
	# }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");
	print R "a <- matrix(c(".join(", ", @_)."), ncol=2, byrow=TRUE)\n";
	print R "a\n" if ($loud == 1);
	print R "f = fisher.test(a)\n";
	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		if (/^ \$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# Contingency tables - Nx2, large values (> 10 each)

sub chisquare
{
	if (scalar(@_) % 2 == 1)
	{
		$loud = pop(@_);
	}
	else
	{
		$loud = 0;
	}

	# Warning if numbers are too small for chi-square
	foreach (@_)
	{
		if ($_ < 10)
		{
			warn("Warning: Numbers are too small for Chi-Square Test") if ($loud == 1);
		}
	}

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");
	print R "a <- matrix(c(".join(", ", @_)."), ncol=2, byrow=TRUE)\n";
	print R "a\n" if ($loud == 1);
	print R "f = chisq.test(a)\n";
	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		if (/^ \$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# Student's t-test (one or two-sample)

sub ttest
{
	my($a, $b, $loud) = @_;

	if (!defined($loud)) { $loud = 0; }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");

	# two-sample test
	if (scalar(@$b) > 1)
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- c(".join(", ", @$b).")\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = t.test(a, b)\n";
	}
	else
	# one-sample test vs. mu
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- $b\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = t.test(a, mu=b)\n";
	}

	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		if (/^ \$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# Wilcoxon rank-sum test (= Mann-Whitney U) (non-parametric, for independent samples, i.e. not the Wilcoxon signed-rank test)

sub wilcoxon
{
	my($a, $b, $loud, $alternative) = @_;

	if (!defined($loud)) { $loud = 0; }
	if (!defined($alternative)) { $alternative = 'two.sided'; }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");

	# two-sample test
	if (scalar(@$b) > 1)
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- c(".join(", ", @$b).")\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = wilcox.test(a, b, alternative='$alternative')\n";
	}
	else
	# one-sample test vs. mu
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- $b\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = wilcox.test(a, mu=b, alternative='$alternative')\n";
	}

	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		if (/^ \$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# KS-test (Kolmogorov-Smirnov)

sub ks
{
	my($a, $b, $loud) = @_;

	if (!defined($loud)) { $loud = 0; }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");

	# two-sample test
	if (scalar(@$b) > 1)
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- c(".join(", ", @$b).")\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = ks.test(a, b)\n";
	}
	else
	# one-sample test vs. distribution
	{
		warn("One-sample KS-test: Testing against normal distribution, check if this makes sense");
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- \"pnorm\"\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = ks.test(a, b)\n";
	}

	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		if (/^ \$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# Bootstrap KS-test (Kolmogorov-Smirnov) (tolerates ties)

sub ksboot
{
	my($a, $b, $loud) = @_;

	if (!defined($loud)) { $loud = 0; }

	$rand = int(rand(1000000000));
	my $rin = "tmp-R-in-$rand.txt";
	my $rout = "tmp-R-out-$rand.txt";

	open(R, ">$rin") or die("\nError: Couldn't open '$rin'\n\n");
	print R "library(\"Matching\")\n";

	# two-sample test
	if (scalar(@$b) > 1)
	{
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- c(".join(", ", @$b).")\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = ks.boot(a, b)\n";
	}
	else
	# one-sample test vs. distribution
	{
		warn("One-sample bootstrap KS-test: Testing against normal distribution, check if this makes sense");
		print R "a <- c(".join(", ", @$a).")\n";
		print R "b <- \"pnorm\"\n";
		print R "a\n" if ($loud == 1);
		print R "b\n" if ($loud == 1);
		print R "f = ks.boot(a, pnorm())\n";
	}

	print R "str(f)\n";
	close(R);

	run("R", "/users/gt/blang/bin/R --vanilla --slave < $rin > $rout", 1);

	open(R, $rout) or die("\nError: Couldn't open '$rout'\n\n");
	my $p = undef;
	print "\n" if ($loud == 1);
	while (<R>)
	{
		print if ($loud == 1);
		chomp;

		# 	  ..$ p.value    : num 1.71e-14
		if (/\$ p\.value\s+: num ([\d\.e\-\+]+)$/)
		{
			$p = $1;
		}
	}
	print "\n" if ($loud == 1);

	if (!defined($p)) { die("Error: Couldn't parse p-value from '$rout'"); }

	run("Remove temporary R files", "rm -f $rin", 1);
	run("Remove temporary R files", "rm -f $rout", 1);

	return $p;
}

# File functions

# Examples:
#
#
# ($outfile) = Open(OUT, ">output-$list.txt");							# > means overwrite, >> means append, no prefix means open for reading
# ($outfile) = Open(OUT, ">output-$list.txt", "Oncomine", 1);			# 1 means be verbose, Oncomine is the file 'type' (>> "opened Oncomine output file ....")
#
# Close(OUT, $outfile);
# Close(OUT, $outfile, "Oncomine", 1);
#
#
# run("Some other script lol :D", "thatotherscript.pl -aples", 1)		# 1 means silent mode!
#

sub Open
{
	my ($handle, $file, $type, $loud) = @_;

	if (!defined($loud)) { $loud = 1; }
	if (!defined($type)) { $type = ''; }

	$file =~ /^(>*)/;
	if (defined($1)) { $mode = $1; } else { $mode = ""; }
	$file =~ s/^>*//;
	if ($mode eq '') { $type = "input ".$type; }							# 		read
	if ($mode eq '>') { $type = "output ".$type; }					# > 	write
	if ($mode eq '>>') { $type = "output (append) ".$type; }		# >>	append

	$type =~ s/\s+$//;

	open($handle, $mode.$file) or die("\nError: Couldn't open $type file '$file'\n\n");

	if ($loud == 1)
	{
		print "Opened $type file '$file'\n";
	}

	return ($file);
}

sub Close
{
	my ($handle, $file, $type, $loud) = @_;

	if (!defined($loud)) { $loud = 1; }
	if (!defined($type)) { my $type = ''; }

	$file =~ /^(>*)/;
	if (defined($1)) { $mode = $1; } else { $mode = ""; }
	$file =~ s/^>*//;
	if ($mode eq '') { $type = "".$type; }							# 		read
	if ($mode eq '>') { $type = "output ".$type; }					# > 	write
	if ($mode eq '>>') { $type = "output (append) ".$type; }		# >>	append

	close($handle);

	if ($loud == 1)
	{
		print "Wrote to $type file '$file'\n";
	}
}

sub run
{
	instakill();

	my ($what, $command, $silent) = @_;

	if (!defined($silent)) { $silent = 0; }

	die("Error: No command in run") if(!defined($command));

	if ($silent != 1)
	{
		print "\n--------------------------------------------------------------------------------\n";
		print "> $what ($command)";
		print "\n--------------------------------------------------------------------------------\n\n";
	}

#	open(COMMAND, $command."|");
#
#	while (<COMMAND>)
#	{
#		if ($silent != 1) { print; }
#	}

	# Replace ~ with /users/gt/blang
	# $command =~ s/^~/\/users\/gt\/blang/;
	$command =~ s/^~/\/users\/gt\/blang/;

	system($command);
	# system(" . /users/gt/blang/.login && $command");

	#print "\n -- Press ENTER to continue -- \n";
	#<STDIN>;#
}

# Current working directory functions

sub cd
{
	my ($directory, $silent) = @_;

	if (!defined($silent)) { $silent = 0; }

	die("Error: No directory in cd") if(!defined($directory));

	if ($silent != 1)
	{
		print "\n( ---------------------------------------------------------------------------- )\n";
		print "( Changing directory to $directory )";
		print "\n( ---------------------------------------------------------------------------- )\n\n";
	}

	# Replace ~ with /users/gt/blang
	$directory =~ s/^~/\/users\/gt\/blang/;

	chdir($directory);

	#print "\n -- Press ENTER to continue -- \n";
	#<STDIN>;#
}

sub Cwd
{
	# current working directory
	return dirname(rel2abs($0));
}

sub Hwd
{
	# hex working directory (to be used when running on pcgn5)
	my $dir = dirname(rel2abs($0));

	# From	/nfs/nfs2/nfs2/
	# To	/lmb/home/

	$dir =~ s/\/nfs\/nfs2\/nfs2\//\/lmb\/home\//;

	return $dir;
}

sub Pwd
{
	# pcgn5 working directory (to be used when running on hex)
	my $dir = dirname(rel2abs($0));

	# From	/lmb/home/
	# To	/nfs/nfs2/nfs2/

	$dir =~ s/\/lmb\/home\//\/nfs\/nfs2\/nfs2\//;

	return $dir;
}

# NCBI Taxon ID functions

# Examples:
#
#
# print species2tax('human');			# prints '9606'!
# print tax2species('9606');			# prints 'human'!
#
#

sub species2tax
{
	my ($species) = @_;

	my $query = Query("SELECT DISTINCT tax FROM tax WHERE species='$species'");
	my ($tax) = Fetch($query);

	die("Error: Couldn't translate species '$species' into an NCBI Taxon ID") if ((Numrows($query) != 1) or (!defined($tax)) or ($tax eq ''));

	return $tax;
}

sub tax2species
{
	my ($tax) = @_;

	my $query = Query("SELECT name FROM tax WHERE tax='$tax' AND type='scientific name'");
	my ($species) = Fetch($query);

	die("Error: Couldn't translate NCBI Taxon ID '$tax' into a species") if ((Numrows($query) != 1) or (!defined($species)) or ($species eq ''));

	return $species;
}


# Amino acid functions

# Examples:
#
#
# print 1to3('K');			# prints 'LYS'
# print 3to1('LYS');		# prints 'K'
#
#

sub onetothree
{
	my ($aa) = @_;
	$aa = uc($aa);

	$res = '';

	if ($aa eq 'A') { $res = 'ALA'; }
	elsif ($aa eq 'R') { $res = 'ARG'; }
	elsif ($aa eq 'N') { $res = 'ASN'; }
	elsif ($aa eq 'D') { $res = 'ASP'; }
	elsif ($aa eq 'C') { $res = 'CYS'; }
	elsif ($aa eq 'E') { $res = 'GLU'; }
	elsif ($aa eq 'Q') { $res = 'GLN'; }
	elsif ($aa eq 'G') { $res = 'GLY'; }
	elsif ($aa eq 'H') { $res = 'HIS'; }
	elsif ($aa eq 'I') { $res = 'ILE'; }
	elsif ($aa eq 'L') { $res = 'LEU'; }
	elsif ($aa eq 'K') { $res = 'LYS'; }
	elsif ($aa eq 'M') { $res = 'MET'; }
	elsif ($aa eq 'F') { $res = 'PHE'; }
	elsif ($aa eq 'P') { $res = 'PRO'; }
	elsif ($aa eq 'U') { $res = 'SEC'; }	# Selenocysteine
	elsif ($aa eq 'S') { $res = 'SER'; }
	elsif ($aa eq 'T') { $res = 'THR'; }
	elsif ($aa eq 'W') { $res = 'TRP'; }
	elsif ($aa eq 'Y') { $res = 'TYR'; }
	elsif ($aa eq 'V') { $res = 'VAL'; }

	# if ($res eq '')
	# {
	# 	die("Error: Couldn't convert one-letter AA code '$aa' (unknown)");
	# }

	return $res;
}

sub threetoone
{
	my ($aa) = @_;
	$aa = uc($aa);

	$res = '';

	if ($aa eq 'ALA') { $res = 'A'; }
	elsif ($aa eq 'ARG') { $res = 'R'; }
	elsif ($aa eq 'ASN') { $res = 'N'; }
	elsif ($aa eq 'ASP') { $res = 'D'; }
	elsif ($aa eq 'CYS') { $res = 'C'; }
	elsif ($aa eq 'GLU') { $res = 'E'; }
	elsif ($aa eq 'GLN') { $res = 'Q'; }
	elsif ($aa eq 'GLY') { $res = 'G'; }
	elsif ($aa eq 'HIS') { $res = 'H'; }
	elsif ($aa eq 'ILE') { $res = 'I'; }
	elsif ($aa eq 'LEU') { $res = 'L'; }
	elsif ($aa eq 'LYS') { $res = 'K'; }
	elsif ($aa eq 'MET') { $res = 'M'; }
	elsif ($aa eq 'PHE') { $res = 'F'; }
	elsif ($aa eq 'PRO') { $res = 'P'; }
	elsif ($aa eq 'SEC') { $res = 'U'; }	# Selenocysteine
	elsif ($aa eq 'SER') { $res = 'S'; }
	elsif ($aa eq 'THR') { $res = 'T'; }
	elsif ($aa eq 'TRP') { $res = 'W'; }
	elsif ($aa eq 'TYR') { $res = 'Y'; }
	elsif ($aa eq 'VAL') { $res = 'V'; }

	# if ($res eq '')
	# {
	# 	die("Error: Couldn't convert three-letter AA code '$aa' (unknown)");
	# }

	return $res;
}


# Array / String functions


sub intersection
{
    # produce intersection array from two arrays
    # case sensitive
    my ($a, $b) = @_;

    my %union = ();
    my %intersection = ();
    foreach (@$a)
    {
        $union{$_} = 1;
    }

    foreach (@$b)
    {
        if (exists($union{$_}))
        {
            $intersection{$_} = 1;
        }
    }

    return keys(%intersection);
}

sub union
{
    # produce union array from two arrays
    # case sensitive
    my ($a, $b) = @_;

    my @union = unique(@$a, @$b);

    return @union;
}

sub lonly
{
	# produce array of values not shared by two arrays
	# case sensitive
	my ($a, $b) = @_;

    # warn("WARNING: This function (lonly) uses List::Compare and I'm not sure this module works correctly");

    $lc = List::Compare->new($a, $b);
	return ($lc->get_Lonly);
}

sub ronly
{
	# produce array of values not shared by two arrays
	# case sensitive
	my ($a, $b) = @_;

    # warn("WARNING: This function (ronly) uses List::Compare and I'm not sure this module works correctly");

    $lc = List::Compare->new($a, $b);
	return ($lc->get_Ronly);
}

sub symdiff
{
	# produce array of values not shared by two arrays
	# case sensitive
	my ($a, $b) = @_;

	warn("WARNING: This function (symdiff) uses List::Compare and I'm not sure this module works correctly");

    $lc = List::Compare->new($a, $b);
	return ($lc->get_symdiff);
}

sub contains
{
	# find string in array
    # "needle, haystack"
	# case sensitive
	my ($s, @a) = @_;

	my $contains = 0;
	foreach (@a)
	{
		if ($s eq $_) { $contains = 1; last; }
	}

	return $contains;
}

sub contains_i
{
	# case INsensitive
	my ($s, @a) = @_;

	my $contains = 0;
	foreach (@a)
	{
		if (lc($s) eq lc($_)) { $contains = 1; last; }
	}

	return $contains;
}

sub positions
{
	# get positions of a substring in a string, returned as an array
	# case sensitive
	my ($s, $l) = @_;

	my @a = ();

	while ($l =~ m/($s)/g)
	{
		# push(@a, pos($l));	# This is broken for substrings longer than 1 (pos() will report the end of the substring, not its beginning)
		# push(@a, pos($l) - (length($s) - 1));	# This is broken when using a pattern (like e.g. [KNSTY] in ~/update/uniprot/unimod_control.pl)

		push(@a, pos($l) - (length($1) - 1));	# This should successfully allow patterns as well as substrings longer than 1
	}

	return @a;
}

sub apositions
{
	# get positions of string in an array, returned as an array
	# case sensitive
	my ($s, @a) = @_;

	my @b = ();

	my $i = 0;
	while ($i < scalar(@a))
	{
		if ($a[$i] eq $s) { push(@b, $i); }
		$i++;
	}

	return @b;
}

sub aindex
{
	# index of string in array
	# case sensitive
	my ($s, @a) = @_;

	my $i = 0;
	foreach (@a)
	{
		if ($s eq $_) { last; }
		$i++;
	}

	return $i;
}

sub aindex_i
{
	# index of string in array
	# case INsensitive
	my ($s, @a) = @_;

	my $i = 0;
	foreach (@a)
	{
		if (lc($s) eq lc($_)) { last; }
		$i++;
	}

	return $i;
}

sub adelete
{
	# delete specific string from an array
	# case sensitive
	my ($s, @a) = @_;

	if (contains($s, @a))
	{
    	splice(@a, aindex($s, @a), 1);
	}

	return @a;
}

sub adelete_i
{
	# delete specific string from an array
	# case INsensitive
	my ($s, @a) = @_;

	if (contains_i($s, @a))
	{
    	splice(@a, aindex_i($s, @a), 1);
	}

	return @a;
}

sub aequal
{
    # from http://perldoc.perl.org/perlfaq4.html#How-do-I-test-whether-two-arrays-or-hashes-are-equal%3f
    my ($first, $second) = @_;
    # no warnings; # silence spurious -w undef complaints
    return 0 unless @$first == @$second;
    for (my $i = 0; $i < @$first; $i++)
    {
        return 0 if $first->[$i] ne $second->[$i];
    }
    return 1;
}

sub Sleep
{
	my ($i) = @_;

	print "Waiting $i seconds...\n";

	sleep($i);
}

sub unique
{
	if (defined($superfastsort))
	{
		# sort alphabetically
		@_ = sort(@_);
	}
	else
	{
		# sort naturally (slooow)
		# state("Unsorted:");
		# show(@_);
		@_ = nsort(@_);
		# state("NatSorted:");
		# show(@_);
	}

	my $prev = "";
	return grep($_ ne $prev && (($prev) = $_), @_);
}

sub uniq
{
	# sort numerically
	@_ = sort {$a <=> $b} @_;

	my $prev = "";
	return grep($_ ne $prev && (($prev) = $_), @_);
}

sub pairs
{
	# make a unique pairs array, tab-separated
	my @a = @_;
	@a = unique(@a);

	my @p = ();
	foreach $a (@a)
	{
		foreach $b (@a)
		{
			next if ($a eq $b);

			# push(@p, "$a\t$b");
			#
			# this would prevent some duplicates, but it would also mix the sets (with this commented out set 1 will always end up in the 'ensp1' column etc.)
			if ($a lt $b) { push(@p, "$a\t$b"); }
			else { push(@p, "$b\t$a"); }
		}
	}

	# unique pairs array!
	@p = unique(@p);

	return @p;
}



# Other functions



sub esc
{
	# escapes things for mysql
	my ($s) = @_;

	# escape \ to \\
	$s =~ s/\\/\\\\/g;

	# escape ' to \'
	$s =~ s/'/\\'/g;

	return $s;
}

sub percent
{
    # formats two numbers into a percentage
	my ($one, $two) = @_;

	my $s = 'N/A';
	if ($two != 0)
	{
    	$s = round(($one/$two * 100), 1)."%";
	}

	return $s;
}

sub log10
{
    my $n = shift;
    return log($n) / log(10);
}

sub decimals
{
    # get number of decimals
	my ($i) = @_;

	my $res = 0;
	if ($i =~ /^-?\d+\.(\d+)$/)
	{
		$res = length($1);
	}

	return $res;
}

sub Continue
{
	my ($text) = @_;
	if ($text eq '') { $text = "Continue?"; }

	print "\n$text (y/n)\n";

	my $input = <STDIN>;
	chomp($input);

	my $res = 0;
	if (lc($input) ne 'y')
	{
		print "\nExiting!\n\n";
		exit;
	}
}

sub Ask
{
	my ($text) = @_;

	print "\n$text (y/n)\n";

	my $input = <STDIN>;
	chomp($input);

	my $res = 0;
	if (lc($input) eq 'y') { $res = 1; }

	return $res;
}

sub chompme
{
    # Chomp a string and return it
    my ($s) = @_;

    if (!defined($s)) { $s = ''; }

    # chomp($s);
    $s =~ s/\n+$//;

    return $s;
}

sub state
{
    instakill();

    # make a statement. (i.e. print something with nice newlines around it for dramatic effect)
    my ($s, $tiny) = @_;

    if (!defined($s)) { $s = ''; }
    if (!defined($tiny)) { $tiny = 0; }

    if ($tiny)
    {
        print "$s\n";
    }
    else
    {
        print "\n$s\n\n";
    }
}

sub nl
{
	instakill();

	print "\n";
}

sub done
{
	my ($tiny) = @_;

    if (!defined($tiny)) { $tiny = 0; }

	if ($tiny)
	{
		print "Done!\n";
	}
	else
	{
		print "\nDone!\n\n";
	}
}

return 1;
