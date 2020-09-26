#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $superreservedmysqldb;

our $usage = "$0 [input file] [table name] [-overwrite] [-nonulls] [-nonas] [-quiet] [-addheader] [-keepcomments] [-notype] [-double] [-moreindices] [-allindices]\n\nImport data from a CSV or tab-separated file into a MySQL table, which gets created automatically.\n -overwrite: Remove MySQL table with this name first, if it already exists.\n -nonulls: Do NOT convert blank fields to NULL.\n -nonas: Do NOT convert NA fields to NULL.\n -quiet: Don't show progress messages.\n -addheader: Artificially add a header of 'field1', 'field2', etc.\n -keepcomments: Keep comment lines (starting with '#' or '%'), e.g. for when the header line starts with '#' or '%'\n -notype: Leave all fields as text and don't add any indices (much faster).\n -moreindices: Use more fields as indices (1% unique, 'bad').\n -allindices: Just use all fields as indices.\n\nNote: The first line of the input file should be a header line with field names.\n\nExample: $0 file.csv table_name";
($infile, $table) = args(2);

our $superloudmysql = 1 if (!switch('quiet'));

if ($table =~ /(.+)\.(.+)/)
{
	Change($1);
	$table = $2;
}

# $infile = "input.txt";
# $outfile = "output.txt";

starttime();
open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

die("Error: Input file is empty") if (!-s $infile);

# Choose linebreak
$/ = "\r\n";
@a = <IN>;
if (scalar(@a) == 1)
{
	$/ = "\r";
	seek(IN, 0, 0);
	@a = <IN>;
	if (scalar(@a) == 1)
	{
		$/ = "\n";
		seek(IN, 0, 0);
		@a = <IN>;
		if (scalar(@a) == 1)
		{
			die("Error: Couldn't determine linebreak type in file '$infile'");
		}
	}
}
seek(IN, 0, 0);

$linebreak = '';
if ($/ eq "\r\n")
{
	$linebreak = '\r\n';
}
elsif ($/ eq "\r")
{
	$linebreak = '\r';
}
elsif ($/ eq "\n")
{
	$linebreak = '\n';
}
else
{
	die("Error: Unknown linebreak type in file '$infile': \n['$/']");
}
state("Linebreak type: $linebreak") if (!switch('quiet'));

# start

$query = Query("SHOW TABLES LIKE '$table'");
if (Numrows($query) != 0)
{
	if (switch('overwrite'))
	{
		Query("DROP TABLE IF EXISTS `$table`");
		state("Dropped table '$table'") if (!switch('quiet'));
	}
	else
	{
		die("Error: Table '$table' already exists (use -overwrite to remove it first)");
	}
}

$headerlines = 1;
if (!switch('keepcomments'))
{
	while ($_ = <IN>)
	{
		chomp;
		if (/^[#%]/)
		{
			state("Skipping comment line $headerlines: $_", 1) if (!switch('quiet'));
			$headerlines++;
			next;
		}
		else
		{
			last;
		}
	}
	state("Skipping $headerlines header lines in SQL import later", 1) if (!switch('quiet'));
}
else
{
	$_ = <IN>; chomp;
}
@fields = split(/\t/, $_);
if (scalar(@fields) == 1)
{
    # @fields = splitcsv($_);
	@fields = split(/,/);
	if (scalar(@fields) == 1)
	{
		die("Error: Format not CSV or tab-separated (can only find one column in table header)");
	}
	state("Assuming format: CSV") if (!switch('quiet'));
	$format = 'CSV';
}
else
{
	state("Assuming format: Tab-separated") if (!switch('quiet'));
	$format = 'TSV';
}
die if (($format ne 'CSV') and ($format ne 'TSV'));
@fields = removetrailingfields(@fields);

die("Error: First line in '$infile' appears empty (tried to get column headers from there)") if (scalar(@fields) == 0);

$doublequoted = 1;
if (switch('addheader'))
{
	$i = 1;
	foreach (@fields)
	{
		# Check whether ALL fields are surrounded by double quotes
		if (($_ !~ /^"/) or ($_ !~ /"$/))
		{
			$doublequoted = 0;
		}

		$_ = "field$i";
		$i++;
	}
	seek(IN, 0, 0);
}
else
{
	foreach (@fields)
	{
		# Check whether ALL fields are surrounded by double quotes
		if (($_ !~ /^"/) or ($_ !~ /"$/))
		{
		    state("Field NOT double-quoted: '$_'") if (!switch('quiet'));
			$doublequoted = 0;
		}

		# Make field names lower case and replace anything non-alphanumeric with underscores, then collapse successive underscores to one and remove leading and trailing underscores
		$_ = lc($_);
		s/[^\w]/_/g;
		s/_+/_/g;
		s/_+$//g;
		s/^_+//g;
		
		# Limit them to 64 characters in length (the MySQL maximum)
		if (length($_) > 64)
		{
			$trim = substr($_, 0, 64);
			state("Trimming column name '$_' to '$trim' (64 character maximum)") if (!switch('quiet'));
			$_ = $trim;
		}

		if ($_ eq '')
		{
			state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
			die("Error: Field name can't be empty in file '$infile' (try -addheader if the first line isn't a header)");
		}
		if ($_ eq 'id')
		{
			state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
			die("Error: Field name can't be 'id' in file '$infile' (it's reserved)");
		}
	}

	state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
}
if ($doublequoted == 1)
{
	state("Stripping double quotes (\"\")") if (!switch('quiet'));
}

state("Creating table '$table'") if (!switch('quiet'));

$q = "CREATE TABLE `$table` (
";

$i = 0;
foreach $field (@fields)
{
	# Rename some problematic fields (for MySQL)
	if ($field eq 'set')
	{
		$field = 'myset';
	    state("Renamed field: 'set' to '$field'") if (!switch('quiet'));
	}
	# if ($field eq 'mode')
	# {
	# 	$field = 'mymode';
	#     state("Renamed field: 'mode' to '$field'") if (!switch('quiet'));
	# }
	
	# Always text to begin with
	$q .= "  `$field` text,\n";

	# Next field number
	$i++;
}

$q =~ s/,\n$//s;

state("", 1) if (!switch('quiet'));

$q .= "
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='$infile'";


# state($q);
Query($q);

# state("Reading '$infile' and inserting into table '$table'...");
# starttime();
# $query = Query("LOAD DATA INFILE '".rel2abs($infile)."' INTO TABLE `$table` FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '$linebreak' IGNORE 1 LINES");
# done();
# state("Inserted ".Numrows($query)." rows");
# stoptime();

nl() if (switch('quiet'));
state("Reading '$infile' and inserting into table '$table'...", switch('quiet'));

# Get path from MySQL server side
$serverpath = `pwd`;
$serverpath =~ s/[\r\n]+$//;
# state("['$serverpath']");
$serverpath =~ s/^\/home\/lang\/home\///;
# $serverpath =~ s/^\/Users\/lang\///;
$serverpath .= "/$infile";
# state("['$serverpath']");
# exit;

seek(IN, 0, 0);
# skip headers
<IN>;
$_ = <IN>;
chomp;

# Construct import query

$q = qq(LOAD DATA LOCAL INFILE '$serverpath' INTO TABLE `$table`);

if ($format eq 'CSV')
{
	$q .= qq( FIELDS TERMINATED BY ',');
}
elsif ($format eq 'TSV')
{
	$q .= qq( FIELDS TERMINATED BY '\\t');
}

if ($doublequoted == 1)
{
	$q .= qq( ENCLOSED BY '"');
}

$q .= qq( LINES TERMINATED BY '$linebreak');

if (!switch('addheader'))
{
	# $q .= q( IGNORE 1 LINES);
	$q .= qq( IGNORE $headerlines LINES);
}

$query = Query($q);

state("Loaded ".commify(Numrows($query))." rows", switch('quiet'));

# Wait for table to appear
state("Waiting for table to appear...") if (!switch('quiet'));
$query = Query("SELECT `TABLE_NAME` FROM `information_schema`.`tables` WHERE `TABLE_SCHEMA`='$superreservedmysqldb' AND `TABLE_NAME`='$table'");
while (Numrows($query) == 0)
{
	sleep(0.5);
}
$query = Query("SELECT * FROM `$table` LIMIT 1");
if (Numrows($query) == 0)
{
	die("Error: Table '$table' is empty");
}
state("OK") if (!switch('quiet'));

if (!switch('nonulls'))
{
	startme("Converting blank fields to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`=''");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field`=''");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (!switch('nonas'))
{
	startme("Converting NA and N/A fields (from R) to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`='NA' OR `$field`='N/A'");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field`='NA' OR `$field`='N/A'");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));


# $q = "CREATE TABLE `$table` (
#   `id` int(10) unsigned NOT NULL auto_increment,
# ";

$q = "ALTER TABLE `$table`
ADD COLUMN `id` INT(10) UNSIGNED NOT NULL AUTO_INCREMENT FIRST,\n";


if (!switch('notype'))
{
	nl() if (switch('quiet'));
	state("Determining best field types and selecting indices...", switch('quiet'));
	%fieldtype = ();
	@indices = ();
	@badindices = ();
	$i = 0;
	foreach $field (@fields)
	{
		# @values = ();

		# $mainquery = Query("SELECT * FROM `$table`");

		# startme("$field [1 of ".scalar(@fields)."]", 1, Numrows($mainquery)) if (!switch('quiet'));
		state(">> $field", 1) if (!switch('quiet'));

		$maxlen = 0;
		$numeric = 1;
		$double = 1;
		$uniq = 0;
		$baduniq = 0;

		# if ((defined($a[$i])) and ($a[$i] =~ /[^\d]/))
		$query = Query("SELECT 1 FROM `$table` WHERE `$field` REGEXP '[^0-9]' LIMIT 1");
		if (Numrows($query) != 0)
		{
			$numeric = 0;
		}

		# if ((defined($a[$i])) and (($a[$i] =~ /[^\d\.\-]/) or ($a[$i] =~ /^.+-/)))
		$query = Query("SELECT 1 FROM `$table` WHERE `$field` REGEXP '^[^0-9]+\$' OR `$field` REGEXP '[^\\-0-9eE.+]' OR `$field` REGEXP '-[^0-9]' OR `$field` REGEXP '^[eE]' OR `$field` REGEXP '[eE]\$' OR `$field` REGEXP '-[^eE]+-' LIMIT 1");
		if (Numrows($query) != 0)
		{
			$double = 0;
		}

		# if ((defined($a[$i])) and (length($a[$i]) > $maxlen))
		$query = Query("SELECT MAX(LENGTH(`$field`)) FROM `$table`");
		# {
		# 	$maxlen = length($a[$i]);
		# }
		($maxlen) = FetchOne($query);
		$maxlen = 0 if (!defined($maxlen));

		# if (defined($a[$i]))
		# {
		# 	push(@values, $a[$i]);
		# }
		# else
		# {
		# 	push(@values, '');
		# }

		# See if field content is fairly unique (# unique values >= 10% of #total values), and therefore informative as an index:
		# if (scalar(unique(@values)) >= (scalar(@values) / 10))
		# {
		# 	$uniq = 1;
		# }
		$query = Query("SELECT COUNT(`$field`), COUNT(DISTINCT `$field`) FROM `$table`");
		($values, $uniquevalues) = FetchOne($query);
		if ($uniquevalues >= ($values / 10))
		{
			$uniq = 1;
		}
		# See if field content is vaguely unique (# unique values >= 1% of #total values), and therefore informative as an index:
		# elsif (scalar(unique(@values)) >= (scalar(@values) / 100))
		# {
		# 	$baduniq = 1;
		# }
		elsif ($uniquevalues >= ($values / 1000))
		{
			$baduniq = 1;
		}

		if ($maxlen == 0)
		{
			# Shouldn't happen
			$maxlen = 1;
		}

		# Choose 'best' field type
		if ($numeric == 1)
		{
			# Int
			# $q .= "  `$field` int(10) default NULL,\n";
			$q .= "CHANGE COLUMN `$field` `$field` INT(10) DEFAULT NULL,\n";
			$fieldtype{$field} = 'numeric';

			# as index
			if ($uniq == 1)
			{
				push(@indices, $field);
			}

			# as bad index
			if ($baduniq == 1)
			{
				push(@badindices, $field);
			}
		}
		elsif ($double == 1)
		{
			# Double
			# $q .= "  `$field` double default NULL,\n";
			# if (switch('double'))
			# {
			# 	$q .= "CHANGE COLUMN `$field` `$field` DOUBLE DEFAULT NULL,\n";
			# }
			# else
			# {
				$q .= "CHANGE COLUMN `$field` `$field` DOUBLE DEFAULT NULL,\n";
			# }
			$fieldtype{$field} = 'double';
		}
		elsif ($maxlen <= 250)
		{
			# Varchar
			# $q .= "  `$field` varchar($maxlen) default NULL,\n";
			# $q .= "  `$field` varchar(250) default NULL,\n";
			$q .= "CHANGE COLUMN `$field` `$field` VARCHAR(250) DEFAULT NULL,\n";
			$fieldtype{$field} = 'varchar';

			# as index
			if ($uniq == 1)
			{
				push(@indices, $field);
			}

			# as bad index
			if ($baduniq == 1)
			{
				push(@badindices, $field);
			}
		}
		else
		{
			# Text
			# $q .= "  `$field` text,\n";
			# $q .= "CHANGE COLUMN `$field` `$field` TEXT,\n";		# already text
			$fieldtype{$field} = 'text';
		}

		state("  >> $fieldtype{$field}", 1) if (!switch('quiet'));

		# Next field number
		$i++;
	}

	# $q .= "
	#   PRIMARY KEY  (`id`)
	# ";

	state("", 1) if (!switch('quiet'));

	# If there are no decent 10% unique indices, get 1% unique indices...or just use all if the switch is set
	if (switch('allindices'))
	{
		state("Using all indices (-allindices)", switch('quiet'));
		@indices = @fields;
	}
	elsif ((@indices == 0) or (switch('moreindices')))
	{
		state("No 10% unique indices -- using 1% unique indices (bad)", switch('quiet'));
		@indices = @badindices;
	}
	else
	{
		state("Using 10% unique indices (good)", switch('quiet'));
	}

	foreach $index (unique(@indices))
	{
		$indextitle = ucfirst($index);

		state(">> $indextitle", 1) if (!switch('quiet'));

		if ($fieldtype{$index} ne 'text')
		{
			# $q .= ",
			#   KEY `$indextitle` (`$index`)";
			$q .= "ADD INDEX `$indextitle` (`$index`),\n";
		}
		else
		{
			# $q .= ",
			#   KEY `$indextitle` USING HASH (`$index`(200))";
			$q .= "ADD INDEX `$indextitle` (`$index`(200)),\n";
		}
	}
}

$q .= q( ADD PRIMARY KEY (`id`));

# $q =~ s/,\n$//s;


state("Changing field types and adding indices...") if (!switch('quiet'));
Query($q);
done() if (!switch('quiet'));

nl() if (!switch('quiet'));
nl() if (!switch('quiet'));
nl() if (switch('quiet'));

stoptime();

Optimize($table);
