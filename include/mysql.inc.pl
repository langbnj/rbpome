#!/users/gt/blang/bin/perl -w
# use Mysql;
use DBD::mysql;
# use DBD::Pg;

# Print out each MySQL query before running it?
# $superloudmysql = 1;

# Connect to default database
# our $superreservedmysqldb = 'blang';
our $superreservedmysqldb = 'blang';
Connect($superreservedmysqldb);

sub Connect
{
	my ($database) = @_;

	instakill();

	$server = 'cluster';

	our $superreservedmysqldb = $database;

	# If connection doesn't work straight away, retry a bit later (probably too many connections open right now)
	while (!defined($superreservedlink))
	{
		eval { our $superreservedlink = DBI -> connect("DBI:mysql:database=$database;host=localhost;mysql_multi_statements=1", 'username', 'password') };
	}
    if ($@)
    {
        warn("Warning: ".$@);
        $superreservedlink = undef;
    }

		if (!defined($superreservedlink))
		{
			sleep(10);
		}
		# else
		# {
		# 	$superreservedlink -> selectdb($database) or $superreservedlink = undef;
		# }
	}
}


sub Connect2
{
	my ($database) = @_;

	instakill();

	if (defined($superreservedlink2))
	{
		# print "\n==================================================================================\n";
		# print   "                ****  WARNING: ALREADY CONNECTED TO MYSQL  ****            ";
		# print "\n==================================================================================\n\n";

		# return 0;

		Disconnect();
	}


	# If connection doesn't work straight away, retry a bit later (probably too many connections open right now)
	while (!defined($superreservedlink2))
	{
		eval { our $superreservedlink2 = DBI -> connect("DBI:mysql:database=$database;host=localhost;mysql_multi_statements=1", 'username', 'password') };
	}
    if ($@)
    {
        warn("Warning: ".$@);
        $superreservedlink2 = undef;
    }

		if (!defined($superreservedlink2))
		{
			sleep(10);
		}
		# else
		# {
		# 	$superreservedlink2 -> selectdb($database) or $superreservedlink2 = undef;
		# }
	}
}

sub Disconnect
{
	our $superreservedlink;

	if (defined($superreservedlink))
	{
		$superreservedlink -> disconnect();
		$superreservedlink = undef;
	}
}

sub Disconnect2
{
	our $superreservedlink2;

	if (defined($superreservedlink2))
	{
		$superreservedlink2 -> disconnect();
		$superreservedlink2 = undef;
	}
}

sub Change
{
	my ($database) = @_;

	instakill();

	Disconnect();

	Connect($database);

	our $superreservedmysqldb = $database;

	# print "\n==================================================================================\n";
	# print   "            ****  WARNING: THIS SCRIPT IS USING AN OLD DATABASE  ****            \n";
	# print   "                              ****  $database  ****                    ";
	# print "\n==================================================================================\n\n";

	warn("\n==================================================================================\n"
	."            ****  WARNING: THIS SCRIPT IS USING AN OLD DATABASE  ****\n"
	."                              ****  $database  ****"
	."\n==================================================================================\n\n");
}

sub Change2
{
	my ($database) = @_;

	instakill();

	Disconnect2();

	Connect2($database);

	our $superreservedmysqldb = $database;

	# print "\n==================================================================================\n";
	# print   "            ****  WARNING: THIS SCRIPT IS USING AN OLD DATABASE  ****            \n";
	# print   "                              ****  $database  ****                    ";
	# print "\n==================================================================================\n\n";

	warn("\n==================================================================================\n"
	."            ****  WARNING: THIS SCRIPT IS USING AN OLD DATABASE  ****\n"
	."                              ****  $database  ****"
	."\n==================================================================================\n\n");
}

sub Query
{
    my ($query) = @_;

    our $superreservedlink;

    instakill();

    print "$query\n" if (defined($superloudmysql));

    # if ($query =~ /^SELECT/)
    # {
        # my $result = $superreservedlink -> query($query);
        $result = $superreservedlink -> prepare($query);

        # Handling for multiple result sets, I hope this works (keep last set)
        # while ($result -> more_results)
        # {
            $result -> execute;
        # }
    # }
    # else
    # {
    #     # This should be faster for UPDATEs and DELETEs
    #     $result = $superreservedlink->do($query);
    # }

    # my $error = $superreservedlink -> errmsg();
    my $error = $superreservedlink -> {'mysql_error'};
    if ($error ne '')
    {
        # die("MySQL error: $error");
        # die("Error: MySQL query failed!\n\n");
        die("Error: MySQL query failed:\n\n$query\n\n");
    }


    return $result;
}

sub Query2
{
    my ($query) = @_;

    our $superreservedlink2;

    instakill();

    print "$query\n" if (defined($superloudmysql));

    # if ($query =~ /^SELECT/)
    # {
        # my $result = $superreservedlink2 -> query($query);
        $result = $superreservedlink2 -> prepare($query);

        # Handling for multiple result sets, I hope this works (keep last set)
        # while ($result -> more_results)
        # {
            $result -> execute;
        # }
    # }
    # else
    # {
    #     # This should be faster for UPDATEs and DELETEs
    #     $result = $superreservedlink2->do($query);
    # }

    # my $error = $superreservedlink2 -> errmsg();
    my $error = $superreservedlink2 -> {'mysql_error'};
    if ($error ne '')
    {
        # die("MySQL error: $error");
        # die("Error: MySQL query failed!\n\n");
        die("Error: MySQL query failed:\n\n$query\n\n");
    }


    return $result;
}

sub Megaquery
{
    # Run a (huge) query without binary logging (so as not to waste disk space)

	my ($query) = @_;

    Query("SET sql_log_bin=0");

    my $result = Query($query);

    Query("SET sql_log_bin=1");

    return $result;
}

sub Megaquery2
{
    # Run a (huge) query without binary logging (so as not to waste disk space)

	my ($query) = @_;

    Query("SET sql_log_bin=0");

    my $result = Query2($query);

    Query("SET sql_log_bin=1");

    return $result;
}

sub Fetch
{
	my ($link) = @_;

	instakill();

	# return $link -> fetchrow;
	return $link -> fetchrow_array;
}

sub FetchOne
{
	my ($link) = @_;

	instakill();

	# if (Numrows($$1)() != 1)
	if ($link -> rows != 1)
	{
		die("Error: Numrows ".$link -> rows." instead of 1\n");
	}

	# return $link -> fetchrow;
	return $link -> fetchrow_array;
}

sub DumpTable
{
	my ($table) = @_;

	state("Dumping table '$table':");

	my $query = Query("SELECT * FROM $table");
	while (my @a = Fetch($query))
	{
		print join("\t", @a)."\n";
	}

	nl();
}

sub DumpTable2
{
	my ($table) = @_;

	state("Dumping table '$table':");

	my $query = Query2("SELECT * FROM $table");
	while (my @a = Fetch($query))
	{
		print join("\t", @a)."\n";
	}

	nl();
}

sub Numrows
{
	my ($link) = @_;

	instakill();

	# return Numrows($$1);
	return $link -> rows;
}

sub Clear
{
	my ($table) = @_;

	instakill();

	# $superreservedlink -> query("TRUNCATE `$table`;");
	$superreservedlink -> do("TRUNCATE $table;");

	print "\nCleared table '$table'\n";
}

sub Clear2
{
	my ($table) = @_;

	instakill();

	# $superreservedlink2 -> query("TRUNCATE `$table`;");
	$superreservedlink2 -> do("TRUNCATE $table;");

	print "\nCleared table '$table'\n";
}


sub Optimize
{
	my ($table, $silent) = @_;

	instakill();

	if (!defined($silent)) { $silent = 0; }

	if ($silent == 0)
	{
		print "Optimizing table '$table'...\n";
	}

	# $superreservedlink -> query("ANALYZE TABLE `$table`;");
	# $superreservedlink -> query("OPTIMIZE TABLE `$table`;");
	$superreservedlink -> do("ANALYZE TABLE $table;");
	$superreservedlink -> do("OPTIMIZE TABLE $table;");

	if ($silent == 0)
	{
		print "Done!\n\n";
	}
}

sub Optimize2
{
	my ($table, $silent) = @_;

	instakill();

	if (!defined($silent)) { $silent = 0; }

	if ($silent == 0)
	{
		print "Optimizing table '$table'...\n";
	}

	# $superreservedlink2 -> query("ANALYZE TABLE `$table`;");
	# $superreservedlink2 -> query("OPTIMIZE TABLE `$table`;");
	$superreservedlink2 -> do("ANALYZE TABLE $table;");
	$superreservedlink2 -> do("OPTIMIZE TABLE $table;");

	if ($silent == 0)
	{
		print "Done!\n\n";
	}
}

sub Exists
{
	my ($table) = @_;
	
	my $query = Query("SHOW TABLES LIKE '$table'");
	
	return Numrows($query);
}

sub Exists2
{
	my ($table) = @_;
	
	my $query = Query2("SHOW TABLES LIKE '$table'");
	
	return Numrows($query);
}




# example loop:
#
# Connect('PTMs');
#
# my $result = Query('SELECT * FROM `expression-human2`');
# while (my %row = $result -> fetchhash)
# {
# 	print $row{name}."\n";
# }


return 1;
