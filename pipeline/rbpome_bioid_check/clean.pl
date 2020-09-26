#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
#require('mysql.inc.pl');


# Safety prompt (input files)
$inclean = 'disabled';
# $inclean = 'some';
if ((-d './input/') and ($inclean ne 'disabled'))
{
	$text = "\nAlso remove input files? (y/n)\n";
	print $text;

	$input = <STDIN>;
	chomp($input);
	if (lc($input) eq 'y') { $inclean = 'all'; }


	if ($inclean eq 'all')
	{
		run("Clean up input files", "rm -f input-*");

		run("Clean up input directory", "mv ./input ./deleting_input");
		run("Clean up input directory", "rm -rf ./deleting_input &");
		run("Recreate input directory", "mkdir ./input");
	}
}

# Safety prompt (output files)
$outclean = 'some';
if (-d './output/')
{
	$text = "\nAlso remove output files? (y/n)\n";
	print $text;

	$input = <STDIN>;
	chomp($input);
	if (lc($input) eq 'y') { $outclean = 'all'; }
}


# Start

run("Clean up temporary files", "rm -f tmp-*");

if (-d './tmp/')
{
	run("Clean up temporary directory", "mv ./tmp ./deleting_tmp");
	run("Clean up temporary directory", "rm -rf ./deleting_tmp &");
	run("Recreate temporary directory", "mkdir ./tmp");
}


if (-d './output/')
{
	if ($outclean eq 'all')
	{
		run("Clean up output files", "rm -f output-*");

		run("Clean up output directory", "mv ./output ./deleting_output");
		run("Clean up output directory", "rm -rf ./deleting_output &");
		run("Recreate output directory", "mkdir ./output");
	}
}
else
{
	run("Clean up output files", "rm -f output-*");
}

done();
