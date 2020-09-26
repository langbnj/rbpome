#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
args(0);


# start

# @fimosources = ('test');
@fimosources = ('attract', 'cisbp', 'dominguez', 'rbpdb', 'rbpmap', 'rnacompete');
@logoformats = ('eps', 'png');
foreach $logoformat (@logoformats)
{
	run("Remove logo directory", "rm -rf output_logos_$logoformat");
	run("Recreate logo directory", "mkdir -p output_logos_$logoformat");

	foreach $fimosource (@fimosources)
	{
		# Plot all input PWMs individually using meme2images
		run("Remove logo directory", "rm -rf output_logos_$logoformat\_$fimosource");
		run("Recreate logo directory", "mkdir -p output_logos_$logoformat\_$fimosource");
		run("meme2images", "meme2images -$logoformat input/$fimosource.meme output_logos_$logoformat\_$fimosource");
		
		cd("output_logos_$logoformat\_$fimosource");
		run("Rename", "rename logo $fimosource\_ *");
		cd("..");
		
		run("Move", "mv output_logos_$logoformat\_$fimosource/* output_logos_$logoformat");

		run("Remove logo directory", "rm -rf output_logos_$logoformat\_$fimosource");
	}
}

done();
