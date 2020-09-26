#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$infile = "input-autoactivators.txt";
$aliasfile = "input-aliases.txt";

open(IN, $infile);
open(ALIAS, $aliasfile) or die("\nError: Couldn't open '$aliasfile'\n\n");







startme("Reading gene symbol alias mapping from '$aliasfile'");
while (<ALIAS>)
{
	# # Alias definitions
	#
	# # I'm going to follow UniProt, not HGNC. HGNC still prefers "C1orf35" over MMTAG2, even though it lists MMTAG2 as a synonym, which is really strange.
	#
	#
	#
	#
	#
	#
	#
	# # HUMAN
	# # Used this type of query to get these: https://www.uniprot.org/uniprot/?query=DBC1+AND+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
	#
	#
	#
	# HUMAN	ATP5B	ATP5F1B
	# HUMAN	ATP5C1	ATP5F1C

	chomp;

	# Skip empty lines
	next if /^$/;

	# Skip comments
	next if /^#/;

	@a = split(/\t/);
	die("Error: Line doesn't contain 3 fields:\n\n$_\n\n") if (scalar(@a) != 3);

	# ($species, $original, $alias) = @a;
	$original = $a[1];
	$alias = $a[2];

	# Convert to upper case
	# $species = uc($species);
	$original = uc($original);
	$alias = uc($alias);

	# $alias{"$species|$original"} = $alias;
	# Discard species
	$alias{$original} = $alias;

	stepme(10);
}
stopme();
close(ALIAS);







@autoactivators = ();
while (<IN>)
{
	chomp;
	
	# Skip comments and empty lines
	next if /^#/;
	next if /^$/;
	
	$symbol1 = $_;
	

	# Make symbol all upper case
	$symbol1 = uc($symbol1);
	
	# # Make species all upper case
	# $species = uc($species);
	
	# # Replace spaces in "RBM10 316F" >> "rbm10_316f", "eMIC SRRM3 mut 25" >> "emic_srrm3_mut_25" etc.
	# $symbol =~ s/ /_/g;
	# Remove spaces in "RBM10 316F" >> "rbm10316f", "eMIC SRRM3 mut 25" >> "emicsrrm3mut25" etc. (as in Sebastian's table)
	$symbol1 =~ s/ //g;

	
	# Apply alias mapping (e.g. ATP5B protein to ATP5F1B, which is its gene)
	if (exists($alias{$symbol1}))
	# if (exists($alias{"$species1|$symbol1"}))
	{
		# addme("aliased from species|symbol", "$species1|$symbol1");
		addme("aliased from symbol1", $symbol1);
		$symbol1 = $alias{$symbol1};
		# $symbol1 = $alias{"$species1|$symbol1"};
		# addme("aliased to species|symbol", "$species1|$symbol1");
		addme("aliased from symbol1", $symbol1);
	}
	

	push(@autoactivators, $symbol1);
	addme("total autoactivators", $symbol1);
}
@autoactivators = unique(@autoactivators);
# @autoactivators = nsort(@autoactivators);

# foreach (@autoactivators) {print "$_\n"};
# exit;



# # Filter out autoactivators as baits (where Protein A is an autoactivator)
# # The 10 AAs from Jae-Seongs figure (top right panel) (Fig1.autoactivator.v1.ai in Gmail):
# # CD2BP2
# # Fez2
# # HNRNPF
# # KPNA1
# # Psmd4
# # Sf3b2
# # TFIP11
# # Tacc3
# # WBP11
# # Zc3h7b
# # The sticky ones, bottom right, will stay in
# @autoactivators = ();
# push(@autoactivators, "CD2BP2");
# push(@autoactivators, "Fez2");
# push(@autoactivators, "HNRNPF");
# push(@autoactivators, "KPNA1");
# push(@autoactivators, "Psmd4");
# push(@autoactivators, "Sf3b2");
# push(@autoactivators, "TFIP11");
# push(@autoactivators, "Tacc3");
# push(@autoactivators, "WBP11");
# push(@autoactivators, "Zc3h7b");
#
# # # Currently just:
# # # Filter out RBM5-NONO
# # # 	RBM5 is an auto-activator (the only one with eCLIP data)
# # # 		RBM5-NONO wasn’t found inversely, so it’s a false positive
# # # 	Filter it out in ~/update/rbpome in an autoactivator removal step that’s late enough so found_inverse is set correctly (after deduplicate.pl)
# # # 		Then later, Jae-Seong will send the final list of autoactivators and I can create a clean master table
# # $query = Query("DELETE FROM rbpome WHERE symbol1='NONO' AND symbol2='RBM5'");
# # state("Deleted ".Numrows($query)." autoactivator interactions");



$deleted = 0;
foreach $symbol1 (@autoactivators)
{
	$query = Query("DELETE FROM rbpome WHERE symbol1='$symbol1'");
	$deleted += Numrows($query);
}
state("Deleted $deleted interactions with autoactivators as baits (= symbol1 = Protein A)");

showmeall(1);

done();
