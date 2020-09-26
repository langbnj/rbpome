#!/users/gt/blang/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = "rbpome";
$scoretype = 'SUM_IS';
$threshold = 7.1;
$h47_threshold = 1.9;

our $usage = "$0";
# ($var) = args(1);
args(0);

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

starttime();


$query = Query("SELECT COUNT(*) FROM rbpome");
($before) = FetchOne($query);

$query = Query("SELECT COUNT(DISTINCT symbol1, symbol2) FROM rbpome");
($uniquepairs_before) = FetchOne($query);

# $scoretypequery = Query("SELECT DISTINCT scoretype FROM rbpome ORDER BY scoretype");
# while (($scoretype) = Fetch($scoretypequery))
# {
	# Filter out duplicate reverse pairs (retaining the configuration with the higher score)
	
	# First, get a list of pairs (flipped so the first is alphanumerically lower)
	$mainquery = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype'");
	@pairs = ();
	startme("$scoretype: Getting all pairs", 0, Numrows($mainquery));
	while (($symbol1, $symbol2) = Fetch($mainquery))
	{
		# # flip to make the pairs alphabetic
		# if ($symbol1 gt $symbol2)
		# {
		# 	push(@pairs, "$symbol2|$symbol1");
		# }
		# else
		# {
		# 	push(@pairs, "$symbol1|$symbol2");
		# }
		#
		
		# Retain original orientation
		push(@pairs, "$symbol1|$symbol2");
		
		stepme(1000);
	}
	stopme();
	@pairs = unique(@pairs);
	$flipped_uniquepairs_before = scalar(@pairs);
	state("$scoretype: Got ".commify(scalar(@pairs))." unique pairs");
	


	# Then, go through the list and check if there are table entries that are the wrong way around. If yes, flip them. If there are reverse duplicates, remove the lower-scoring one.
	startme("$scoretype: Merging duplicate pairs (retaining highest scores and SUMing or MAXing information)...", 0, scalar(@pairs));
	foreach $pair (@pairs)
	{
		($symbol1, $symbol2) = split(/\|/, $pair);
	
		addme("total unique 'correct' pairs before", $pair);
	
		# Initial query
		$query = Query("SELECT id, symbol1, symbol2, fullsymbol1, fullsymbol2, species1, species2, accs1, accs2, names1, names2, same_species, times_detected, bg_alldirect, bg_all, found_inverse, eclip1, eclip2, hc, unique_ab_avg_is, avg_is, avg_rs, avg_ris, h47_avg_is, h47_sum_is, h47_times_detected, h47_found_in_both_orientations, h47_avg_rs, h47_avg_ris, hippie, huri, nanobret, nanobret_mbu_avg, nanobret_mbu_sd FROM $table WHERE ((symbol1='$symbol1' AND symbol2='$symbol2') OR (symbol1='$symbol2' AND symbol2='$symbol1')) AND scoretype='$scoretype'");
		if (Numrows($query) == 0)
		{
			die;
		}
		# elsif (Numrows($query) == 1)
		# {
		# 	addme("exactly one row existed for pair (unchanged)", $pair);
		# }
		# elsif (Numrows($query) > 1)
		else
		{
			addme("more than one row existed for pair (merged)", $pair);
			
			@our_ids = ();
			$our_max_id = undef;
			$our_symbol1 = undef;
			$our_symbol2 = undef;
			$our_fullsymbol1 = undef;
			$our_fullsymbol2 = undef;
			$our_species1 = undef;
			$our_species2 = undef;
			$our_accs1 = undef;
			$our_accs2 = undef;
			$our_names1 = undef;
			$our_names2 = undef;
			$our_avg_is = undef;
			$our_unique_ab_avg_is = undef;
			$our_avg_rs = undef;
			$our_avg_ris = undef;
			$our_h47_avg_is = undef;
			$our_h47_sum_is = undef;
			$our_h47_avg_rs = undef;
			$our_h47_avg_ris = undef;
			$our_correctpair = 0;
			$our_inversepair = 0;
			$our_found_inverse = 0;
			$our_found_inverse_max = 0;
			$our_h47_correctpair = '';
			$our_h47_inversepair = '';
			$our_same_species = 0;
			$our_times_detected = 0;
			$our_bg_alldirect = 0;
			$our_bg_all = 0;
			$our_homodimer = 0;
			$our_eclip1 = 0;
			$our_eclip2 = 0;
			$our_hc = 0;
			$our_h47_found_in_both_orientations = '';
			$our_h47_found_in_both_orientations_max = '';
			$our_h47_times_detected = '';
			$our_hippie = 0;
			$our_huri = 0;
			$our_nanobret = '';
			$our_nanobret_mbu_avg = '';
			$our_nanobret_mbu_sd = '';
			while (($this_id, $this_symbol1, $this_symbol2, $this_fullsymbol1, $this_fullsymbol2, $this_species1, $this_species2, $this_accs1, $this_accs2, $this_names1, $this_names2, $this_same_species, $this_times_detected, $this_bg_alldirect, $this_bg_all, $this_found_inverse, $this_eclip1, $this_eclip2, $this_hc, $this_unique_ab_avg_is, $this_avg_is, $this_avg_rs, $this_avg_ris, $this_h47_avg_is, $this_h47_sum_is, $this_h47_times_detected, $this_h47_found_in_both_orientations, $this_h47_avg_rs, $this_h47_avg_ris, $this_hippie, $this_huri, $this_nanobret, $this_nanobret_mbu_avg, $this_nanobret_mbu_sd) = Fetch($query))
			{
				push(@our_ids, $this_id);
			
				# Pair flipping
				if ($this_symbol1 gt $this_symbol2)
				{
					if ($this_avg_is >= $threshold)
					{
						# "inverse pair" is present (for found_inverse)
						$our_inversepair = 1;
					}
				
					if (defined($this_h47_sum_is))
					{
						# "inverse pair" was tested for h47 (for h47_found_inverse)
						$our_h47_inversepair = 0;
						if ($this_h47_sum_is >= $h47_threshold)
						{
							# "inverse pair" is present for h47 (for h47_found_inverse)
							$our_h47_inversepair = 1;
						}
					}

					# # Flip
					# ($this_symbol1, $this_symbol2) = ($this_symbol2, $this_symbol1);
					# ($this_fullsymbol1, $this_fullsymbol2) = ($this_fullsymbol2, $this_fullsymbol1);
					# ($this_species1, $this_species2) = ($this_species2, $this_species1);
					# ($this_accs1, $this_accs2) = ($this_accs2, $this_accs1);
					# ($this_names1, $this_names2) = ($this_names2, $this_names1);
					# ($this_eclip1, $this_eclip2) = ($this_eclip2, $this_eclip1);
				}
				else
				{
					if ($this_avg_is >= $threshold)
					{
						# "correct pair" is present (for $found_inverse)
						$our_correctpair = 1;
					}

					if (defined($this_h47_sum_is))
					{
						# "correct pair" was tested for h47 (for h47_found_inverse)
						$our_h47_correctpair = 0;
						if ($this_h47_sum_is >= $h47_threshold)
						{
							# "correct pair" is present for h47 (for h47_found_inverse)
							$our_h47_correctpair = 1;
						}
					}
				}
			
			
			
				# Get MAX scores
				if (!defined($our_avg_is) or ($our_avg_is < $this_avg_is))
				{
					$our_avg_is = $this_avg_is;
				
					# This is the row to retain (it'll be completely updated, though)
					$our_max_id = $this_id;
				
					$our_symbol1 = $this_symbol1;
					$our_symbol2 = $this_symbol2;
					$our_fullsymbol1 = $this_fullsymbol1;
					$our_fullsymbol2 = $this_fullsymbol2;
					$our_species1 = $this_species1;
					$our_species2 = $this_species2;
					$our_accs1 = $this_accs1;
					$our_accs2 = $this_accs2;
					$our_names1 = $this_names1;
					$our_names2 = $this_names2;
					$our_eclip1 = $this_eclip1;
					$our_eclip2 = $this_eclip2;
				}
				$our_unique_ab_avg_is = $this_unique_ab_avg_is if (!defined($our_unique_ab_avg_is) or ($our_unique_ab_avg_is < $this_unique_ab_avg_is));
				$our_avg_rs = $this_avg_rs if (!defined($our_avg_rs) or ($our_avg_rs < $this_avg_rs));
				$our_avg_ris = $this_avg_ris if (!defined($our_avg_ris) or ($our_avg_ris < $this_avg_ris));
				$our_h47_avg_is = $this_h47_avg_is if (defined($this_h47_avg_is) and (!defined($our_h47_avg_is) or ($our_h47_avg_is < $this_h47_avg_is)));
				$our_h47_sum_is = $this_h47_sum_is if (defined($this_h47_sum_is) and (!defined($our_h47_sum_is) or ($our_h47_sum_is < $this_h47_sum_is)));
				$our_h47_avg_rs = $this_h47_avg_rs if (defined($this_h47_avg_rs) and (!defined($our_h47_avg_rs) or ($our_h47_avg_rs < $this_h47_avg_rs)));
				$our_h47_avg_ris = $this_h47_avg_ris if (defined($this_h47_avg_ris) and (!defined($our_h47_avg_ris) or ($our_h47_avg_ris < $this_h47_avg_ris)));
			
			

				# # Deal with NULL values
				# if (!defined($this_avg_is)) { $this_avg_is = ''; }
				# if (!defined($this_avg_ris)) { $this_avg_ris = ''; }
				# if (!defined($this_avg_rs)) { $this_avg_rs = ''; }
				# if (!defined($this_h47_avg_is)) { $this_h47_avg_is = ''; }
				# if (!defined($this_h47_avg_ris)) { $this_h47_avg_ris = ''; }
				# if (!defined($this_h47_avg_rs)) { $this_h47_avg_rs = ''; }
				# if (!defined($this_h47_sum_is)) { $this_h47_sum_is = ''; }
				# if (!defined($this_h47_found_in_both_orientations)) { $this_h47_found_in_both_orientations = ''; }
				# if (!defined($this_h47_times_detected)) { $this_h47_times_detected = ''; }
				# if (!defined($this_nanobret)) { $this_nanobret = ''; }
				# if (!defined($this_nanobret_mbu_avg)) { $this_nanobret_mbu_avg = ''; }
				# if (!defined($this_nanobret_mbu_sd)) { $this_nanobret_mbu_sd = ''; }
				# if (!defined($this_unique_ab_avg_is)) { $this_unique_ab_avg_is = ''; }
			
			
			
				# Do MAX/SUM integration of values across pair entries
				$our_same_species = 1 if ($this_same_species == 1);	# MAX
				$our_times_detected += $this_times_detected;	# SUM
				$our_bg_alldirect = 1 if ($this_bg_alldirect == 1);	# MAX
				$our_bg_all = 1 if ($this_bg_all == 1);	# MAX
				# $our_homodimer = 1 if (($this_symbol1 eq $this_symbol2) or ($this_homodimer == 1)); # MAX
				# die("Error: 'homodimer' in table 'rbpome' is '$this_homodimer', but gene symbols are '$this_symbol1' and '$this_symbol2'") if (($this_symbol1 eq $this_symbol2) and ($this_homodimer != 1));
				# die("Error: 'homodimer' in table 'rbpome' is '$this_homodimer', but gene symbols are '$this_symbol1' and '$this_symbol2'") if (($this_symbol1 ne $this_symbol2) and ($this_homodimer == 1));
				$our_homodimer = 1 if ($this_symbol1 eq $this_symbol2); # MAX
				$our_found_inverse_max = 1 if ($this_found_inverse == 1);	# MAX	 # There doesn't seem to be any threshold checking for this column initially, don't trust it. Just keeping it for diagnosis
				# $our_eclip1 = 1 if ($this_eclip1 == 1);
				# $our_eclip2 = 1 if ($this_eclip2 == 1);
				$our_hc = 1 if ($this_hc == 1);
				die if (($this_hc == 1) and ($this_avg_is < $threshold));
				# $our_$h47_found_in_both_orientations = 1 if ($this_h47_found_in_both_orientations == 1);	# MAX	 # No, there's no threshold checking for this column initially. Don't trust it
				if (defined($this_h47_times_detected))
				{
					if ($our_h47_times_detected ne '')
					{
						$our_h47_times_detected += $this_h47_times_detected; # SUM
					}
					else
					{
						$our_h47_times_detected = $this_h47_times_detected; # SUM
					}
				}
				if (defined($this_h47_found_in_both_orientations))
				{
					if ($our_h47_found_in_both_orientations_max eq '')
					{
						$our_h47_found_in_both_orientations_max = $this_h47_found_in_both_orientations;	# MAX	 # There doesn't seem to be any threshold checking for this column initially, don't trust it. Just keeping it for diagnosis
					}
					else
					{
						$our_h47_found_in_both_orientations_max = $this_h47_found_in_both_orientations if ($this_h47_found_in_both_orientations == 1);
					}
				}
				$our_hippie = 1 if ($this_hippie == 1);
				$our_huri = 1 if ($this_huri == 1);
			
				if (defined($this_nanobret))
				{
					$our_nanobret = 'pos' if ($this_nanobret eq 'pos');	# pos overrides neg
					$our_nanobret = 'neg' if (($this_nanobret eq 'neg') and ($our_nanobret eq ''));	# neg only overrides ''
					if (($our_nanobret_mbu_avg eq '') or ($this_nanobret_mbu_avg > $our_nanobret_mbu_avg))
					{
						$our_nanobret_mbu_avg = $this_nanobret_mbu_avg; # MAX
						$our_nanobret_mbu_sd = $this_nanobret_mbu_sd; # use the sd value from the maximum
					}
				}
			}
		
			# Set found_inverse
			if (($our_correctpair == 1) and ($our_inversepair == 1))
			{
				$our_found_inverse = 1;
			}
			# Set h47_found_inverse
			if (($our_h47_correctpair eq '1') and ($our_h47_inversepair eq '1'))
			{
				# If both the correct and the inverse pair were found, set it to 1
				$our_h47_found_in_both_orientations = 1;
			}
			elsif (($our_h47_correctpair ne '') or ($our_h47_inversepair ne ''))
			{
				# If only one was found, set it to 0. Otherwise it remains blank ('')
				$our_h47_found_in_both_orientations = 0;
			}
			
			

			# Check found_inverse_max (diagnostic only)
			if ($our_found_inverse ne $our_found_inverse_max)
			{
				addme("mismatch between found_inverse and found_inverse_max for pair (kept my calculated found_inverse)", $pair);
				addme("mismatch between found_inverse and found_inverse_max for found_inverse|found_inverse_max (kept my calculated found_inverse)", "$our_found_inverse|$our_found_inverse_max");
			}
			# Use max of the two (sometimes mine was 0 while the previous one was 1, so I was losing found_inverse pairs)
			# print "\n\nOUR_FOUND_INVERSE / OUR_FOUND_INVERSE_MAX / OUR_FOUND_INVERSE NEW = $our_found_inverse / $our_found_inverse_max / ";
			$our_found_inverse = max($our_found_inverse, $our_found_inverse_max);
			# print "$our_found_inverse <<< \n\n";

			# Check h47_found_in_both_orientations_max (diagnostic only)
			# >> This shows that whenever they differ, I have 0 while _max (the previous data in the table) has 1. Never the other way around. I'll use $our_h47_found_in_both_orientations_max below since it has information I don't have in the table.
			if ($our_h47_found_in_both_orientations ne $our_h47_found_in_both_orientations_max)
			{
				addme("mismatch between h47_found_in_both_orientations and h47_found_in_both_orientations_max for pair (used h47_found_in_both_orientations_max)", $pair);
				addme("mismatch between h47_found_in_both_orientations and h47_found_in_both_orientations_max for h47_found_in_both_orientations|h47_found_in_both_orientations_max (used h47_found_in_both_orientations_max)", "$our_h47_found_in_both_orientations|$our_h47_found_in_both_orientations_max");
			}
			# Use $our_h47_found_in_both_orientations_max
			$our_h47_found_in_both_orientations = $our_h47_found_in_both_orientations_max;
		
		
		
			# Update table
		
			# DELETE everything but our_max_id (it'll be updated to contain the merged values)
			foreach $this_id (@our_ids)
			{
				if ($this_id == $our_max_id)
				{
					if ($our_homodimer == 0) 
					{
						# Keep heterodimers
						addme("kept row for id", $this_id);
						addme("kept one row for pair", $pair);
						next;
					}
					else
					{
						# Delete homodimers
						addme("deleted homodimer row for id", $this_id);
						addme("deleted homodimer for pair", $pair);
					}
				}
			
				$q = "DELETE FROM $table WHERE id=$this_id";
				addme("deleted row for id", $this_id);
				addme("deleted at least one row for pair", $pair);
			
				if (switch('debug'))
				{
					state($q);
				}
				else
				{
					$query = Query($q);
					die("Error: Numrows!=1 for query '$q'") if (Numrows($query) != 1);
				}
			}
			
			# Handle NULLs
			$our_accs1 = '' if (!defined($our_accs1));
			$our_accs2 = '' if (!defined($our_accs2));
			$our_names1 = '' if (!defined($our_names1));
			$our_names2 = '' if (!defined($our_names2));
			$our_h47_avg_is = '' if (!defined($our_h47_avg_is));
			$our_h47_sum_is = '' if (!defined($our_h47_sum_is));
			$our_h47_avg_rs = '' if (!defined($our_h47_avg_rs));
			$our_h47_avg_ris = '' if (!defined($our_h47_avg_ris));
			


			# die("Error: Pair is in incorrect orientation: $our_symbol1|$our_symbol2 (shouldn't happen at this point)") if ($our_symbol1 gt $our_symbol2);
			

			
			# Update heterodimers
			# If it's a homodimer: already DELETEd above, nothing to UPDATE here
			if ($our_homodimer == 0)
			{
				# UPDATE our_max_id (to contain the merged values)
				$q = "UPDATE $table SET
				symbol1='$our_symbol1', 
				symbol2='$our_symbol2', 
				fullsymbol1='$our_fullsymbol1', 
				fullsymbol2='$our_fullsymbol2', 
				species1='$our_species1', 
				species2='$our_species2', 
				accs1='$our_accs1', 
				accs2='$our_accs2', 
				names1='$our_names1', 
				names2='$our_names2', 
				avg_is='$our_avg_is', 
				unique_ab_avg_is='$our_unique_ab_avg_is', 
				avg_rs='$our_avg_rs', 
				avg_ris='$our_avg_ris', 
				h47_avg_is='$our_h47_avg_is', 
				h47_sum_is='$our_h47_sum_is', 
				h47_avg_rs='$our_h47_avg_rs', 
				h47_avg_ris='$our_h47_avg_ris', 
				found_inverse='$our_found_inverse', 
				same_species='$our_same_species', 
				times_detected='$our_times_detected', 
				bg_alldirect='$our_bg_alldirect', 
				bg_all='$our_bg_all', 
				homodimer='$our_homodimer', 
				eclip1='$our_eclip1',
				eclip2='$our_eclip2',
				hc='$our_hc', 
				h47_found_in_both_orientations='$our_h47_found_in_both_orientations', 
				h47_times_detected='$our_h47_times_detected', 
				hippie='$our_hippie', 
				huri='$our_huri', 
				nanobret='$our_nanobret', 
				nanobret_mbu_avg='$our_nanobret_mbu_avg', 
				nanobret_mbu_sd='$our_nanobret_mbu_sd'
				WHERE id=$our_max_id";

				$q =~ s/=''/=NULL/g;
		
				addme("updated row for id", $our_max_id);
				addme("updated one row for pair", $pair);

				if (switch('debug'))
				{
					# state($q);
					$query = Query("SELECT found_inverse FROM rbpome WHERE id=$our_max_id");
					($tmp_found_inverse) = FetchOne($query);
					die("Weird: found_inverse is being set from 1 to 0 for id '$our_max_id' $tmp_found_inverse -> $our_found_inverse") if ($tmp_found_inverse > $our_found_inverse);
				}
				else
				{
					$query = Query($q);
					die("Error: Numrows!=1 for query '$q'") if (Numrows($query) != 1);
				}
			}
		}
		# else
		# {
		# 	die;
		# }
		
		stepme(100);
	}
# }
stopme();
stoptime();

# showmeall(1);
showmesome(10);

$query = Query("SELECT COUNT(*) FROM rbpome");
($after) = FetchOne($query);

$query = Query("SELECT COUNT(DISTINCT symbol1, symbol2) FROM rbpome");
($uniquepairs_after) = FetchOne($query);

# Get list of pairs (flipped so the first is alphanumerically lower)
$mainquery = Query("SELECT DISTINCT symbol1, symbol2 FROM $table WHERE scoretype='$scoretype'");
@tmppairs = ();
while (($symbol1, $symbol2) = Fetch($mainquery))
{
	if ($symbol1 gt $symbol2)
	{
		push(@tmppairs, "$symbol2|$symbol1");
	}
	else
	{
		push(@tmppairs, "$symbol1|$symbol2");
	}
}
@tmppairs = unique(@tmppairs);
$flipped_uniquepairs_after = scalar(@tmppairs);

nl();
state("Total rows in the table before: ".commify($before), 1);
state("Total rows in the table now:    ".commify($after), 1);
state("Total rows difference:          ".commify($after - $before), 1);
nl();
state("Total 'unflipped' unique pairs in the table before: ".commify($uniquepairs_before), 1);
state("Total 'unflipped' unique pairs in the table now:    ".commify($uniquepairs_after), 1);
state("Total 'unflipped' unique pairs difference:          ".commify($uniquepairs_after - $uniquepairs_before), 1);
nl();
state("Total 'flipped' (correct) unique pairs in the table before:                            ".commify($flipped_uniquepairs_before), 1);
state("Total 'flipped' (correct) unique pairs in the table now:                               ".commify($flipped_uniquepairs_after), 1);
state("Total 'flipped' (correct) unique pairs difference (probably list homodimers):          ".commify($flipped_uniquepairs_after - $flipped_uniquepairs_before), 1);
nl();

done();

Optimize($table);
