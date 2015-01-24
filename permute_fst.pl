#!/usr/bin/env perl

# ----
print "\n","-"x60,"\n";
print "permute_fst.pl v 1.0 R Capper\n";
print "Last modified: 25 March 2014\n";
print "-"x60,"\n";

unless ($#ARGV == 2)
	{print "this script takes an input file of RAW Fst scores (before smoothing)\n";
	print "produced by vcf2fst_per_locus.pl, which looks like this: \ncontig\tpos\t#pop1_gt\t#pop2_gt\tfst\n";
	print "then permutes the raw values and performs smoothing however many times you'd like.\n";
	print "there is a single output file, '.raw', which has one column of Fsts per permutation performed.\n";
	print "This script should be followed by 'smooth_permuted_fsts.pl' which will then smooth each column\n";
	print "by overlapping window frames.\n";
	print "This script outputs a log file to 'permute.log' so you can keep track of how fast it's going.\n";
	print "This script works like this: for each position, it draws a random Fst value from all available\n";
	print "data (i.e., across the genome) and fills it into the .raw Permute_1 column.  Then again, for the \n";
	print "same contig/position, and again (so, it's writing a row) for as many times as you'd like.  Then,\n";
	print "it moves to the next position and writes that row of 1000 or however many permutations you want.\n";
	print "Compare this method to permuting Fsts and filling in as columns for every position, then moving to\n";
	print "the next permutation.  Same idea, slightly different code.\n";
	print "\nusage: script in.fst number_of_permutations out.raw\n";
	print "example: permute_fst.pl in.fst.tab 1000 KxN_permuted.raw\n";
	exit;}

open(FST,$ARGV[0]);
$reps = $ARGV[1];
open(RAW,">$ARGV[2]");

@name = split(/\./,$ARGV[0]);
$new = $name[0].".permute.log";
open(LOG,">$new");

$contig_count=0;
$total_positions=0;
$pos_count = 0;
$time = localtime();
print LOG "$time: reading in the fst file!\n";

while(<FST>) {chomp; $line = $_; if ($line =~/^contig/) {next;}
	($ctg,$pos,$pop1,$pop2,$fst) = split(/\t/,$line);
	$hash{$ctg}{$pos}++;
	$total_positions++;
	push(@fsts,$fst);}		
$time = localtime(); print LOG "$time: finished hashing the raw Fst scores!\n";

##### ---- making the table of raw permuted values
print RAW "contig\tpos\t";
for ($i = 1; $i <= $reps; $i++)
	{$note = join("_","permute",$i); print RAW "$note\t";}
print RAW "\n";

foreach $c (sort keys %hash)
	{
	$contig_count++;
	foreach $p (sort {$a<=>$b} keys %{$hash{$c}})
		{
		print RAW "$c\t$p\t";
		$pos_count++; $pos_total++;
		$permute_counter = 0;
		until ($permute_counter == $reps)
			{
			#draw a random index of the fst array and print that out (with replacement).
			$random = int(rand($total_positions));	# $random is random value from 0 to $total_positions-1 inclusive, perfect for 
								# selecting random index from fst array
			print RAW "$fsts[$random]\t";
			$permute_counter++;
			}	
		print RAW "\n";
		}
	if ($contig_count % 100 == 0) 
		{$time = localtime();
		print LOG "$time - 100 more contigs processed\n";
#		$total_perms = $pos_count * $reps * 10;
#		print "$total_perms permutations performed as of $time; 10 contigs, $pos_count total positions, $reps permutations each.\n";
#		$pos_count = 0;
		}
	}
$time = localtime();
print LOG "Finished making the RAW FST table!  $pos_total total positions across $contig_count contigs, permuted $reps times; finished at $time\n";
close(RAW);

