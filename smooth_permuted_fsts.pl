#!/usr/bin/env perl

use Bio::Perl;
use Bio::SeqIO;

# -----
print "\n","-"x60, "\n";
print "smooth_permuted_fsts.pl v 1.0\m";
print "Roxana Capper, Univ of Texas at Austin\n";
print "Last modified: 26 March 2014\n";
print "-"x60,"\n";

unless ($#ARGV == 4)
	{print "This script is meant to follow 'permute_fst.pl', which generates a giant table of\n";
	print "contig\tpos\tfst_permute1\tfst_permute2 ...\n";
	print "This script reads in that raw file contig by contig and calculates gaussian-smoothed\n";
	print "Fst values using overlapping windows with the same weighting scheme as used in Hohenlohe (2010)\n";
	print "Smoothing occurs within each column.  \nNote: make sure to use the same sigma and step size that you used in the 'fst2gauss_v2.pl' file you used to calculate the initial Fst values.\n";
	print "Note2: yeah, this requires the reference genome just to calculate the length of each contig.  Lame, I know.\n";
	print "Note3: updates are sent to an automatically generated LOG file.  Check it out: file 'prefix.smooth.log'\n\n";
	print "usage: script in.permuted.fst.tab ref_genome.fasta sigma step_size smoothed.fst.tab\n";
	exit;}

open(RAW,$ARGV[0]);
$sigma = $ARGV[2];
$step = $ARGV[3];
open(SMOOTH,">$ARGV[4]");
@name = split(/\./,$ARGV[0]); $new = $name[0].".smooth.log";
open(LOG,">$new");
$time = localtime();
$contig_count=0;

print "$time: reading in the ref genome!\n";
$in = Bio::SeqIO->new(-file=>"$ARGV[1]",-format=>"fasta");
while($genome = $in->next_seq)
	{$name = $genome->id;
	@ct = split(/_/,$name);
	$cc = join("","c",$ct[1]);
	$len = $genome -> length;
	$ref{$cc} = $len;}

undef(%contig);
undef(%range);

$time = localtime();
print "$time: starting smoothing within permutation columns.\n";

while(<RAW>) 
	{
	chomp; $line = $_; 
	@line = split(/\t/,$line);
	if ($line[0] =~/^contig/) {		### Make the header line for arbitrary length of permutations...
		print SMOOTH "contig\tpos\t";
		foreach $number (1..$#line-1) {$note = join("_","smoothed",$number); print SMOOTH "$note\t";} print SMOOTH "\n";
		next;
		}	
	$ctg = $line[0]; $pos = $line[1];

	if (! exists $contig{$ctg})		### If it is a new contig
		{ 
		$contig_count++;
		if (! defined %range)		### If it is the first record in the file; only then will NEITHER hash be initialized
			{
			print LOG "*** first record in the vcf is for contig $ctg and position $pos.\n";
			}

		if (defined %range)		### If it is just a new contig; calculate windows for OLD contig
			{
#			print LOG "1. starting a new contig.  ";
			foreach $ctr (sort {$a<=>$b} keys %range)
				{
				$start = $range{$ctr} - 3*$sigma;
				$end = $range{$ctr} + 3*$sigma - 1;
				if ($start < 1) {$start = 1;}
				if ($end > $ref{$old_contig}) {$end = $ref{$old_contig};}
				if ($range{$ctr} > $ref{$old_contig}) {$range{$ctr} = $ref{$old_contig};}	### if center + step is longer than the contig
				print SMOOTH "$old_contig\t$range{$ctr}\t";

				foreach $permut (sort {$a<=>$b} keys %{$hash{$old_contig}})
					{
					undef(@window_fst); undef(@denominator);$numer=0;$denom=0;
					foreach $snp (sort {$a<=>$b} keys %{$hash{$old_contig}{$permut}})
						{
						if (($snp >= $start) && ($snp <= $end))
							{
#							print "$ctg: window $ctr: $start - $range{$ctr} - $end\n";
							$weight = exp(-($snp-$range{$ctr})**2/(2*$sigma**2));
							$weighted_fst = $hash{$old_contig}{$permut}{$snp}*$weight;
							push(@window_fst,$weighted_fst);
							push(@denominator,$weight);
							}
						}
					if (! @window_fst) {print SMOOTH "NA\t";next;}
					$numer += $_ for @window_fst;
					$denom += $_ for @denominator;
#					print "$numer over $denom\n";
					$avg = $numer/$denom;
					print SMOOTH "$avg\t";
					}
				print SMOOTH "\n";
				}			
			$time = localtime();
#			print LOG "The last contig processed was $old_contig; finished at $time.\n";
			}

		### Reset arrays:
		undef(%contig);		### clear out all the old information to keep the hash "small"
		$contig{$ctg}++;	### initialize with new contig 
		$old_contig = $ctg;	### keep track of the new contig such that when you hit the NEXT new contig you can recall this
		undef(%range);		### clear all the old windows
		undef(%hash);		### clear out all the permutations for prev contig

		$window_count = 1;	### Set up all the new windows
		$center = 1;
		$range{$window_count}=$center;
		until ($center >= $ref{$ctg})
			{$window_count++; $center=$center+$step; $range{$window_count}=$center;}

		foreach $permutation (1..$#line-1)
			{
			$index = $permutation+1;
			$hash{$ctg}{$permutation}{$pos} = $line[$index];
#			print "$ctg\t$permutation\t$pos: $line[$index]\n";
			}
		if ($contig_count % 100 == 0) {$time = localtime(); print LOG "$time = 100 more contigs processed\n";}
		next;
		}

	if (exists $contig{$ctg})
		{
		$contig{$ctg}++;
                foreach $permutation (1..$#line-1)
			{
			$index = $permutation+1;
			$hash{$ctg}{$permutation}{$pos} = $line[$index];
#			print "$ctg\t$permutation\t$pos: $line[$index]\n";
			}
		}
	}
print LOG "last record: $ctg\t$pos\n";
foreach $ctr (sort {$a<=>$b} keys %range)
	{
	$start = $range{$ctr} - 3*$sigma;
	$end = $range{$ctr} + 3*$sigma - 1;
	if ($start < 1) {$start = 1;}
	if ($end > $ref{$old_contig}) {$end = $ref{$old_contig};}
	if ($range{$ctr} > $ref{$old_contig}) {$range{$ctr} = $ref{$old_contig};}
	print SMOOTH "$old_contig\t$range{$ctr}\t";
	foreach $permut (sort {$a<=>$b} keys %{$hash{$old_contig}})
		{
		undef(@window_fst); undef(@denominator);$numer=0;$denom=0;
		foreach $snp (sort {$a<=>$b} keys %{$hash{$old_contig}{$permut}})
			{
			if (($snp >= $start) && ($snp <= $end))
				{
				$weight = exp(-($snp-$range{$ctr})**2/(2*$sigma**2));
				$weighted_fst = $hash{$old_contig}{$permut}{$snp}*$weight;
				push(@window_fst,$weighted_fst);
				push(@denominator,$weight);
				}
			}
		if (! @window_fst) {print SMOOTH "NA\t";}
		$numer += $_ for @window_fst;
		$denom += $_ for @denominator;
		$avg = $numer/$denom;
		print SMOOTH "$avg\t";                         
		}
	}
$time = localtime();
print "$time: finished all contigs, all windows, all permutations.\n";	




