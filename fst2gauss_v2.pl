#!/usr/bin/env perl

# -----
print "\n","-"x60, "\n";
print "fst2gauss.pl v 1.0 R Capper\n";
print "Last modified: 28 March 2014\n";
print "-"x60,"\n";

use Bio::Perl;
use Bio::SeqIO;

unless ($#ARGV == 4) {
	print "this script takes the output of vcf2fst_per_locus.pl v 1.0 R. Capper\n";
	print "which looks like this: \ncontig\tpos\t#pop1_gt\t#pop2_gt\tfst\n";
	print "and applies a Gaussian-weighted moving average applied via distance from\n";
	print "some center point.  See Hohenlohe et al. 2010 for the template I'm using,\n";
	print "but basically I multiply each Fst value by exp((-(p-c)^2)/(2*sigma^2))\n";
	print "where p is position from center value, c = center, and sigma = supplied.\n\n";
	print "This script requires the reference genome such that it can find the total length\n";
	print "of the contig being smoothed.  The last window of the contig is moved inwards\n";
	print "such that the step size may be smaller than every other step.\n\n";
	print "Briefly, this script 1) reads in the length of each contig from the ref genome,\n";
	print "2) reads the whole Fst file into memory at once (see note below), then 3)\n";
	print "iteratres over each contig in OVERLAPPING windows.  If there are no SNPs in a given\n";
	print "window, the script outputs 'NA' instead.\n";
	print "Note: A better way to have coded this script would be to read in the fst file\n";
	print "line-by-line and calculate weights on the fly, but my goal was to simply\n";
	print "finish this script quickly and I don't have a huge number of SNPs.  Be careful\n";
	print "if you have many SNPs or many contigs, or many SNPs per contig.  Updating this\n";
	print "script to be less memory intensive is the next version.\n\n";
	print "\nusage: script ref_genome.fasta in.fst.tab sigma step_size out.tab\n\n";
	exit;}

open(FST,$ARGV[1]);
$sigma = $ARGV[2];
$step = $ARGV[3];
open(OUT,">$ARGV[4]");
$contig_count=0;
$time = localtime();

print "$time: reading in the Fst file!\n";
while(<FST>) {chomp; $line=$_; if ($line =~ /^contig/) {next;}
	($ctg,$pos,$pop1,$pop2,$fst) = split(/\t/,$line);
	$hash{$ctg}{$pos}=$fst;}
$time = localtime(); print "$time: finished hashing Fst file!\n";	

$in = Bio::SeqIO->new(-file=>"$ARGV[0]",-format=>"fasta");
while($genome = $in -> next_seq)
	{$name = $genome->id;
	@ct = split(/_/,$name);
	$cc = join("","c",$ct[1]);
	if (exists $hash{$cc}) {$len = $genome -> length; $ref{$cc}=$len;}
	}
$time = localtime(); print "$time: finished hashing genome!\n";



print OUT "contig\tcenter_pos\tweighted_fst\n";
foreach $c (sort keys %hash) {
	$window_count=1;undef(%range);$contig_count++;	

	#define ranges:
	$center = 1; 
	$range{$window_count} = $center;
	until ($center >= $ref{$c}) 
		{$window_count++; $center=$center+$step; $range{$window_count}=$center;}

	foreach $win (sort {$a<=>$b} keys %range) 
		{
		$numer=0;$denom=0;
		undef(@window_fst);undef(@denominator);
		$start = $range{$win} - 3*$sigma; 
		$end = $range{$win} + 3*$sigma -1;
		if ($start < 1) {$start = 1;}
		if ($end > $ref{$c}) {$end = $ref{$c};}
		if ($range{$win} > $ref{$c}) {$range{$win} = $ref{$c};}
#		print "contig $c, window $win, center $range{$win}: range:\t$start\t$end\n";

		foreach $snp (sort {$a<=>$b} keys %{$hash{$c}})
			{ 
			if (($snp >= $start) && ($snp <= $end))
				{
				$weight = exp(-($snp-$range{$win})**2/(2*$sigma**2));
#				print "\t$snp is between $start and $end; fst = $hash{$c}{$snp}\tweight=$weight\n";
				$weighted_fst = $hash{$c}{$snp}*$weight;
				push(@window_fst,$weighted_fst);
				push(@denominator,$weight);			
				}
			}
		if (! @denominator) {print OUT "$c\t$range{$win}\tNA\n"; next;}
		$numer += $_ for @window_fst;
		$denom += $_ for @denominator;
		$avg = $numer/$denom;
#		print "\tcenter = $range{$win}, weighted fst = $avg\n";
		print OUT "$c\t$range{$win}\t$avg\n";

		}
	}	

$time_end = localtime();

print "$time_end: finished calculating Gaussian-weighted Fst values for all contigs!\n";

