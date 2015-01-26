#!/usr/bin/env perl

# -- program name, header
print "-"x60, "\n";
print "fst_pvalues.pl v 1.0 R Capper\n";
print "Last modified: 28 March 2012\n";
print "-"x60,"\n";

unless ($#ARGV == 2)
        {
	print "this script takes the output from fst2gauss_v2.pl and from smooth_permuted_fsts.pl,\n";
	print "ranks each of the permutation values by position and compares to the initial fst values\n";
	print "at that position to find a p-value.  All positions, all pvalues are written to an out file\n";
	print "and the pval < 0.05 are printed to standard out as well.\n";
	print "don't forget to follow this script with multiple test corrections in R!  {p.adjust}\n";
	print "\nusage: script initial.fst.tab smoothed_permuted.fst.tab out.tab\n";
	exit;
	}

open(FST, $ARGV[0]);
open(P_S, $ARGV[1]);	# stands for 'permuted and smoothed'
open(OUT, ">$ARGV[2]");
$pval = 0.05;
$time = localtime();
undef(@pvalues);

print "$time: starting to read in initial_Fst_values file!\n";
while(<FST>)
	{chomp; $row = $_; if ($row=~/^contig/) {next;}
	($ctg, $pos, $fst) = split(/\t/,$row);
	$hash{$ctg}{$pos} = $fst;}

while(<P_S>)
	{chomp; @line = split(/\t/,$_); 
	if ($line[0]=~/^contig/) 
		{$length = scalar(@line) - 2; $time = localtime(); print "$time: starting to read in the permuted and smoothed file - has $length columns\n";next;}
	$contig = $line[0]; $snp = $line[1];
	$fst_gt_initial=0;$fst_lt_initial=0;
	@sorted = sort {$a<=>$b} @line[2..$#line];

	foreach $idx (0..$#sorted)
		{

		### how many times is the hashed/initial value greater than the distribution of permuted values?

		if ($sorted[$idx] > $hash{$contig}{$snp}) {$fst_gt_initial++;}
		if ($sorted[$idx] <= $hash{$contig}{$snp}) {$fst_lt_initial++;}	
		}
	$total = $fst_gt_initial + $fst_lt_initial;
	$pv = $fst_gt_initial/$total;

	if ($pv < 0.05) {$sig{$contig}{$snp}=$pv;}
	print OUT "$contig\t$snp\t$pv\n";
	}
	
print "Significant outliers (pval < 0.05) BEFORE multiple test corrections:\n";
if (! %sig) {print "No significant outliers.\n";}
foreach $c (sort keys %sig)
	{
	foreach $p (sort {$a<=>$b} keys %{$sig{$c}})
		{
		print "$c\t$p\t$sig{$c}{$p}\n";
		}
	}
