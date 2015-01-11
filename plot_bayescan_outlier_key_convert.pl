#!/usr/bin/env perl

# -- script description and required arguments
unless ($#ARGV == 3)
	{
	print "\nplot_bayescan_outlier_key_convert.pl\n";
	print "v. 2.0, Roxana Capper\n";
	print "last modified 11 January 2015\n";
	print "\n\ntakes plot_bayescan-identified outliers and matches them to the\n";
	print "key made when the files were split from the original genepop file.\n";
	print "the key was generated through the script 'match_genepop_locus_names_to_bayescan_index.pl'\n";
	print "The key input is simply a single column populated by locus names in order.\n";
	print "The output of this script is a tab file that looks like this:\n";
	print "chrom\tpos\tFst\n";
	print "\nCaveat: this script was written for loci with the naming convention\n";
	print "contignumberXXXX_position; i.e., it splits the key file on the underscore etc.\n";
	print "check to see if you have to modify it for your purposes.\n";
	print "\nusage: script plot_bayescan.outliers file.key in.bayes_fst.txt out.converted.outliers\n\n";
	exit;
	}

open(OUTLIERS,$ARGV[0]);
open(KEY,$ARGV[1]);
open(FST,$ARGV[2]);
open(OUT,">$ARGV[3]");

while(<KEY>)
	{
	chomp;
	if ($_ =~ /^BayeScanIndex/) {next;}
	($number,$name) = split(/\t/,$_);
	($chr,$pos) = split(/_/,$name);
	$line = join("\t",$chr,$pos);
	$hash{$number}=$line;
	}

while(<FST>)
	{chomp;
	if ($_ =~ /^\s/) {next;}
	@line = split(/\s+/,$_);
	$fst{$line[0]}=$line[5];
	}

print OUT "contig\tpos\tBayeScanFst\n";
while(<OUTLIERS>)
	{
	chomp;
	print OUT "$hash{$_}\t$fst{$_}\n";
	}

