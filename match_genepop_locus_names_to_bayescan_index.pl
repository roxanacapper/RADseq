#!/usr/bin/env perl

# -- program name, header
# match_genepop_locus_names_to_bayescan_index.pl
# v. 2.0,  Roxana Capper
# Last modified 11 Jan 2015

unless ($#ARGV == 1) {
        print "\nthis script extracts the locus names from line 2 of a \n";
        print "genepop file, which should be the list of loci separated by\n";
        print "a comma and a space, then transposes it and prints an index to the right\n";
        print "of it, such that one can take the produced file and use it as \n";
        print "a key to rename the unlabeled locus numbers from BayeScan\n";
        print "output .fst files.\n";
        print "Can be easily modified to extract the first line instead if your\n";
        print "file is set up a little differently.\n"
        print "\nusage: script in.genepop out.tab\n";
        exit;}

open(IN, $ARGV[0]);
open(OUT, ">$ARGV[1]");
$line_count = 0;

while(<IN>)
	{
	$line_count++;
	chomp;
	$line = $_;
  if ($line_count==2)   #Change this number to whatever line is the locus list in your genepop; should be 1 or 2, depending on headers
    {@contigs = split(/\, /,$line);}
    else {next;}
  }
  
print OUT "BayeScanIndex\tChrom_posID\n";  #modify this line to print out whatever format header that is relevant for your data
$index = 1;

foreach $name (@contigs)
       	{
	print OUT "$index\t$name\n";
       	$index++;
       	}
