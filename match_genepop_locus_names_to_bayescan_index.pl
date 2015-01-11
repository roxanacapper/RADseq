``` perl

#!/usr/bin/env perl

# match_genepop_locus_names_to_bayescan_index.pl
# v. 2.0,  Roxana Capper
# Last modified 11 Jan 2015

unless ($#ARGV == 1) {
        print "\nthis script extracts the locus names from line 2 of a \n";
        print "genepop file, transposes it and prints an index to the right\n";
        print "of it, such that one can take the produced file and use it as \n";
        print "a key to rename the unlabeled locus numbers from BayeScan\n";
        print "output .fst files.\n";
        print "\nusage: script in.genepop out.tab\n";
        exit;}
```
