#!/usr/bin/env perl

# -----
print "\n","-"x60, "\n";
print "vcf2fst_per_locus.pl v 1.0\n";
print "Roxana Capper, Univ of Texas at Austin\n";
print "Last modified: 12 Feb 2014\n";
print "-"x60,"\n";

unless ($#ARGV == 3) {
	print "this script simply calculates UNWEIGHTED Fsts per LOCUS according to\n";
	print "Mark Kirkpatrick's formula Fst = (p1 – p2)^2 / (4 pBar qBar)\n";
	print "Output of this script can then be plugged into the next script\n";
	print "(to be written after this one, later than 12 Feb 2014) which will\n";
	print "take a supplied sigma and step values to calculate Gaussian-weighted\n";
	print "Fsts per sliding window.  But, that's next.  \n\n";
	print "output file looks like this:\ncontig\tpos\t#pop1_genotyped_alleles\t#pop2_genotyped_alleles\tFst\n\n";
	print "usage: script in.vcf pop1 pop2 out.tab\n";
	print "example: vcf2fst_per_locus.pl K N KxN.fst.tab\n";
	exit
	}
	
open(VCF, $ARGV[0]);
$pop1 = $ARGV[1];
$pop2 = $ARGV[2];
open(OUT, ">$ARGV[3]");
$no_data_for_one_pop=0;$monomorph=0;$triallelic=0;$biallelic=0;
undef(@all_fst);

print "Reading in vcf file line by line!\n";
print OUT "contig\tpos\t#genotyped_pop1\t#genotyped_pop2\tFst\n";

while(<VCF>)
	{chomp;
	$line = $_;	
	if ($line =~ /^##/) {next;}
	if ($line =~ /^#CHROM/)  	###  Split up the populations by noting their indices
		{@line = split(/\t/,$line);
		for ($index = 9; $index <= $#line; $index++) 
			{if ($line[$index] =~ (/^$pop1/)) {push (@pop1_index, $index);}
			if ($line[$index] =~ (/^$pop2/)) {push (@pop2_index, $index);}}
		next;}
   ### Start pulling out genotypes
	@pop1_genot=();@pop2_genot=();@combo=();
	@spl = split(/\t/,$line); $ctg = $spl[0]; $pos = $spl[1];
	foreach $idx_pop1 (@pop1_index)
		{@pop1_info = split(/:/,$spl[$idx_pop1]); 
		@pop1_gts = split(/\//,$pop1_info[0]);
		if ($pop1_gts[0] eq ".") {next;}
		push(@pop1_genot, @pop1_gts);}
	foreach $idx_pop2 (@pop2_index)
		{@pop2_info = split(/:/,$spl[$idx_pop2]); 
		@pop2_gts = split(/\//,@pop2_info[0]);
		if ($pop2_gts[0] eq ".") {next;}
		push(@pop2_genot, @pop2_gts);}
	@combo = @pop1_genot; push(@combo, @pop2_genot);	
	
	$test = scalar(@combo);


	### Pop1 genotypes: in @pop1_genot; # alleles sampled = $pop1_count
	### Pop2 genotypes: in @pop2_genot; # alleles sampled = $pop2_count
	### Both pops together: in @combo; # alleles sampled = $total_count

	### Skip loci that are only genotyped in a single population:
	if ((scalar @pop1_genot == 0)|(scalar @pop2_genot == 0)) {$no_data_for_one_pop++; next;}

   	### note invariants, skip triallelics+ (for now)	
   	$var_ct=0;
   	undef(@unique); my %seen=();
   	my @unique = grep { ! $seen{$_}++} @combo;
   	foreach $elem (@unique) {$var_ct++;}
   	if ($var_ct == 1) {$monomorph++;next;}
   	if ($var_ct > 2) {$triallelic++;next;}  ## <--- Here is where to address calculations of triallelics
	$biallelic++;  #every SNP that made it to this point is biallelic
	
	### testing
	
	
### BUG ALERT:
### what about cases where one pop comparison is only 0's and 2's, not 0's and 1's?
### within that pop it's bilallelic but you have to explicitely call $p1_ct{2} not 
### WAIT NO IT's OKAY; unless every ind has the alt allele, in which case there will 
### be ZERO $p1_ct{0}, etc etc.  (i.e. those two pops have 1 and 2 alleles, not 0 and 1.	
	
	
  ### Calculate p1 and p2:  Fst = (p1 – p2)^2 / (4 pBar qBar)\n";
  ### p = freq of "0" allele, q = freq of "1" allele; arbitrarily defined.  	
  ### p1 = freq of "0" in pop1, p2 = freq of "0" allele in pop2
	my %p1_ct=();my %p2_ct=();
	foreach $gt (@pop1_genot) {$p1_ct{$gt}++;}
	foreach $gt (@pop2_genot) {$p2_ct{$gt}++;}
	$pop1_totalct = scalar(@pop1_genot);
	$pop2_totalct = scalar(@pop2_genot);
	if ((exists $p1_ct{0}) | (exists $p2_ct{0})) {
		$p1 = $p1_ct{0} / $pop1_totalct;
		$p2 = $p2_ct{0} / $pop2_totalct;}
	else {$p1 = $p1_ct{1}/$pop1_totalct;$p2 = $p2_ct{1}/$pop2_totalct;	}
	$pbar = ($p1+$p2)/2;
	$qbar = 1 - $pbar;
#	$numerator = ($p1-$p2)**2;
#	$denom = 4*$pbar*$qbar;
#	$fst = $numerator/$denom;
	$fst = (($p1-$p2)**2) / (4 * $pbar * $qbar);
	push(@all_fst,$fst);
	print OUT "$ctg\t$pos\t$pop1_totalct\t$pop2_totalct\t$fst\n";	
#	print "$ctg\t$pos\t$pop1_totalct\t$pop2_totalct\t$p1\t$p2\t$pbar\t$qbar\t$fst\n";
	}   
	
$sum += $_ for @all_fst;
$global_fst = $sum/$biallelic;
print "Finished!\n\nThe global Fst value, averaged over $biallelic biallelic SNPs, is $global_fst.\n\n";	   
print "there are $no_data_for_one_pop loci that are missing data for one pop or the other and were therefore skipped over.\n";
print "There were $monomorph loci that are totally monomorphic between these two pops.\n";
print "There were $triallelic loci that had more than two alleles and were therefore skipped (in this version of the script).\n";   
   
   
   
   
   
