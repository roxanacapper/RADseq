#!/usr/bin/env perl
use Bio::Perl;
use Bio::SeqIO;

print "\n","-"x60, "\n";
print "reorder_vcf_to_match_genome.pl v 1.1 Roxana Capper\n";
print "Last modified: 22 May 2013\n";
print "-"x60,"\n";

unless ($#ARGV == 2)
	{
	print "takes an input vcf file and reorders the contigs\n";
	print "to match the order of an input genome.fasta\n";
	print "usage: script input.genome.fasta input.vcf out.vcf\n\n";
	exit;
	}

$genome = new Bio::SeqIO(-file=>$ARGV[0],-format=>"fasta");
open(VCF, $ARGV[1]);
open(OUT, ">$ARGV[2]");

print "starting to read in the genome\n";
while($seq = $genome->next_seq)
	{
	push(@genome_order,$seq->display_id);
	}
print "finished reading in the genome\n";
while(<VCF>)
	{
	chomp;
	$whole_line = $_;
	if ($whole_line=~/^#/) {print OUT "$_\n"; next;}
	@line = split(/\t/,$whole_line);
	$node = $line[0];
	$pos = $line[1];
	$hash{$node}{$pos}=$whole_line;
	}
print "finished reading in the .vcf\n";
foreach $contig (@genome_order)
	{
	if (exists $hash{$contig})
		{
		foreach $position (sort {$a <=> $b} keys %{$hash{$contig}})
			{
			print OUT "$hash{$contig}{$position}\n";
			}
		}
	}

print "done!\n";
