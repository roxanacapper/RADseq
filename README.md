RADseq
======
RADseq tutorial for non-model systems

Overview:
======
0.  **TypeIIB RAD Library Prep**
    1.  Digest genomic DNA with BcgI (decision point: which enzyme to use, AlfI or BcgI?  Things to consider: heat-inactivation ability, methylation sensitivity.)
1.  **Processing Raw Reads**: 
    1.  Quantify number of tags mapping, read depth per individual, etc (Hohenlohe 2010)
    3.  Meaning of "heterozygosity" and its usage in excluding Night dudes, other "related" guys; how to explain this in the methods section?   
    8.  Ascertainment bias and its effect on pop gen statistics (see [RADseq underestimates diversity and introduces genealogical biases due to nonrandom haplotype sampling](http://onlinelibrary.wiley.com/store/10.1111/mec.12276/asset/mec12276.pdf?v=1&t=hoc0vd7e&s=f99b10dd0191894fb44d066b189f18c2b4b595eb), B. ARNOLD,1 R. B. CORBETT-DETIG,1 D. HARTL and K. BOMBL IE S 2013; and [Emily McTavish's dissertation](http://repositories.lib.utexas.edu/bitstream/handle/2152/21925/MCTAVISH-DISSERTATION-2013.pdf?sequence=1))
2. **Connectivity (demographics) among populations**
      1. Fst/distance (IBD) - SICB 2013
      2. STRUCTURE (neutral/non-outlier SNPs only)
      3. Global population statistics
      4. Migration rates between reefs
      5. HWE (maybe)
      6. Landscape Genomics: PCA nad geography (see Novembre 2008; Helyar 2011 pg 128, bottom right; Paschou et al. 2007); PCA may outperform STRUCTURE in some cases (Helyar 2011); effectively selects the subset of SNPs that describe the differences in the data
      7. Population Size (Ne) (see the program [LDNe](http://conserver.iugo-cafe.org/user/Robin%20Waples/LDNe))
      8. MDS (multidimensional scaling plot of "neutral" loci)
3.  **History (past demographics) among populations**
      1.  Coalescent modeling (can this be done?)    
4.  **Local Adaptation**
      1. Identify regions with molecular signatures of selection, etc.
        1. Define signatures;
        2. Calculate continuous pop gen stats along the gemome:
          1. allele frequency spectrum (MAF, TajD, pi, Heterozygosity _H_)
          2. LD (maybe not)
          3. Fst
      2.  Annotate those regions (Maker 2)
      3.  Logistic regression mapping outlier presence/absence onto continuous environmental variables/gradients (assosiate genotypes to clines)
5.  **Zoox connectivity** (not going to work; 91% of zoox-mapped reads (~20% of all reads) also map to the coral genome -- i.e., likely very nonspecific)
      1.  Patterns of zoox signatures
      2.  Connectivity among populations

Type II B RAD Library Prep and background information
======

[Type II B RAD](http://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2023.html) differs from other RAD flavors because it uses a special restriction enzyme that cuts outside of its recognition site.  For the enzyme BcgI, the excised tag is  ...∇(N10)CGA(N6)TGC(N12)∇...; the recognition site is six bases long but the rest of the tag has 30 arbitrary bases.  This means that the enzyme can excise the same tag many times from different copies of the genome, yielding high coverage, and that each tag is effectively randomly scattered throughout the genome with respect to other tags.  Other enzymes can be used such as AlfI, but keep in mind that some enzymes are not heat-inactivatable and some enzymes, such as BcgI, are methylation sensitive.  We're not sure what kind of bias, if any, is introduced via methylation sensitivity, but it's definitely something to consider.

Here, we use TypeIIB RAD with the enzyme BcgI to digest the genome of the scleractinian coral __Acropora millepora__.  We followed the protocol outlined on the [Matz Lab website -> methods -> 2bRAD protocol](http://www.bio.utexas.edu/research/matz_lab/matzlab/Methods_files/2bRAD_protocol-1.pdf) (last updated 9 July 2014.  Briefly, 
1.  DNA is digested with BcgI to excise 36-base tags distributed throughout the genome.  
2.  Next, adaptors are ligated to the tags and the tags are amplified via PCR.  
3.  The tags are then purified out of the digestion mixture by gel electrophoresis and gel excision of the short, 170-base fragments (36-base tags + adaptors, etc).  
4.  Individual samples are then quantified and sent for Illumina sequencing.

Processing Raw Reads: 
======

So you have a bunch of raw sequence data.  Now what?

The first step in any sequencing project analysis is to remove non-informative sequences, such as adaptor sequences and low quality bases.  But, after that, there are quite a few decision points that can affect the output of your analysis.  The effective magnitudes of those decisions is totally unknown and data-dependent.

* __Read filtering__: more than trimming quality and adaptors -- maybe.
    - do you require each retained read to contain the canonical restriction site? 
    - do you trim the overhanging bases from the digestion?  The enzyme leaves sticky ends when it cuts.  The question here is if we should analyze the length of the tag or if we should remove the two overhanging bases from each end of the reads.  Using some combinations of things we found that there is an overabundance of SNPs on bases 2 and 35 of the reads (due to ligation? degradation?); this goes away under other conditions (namely, GATK's Unified Genotyper works fine).  Might be most conservative to get rid of them, but on the other hand you could be throwing out a ton of data...

Reference Selection
* __Reference decision__:  
    -  Could __in silico__ extract each potential, canonical RAD tag from the genome and map to those; 
    -  Could __in silico__ extract each potential RAD tag allowing for up to one mismatch in each restriction site
    * GATK can't call SNPs from large ref databases, so _extracting all the potential genomic RAD tags_ is not technically feasible.  
    * Also, some colleages are concerned that by forcing mapping to a reference genome's digest you may force a bad read to map to a bad reference.  However, those bad reads should be removed via filtering for quality, so it's likely a moot point.
    -  Mapping to the 12.5k contigs (_whole genome reference_) is much faster.  However, there is some thoguht about shifty restriction site mapping (no data to support this either way).  

* __Mapping__: what reference resource do you map to?  Whole-genome?  Extracted tags?  Tags extracted but allowing 1-mismatch?  32-base tags or 36-base tags?  (i.e., are there more errors at the ends of the reads?  Some evidence says yes.)
    0) De novo mapping?  
    1) extract all possible tags, including 1 mismatch away from the canonical restriction site.  (GATK can't handle so many "chromosomes"/reference contigs so this isn't practical, unless you only give it the reference tags that are mapped to at least 5x or whatever)
    2) extract all exact tags, matching the restriction site exactly (GATK can't handle so many "chromosomes"/reference contigs so this isn't practical, unless you only give it the reference tags that are mapped to at least 5x or whatever)
    3) map to the whole genome    
    4) concatenate all potential RAD tags into one mega-gene; either map to this, or map to tags separately and then glue them together for GATK (totally unnecessary)
    5) Unified Genotyper.  This may not work easily, though it will make GATK run super fast.)

* __Mappers__: what algorithm to use?  Does it make a difference? (yes)  Is there any way to validate them? (no)
    *  Bowtie1 doesn't play nicely with GATK; it sets every MAPQ value to 255 for "yes mapped" and 0 for "not mapped".  See [seqanswers1](http://seqanswers.com/forums/showthread.php?t=3142) and [seqanswers2](http://seqanswers.com/forums/showthread.php?t=10594)
    *  Bowtie2
    *  SHRiMP
    *  BWA
    *  Stampy?

* __SNP callers__: 
    - home-made scripts, 
    - STACKS, 
    - frequency-based SNP calling, 
    - samtools mpileup, 
    - GATK

* __GATK__: 
    - to use local indel realignment or not?  Some colleagues think that ignoring this kind of biology is a Very Bad Thing.
    - Unified Genotyper or HaplotypeCaller?
    - Do you follow the Best Practices including VQSR?  Do you train VQSR on the top 5% or top 10% of SNPs by quality?  by another metric?  Yes.  But, are these all artifacts??

*__VCF filtering__:
    - filter to retain 90% of "true variants"
    - Keep only loci that are genotyped in >80% of individuals
    - Keep only individuals in which >60% of all loci remaining are genotyped.
    - Remove clonal individuals (those with extremely high relatedness; >95% the same genotypes, which is totally unlikely in natural populations.  Still not sure why this is a problem; multiple samplings of the same individual?  Tube confusion during library prep?)
    - Remove individuals with low heterozygosity (`vcftools --vcf final.vcf --het`). 
    

Anyway, you've made all these decisions and now you have a vcf file that is really pretty.
 
    
