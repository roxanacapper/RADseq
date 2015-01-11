RADseq
======
RADseq tutorial for non-model systems

Overview:
======
0.  **TypeIIB RAD Library Prep**
    1.  Digest genomic DNA with BcgI 
    2.  Make libraries
    3.  Sequence (Illumina)
1.  **Processing Raw Reads**: 
    1.  Quantify number of tags mapping, read depth per individual, etc (Hohenlohe 2010)
    2.  Meaning of "heterozygosity" and its usage in excluding Night dudes, other "related" guys; how to explain this in the methods section?   
    3.  Ascertainment bias and its effect on pop gen statistics (see [RADseq underestimates diversity and introduces genealogical biases due to nonrandom haplotype sampling](http://onlinelibrary.wiley.com/store/10.1111/mec.12276/asset/mec12276.pdf?v=1&t=hoc0vd7e&s=f99b10dd0191894fb44d066b189f18c2b4b595eb), B. ARNOLD,1 R. B. CORBETT-DETIG,1 D. HARTL and K. BOMBL IE S 2013; and [Emily McTavish's dissertation](http://repositories.lib.utexas.edu/bitstream/handle/2152/21925/MCTAVISH-DISSERTATION-2013.pdf?sequence=1))
2. **Connectivity (demographics) among populations**
      2. STRUCTURE (neutral/non-outlier SNPs only)
      3. Fst/distance (IBD) - SICB 2014
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

Here, we use TypeIIB RAD with the enzyme BcgI to digest the genome of the scleractinian coral __Acropora millepora__.  We followed the protocol outlined on the [Matz Lab website -> methods -> 2bRAD protocol](http://www.bio.utexas.edu/research/matz_lab/matzlab/Methods_files/2bRAD_protocol-1.pdf) (last updated 9 July 2014). 

Briefly, 
1.  DNA is digested with BcgI to excise 36-base tags distributed throughout the genome.  
2.  Next, adaptors are ligated to the tags and the tags are amplified via PCR.  
3.  The tags are then purified out of the digestion mixture by gel electrophoresis and gel excision of the short, 170-base fragments (36-base tags + adaptors, etc).  
4.  Individual samples are then quantified and sent for Illumina sequencing.

---

Processing Raw Reads: 
======

So you have a bunch of raw sequence data.  Now what?

The first step in any sequencing project analysis is to remove non-informative sequences, such as adaptor sequences and low quality bases.  But, after that, there are quite a few decision points that can affect the output of your analysis.  The effective magnitudes of those decisions is totally unknown and data-dependent.

* __Read filtering__: more than trimming quality and adaptors -- maybe.
    - do you require each retained read to contain the canonical restriction site? 
    - do you trim the overhanging bases from the digestion?  The enzyme leaves sticky ends when it cuts.  The question here is if we should analyze the length of the tag or if we should remove the two overhanging bases from each end of the reads.  Using some combinations of things we found that there is an overabundance of SNPs on bases 2 and 35 of the reads (due to ligation? degradation?); this goes away under other conditions (namely, GATK's Unified Genotyper works fine).  Might be most conservative to get rid of them, but on the other hand you could be throwing out a ton of data...

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

* __VCF filtering__:
    * filter to retain 90% of "true variants"
    * Keep only loci that are genotyped in >80% of individuals
    * Keep only individuals in which >60% of all loci remaining are genotyped.
    * Remove clonal individuals (those with extremely high relatedness; >95% the same genotypes, which is totally unlikely in natural populations.  Still not sure why this is a problem; multiple samplings of the same individual?  Tube confusion during library prep?)
    * Remove individuals with low heterozygosity (`vcftools --vcf final.vcf --het`). 
    

Anyway, you've made all these decisions and now you have a vcf file that is really pretty.
 
---

Connectivity (demographics) among populations
=====

**STRUCTURE** with subset of all SNPs: 
===

Visualizing population structure within the metapopulation with a subset of all SNPs, neutral plus potential outliers.
Note: there are a lot of fussy reformatting steps here.  There are, for sure, better or more efficient ways to do this same thing, but this is my method cobbled together with duct tape and a little hope.  It should get you by but definitely there is room for improvement.

1.  Convert .vcf into .genepop format: `vcf2genepop.pl vcf=variants.vcf pops=K,M,N,O,S > variants.genepop`
        *  Script notes: pops = single-letter prefix of each individual, assigning them to whichever population.  You might need to reformat your individuals' names to fit this, or edit the script to handle your individual names.
2.  Split genepop into pieces.  In my experience, STRUCTURE gets upset when you throw more than 5,000 SNPs at it (takes too long to run).  Obviously, this is data-dependent on the number of SNPs, number of individuals, k value, and however many million runs you want to do.  `genepop_split.pl variants.genepop 18 r`
    - '18' = number of pieces to split your file into
    - 'r' = split **r**andomly or sequentially ('s', I think) 
3.  Rename the first subset something easy.  We only need to work with one subset of all SNPs, though for sure you can (and probably should) run a few additional SNP samples just to check that everything looks the same.  `cp sub1_variants.genepop sub_var.genepop`
4.  Convert `sub_var.genepop` into `structure` format via PGDSpider2: `java -jar PGDSpider2-cli.jar -inputfile sub_var.genepop -inputformat GENEPOP -outputfile sub_var.structure -outputformat STRUCTURE -spid genepop2structure.spid`
    - Generate the -spid file via the [PGDSpider manual](http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider%20manual_vers%202-0-7-2.pdf)
5.  Run STRUCTURE on the command line (hopefully using a cluster): `structure -K 2 -i sub_var.structure -o sub_var_k2_r1.out`
    - Run each K value 5x (r1, r2, r3, ...)
    - Run K values for 1..6 or whatever range fits your data best.
6.  Run CLUMPP, DISTRUCT, etc;  See http://bodegaphylo.wikispot.org/Structure

**STRUCTURE** with all previously-identified OUTLIER SNPs
===

I'm skipping ahead a little.  To do this, obviously you have to identify your outlier SNPs beforehand (ex., Fst analysis).
The idea is to see if you can identify gene trees, or use gene clusters (linked loci) (yes, even in STRUCTURE which assumes zero linked loci), to reveal population history, unusual patterns, or anything else.

1.  Convert `variants_only.vcf` into `variants_only.genepop` (vcf2genepop.pl)
2.  No need to split outliers list into subpieces; there probably aren't that many on the list.  For example, I only had 90 identified with BayeScan.
3.  Convert `outliers_only.genepop` into `structure` format (PDGSpider)
4.  Run STRUCTURE for 5 reps, K = 2..8
  


**Model-based Fst with BayeScan**
===

---> NOTES ABOUT model-based Fst calcs <----
---> NOTES ABOUT model-based vs. model free calcs <----
---> NOTES ABOUT BAYESCAN SPECIFICALLY <----
---> NOTES ABOUT running Bayescan, particuarly about the number/type of pop comparisons  <----


- Run BayeScan once per population pair
- Run BayeScan once for all five pops (metapopulation)
- Run BayeScan twice for a single population pair (check for convergence)
 
1.  Convert .vcf to .genepop: `vcf2genepop.pl vcf='KxO.vcf' pops=K,O > KxO.genepop`
2.  If `variants.genepop` is too large to run in the 24-hour cluster window (particularly for the metapop analysis), split it into two pieces or more and run them as two jobs:  
    - `genepop_split.pl variants.genepop 2 r`
3.  Convert genepop to bayescan format: `java -jar PGDSpider2-cli.jar -inputfile KxO.genepop -inputformat GENEPOP -outputfile KxO.bayescan -outputformat GEST_BAYE_SCAN -spid genepop2bayescan.spid`
4.  Run BayeScan: `bayescan_2.1 KxO.bayescan -threads 24`
5.  Locus names are overwritten in BayeScan, so you have to make a key to match them back to their original names, then replace numbers with names.  Sigh. 
    1.  Identify outliers by locus NUMBER (.)
        - Run full `plot_R.r`in [R] to source (script is packaged with BayeScan)
        - in [R], `kxo <- plot_bayescan('KxO.baye_fst.txt',0,FDR=0.05)`
        - in [R], `write(kxo$outliers,file="KxO.outliers",ncolumns=1)`
    2.  Make key of all named:numbered loci for use in matching outlier numbers to locus names:
        - `match_genepop_locus_names_to_bayescan_index.pl KxO.genepop KxO.genepop_key`
    3.  Extract outlier loci (numbered, not renamed) and match them to their correponding "real" names AND bayescan Fst values (which are missing from the BayeScan outliers list): 
        - `plot_bayescan_outlier_key_convert.pl KxO.outliers KxO.genepop_key KxO.baye_fst.txt KxO.renamed.outliers_w_fst`
    4.  Sanity check: check convergence of BayeScan output: run the same input.bayescan file a second time. 
        - `cp KxO.bayescan KxO_x2.bayescan`
        - `bayescan_2.1 KxO_x2.bayescan -threads 24`
        - compare outlier lists.  They should be mostly the same.




**Global Fst**
===

Average BayeScan Fst values over all SNPs.  This is the same data (per locus) that is used to graph model-based Fst along the genome.

**Graphing Fst along the chromosomes** (or reference contigs)


**Signatures of Selection and Local Adaptation**
