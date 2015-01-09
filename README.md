RADseq
======
RADseq tutorial for non-model systems

Overview:
======
0.  **TypeIIB RAD Library Prep**
    1.  Digest genomic DNA with BcgI (decision point: which enzyme to use, AlfI or BcgI?  Things to consider: heat-inactivation ability, methylation sensitivity.)
1.  **Processing Raw Reads**: 
    1.  Quantify number of tags mapping, read depth per individual, etc (Hohenlohe 2010)
    2.  Why to exclude "A" reef, Night dudes
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

Type II B RAD
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

The first step in any sequencing project analysis is to remove non-informative sequences, such as adaptor sequences and low quality bases.
