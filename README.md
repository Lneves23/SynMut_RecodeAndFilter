Script built off the SynMut package that uses weighted codon values to:

	1. Assess an input sequence
	2. Generate candidates with recoded codon sequences with synonymous mutations
	3. Filter candidates based on the average codon weight

This code draws heavily on the SynMut package which can be found here: 
  
  	Gu H (2021). SynMut: SynMut: Designing Synonymously Mutated Sequences with Different Genomic
  	Signatures. R package version 1.10.0, https://github.com/Koohoko/SynMut.

Required input files:

	1. Original coding sequence in FASTA file format
	2. Weighted codon value table
		Examples of weighted codon values are tRNA adaptive index and per-codon translation rates
		Script can be used with any type of codon value as long as user provide a codon-to-value table. An example of translation rate is included. 

Example files are included:

	1. ExampleSequence.fasta - CDS for the E. Coli gene ynfD
	2. Escherichia_coli_codonSpecificElongation.csv - codon sequence vs per-codon translation rate
		Translation rates are sourced Supplementary Figure 7 of: 
			Trösemeier, JH., Rudorf, S., Loessner, H. et al. Optimizing the dynamics of protein expression. 
			Sci Rep 9, 7511 (2019). https://doi.org/10.1038/s41598-019-43857-5
