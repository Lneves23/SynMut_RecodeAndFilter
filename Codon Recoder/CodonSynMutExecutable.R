#set working directory and load in dependencies
wd <- getwd()
setwd(wd)
source("CodonSynMutFunctions.R")


#adjustable commands using elongation rates

#e.coli translation rates:
#Paper: https://www.nature.com/articles/s41598-019-43857-5#Sec29
#Supplementary # 7: https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-43857-5/MediaObjects/41598_2019_43857_MOESM1_ESM.pdf
codon_rates_raw <-
  as.data.frame(read.csv(
    "Escherichia_coli_codonSpecificElongation.csv"))

#reformats list of codon values for downstream use
codon_rates <- get_CodonIndexValues(codon_rates_raw)



#read in and get translation rate for original full sequence
#get_totalSeqIndexValue required: (DNAsequence, codon_index values)
#get_totalSeqIndexValue can be used for a single value (as in DNAoriginal) or a dataframe of values(as in SynMut_candidates_rates below)
DNAoriginal <- readDNAStringSet("ExampleSequence.fasta", format = "fasta")

DNAoriginal_rate <- 
  get_totalSeqIndexValue(
    DNAsequence = DNAoriginal, 
     codon_indexValues = codon_rates)


#can also trim the original sequence if desired, in this case only targeting the first 20 AAs after the start codon
DNAoriginal_5prime <- 
  subseq(x = DNAoriginal, start = 4, end = 63)

DNAoriginal_5prime_rate <-
  get_totalSeqIndexValue(DNAoriginal_5prime, codon_rates)


#get list of mutated sequences
#generate_SynMut requires: (seq, repetitions, repetition_percent, iterations, iteration_percent)
#repetitions = number of times to restart from the original sequence
#repetition_percent = percent of sequence to recode
#iterations = number of times to iterate off a repetition
#repetition_percent = percent of sequence to recode during an iteration
SynMut_candidates <-
  generate_SynMutCandidates(
    DNAoriginal_5prime,
    repetitions = 10,
    repetition_percent = .3,
    iterations = 10,
    iteration_percent = .3
  )



#calculate rate values for mutated sequences
SynMut_candidates_rates <-
  get_totalSeqIndexValue(
    DNAsequence = SynMut_candidates, 
    codon_indexValues = codon_rates)



#filter based on upper tAI cutoff value
#get_filteredCandidates requires: (candidatelist_seq, candidatelist_indexValue, original_indexValue, indexValue_cutoff, is.upperbound = FALSE)
#Default is to treat cutoff value as the the lower limit. is.upperbound = TRUE sets cutoff as upper limit.
filtered_candidates <-
  get_filteredCandidates(
    candidatelist_seq = SynMut_candidates,
    candidatelist_indexValue = SynMut_candidates_rates,
    original_indexValue = DNAoriginal_5prime_rate,
    indexValue_cutoff = 4,
    is.upperbound = FALSE
  )



#export
write.csv(filtered_candidates,
          "FilteredSynonymousMutantSequences_translationRate.csv")



