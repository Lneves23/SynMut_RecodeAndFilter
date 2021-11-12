#installing and loading required dependencies
library(BiocManager)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("SynMut")
if (!requireNamespace("SynMut", quietly = TRUE))
  BiocManager::install("SynMut")
library("Biostrings")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
library("coRdon")
if (!requireNamespace("coRdon", quietly = TRUE))
  BiocManager::install("coRdon")
library("dplyr")
if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")
library("pracma")
if (!requireNamespace("pracma", quietly = TRUE))
  BiocManager::install("pracma")



  #Functions required for adjustable commands
  
  #function: get_CodonIndexValue  -- Generates total sequence value from input table of individual codons (such as codon weights from raw STADIUM table or translation rates)}
  get_CodonIndexValues <- function(codon_values_raw){
    
    #re-order value table to be alphabetical (to match output from coRdon)
    codon_values_raw <- codon_values_raw[,-1]
    rownames(codon_values_raw) <- codon_values_raw[,1]
    codon_values <- as.data.frame(codon_values_raw)
    codon_values <- as.data.frame(codon_values[order(row.names(codon_values)),])
    codon_values <- codon_values[,-1, drop=FALSE]
    
    
    return(codon_values)
  }



#function: get_totalCodonValue  -- Counts codons of input sequence(s) and determine value (tAI or rate) for whole sequence(s)}
get_totalSeqIndexValue <- function(DNAsequence, codon_indexValues){
  
  #structure fasta file as codon table
  DNA_codontable <- codonTable(DNAsequence)
  #removes metadata and structures as data frame instead of formal class
  DNA_codoncounts <- as.data.frame(codonCounts(DNA_codontable))
  #replaces generic line name with line names from Fasta file (ex. >gene 1)
  rownames(DNA_codoncounts)<- getID(DNA_codontable)
  #remove stop codon count
  DNA_codoncounts <- DNA_codoncounts[,!names(DNA_codoncounts)%in%c("TAA","TAG", "TGA")]
  
  #count and average tAI for input sequence
  weight <- codon_indexValues[,0]
  finalValue <- DNA_codoncounts[,0]
  finalValue$indexValue <- 0
  
  total <- nrow(DNA_codoncounts) * 60
  start_time <- Sys.time()
  for (j in c(1: nrow(DNA_codoncounts))){
    row1 <- DNA_codoncounts[j,]
    for (i in c(1:61)){
      weight[i,1] = (row1[i]*log(codon_indexValues[i,]))
      finalValue$indexValue[j] <- exp(sum((weight))/sum(row1))
    }
    counter <- 100 * ((i * j) / (total))
    if(mod(counter,10) < 0.01){
      print(paste(round(counter,0),"% complete", sep = ""))
      end_time <- Sys.time()
      time_elapsed <- end_time - start_time
      print(time_elapsed)
    }
  }
  
  return(finalValue)
}



#More details on SynMut here: https://uclouvain-cbio.github.io/WSBIM1322/sec-biostrings.html
#Includes both an iterative and repetitive loop to introduce synonymous changes to codon sequence
#function: generate_SynMutCandidates  -- Generates candidates with synonymous mutations}
generate_SynMutCandidates <- function(seq, repetitions, repetition_percent, iterations, iteration_percent){
  
  if(repetitions < 1){
    print("Repetition value must be an integer 1 or greater")
  } else if(iterations < 1){
    print("Iteration value must be an integer 1 or greater")
  } else if(repetition_percent > 1 || repetition_percent < 0){
    print("Repetition percent must be a decimal value between 0 and 1")
  } else if(repetition_percent > 1 || repetition_percent < 0){
    print("Iteration percent must be a decimal value between 0 and 1")
  } else {
    
    #Define DNAStringSet as normal string for manipulation
    origseq <- input_seq(seq)
    #define optimized DNA sequence as new object to maintain DNA_opt as same sequence
    mutseq <- origseq
    
    #create new data frame, assign first row as original sequence
    candidate_mutseqs <- origseq@dnaseq
    candidate_mutseqs@ranges@NAMES[1] <-"original"
    
    
    #Repetitive and iterative for loop- iteratively introduce X% random mutations to codons from SynMut. Add each iteration to candidate_mutseqs dataframe and name as mutant. Once inside loop is complete, returns to fully optimized DNA sequence to start new set of random mutations.
    
    total <- repetitions * iterations
    start_time <- Sys.time()
    for (j in c(1:repetitions)){
      mutseq <- codon_random(origseq, n = repetition_percent)
      for (i in c(1:iterations)){
        var <- ((j-1)*iterations)+i
        mutseq <- codon_random(mutseq, n = iteration_percent)
        candidate_mutseqs[var+1] <- mutseq@dnaseq[1]
        candidate_mutseqs@ranges@NAMES[var+1] <-paste("mut",var,sep="")
      }
      counter <- 100 * ((i * j) / (total))
      if(mod(counter,10) == 0){
        print(paste(counter,"% complete", sep = ""))
        end_time <- Sys.time()
        time_elapsed <- end_time - start_time
        print(time_elapsed)
      }
    }
    
    return(candidate_mutseqs)
  }
}



#function: get_filteredCandidates  -- Sets tAI or rate bounds and filters candidate list}
get_filteredCandidates <- function(candidatelist_seq, candidatelist_indexValue, original_indexValue, indexValue_cutoff, is.upperbound = FALSE){
  
  #filter list of candidates based on upper and lower tAI bounds of original sequence 
  indexValue_filteredList <- candidatelist_indexValue[0,0]
    
    #If is.upperbound is set to true, the cutoff value becomes the upper cutoff for tAI values. So can choose sequences with lower tAI than original
    if(is.upperbound == TRUE){
      indexValuebound <- (as.numeric(original_indexValue) - indexValue_cutoff)
      boundaryconditions <- candidatelist_indexValue$indexValue < indexValuebound
    } else {
      indexValuebound <- (as.numeric(original_indexValue) + indexValue_cutoff)
      boundaryconditions <- candidatelist_indexValue$indexValue > indexValuebound
    }
    
    indexValue_filteredList <- candidatelist_indexValue %>% filter(boundaryconditions | match(rownames(candidatelist_indexValue), "original"))
    
    #calculate ∆tAI values and add in column
    indexValue_filteredList$delta <- indexValue_filteredList$indexValue - indexValue_filteredList$indexValue[1]
    
    #filter candidate sequence list by tAI list and restructure as dataframe
    indexValuehits_seqs <- as.data.frame(candidatelist_seq[rownames(indexValue_filteredList)])
    
    
    #create new data frame with final tAI and sequence
    indexValue_compiledlist <- indexValue_filteredList
    indexValue_compiledlist$sequence <- indexValuehits_seqs[,1]
    
    #filter for unique hits
    indexValue_compiledlist <- unique(indexValue_compiledlist)
    
    
    return(indexValue_compiledlist)
}


