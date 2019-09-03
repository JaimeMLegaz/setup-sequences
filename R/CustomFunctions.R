length_longest_sequence <- function(sequences){
  return (max(nchar(sequences)))
}

padding_sequence <- function(sequence, len) {

  newseq <- sequence

  if (nchar(sequence) < len){
    to_pad <- len - nchar(sequence)
    for (i in 1:to_pad){
      newseq <- paste0(newseq,"X")
    }
  }
  else if (nchar(sequence) > len){
    newseq <- substr(sequence,1,len)
  }

  return (newseq)

}

padding_sequence_left <- function(sequence, len) {

  newseq <- sequence

  if (nchar(sequence) < len){
    to_pad <- len - nchar(sequence)
    for (i in 1:to_pad){
      newseq <- paste0("X",newseq)
    }
  }
  else if (nchar(sequence) > len){
    newseq <- substr(sequence,1,len)
  }

  return (newseq)

}

to_onehot <- function(sequences){

  newsequences <- sapply(sequences,gsub,pattern="A",replacement="10000", USE.NAMES = FALSE)
  newsequences <- sapply(newsequences,gsub,pattern="C",replacement="01000", USE.NAMES = FALSE)
  newsequences <- sapply(newsequences,gsub,pattern="G",replacement="00100", USE.NAMES = FALSE)
  newsequences <- sapply(newsequences,gsub,pattern="T",replacement="00010", USE.NAMES = FALSE)
  newsequences <- sapply(newsequences,gsub,pattern="X",replacement="00001", USE.NAMES = FALSE)
  newsequences <- sapply(newsequences,gsub,pattern="[^01]",replacement="W", USE.NAMES = FALSE)

  return(newsequences)

}

padding_sequences <- function(sequences, length=250) {
  return(sapply(sequences,padding_sequence,len=length, USE.NAMES = FALSE))
}

padding_sequences_left <- function(sequences, length=250) {
  return(sapply(sequences,padding_sequence_left,len=length, USE.NAMES = FALSE))
}

setup_sequences <- function(sequences) {
  return (pad_sequences(to_onehot(sequences)))
}

mirror_sequences <- function(sequences){

  return (sapply(sequences,mirror_sequence, USE.NAMES = FALSE))
}

mirror_sequence <- function(sequence){
  return (intToUtf8(rev(utf8ToInt(sequence))))
}

strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste,
                                 collapse="")

prefix_statistics <- function(dataframe){ # This dataframe has to contain all sequences we are working with
  members <- 0
  totals <- nrow(dataframe)

  df <- dataframe

  longitudes <- NULL
  compartidos <- NULL
  secuenciasPorGen <- NULL

  # Repeating this while there are more than one sequence in the dataframe
  while (nrow(df) > 1)
  {
    # Selecting the sequences of the first gene
    select <- df$ensembl_gene_id == df$ensembl_gene_id[1]
    actual <- df[select,]
    df <- df[!(select),]

    secuenciasPorGen <- c(secuenciasPorGen,nrow(actual))

    # This will only be done if there is more than one sequence
    while (nrow(actual) > 1){
      # Extracting the first sequence and moving it to another dataframe
      temp <- actual[1,]
      actual <- actual[-1,]

      # Doing the same with the rest of sequences that start with the same character
      charac <- substr(temp$`5utr`,1,1)

      samePrefix <- actual[substr(actual$`5utr`,1,1) == charac,]
      actual <- actual[!(substr(actual$`5utr`,1,1) == charac),]
      temp <- rbind(temp,samePrefix)

      if (nrow(temp) > 1){ # The following will only be done if there are more than one sequence sharing the prefix

        # Checking the length of sharing to see if it is the same as the prefix
        maxlen <- min(nchar(temp$`5utr`))

        j <- 1

        while (j <= maxlen && nrow(temp) > 1){
          charac <- substr(temp$`5utr`[1],j,j)
          comp <- substr(temp$`5utr`,j,j) == charac
          temp <- temp[comp,]
          j <- j + 1
        }

        j <- j - 1

        if (j >= maxlen && nrow(temp) > 1) { # Checking if we went through a whole sequence, which will be our prefix
          members <- members + nrow(temp)
          compartidos <- c(compartidos, nrow(temp))
          longitudes <- c(longitudes, j)
        }

        }
      }
  }

  # Percentage of sequences sharing a prefix
  prefixSharePerc <- sum(compartidos) / totals

  # Longest prefix
  maxPref <- max(longitudes)

  # Shortest prefix
  minPref <- min(longitudes)

  # Prefix length mean
  meanPref <- mean(longitudes)

  # Most shared prefix
  mostShared <- max(compartidos)

  # Less shared prefix (from the shared prefixes of course)
  leastShared <- min(compartidos)

  # Mean of sharing rate of the prefixes
  meanShared <- mean(compartidos)

  # Mean number of sequences a gene has
  meanSeqPerGen <- mean(secuenciasPorGen)


  result <- list(prefixSharePerc,maxPref,minPref,meanPref,mostShared, leastShared, meanShared, meanSeqPerGen)
  names(result) <- c("prefixSharePerc","maxLen","minLen","meanLen","mostShared","leastShared","meanShared","meanSeqPerGen")
  return (result)

}
