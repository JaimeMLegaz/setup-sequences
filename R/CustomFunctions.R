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

remove_unavailables <- function(dataframe, col) { # Aparentemente no funciona bien
  df <- dataframe[!(dataframe[,col]== "Sequence unavailable"),]
  return (df)
}

mirror_sequences <- function(sequences){

  return (sapply(sequences,mirror_sequence, USE.NAMES = FALSE))
}

mirror_sequence <- function(sequence){
  return (intToUtf8(rev(utf8ToInt(sequence))))
}

strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste,
                                 collapse="")

prefix_statistics <- function(dataframe){ # Tener en cuenta que este dataframe que pasaremos será el dataframe de todas las secuencias. Todas todas.

  # Después almacenar todas las j resultantes en un array, trabajar las estadísticas correspondientes al final, y devolver esos valores

  members <- 0
  totals <- nrow(dataframe)

  df <- dataframe

  longitudes <- NULL
  compartidos <- NULL
  secuenciasPorGen <- NULL

  # Repetimos esto mientras haya más de una secuencia en el dataframe
  while (nrow(df) > 1)
  {
    # Seleccionar las secuencias del primer gen
    select <- df$ensembl_gene_id == df$ensembl_gene_id[1]
    actual <- df[select,]
    df <- df[!(select),]

    secuenciasPorGen <- c(secuenciasPorGen,nrow(actual))

    # Lo de ahora en adelante solo lo voy a hacer si hay más de una secuencia
    while (nrow(actual) > 1){
      # Me saco la primera secuencia y me la llevo a otro dataframe más
      temp <- actual[1,]
      actual <- actual[-1,]

      # Hago lo mismo con aquellas secuencias que empiecen por el mismo caracter
      charac <- substr(temp$`5utr`,1,1)

      samePrefix <- actual[substr(actual$`5utr`,1,1) == charac,]
      actual <- actual[!(substr(actual$`5utr`,1,1) == charac),]
      temp <- rbind(temp,samePrefix)

      if (nrow(temp) > 1){ # Lo de aqui en adelante solo tiene sentido si hay más de una cadena que comparta ese comienzo de cadena

        # Ahora vamos a ver la longitud de su coincidencia para ver si de verdad comparten el prefijo
        maxlen <- min(nchar(temp$`5utr`))

        j <- 1

        while (j <= maxlen && nrow(temp) > 1){
          charac <- substr(temp$`5utr`[1],j,j)
          comp <- substr(temp$`5utr`,j,j) == charac
          temp <- temp[comp,]
          j <- j + 1
        }

        j <- j - 1

        if (j >= maxlen && nrow(temp) > 1) { # Es decir, si hemos llegado a completar alguna cadena entera, que sería considerada nuestro prefijo
          members <- members + nrow(temp)
          compartidos <- c(compartidos, nrow(temp))
          longitudes <- c(longitudes, j)
        }

        }
      }
  }

  # Porcentaje de secuencias que comparten prefijo
  prefixSharePerc <- sum(compartidos) / totals

  # Prefijo más largo
  maxPref <- max(longitudes)

  # Prefijo más corto
  minPref <- min(longitudes)

  # Media de longitud de los prefijos
  meanPref <- mean(longitudes)

  # Prefijo más compartido
  mostShared <- max(compartidos)

  # Prefijo menos compartido (de entre los compartidos)
  leastShared <- min(compartidos)

  # Media de compartición de los prefijos
  meanShared <- mean(compartidos)

  # Media de secuencias que tiene un gen
  meanSeqPerGen <- mean(secuenciasPorGen)


  result <- list(prefixSharePerc,maxPref,minPref,meanPref,mostShared, leastShared, meanShared, meanSeqPerGen)
  names(result) <- c("prefixSharePerc","maxLen","minLen","meanLen","mostShared","leastShared","meanShared","meanSeqPerGen")
  return (result)

} # Fin funcion
