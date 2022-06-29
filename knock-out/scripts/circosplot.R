library(dplyr)
library(ggplot2)
library(BioCircos)
library(stringr)

true_print <- function(txt){
  cat(txt, "\n")
}

split_col <- function(txt){
  strsplit(txt, ";")
} 

parse_point <- function(alt_str){
  split_ary <- (strsplit(as.character(alt_str), "[]:[]"))
  for (word in split_ary[[1]]){
    if (!is.na(as.numeric(word))){
      return(as.numeric(word))
    }
  }
}


parse_chr <- function(alt_str){
  print(alt_str)
  split_ary <- (strsplit(as.character(alt_str), "[]:[]"))
  
  for (word in split_ary[[1]]){
    if (grepl("chr", word)){
      return((substr(word,4,nchar(word))))
    }
  }
}

produceCircosPlot <- function(table) {
  res_pos <- unlist(lapply(table$POS2, parse_point))
  res_chr <- lapply(table$CHROM2, parse_chr)
  st_pos <- unlist(lapply(table$POS1, parse_point))
  st_chr <- lapply(table$CHROM1, parse_chr)
  link_names = 1:length(res_pos)
  tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 1,
                                       borderSize = 0, fillColors = "#EEFFEE")
  tracklist = tracklist + BioCircosLinkTrack('Links', st_chr, st_pos,(st_pos+1), res_chr, res_pos, (res_pos+1), maxRadius = 1, displayLabel = FALSE, color = "#CC0000")
  return(BioCircos(tracklist, genomeFillColor = "Spectral",
            chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  TRUE,
            genomeTicksDisplay = FALSE,  genomeLabelTextSize = "8pt", genomeLabelDy = 0))
}
