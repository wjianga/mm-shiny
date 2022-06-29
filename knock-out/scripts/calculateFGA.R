# sampleData is a csv file used to generate the FGA

calc_FGA <- function(sampleData) {
  exomeSize <- 3099734149
  numChanges <- 0
  for(row in 1:nrow(sampleData)) {
    # this line excludes chromosomes Y and M (chromosome X shouldn't be in the dataset)
    if (sampleData[row, 7] == "chrY" | sampleData[row, 7] == "chrM") next
    numChanges <- numChanges + (sampleData[row, 6] - sampleData[row, 5])
  } 
  fga <- numChanges / exomeSize  
  return(fga)
}