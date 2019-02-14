# Demo consts
exp.nr <- 1
gender <- "Male"

#' Run demo using expr number 1 and male = gender
#' Asks for file to be selected
#' @export
#'
run.demo <- function() {
  # Run demo script
  cat("Select snp data","\n")
  snp.file <- file.choose()
  cat("Read snp data","\n")
  snpm.data <- read.table(snp.file, skip = 9, header = T, sep = "\t")[,c(3,4,6)]

  cat("Select snp data","\n")
  deviations.file <- file.choose()
  cat("Read events","\n")
  deviations <- read.table(deviations.file, header = T, sep = "\t")

  cat("Run analysis","\n")
  MosaicCalculator(exp.nr, gender, snpm.data, deviations)

  cat("Done","\n")
}

