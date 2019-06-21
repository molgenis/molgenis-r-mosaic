#' Run demo using expr number 1 and male = gender
#' Asks for file to be selected
#' @export
#'
run.demo <- function() {
  # Demo consts
  exp.nr <- 1
  gender <- "Male"
  # Run demo script
  cat("Select snp data","\n")
  snp.file <- file.choose()
  cat("Read snp data","\n")

  #In the SNP-array file only columns 'Chr', 'Position and 'B Allele Freq are of interest
  #First find row with headers (following) row with '[Data]' in column 1.
  snp.headerlines <- readLines( snp.file, n=20)
  snp.startdata <- which(grepl(pattern = "Data", x= snp.headerlines))
  #The column headers are in the following line
  snp.colnamesrow <- snp.startdata + 1
  chrom <- which(strsplit(snp.headerlines, "\t")[[snp.colnamesrow]]=="Chr")
  chrompos <- which(strsplit(snp.headerlines, "\t")[[snp.colnamesrow]]=="Position" )
  BAF <- which(strsplit(snp.headerlines, "\t")[[snp.colnamesrow]]=="B Allele Freq" )
  snpm.data <- read.table(snp.file, skip = 9, header = T, sep = "\t")[,c(chrom,chrompos,BAF)]


  cat("Select snp data","\n")
  deviations.file <- file.choose()
  cat("Read events","\n")
  deviations <- read.table(deviations.file, header = T, sep = "\t")

  cat("Run analysis","\n")
  events.filter <- MosaicCalculator(exp.nr, gender, snpm.data, deviations, "~/mosaic-result.pdf")
  print(events.filter)

  cat("Done","\n")
}

