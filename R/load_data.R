#' Runner script to try out MosaicCalculator
#' 
#' Run demo defaults to using using experiment_number 1, gender "Male" and results_path "./mosaic-result.pdf"
#' Asks for file to be selected
#'
#' @param experiment_number number used for tracking experiment.
#' @param gender string one of (Male, Female or Unknown).
#' @param results_path string path and file name of pdf file contaiing results
#' @export
#'
run.demo <- function(experiment_number=1, gender="Male", results_path="./mosaic-result.pdf") {
  # Run demo script
  cat("Select snp data","\n")
  snp_file <- file.choose()
  cat("Read snp data","\n")

  #In the SNP-array file only columns 'Chr', 'Position and 'B Allele Freq are of interest
  #First find row with headers (following) row with '[Data]' in column 1.
  snp_headerlines <- readLines( snp_file, n=20)
  snp_startdata <- which(grepl(pattern = "Data", x= snp_headerlines))
  #The column headers are in the following line
  snp_colnames_row <- snp_startdata + 1
  chrom <- which(strsplit(snp_headerlines, "\t")[[snp_colnames_row]]=="Chr")
  chrom_pos <- which(strsplit(snp_headerlines, "\t")[[snp_colnames_row]]=="Position" )
  BAF <- which(strsplit(snp_headerlines, "\t")[[snp_colnames_row]]=="B Allele Freq" )
  snpm_data <- read.table(snp_file, skip = 9, header = T, sep = "\t")[,c(chrom, chrom_pos, BAF)]

  cat("Select CNV data","\n")
  deviations_file <- file.choose()
  
  cat("Read events","\n")
  deviations <- read.table(deviations_file, header = T, sep = "\t")

  cat("Run analysis","\n")
  events.filter <- MosaicCalculator(experiment_number, gender, snpm_data, deviations, results_path)

  cat("Done","\n")
}
