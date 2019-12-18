#' Calculate mosaic percentage
#'
#' Takes a experiment number, gender, snp dataframe and events dataframe
#' and calculates the percentage of mosaic dna
#'
#' @param exp.nr numerical experiment number
#' @param gender string denoting gender
#'  ( "Male", "Female" or "Unknown")
#' @param snpm.data dataframe containing snp data
#'  ( "Chr", "Position", "B Allele Freq")
#' @param deviations dataframe containing deviations
#' ( "Chromosome Region",	"Event",	"Length",	"Cytoband",	precentage of CNV Overlap",
#' 	"Probe Median",	"precentage Heterozygous",	"Probes",	"Count of Gene Symbols")
#' @param outfile optional path to direct pdf output to, if not set './Rplots.pdf' is used
#' @import e1071
#' @import stats
#' @import grDevices
#' @import grid
#' @import gridExtra
#' @importFrom utils read.table setTxtProgressBar txtProgressBar
#' @importFrom graphics plot text
#' @export
#'
#Versionnumber: 0.6.0.0
MosaicCalculator <- function(exp.nr, gender, snpm.data, deviations, outfile){
  #Set the option of the output decimal to .
  options(OutDec = ".")

  cat("Conversion from deviations to usable parameters.","\n")
  eventsoutput <- process.deviation.line(deviations, gender)

  cat("Calculate quality aspects array.","\n")

  #Calculating average BAF outside of any event (within 0.2-0.8 BAF range)
  averageBAF <- snpm.data[with(snpm.data, B.Allele.Freq >= 0.2 & B.Allele.Freq <= 0.8),]
  #Calculating over all heterozygous BAF areas
  if(gender == "Female"){
    averageBAF[averageBAF$Chr %in% c("Y","0","MT"),] <- NA
    averageBAF <- averageBAF[complete.cases(averageBAF),]
  }
  if(gender == "Male"){
    averageBAF[averageBAF$Chr %in% c("X","Y","XY","0","MT"),] <- NA
    averageBAF <- averageBAF[complete.cases(averageBAF),]
  }
  if(gender == "Unknown"){
    averageBAF[averageBAF$Chr %in% c("X","Y","XY","0","MT"),] <- NA
    averageBAF <- averageBAF[complete.cases(averageBAF),]
  }
  #Overwriting averageBAF dataframe to exclude SNPs that are within deviations  ##COMMENT: CAN THIS BE PROGRAMMED MORE EFFICIENTLY????
  #WARNING: This apply takes up a lot of memory unless gc() is used!
  #Including a progress bar to show that the tool isn't freezing
  pb <- txtProgressBar(min = 0, max = nrow(eventsoutput), style = 3)
  rownumber <- 0
  apply(eventsoutput,1, function(x){
    rownumber <<- rownumber + 1
    setTxtProgressBar(pb,rownumber)
    chromosome <- as.character(x[1])
    startingpoint <- as.integer(x[2])
    stoppoint <- as.integer(x[3])
    averageBAF[with(averageBAF, Chr == chromosome & Position >= startingpoint & Position <= stoppoint),] <- NA
    averageBAF <<- averageBAF[complete.cases(averageBAF),]
    #Clearing the memory of junk
    suppressMessages(gc())
  })

  cat("\n")

  #The making of a correctionfactor (0.5/average BAF (0.2-0.8))
  correctionfactor <- 0.5/mean(averageBAF$B.Allele.Freq)

  #Making of array specific quality parameters
  percentageofSNPs <- round((nrow(averageBAF)/nrow(snpm.data))*100, digits = 2)
  meanaverageBAF <- paste(round(mean(averageBAF$B.Allele.Freq), digits = 4), sep = "")
  SNPs_used_in_averageBAF <- paste("[",formatC(nrow(averageBAF), big.mark = ","),"/",formatC(nrow(snpm.data), big.mark = ","),"][",percentageofSNPs,"%]", sep = "")
  sdaverageBAF <- round(sd(averageBAF$B.Allele.Freq), digits = 4)
  MADaverageBAF <- round(mad(x = averageBAF$B.Allele.Freq, center = 0.5), digits = 4)
  remove(averageBAF)

  ###For now: fully remove all allelic imbalance calls
  eventsoutput[eventsoutput$Event %in% c("LOH","Allelic Imbalance"),] <- NA
  eventsoutput <- eventsoutput[complete.cases(eventsoutput),]

  #Applying of minimum length filters:
  ##For number of probes: 10
  ##For length of event: 100kb
  ##For length of an allelic imbalance event: 5Mb
  #Deleting the Length and Probes columns afterwards
  eventsoutput[eventsoutput$Probes <= 10,] <- NA
  eventsoutput <- eventsoutput[complete.cases(eventsoutput),]
  eventsoutput[eventsoutput$Length <= 100000,] <- NA
  eventsoutput <- eventsoutput[complete.cases(eventsoutput),]
  eventsoutput[with(eventsoutput, Event %in% c("LOH", "Allelic Imbalance") & Length <= 5000000),] <- NA
  eventsoutput <- eventsoutput[complete.cases(eventsoutput),]
  eventsoutput$Length <- NULL
  eventsoutput$Probes <- NULL

  #Returning the number of rows that are about to be used in the apply calculation
  numberevents <- nrow(eventsoutput)
  line1 <- paste(numberevents, " events that are larger than 150kb, have more than 10 probes and LOH / AO areas larger than 5Mb", sep = "")
  cat(line1, "\n")
  cat("Busy filtering and calculating the remaining events, please wait ...","\n")

  #Turn on the graphical PDF output
  if ( missing(outfile)) {
    pdf("~/result.pdf")
  } else {
    pdf(outfile)
  }

  #Function that uses the remaining events to manipulate the SNP-array data
  eventsfilter <- apply(eventsoutput, 1, function(x) {
    process.deviation.event(x, snpm.data, correctionfactor, gender)
  })
  remove(eventsoutput)

  #Reset eventsfilter to vertical dataframe
  eventsfilter <- data.frame(matrix(unlist(eventsfilter), nrow = length(eventsfilter)/10, byrow = T), stringsAsFactors = F)

  #Returning the common array values
  result.lines <- NULL
  result.lines <- c(result.lines, "RESULTS:")
  result.lines <- c(result.lines, "")

  #Returning the number of rows that are about to be used in the apply calculation
  result.lines <- c(result.lines, paste(numberevents, "events that are larger than 150kb,"))
  result.lines <- c(result.lines, "have more than 10 probes and LOH / AO areas larger than 5Mb")
  result.lines <- c(result.lines, "")
  result.lines <- c(result.lines, c("SNPs outside of deviations:", SNPs_used_in_averageBAF))
  result.lines <- c(result.lines, "")
  result.lines <- c(result.lines, c("Average BAF: ", meanaverageBAF))
  result.lines <- c(result.lines, "")
  result.lines <- c(result.lines, c("STDev BAF: ", sdaverageBAF))
  result.lines <- c(result.lines, "")
  result.lines <- c(result.lines, c("Average deviation from BAF 0.5:", MADaverageBAF))
  result.lines <- c(result.lines, "")
  result.lines <- c(result.lines, c("Correction factor: ", correctionfactor))
  result.lines <- c(result.lines, "")
  OutputToPdf(cat(result.lines, sep = "\n" ))

  # Scale down the font size to make the table fix A4
  eventsfilter.theme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.4)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))

  cols <- c("Chromosome Region","CNV","Mean BAF","NormDist BAF","P>0.5","P<0.5","Skew>0.5","Skew<0.5","NormDist >0.5","NormDist <0.5")
  grb <- tableGrob(eventsfilter, cols = cols, rows = NULL, theme = eventsfilter.theme)
  grid.arrange(grb)

  #Turn off the graphical PDF output
  dev.off()

  #Clearing the memory of junk
  suppressMessages(gc())

  #Returning the final dataframe
  return(eventsfilter)
}
