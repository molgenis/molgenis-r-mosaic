#' Calculate mosaic percentage
#' 
#' Takes a experiment number, gender, snp dataframe and events dataframe
#' and calculates the percentage of mosaic dna
#' 
#' @param exp.nr numbarical experiment number
#' @param gender string denoting gender
#'  ( "Male", "Female" or "Unknown")
#' @param snpm.data dataframe containing snp data
#'  ( "Chr", "Position", "B Allele Freq")
#' @param deviations dataframe containing deviations
#' ( "Chromosome Region",	"Event",	"Length",	"Cytoband",	precentage of CNV Overlap",
#' 	"Probe Median",	"precentage Heterozygous",	"Probes",	"Count of Gene Symbols")
#' @param outfile optional path to direct pdf output to, if not set './Rplots.pdf' is used
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
  print("Dank u voor het kiezen van Mozaiek Bereken Tool, ontwikkeld door Robert Sietsma")
  #Set the option of the output decimal to comma
  options(OutDec = ",")
  
  cat("Ombouwen van deviations naar bruikbare parameters.","\n")
  eventsoutput <- process.deviation.line(deviations, gender)

  cat("Berekenen van kwaliteit aspecten array.","\n")

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
  #Overwriting averageBAF dataframe to exclude SNPs that are within deviations
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
  SNPs_used_in_averageBAF <- paste("[",formatC(nrow(averageBAF), big.mark = "."),"/",formatC(nrow(snpm.data), big.mark = "."),"][",percentageofSNPs,"%]", sep = "")
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
  aantalevents <- nrow(eventsoutput)
  line1 <- paste(aantalevents, " events die groter zijn dan 150kb, meer dan 10 probes hebben en LOH/AO gebieden groter dan 5Mb.", sep = "")
  cat(line1, "\n")
  cat("Bezig met het filteren en berekenen van de overgebleven events, een ogenblik geduld alstublieft...","\n")

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
  # colnames(eventsfilter) <- c("Chromosoom Regio","CNV","Gem. BAF","N.V. BAF","P>0.5","P<0.5","Skew>0.5","Skew<0.5","N.V. >0.5","N.V. <0.5")

  #Returning the common array values
  line2 <- paste0("SNPs buiten afwijkingen: ",SNPs_used_in_averageBAF)
  line3 <- paste0("Gemiddelde BAF: ",meanaverageBAF)
  line4 <- paste0("STDev BAF: ", sdaverageBAF)
  line5 <- paste0("Gemiddelde afwijking van BAF 0.5: ", MADaverageBAF)
  line6 <- paste0("Correctiefactor: ", correctionfactor)
  
  print(line1)
  print(line2)
  print(line3)
  print(line4)
  print(line5)
  print(line6)

  plot(NA, xlim = c(0, 6), ylim = c(0, 6), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  text(0.1, 6, line1, pos = 4, cex = 0.5)
  text(0.1, 5, line2, pos = 4, cex = 0.5)
  text(0.1, 4, line3, pos = 4, cex = 0.5)
  text(0.1, 3, line4, pos = 4, cex = 0.5)
  text(0.1, 2, line5, pos = 4, cex = 0.5)
  text(0.1, 1, line6, pos = 4, cex = 0.5)
  
  # Scale down the font size to make the table fix A4
  eventsfilter.theme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.4)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))
  
  cols <- c("Chromosoom Regio","CNV","Gem. BAF","N.V. BAF","P>0.5","P<0.5","Skew>0.5","Skew<0.5","N.V. >0.5","N.V. <0.5")
  grb <- tableGrob(eventsfilter, cols = cols, rows = NULL, theme = eventsfilter.theme)
  grid.arrange(grb)
  
  #Turn off the graphical PDF output
  dev.off()
  
  #Clearing the memory of junk
  suppressMessages(gc())
  
  #Returning the final dataframe
  return(eventsfilter)
}
