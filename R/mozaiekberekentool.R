#' @param exp.nr numbarical experiment number
#' @param gender string denoting gender ( "Male", "Female" or "Unknown")
#' @param snpm.data dataframe containing snp data ( "Chr", "Position", "B Allele Freq")
#' @param deviations dataframe containing deviations( "Chromosome Region"	"Event"	"Length"	"Cytoband"	"% of CNV Overlap"	"Probe Median"	"% Heterozygous"	"Probes"	"Count of Gene Symbols")
#'
#Versionnumber: 0.6.0.0
MosaicCalculator <- function(exp.nr, gender, snpm.data, deviations){
  #Set the option of the output decimal to comma
  options(OutDec = ",")

  #Building of EXP-n number and path to SNPanalysis Experiments
  EXPnumber <- paste("EXP-",exp.nr,sep = "")
  #exp.files.path <- "I:/SNPanalyse/NexusAnalyses/Experiments/"
  cat("Dank u voor het kiezen van Mozaiek Bereken Tool, ontwikkeld door Robert Sietsma","\n")
  #cat("Zoeken naar bestanden.","\n")

  #Search EXPnumber in exp.files.path location
  #MatchEXP <- list.files(exp.files.path, EXPnumber)

  #Saving of folder name for critical information
  #TXTName <- strsplit(MatchEXP, "_")[[1]][1:4]

  #Making of the PDF output directory and filename (for the histograms of all deviations)
  outputdir <- "~/rout"
  outputname <- "result.pdf"

  #Make the outputdir directory if it doesn't exist
  if(!dir.exists(outputdir)){
    dir.create(outputdir)
  }


  #Combinating of critical information for the SNP-array file
  #TXTName <- paste(TXTName[1], "_", TXTName[2],"_",TXTName[3], "_", TXTName[4], ".txt", sep = "")

  #Search for "events.txt" in MatchEXP path, if 2 hits occur, search for "copy_events.txt"
  #eventstxt <- list.files(paste(exp.files.path, MatchEXP, sep = ""), "events.txt")
  #if(length(eventstxt) > 1){
    #Save the very first copy_events hit (in case of copy_events and copy copy_events, copy copy_events gets used)
   # eventstxt <- list.files(paste(exp.files.path, MatchEXP, sep = ""), "copy_events.txt")[1]
  #}

  #The making of the path to: (copy_)events.txt and the SNP-array data
  #path <- paste(exp.files.path, MatchEXP, "/", TXTName, sep = "")
 # nexus <- paste(exp.files.path, MatchEXP, "/", eventstxt, sep = "")
  cat("Bestanden gevonden.","\n")
  cat("Inlezen van bestanden, een ogenblik geduld alstublieft...","\n")

  #Reading of SNP-array file and (copy_)events.txt
  #In the SNP-array file only column 3,4 and 6 are of interest
  #snpm.data <- read.table(path, skip = 9, header = T, sep = "\t")[,c(3,4,6)]
  #deviations <- read.table(nexus, header = T, sep = "\t")

  #The removal of variables that are not used anymore
  #remove( path, nexus, TXTName, MatchEXP, exp.files.path, EXPnumber)

  cat("Bestanden ingeladen.","\n")
  #cat("Ombouwen van ",eventstxt," naar bruikbare parameters.","\n")

  #Function to rebuild (copy_)events.txt to usefull information
  events <- function(x){
    #Selection of usefull information
    events <- cbind(x[,1:3],x$Probes)
    colnames(events) <- c("Chromosome.Region","Event","Length","Probes")

    #Removing "chr" in front of the chromosome and the comma in the start and stop
    events$Chromosome.Region <- gsub("chr","",events$Chromosome.Region)
    events$Chromosome.Region <- gsub(",","",events$Chromosome.Region)

    #Splitsing of the chromosome and start-stop
    temp <- apply(events, 1, function(x){strsplit(x[1], ":")[[1]]})
    temp <- data.frame(matrix(unlist(temp), nrow = length(temp) / 2, byrow = T), stringsAsFactors = F)

    #Splitsing of the start and stop
    temp2 <- apply(temp,1,function(x){strsplit(x[2], "-")[[1]]})
    temp2 <- data.frame(matrix(unlist(temp2), nrow = length(temp2) / 2, byrow = T), stringsAsFactors = F)

    #Combining the chromosome, start, stop and event
    events <- cbind(temp[,1],temp2,events[,2:4])
    colnames(events) <- c("Chr","Start","Stop","Event","Length","Probes")
    remove(temp,temp2)

    #Alterating events on sex chromosomes in Male patients
    if(gender == "Male"){
      #Only the PAR areas are heterozygous, which are marked XY
      events$Chr <- gsub("Y|X","XY",events$Chr)
      #Removing allelic imbalance calls
      events[with(events, Chr == "XY" & Event %in% c("LOH","Allelic Imbalance")),] <- NA
      events <- events[complete.cases(events),]
      #If more than 1 event remains, only the first one will be used specifically
      if(nrow(events[events$Chr == "XY",]) > 1){
        events <- rbind(events[events$Chr != "XY",],events[events$Chr == "XY",][1,])
      }
      if(nrow(events[events$Chr == "XY",]) == 1){
        events[events$Chr == "XY",][,5] <- 100001
        events[events$Chr == "XY",][,6] <- 11
      }
    }
    if(gender == "Female"){
      gsub("XY","X",events$Chr)
    }

    #Setting Start and Stop as integer
    events$Start <- as.integer(events$Start)
    events$Stop <- as.integer(events$Stop)

    return(events)
  }

  #Calling events function
  eventsoutput <- events(deviations)

  cat("Berekenen van kwaliteit aspecten array.","\n")


  #The build of an outlier removal function through quantiles
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }

  #Saving parameters for homozygous BAF
  minbaf <- c(0,0.1)
  maxbaf <- c(1,0.9)

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
  cat(paste(aantalevents, " events die groter zijn dan 150kb, meer dan 10 probes hebben en LOH/AO gebieden groter dan 5Mb.", sep = ""),"\n")
  cat("Bezig met het filteren en berekenen van de overgebleven events, een ogenblik geduld alstublieft...","\n")

  #Turn on the graphical PDF output
  pdf(file = paste(outputdir, outputname, sep = ""))

  #Function that uses the remaining events to manipulate the SNP-array data
  eventsfilter <- apply(eventsoutput,1,function(x){
    #Making of used variables
    chr <- x[1]
    start <- as.integer(x[2])
    stop <- as.integer(x[3])
    CNVType <- x[4]

    #Manipulating Allelic Imbalance to be Loss Of Heterozygosity
    ##Only called when no matching copy number variation was called
    if(CNVType == "Allelic Imbalance"){
      CNVType <- "LOH"
    }

    #If event matches XY (male), reapply the start and stop to fully cover the XY area
    #(Ref start/stop: HumanCytoSNP 850K; Illumina)
    if(all(chr == "XY", gender == "Male")){
      start <- 0
      stop <- as.integer(155270561)
      if(CNVType == "Homozygous Copy Loss"){
        CNVType <- "CN Loss"
      }
    }

    #Making of Output Chromosome Region (OCR)
    OCR <- paste("chr",chr,":",formatC(start, big.mark = "."),"-",formatC(stop, big.mark = "."), sep = "")

    #Homozygous copy losses disturb the BAF signal
    if(CNVType == "Homozygous Copy Loss"){
      return(c(OCR,"Homozygote deletie","N.V.T.","N.V.T.","N.V.T.","N.V.T.","N.V.T.","N.V.T.","N.V.T.","N.V.T."))
    }
    else{
      #Manipulating the SNP-array file to only contain the matching chromosome.
      df1 <- snpm.data[snpm.data$Chr == chr,]
      df1 <- df1[complete.cases(df1),]

      #Filtering the SNP-array file to contain the SNPs that are within range of start and stop
      region_of_deviation <- c(as.integer(start), as.integer(stop))
      df1 <- df1[with(df1, Position <= max(region_of_deviation) & Position >= min(region_of_deviation)),]
      df1 <- df1[complete.cases(df1),]
      df1 <- df1[with(df1, B.Allele.Freq >= max(minbaf) & B.Allele.Freq <= min(maxbaf)),]
      df1 <- df1[complete.cases(df1),]

      #Check to control the amount of rows that remain
      if(any(nrow(df1) < 1,nrow(df1[df1$B.Allele.Freq > 0.5,]) <= 3, nrow(df1[df1$B.Allele.Freq < 0.5,]) <= 3)){
        #With deletions, its highly likely that all SNPs are filtered out because of a high percentage mosaicism
        if(CNVType %in% c("CN Loss", "LOH")){
          #Building of a percentage output: Number of Rows Top (NRT) and Number of Rows Bottom (NRB)
          NRT <- paste("[>90/T.W.B.S.][",nrow(df1[df1$B.Allele.Freq > 0.5,]),"]", sep = "")
          NRB <- paste("[>90/T.W.B.S.][",nrow(df1[df1$B.Allele.Freq < 0.5,]),"]", sep = "")

          #Make CNVType more clear
          if(CNVType == "CN Loss"){
            CNVType <- "Deletie"
          }
          return(c(OCR,CNVType,"N.V.T.","N.V.T.",NRT,NRB,"N.V.T.","N.V.T.","N.V.T.","N.V.T."))
        }
        else{
          #Else: there are not enough SNPs to accuratly calculate a percentage range
          NRT <- paste("[T.W.B.S.][",nrow(df1[df1$B.Allele.Freq > 0.5,]),"]", sep = "")
          NRB <- paste("[T.W.B.S.][",nrow(df1[df1$B.Allele.Freq < 0.5,]),"]", sep = "")

          #Make CNVType more clear
          if(CNVType == "CN Gain"){
            CNVType <- "Duplicatie"
          }
          if(CNVType == "High Copy Gain"){
            CNVType <- "Tetrasomie"
          }
          return(c(OCR,CNVType,"N.V.T.","N.V.T.",NRT,NRB,"N.V.T.","N.V.T.","N.V.T.","N.V.T."))
        }
      }
      else{
        #Splitsing the dataframe in 2 dataframes, one containing SNP with a BAF >0.5 and one <0.5 for outlier removing
        df1 <- df1[complete.cases(df1),]
        df2 <- df1[df1$B.Allele.Freq > 0.5,]
        df1 <- df1[df1$B.Allele.Freq < 0.5,]

        #Apply the quantile outlier remove filter
        df1$B.Allele.Freq <- remove_outliers(df1$B.Allele.Freq)
        df1 <- df1[complete.cases(df1),]
        df2$B.Allele.Freq <- remove_outliers(df2$B.Allele.Freq)
        df2 <- df2[complete.cases(df2),]

        #Combine the 2 dataframes to correct the BAF
        equaldf <- rbind(df1, df2)

        #Applying the correction factor to always correct towards 0.5
        MeanOfMutation <- c(mean(df2$B.Allele.Freq),mean(df1$B.Allele.Freq))
        MeanOfMutation <- mean(MeanOfMutation)

        #Remove df1 and df2 to make sure no numbers get mixed up
        remove(df1,df2)

        if(all(MeanOfMutation > 0.5, correctionfactor > 1)){
          equaldf$B.Allele.Freq <- equaldf$B.Allele.Freq / correctionfactor
        }
        if(all(MeanOfMutation < 0.5,correctionfactor < 1)){
          equaldf$B.Allele.Freq <- equaldf$B.Allele.Freq / correctionfactor
        }
        if(all(MeanOfMutation > 0.5,correctionfactor < 1)){
          equaldf$B.Allele.Freq <- equaldf$B.Allele.Freq * correctionfactor
        }
        if(all(MeanOfMutation < 0.5,correctionfactor > 1)){
          equaldf$B.Allele.Freq <- equaldf$B.Allele.Freq * correctionfactor
        }

        #Reapply the above and below 0.5 BAF dataframes
        df2 <- equaldf[equaldf$B.Allele.Freq > 0.5,]
        df1 <- equaldf[equaldf$B.Allele.Freq < 0.5,]

        #Create histogramms of the deviation
        #Create histogram CNVtype
        if(CNVType == "CN Gain"){
          hCNVType <- "Duplicatie"
        }
        if(CNVType == "CN Loss"){
          hCNVType <- "Deletie"
        }
        if(CNVType == "High Copy Gain"){
          hCNVType <- "Tetrasomie"
        }
        hist(equaldf$B.Allele.Freq, main = paste("[",chr,":",formatC(start, big.mark = "."),"-",formatC(stop, big.mark = ".")," ", hCNVType,"] [Gecombineerde BAF]", sep = ""), xlab = "BAF")
        hist(df2$B.Allele.Freq, main = paste("[",chr,":",formatC(start, big.mark = "."),"-",formatC(stop, big.mark = ".")," ", hCNVType,"] [BAF > 0.5]", sep = ""), xlab = "BAF")
        hist(df1$B.Allele.Freq, main = paste("[",chr,":",formatC(start, big.mark = "."),"-",formatC(stop, big.mark = ".")," ", hCNVType,"] [BAF < 0.5]", sep = ""), xlab = "BAF")


        #A high chanse occurs for memory leaks, apply gc()
        suppressMessages(gc())

        #Building of an equaltest for low percentage mosaicisms
        #Note: The statistical P values will always be rounded down!

        #More than 2000 rows: Kolmogorov-Smirnov
        if(nrow(equaldf) > 2000){
          equaltest <- suppressWarnings(ks.test(x = equaldf$B.Allele.Freq, "pnorm", 0.5, sd(equaldf$B.Allele.Freq))[1:2][2][[1]])
          equaltest <- floor(equaltest * 1000) / 1000
          equaltest <- c(equaltest, "[K]")
        }
        #Between 30 and 2000 rows: Shapiro-Wilk
        if(all(nrow(equaldf) <= 2000, nrow(equaldf[equaldf$B.Allele.Freq > 0.5,]) > 30, nrow(equaldf[equaldf$B.Allele.Freq < 0.5,]) > 30)){
          equaltest <- shapiro.test(equaldf$B.Allele.Freq)[1:2][2][[1]]
          equaltest <- floor(equaltest * 1000) / 1000
          equaltest <- c(equaltest, "[S]")
        }
        #Lower than 30: No accurate test possible
        if(all(nrow(equaldf) <= 2000, any(nrow(equaldf[equaldf$B.Allele.Freq > 0.5,]) <= 30, nrow(equaldf[equaldf$B.Allele.Freq <0.5,]) <= 30))){
          equaltest <- "N.V.T."
        }

        #Apply tag "Goed" (N.V. < 0.05) or "M.G." (N.V. >= 0.05) for the combined BAF of mutation
        if(equaltest[1] >= 0.05){
          outputequaltest <- paste("[M.G]",equaltest[2],"[",equaltest[1],"]", sep = "")
        }
        if(equaltest[1] < 0.05){
          outputequaltest <- paste("[Goed]",equaltest[2],"[",equaltest[1],"]", sep = "")
        }
        if(equaltest[1] == "N.V.T."){
          outputequaltest <- "N.V.T."
        }

        #Remove MeanOfMutation so numbers don't get mixed up
        remove(MeanOfMutation)

        #Saving the mean of the mean, with applied correction factor
        ##Mean of mean, because now the mean of 2 points (both means) are taken
        MeanOfMutation <- c(mean(df1$B.Allele.Freq),mean(df2$B.Allele.Freq))
        MeanOfMutation <- round(mean(MeanOfMutation),digits = 4)

        #Applying the Skewness, Shapiro-Wilk/Kolmogorov-Smirnov over the split dataframes
        #Same parameters apply here as in equaltest
        #Note: The statistical P values will always be rounded down!
        suppressWarnings(if(!require(e1071)){
          install.packages("e1071")
          library("e1071")
        })
        SkewnessPValueT <- round(skewness(df2$B.Allele.Freq), digits = 4)
        SkewnessPValueB <- round(skewness(df1$B.Allele.Freq), digits = 4)
        if(nrow(equaldf) <= 2000){
          if(any(nrow(df1) <= 30, nrow(df2) <= 30)){
            shapiro1 <- "N.V.T."
            shapiro2 <- "N.V.T."
          }
          else{
            shapiro1 <- shapiro.test(df2$B.Allele.Freq)[1:2][2][[1]]
            shapiro1 <- floor(shapiro1 * 1000) / 1000
            shapiro1 <- c(shapiro1, "[S]")
            shapiro2 <- shapiro.test(df1$B.Allele.Freq)[1:2][2][[1]]
            shapiro2 <- floor(shapiro2 * 1000) / 1000
            shapiro2 <- c(shapiro2, "[S]")
          }
        }
        if(nrow(equaldf) > 2000){
          shapiro1 <- suppressWarnings(ks.test(x = df2$B.Allele.Freq, "pnorm",mean(df2$B.Allele.Freq),sd(df2$B.Allele.Freq))[1:2][2][[1]])
          shapiro1 <- floor(shapiro1 * 1000) / 1000
          shapiro1 <- c(shapiro1, "[K]")
          shapiro2 <- suppressWarnings(ks.test(x = df1$B.Allele.Freq, "pnorm",mean(df1$B.Allele.Freq),sd(df1$B.Allele.Freq))[1:2][2][[1]])
          shapiro2 <- floor(shapiro2 * 1000) / 1000
          shapiro2 <- c(shapiro2, "[K]")
        }

        #Calculating the standard deviation over both dataframes
        stdevH <- sd(df2$B.Allele.Freq)
        stdevL <- sd(df1$B.Allele.Freq)

        #Calculating a 95% confidence interval over both dataframes
        #BI95N: confidence interval, percentage, H for High (>0.5) and L for Low (<0.5)
        BI95H <- qnorm(0.975)*stdevH/sqrt(nrow(df2))
        BI95L <- qnorm(0.975)*stdevL/sqrt(nrow(df1))
        StatsHigh95 <- c(mean(df2$B.Allele.Freq)-BI95H, mean(df2$B.Allele.Freq)+BI95H)
        StatsLow95 <- c(mean(df1$B.Allele.Freq)-BI95L, mean(df1$B.Allele.Freq)+BI95L)

        #Applying the correct formula for the current event
        if(CNVType == "CN Loss"){
          #Loading calculation for the BAF >0.5 for a deletion
          calcdel2 <- function(x){-279.71*x^2+609.65*x-232.73}
          #Loading calculation for the BAF <0.5 for a deletion
          calcdel <- function(x){-279.71*x^2-50.244*x+97.219}
          #Calculating percentage range
          P95H <- c(round(calcdel2(min(StatsHigh95))), round(calcdel2(max(StatsHigh95))))
          P95L <- c(round(calcdel(min(StatsLow95))), round(calcdel(max(StatsLow95))))
        }
        if(CNVType == "CN Gain"){
          #Loading calculation for the BAF >0.5 for a duplication
          calcdup2 <- function(x){1464*x^2-1118.2*x+193.99}
          #Loading calculation for the BAF <0.5 for a duplication
          calcdup <- function(x){1464*x^2-1809.8*x+539.83}
          #Calculating percentage range
          P95H <- c(round(calcdup2(min(StatsHigh95))), round(calcdup2(max(StatsHigh95))))
          P95L <- c(round(calcdup(min(StatsLow95))), round(calcdup(max(StatsLow95))))
        }
        if(CNVType == "High Copy Gain"){
          #Loading calculation for the BAF >0.5 for a tetrasomy
          calctetra2 <- function(x){1118.8*x^2-1018.3*x+232.24}
          #Loading calculation for the BAF <0.5 for a tetrasomy
          calctetra <- function(x){1118.8*x^2-1219.3*x+332.73}
          #Calculating percentage range
          P95H <- c(round(calctetra2(min(StatsHigh95))), round(calctetra2(max(StatsHigh95))))
          P95L <- c(round(calctetra(min(StatsLow95))), round(calctetra(max(StatsLow95))))
        }
        if(CNVType == "LOH"){
          #Loading calculation for the BAF >0.5 for a uniparental disomy event
          calchomo2 <- function(x){200*x-100}
          #Loading calculation for the BAF <0.5 for a uniparental disomy event
          calchomo <- function(x){-200*x+100}
          #Calculating percentage range
          P95H <- c(round(calchomo2(min(StatsHigh95))), round(calchomo2(max(StatsHigh95))))
          P95L <- c(round(calchomo(min(StatsLow95))), round(calchomo(max(StatsLow95))))
        }

        #If the minimum percentage is higher than 100 (which is only possible with a duplication or tetrasomy):
        #Calculate for tetrasomy
        if(min(c(P95H,P95L)) > 100){
          #Loading of calculations for a tetrasomy
          calctetra2 <- function(x){1118.8*x^2-1018.3*x+232.24}
          calctetra <- function(x){1118.8*x^2-1219.3*x+332.73}
          #Calculating percentage range
          P95H1 <- c(round(calctetra2(min(StatsHigh95))), round(calctetra2(max(StatsHigh95))))
          P95L1 <- c(round(calctetra(min(StatsLow95))), round(calctetra(max(StatsLow95))))

          #Reset percentages higher than 100 to 100
          P95H1[P95H1 > 100] <- 100
          P95L1[P95L1 > 100] <- 100

          #Check whenever the N.V. > 0.05 has a P-value lower than 0.05
          #If so, return with N.B tag
          if(shapiro1[1] < 0.05){
            #Check whenever percentages are the same
            if(min(P95H1) == max(P95H1)){
              #If they are, return only 1
              OPT <- paste("[N.B][",min(P95H1),"][",nrow(df2),"]", sep = "")
            }
            else{
              OPT <- paste("[N.B][",min(P95H1),"-",max(P95H1),"][",nrow(df2),"]",sep = "")
            }
          }
          #If not, return normally
          else{
            if(min(P95H1) == max(P95H1)){
              OPT <- paste("[",min(P95H1),"][",nrow(df2),"]", sep = "")
            }
            else{
              OPT <- paste("[",min(P95H1),"-",max(P95H1),"][",nrow(df2),"]",sep = "")
            }
          }

          #Check whenever the N.V. < 0.05 has a P-value lower than 0.05
          #If so, return with N.B tag
          if(shapiro2[1] < 0.05){
            #Check whenever percentages are the same
            if(min(P95L1) == max(P95L1)){
              #If they are, return only 1
              OPB <- paste("[N.B][",min(P95L1),"][",nrow(df1),"]", sep = "")
            }
            else{
              OPB <- paste("[N.B][",min(P95L1),"-",max(P95L1),"][",nrow(df1),"]",sep = "")
            }
          }
          #If not, return normally
          else{
            if(min(P95L1) == max(P95L1)){
              OPB <- paste("[",min(P95L1),"][",nrow(df1),"]", sep = "")
            }
            else{
              OPB <- paste("[",min(P95L1),"-",max(P95L1),"][",nrow(df1),"]",sep = "")
            }
          }

          #Combining the 2 shapiro variables to 1 separately
          if(length(shapiro1) > 1){
            #Combining the 2 shapiro variables to 1 separately
            shapiro1 <- paste(shapiro1[1],shapiro1[2], sep = "")
            shapiro2 <- paste(shapiro2[1],shapiro2[2], sep = "")
          }

          #Return output, with event altered to !High Copy Gain!
          return(c(OCR,"!Tetrasomie!",MeanOfMutation,outputequaltest,OPT,OPB,SkewnessPValueT,SkewnessPValueB,shapiro1,shapiro2))
        }
        else{

          #Reset percentages higher than 100 to 100
          P95H[P95H > 100] <- 100
          P95L[P95L > 100] <- 100

          #Check whenever P > 0.5 has a normal distribution
          if(shapiro1[1] < 0.05){
            #Check if it's a deletion with a minimum percentage of 80
            if(all(CNVType %in% c("CN Loss","LOH"),min(P95H) > 80)){
              OPT <- paste("[>",min(P95H),"][",nrow(df2),"]", sep = "")
            }
            else{
              #Check whenever the percentages are the same
              if(min(P95H) == max(P95H)){
                #If they are, return only 1
                P95H <- P95H[1]
                OPT <- paste("[N.B][",P95H,"][",nrow(df2),"]",sep = "")
              }
              else{
                #If they aren't, return minimum to maximum
                OPT <- paste("[N.B][",min(P95H),"-",max(P95H),"][",nrow(df2),"]",sep = "")
              }
            }
          }
          #Else it is normally distributed
          else{
            #Check whenever the percentages are the same
            if(min(P95H) == max(P95H)){
              #If they are, return only 1
              P95H <- P95H[1]
              OPT <- paste("[",P95H,"][",nrow(df2),"]",sep = "")
            }
            else{
              #If they aren't, return minimum to maximum
              OPT <- paste("[",min(P95H),"-",max(P95H),"][",nrow(df2),"]",sep = "")
            }
          }
          #Apply the same check for P < 0.5
          if(shapiro2[1] < 0.05){
            #Check whenever it's a deletion with a minimum percentage of 80
            if(all(min(P95L) > 80, CNVType %in% c("CN Loss","LOH"))){
              OPB <- paste("[>",min(P95L),"][",nrow(df1),"]", sep = "")
            }
            else{
              #Check if minimum and maximum percentage are the same
              if(min(P95L) == max(P95L)){
                P95L <- P95L[1]
                OPB <- paste("[N.B][",P95L,"][",nrow(df1),"]", sep = "")
              }
              else{
                OPB <- paste("[N.B][",min(P95L),"-",max(P95L),"][", nrow(df1),"]", sep = "")
              }
            }
          }
          #Else it is normally distributed
          else{
            if(min(P95L) == max(P95L)){
              P95L <- P95L[1]
              OPB <- paste("[",P95L,"][",nrow(df1),"]",sep = "")
            }
            else{
              OPB <- paste("[",min(P95L),"-",max(P95L),"][",nrow(df1),"]",sep = "")
            }
          }

          #Combining the 2 shapiro variables to 1 separately
          if(length(shapiro1) > 1){
            #Combining the 2 shapiro variables to 1 separately
            shapiro1 <- paste(shapiro1[1],shapiro1[2], sep = "")
            shapiro2 <- paste(shapiro2[1],shapiro2[2], sep = "")
          }

          #Reapply chromosome Y or X with analysis of Y or X with male patients
          if(all(gender == "Male", chr == "XY")){
            OCR <- paste("PAR(X/Y)", sep = "")
          }

          #Make CNV Type more clear
          if(CNVType == "CN Gain"){
            CNVType <- "Duplicatie"
          }
          if(CNVType == "CN Loss"){
            CNVType <- "Deletie"
          }
          if(CNVType == "High Copy Gain"){
            CNVType <- "Tetrasomie"
          }

          #Return output
          return(c(OCR,CNVType,MeanOfMutation,outputequaltest,OPT,OPB,SkewnessPValueT,SkewnessPValueB,shapiro1,shapiro2))
        }

      }

    }})
  remove(eventsoutput)

  #Turn off the graphical PDF output
  dev.off()

  #Reset eventsfilter to vertical dataframe
  eventsfilter <- data.frame(matrix(unlist(eventsfilter), nrow = length(eventsfilter)/10, byrow = T), stringsAsFactors = F)
  colnames(eventsfilter) <- c("Chromosoom Regio","CNV","Gem. BAF","N.V. BAF","P>0.5","P<0.5","Skew>0.5","Skew<0.5","N.V. >0.5","N.V. <0.5")

  #Returning the common array values
  cat("SNPs buiten afwijkingen:",SNPs_used_in_averageBAF,"\n")
  cat("Gemiddelde BAF:",meanaverageBAF,"\n")
  cat("STDev BAF: ", sdaverageBAF,"\n")
  cat("Gemiddelde afwijking van BAF 0.5: ", MADaverageBAF,"\n")
  cat("Correctiefactor: ", correctionfactor, "\n")

  #Clearing the memory of junk
  suppressMessages(gc())

  #Returning the final dataframe
  return(eventsfilter)
}
suppressMessages(  while (!is.null(dev.list()))  dev.off())
