#' OutputToPdf
#' 
#' Helper function to print the text output to a file
#' Takes a output object and prints it to the current output file.
#' 
#' @param output out put object
#' @param cex character expansion factor, default to 0.7
#' 
#' @importFrom graphics box plot.new
#' @importFrom utils capture.output
#'    
OutputToPdf <- function(output, cex = 0.7) {
  tmp <- capture.output(output)
  plot.new()
  text(0, 1, paste(tmp, collapse='\n'), adj = c(0,1), family = 'mono', cex = cex)
  box()
}