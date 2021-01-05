#' clone_signature
#'
#' Get the mutational signature of each clone in a canopy run.
#' @param clone_mut A list containing the variants observed in each variant of interest
#' @param maf The corresponding maf file to the samples used in canopy.
#' @param sample.names Sample names to consider in the maf file.
#' @param num_sign Number of signatures to look for.
#' @param projectname Name of the project to save results.
#' @param path Relative path to where results should be saved. Default is current directory.
#' @param num_parallelCores number of cores to use
#' @return Results pdf
#' @export
#' @examples
#' @import
#' mutSignatures
#' BSgenome.Hsapiens.UCSC.hg19
#' ggplot2
#' dplyr
#' dtplyr
#' reshape2
#' gridExtra


clone_signature <- function(clone_mut, maf, sample.names, num_sign = 5, projectname = "", path = ".", num_parallelCores = 1){

}
