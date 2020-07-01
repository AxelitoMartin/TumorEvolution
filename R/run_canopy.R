#' Run Canopy
#' @param sna_obj An output object from the make_sna_mat() function.
#' @param cna_obj An output object from the make_cna_mat() function.
#' @param Y An output object from the make_Y_mat() function.
#' @param projectname Name of the project to save results.
#' @param K Vector of the number of clones to be used.
#' @param numchain Number of chains for the MCMC. Default is 15.
#' @param max.simrun Maximum number of simulations to perform before convergence. Default is 100000.
#' @param min.simrun Minimum number of simulations to perform before convergence. Default is 10000.
#' @param burnin Burn-in period for the MCMC. Default is 10.
#' @param thin Thining parameter for the MCMC. Default is 5.
#' @param parallel Should this be run in parallel. Default is FALSE.
#' @return Y
#' @export
#' @examples
#' @import
#' dplyr
#' dtplyr
#' Canopy


run_canopy <- function(sna_obj, cna_obj, Y, projectname, K, numchain = 15, max.simrun = 100000, min.simrun = 10000,
                       burnin = 10, thin = 5,parallel = FALSE){

  R <- sna_obj$R
  X <- sna_obj$X
  WM <- cna_obj$WM
  Wm <- cna_obj$Wm
  epsM <- cna_obj$epsM
  epsm <- cna_obj$epsm

  sampchain = canopy.sample(R = as.matrix(R), X = as.matrix(X), WM = as.matrix(WM),
                            Wm = as.matrix(Wm), epsilonM = as.matrix(epsM),
                            epsilonm = as.matrix(epsm),
                            C = NULL, Y = as.matrix(Y), K = K,
                            numchain = numchain, max.simrun = max.simrun,
                            min.simrun = min.simrun, writeskip = 200,
                            projectname = projectname, cell.line = FALSE,
                            plot.likelihood = TRUE)

  return(sampchain)

}
