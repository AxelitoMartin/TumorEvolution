#' Run Canopy
#' @param sna_obj An output object from the make_sna_mat() function.
#' @param cna_obj An output object from the make_cna_mat() function.
#' @param Y An output object from the make_Y_mat() function.
#' @param projectname Name of the project to save results.
#' @param K Vector of the number of clones to be used.
#' @param numchain Number of chains for the MCMC. Default is 15.
#' @param max.simrun Maximum number of simulations to perform before convergence. Default is 15000.
#' @param min.simrun Minimum number of simulations to perform before convergence. Default is 5000.
#' @param path Relative path to where results should be saved. Default is current directory.
#' @param burnin Burn-in period for the MCMC. Default is 100.
#' @param thin Thining parameter for the MCMC. Default is 5.
#' @return Y
#' @export
#' @examples
#' @import
#' dplyr
#' dtplyr


run_canopy <- function(sna_obj, cna_obj, Y, projectname, K, numchain = 15, max.simrun = 15000, min.simrun = 5000,
                       path = ".", burnin = 100, thin = 5,parallel = FALSE){

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

  bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, pdf = FALSE)
  optK = K[which.max(bic)]

  post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                     numchain = numchain, burnin = burnin, thin = thin,
                     optK = optK, post.config.cutoff = 0.05)
  samptreethin = post[[1]]   # list of all post-burnin and thinning trees
  samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
  config = post[[3]]
  config.summary = post[[4]]

  # choose the configuration with the highest posterior likelihood
  config.i = config.summary[which.max(config.summary[,3]),1]
  output.tree = canopy.output(post, config.i, C=NULL)
  # pdf.name = paste(projectname, '_config_highest_likelihood.pdf', sep='')
  # canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
  # canopy.plottree(output.tree, pdf = FALSE,txt = TRUE,txt.name = "results/MSK0005_example_muts.txt")
  canopy.plottree.mod(output.tree, pdf = FALSE,save = TRUE,rdata.name = paste0(path,"/",projectname,"_muts.Rdata"))
  save(output.tree,file = paste0(path,"/",projectname,".Rdata"))

}
