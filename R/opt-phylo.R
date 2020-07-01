#' Get optimal phylogeny
#' @param sampchain An output object from the run_canpy() function.
#' @param projectname Name of the project to save results.
#' @param K Vector of the number of clones to be used.
#' @param numchain Number of chains for the MCMC. Default is 15.
#' @param path Relative path to where results should be saved. Default is current directory.
#' @param burnin Burn-in period for the MCMC. Default is 10.
#' @param thin Thining parameter for the MCMC. Default is 5.
#' @param cutoff post configuration cutoff. Default is 0.05
#' @return Optimal phylogeny.
#' @export
#' @examples
#' @import
#' dplyr
#' dtplyr
#' Canopy


opt_phylo <- function(sampchain, projectname, K, numchain = 15,
                      path = ".", burnin = 10, thin = 5, cutoff = 0.05){
  bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, pdf = FALSE)
  optK = K[which.max(bic)]

  post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                     numchain = numchain, burnin = burnin, thin = thin,
                     optK = optK, post.config.cutoff = cutoff)
  samptreethin = post[[1]]   # list of all post-burnin and thinning trees
  samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
  config = post[[3]]
  config.summary = post[[4]]

  # choose the configuration with the highest posterior likelihood
  config.i = config.summary[which.max(config.summary[,3]),1]
  output.tree = canopy.output(post, config.i, C=NULL)
  return(output.tree)
}
