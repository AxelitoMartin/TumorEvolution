#' get_results
#'
#' Get the results from a canopy optimized run
#' @param cna.obj
#' @param tree
#' @param projectname Name of the project to save results.
#' @param path Relative path to where results should be saved. Default is current directory.
#' @return Results pdf
#' @export
#' @examples
#' @import
#' pheatmap
#' gnomeR

get_results_canopy <- function(cna.obj, tree, projectname, path = "."){

  pdf(paste0(path,"/",projectname,"_results.pdf"),width = 12)

  dat_facets <- cna$dat_facets
  outcome = dat_facets$ID
  names(outcome) <- dat_facets$ID
  out <- facets.heatmap(seg = dat_facets,epsilon = 0,outcome = outcome,patients=dat_facets$ID)
  out$p

  out <- pheatmap(t(tree$CCF),fontsize_col = 4.2)
  out
  canopy_plottree_mod(tree = tree, save = T, rdata.name = paste0(path,"/",projectname,"_muts.Rdata"))

  load(paste0(path,"/",projectname,"_muts.Rdata"))

  sna.mut <- lapply(mut.list,function(x){
    temp <- gsub(" ","",strsplit(x,split = ",")[[1]])
    temp <- temp[-grep("chr",temp)]
  })

  temp <- t(tree$CCF)

  replace.mut <- do.call('cbind',lapply(
    lapply(colnames(temp), function(x){
      ind <- c()
      for(i in 1:length(sna.mut)){
        ii <- grep(x,sna.mut[[i]])
        if(length(ii) != 0) ind[i] <- i
        else ind[i] <- NA
      }
      return(ind)
    }), function(y){
      return(rep(y[!is.na(y)],length(mut.files)))
    }))
  colnames(replace.mut) <- colnames(temp)
  rownames(replace.mut) <- rownames(temp)

  out.mut <- pheatmap(replace.mut[out$tree_row$order,out$tree_col$order],fontsize_col = 4.2,cluster_rows = F,cluster_cols = F)
  out.mut

  dev.off()
}
