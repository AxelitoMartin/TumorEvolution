#' get_results
#'
#' Get the results from a canopy optimized run
#' @param cna.obj
#' @param tree
#' @param projectname Name of the project to save results.
#' @param path Relative path to where results should be saved. Default is current directory.
#' @param path_mut Relative path to the mutation files used to run canopy
#' @param mut.files Mutation file names (same order as in the canopy run)
#' @param sample.names Sample names for the files
#' @param cex_cna size of font in CNA heatmap
#' @param epsilon Segment aggregation
#' @return Results pdf
#' @export
#' @examples
#' @import
#' pheatmap
#' gnomeR
#' RColorBrewer

get_results_canopy <- function(cna.obj, tree, projectname, path = ".",
                               path_mut = ".", mut.files,sample.names = NULL,
                               cex_cna = 5,epsilon = 0){

  pdf(paste0(path,"/",projectname,"_results.pdf"),width = 12)

  # dat_facets <- cna.obj$dat_facets
  # outcome = dat_facets$ID
  # names(outcome) <- dat_facets$ID
  # out <- facets.heatmap(seg = dat_facets,epsilon = 0,outcome = outcome,patients=dat_facets$ID)

  # original #
  temp <- facets.dat(seg = cna.obj$dat_facets,epsilon = epsilon)
  print(plot_cna(dat = temp$out.cn, outcome = rownames(temp$out.cn)))

  # filtered #
  temp <- facets.dat(seg = cna.obj$cncfs %>%
                       dplyr::rename(ID = sample, loc.start = start, loc.end = end),
                     epsilon = epsilon)
  out.cn <- as.data.frame(temp$out.cn) %>% dplyr::select(rownames(cna.obj$WM))
  print(plot_cna(dat = out.cn, outcome = rownames(out.cn)))
  # print(out$p)

  for(i in 1:length(mut.files)){
    assign(paste0("s",i,"_sna"),read.delim(paste0(path_mut,"/",mut.files[i])))
    # computeCCF(vaf = get(paste0("s",i,"_sna"))$alt/get(paste0("s",i,"_sna"))$total,
    #            tt = get(paste0("s",i,"_sna"))$tcn, minor = get(paste0("s",i,"_sna"))$lcn,
    #            purity = cna.obj$purity[i])

    assign(paste0("s",i,"_sna") , cbind(get(paste0("s",i,"_sna")),apply(get(paste0("s",i,"_sna")),1,function(r){
      get_ccf(alt = r["alt"], total = r["total"], tcn = r["tcn"], lcn = r["lcn"], purity = cna.obj$purity[i])$ccf
    }))
    )

    assign(paste0("s",i,"_sna"), get(paste0("s",i,"_sna"))%>%
             dplyr::rename(ccf = "apply(get(paste0(\"s\", i, \"_sna\")), 1, function(r) {") %>%
             dplyr::filter(!is.na(ccf)))
  }

  pos <- c()
  for(i in 1:length(mut.files)){
    pos <- c(pos,get(paste0("s",i,"_sna"))$pos)
  }
  pos <- unique(pos)
  CCF <- as.data.frame(matrix(ncol = length(mut.files), nrow = length(pos)))
  rownames(CCF) <- pos
  if(is.null(sample.names))
    colnames(CCF) <- mut.files
  else
    colnames(CCF) <- sample.names

  for(i in 1:nrow(CCF)){
    j <- rownames(CCF)[i]
    for(f in 1:length(mut.files)){
      temp <- get(paste0("s",f,"_sna")) %>% filter(pos == j)
      if(nrow(temp)>0){
        CCF[i,f] <- as.numeric(temp$ccf)
      }
    }
  }
  CCF[is.na(CCF)] <- 0

  pos <- c()
  for(i in 1:length(mut.files)){
    pos <- c(pos,paste0(get(paste0("s",i,"_sna"))$chr,"_",
                        get(paste0("s",i,"_sna"))$pos))

  }
  pos <- unique(pos)
  rownames(CCF) <- pos
  out <- pheatmap(t(CCF[order(apply(CCF,1,mean),decreasing = T),]),
                  fontsize_col = 4.2,cluster_cols = F,main = "CCF Heatmap")
  out
  temp.order <- rownames(CCF[order(apply(CCF,1,mean),decreasing = T),])

  canopy_plottree_mod(tree = tree, save = T, rdata.name = paste0(path,"/",projectname,"_muts.Rdata"))

  load(paste0(path,"/",projectname,"_muts.Rdata"))

  sna.mut <- lapply(mut.list,function(x){
    temp <- gsub(" ","",strsplit(x,split = ",")[[1]])
    if(length(grep("chr",temp))>0)
      temp <- temp[-grep("chr",temp)]
    return(temp)
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
      if(!all(is.na(y))) return(rep(y[!is.na(y)],length(mut.files)))
      else return(rep(length(mut.list)+1,length(mut.files)))
    }))
  colnames(replace.mut) <- colnames(temp)
  rownames(replace.mut) <- rownames(temp)

  temp.color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

  out.mut <- pheatmap(replace.mut[out$tree_row$order,match(temp.order, colnames(replace.mut))],
                      fontsize_col = 4.2,cluster_rows = F,cluster_cols = F,main = "Branch Heatmap",
                      color = temp.color[length(temp.color):1])
  out.mut

  cna.mat <- lapply(1:length(mut.list), function(x){
    temp <- gsub(" ","",strsplit(mut.list[[x]],split = ",")[[1]])
    temp <- temp[grep("chr",temp)]
    final <- as.data.frame(matrix(rep(x, length(sample.names)*length(temp)),
                                  nrow = length(sample.names), ncol = length(temp)))
    rownames(final) <- sample.names
    colnames(final) <- temp
    return(final)
  })

  cna.mat <- do.call('cbind',cna.mat)
  if(ncol(cna.mat) > 0)
    pheatmap(cna.mat,color = temp.color[length(temp.color):1],
             fontsize_col = cex_cna,main = "CNA Split Heatmap",cluster_rows = F,cluster_cols = F)

  out <- pheatmap(t(tree$CCF)[out$tree_row$order,order(apply(CCF,1,mean),decreasing = T)],fontsize_col = 4.2,main = "Tree CCF Heatmap",cluster_rows = F,cluster_cols = F)
  out

  dev.off()
}
