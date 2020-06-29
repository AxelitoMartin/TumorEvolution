#' Plot expected germline and somatic VAF along with observed somatic VAF.
#' And retrieve the the mutational status infered by the method.
#' Variants found to have an observed VAF below 2 standard deviations below the expected value are considered "subclonal".
#' Those falling in the +/- 2 standard deviation range from the expected VAF are annoatated "neutral".
#' Those with VAF greater than 2 starndard deviation above the expected value are annotated as "artifacts".
#' @param mat
#' @param purity
#' @param alpha
#' @param plot
#' @param sample
#' @return mat.out
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

germline_filter<-function(mat, purity, alpha=.05, plot=TRUE,sample){
  all<-list()
  for(i in 1:nrow(mat)){

    ccf = get_ccf(mat[i,"alt"], mat[i,"total"], mat[i,"tcn"],mat[i,"lcn"] ,purity)
    evaf = get_expectedVAF(ccf$multiplicity, mat[i,"tcn"], purity, mat[i,"alt"], mat[i,"total"], alpha)
    all[[i]]<-c(ccf, evaf)
  }

  mat.out<-ldply(all, data.frame)

  if(plot){
    n = nrow(mat.out)
    uppers = as.numeric(sapply(mat.out$vafsCI, function(x) get_confBound(x)$xx.u))

    lowers = as.numeric(sapply(mat.out$vafsCI, function(x) get_confBound(x)$xx.l))

    upperg = as.numeric(sapply(mat.out$vafgCI, function(x) get_confBound(x)$xx.u))

    lowerg = as.numeric(sapply(mat.out$vafgCI, function(x) get_confBound(x)$xx.l))


    plot(1:n,mat.out$vafg, type="o",ylim=c(0,1), bty="l", xlab="mutations index",
         ylab="vaf", cex.axis =1.2, cex.lab=1.2,main = paste0(sample, "(", purity,")"))
    lines(1:n,upperg, type="l",lty=2)
    lines(1:n,lowerg, type="l", lty=2)

    lines(1:n,mat.out$vafs, type="o", col="cadetblue")
    lines(1:n,uppers, type="l",lty=2, col="cadetblue")
    lines(1:n,lowers, type="l", lty=2, col="cadetblue")

    lines(1:n, mat[,"alt"]/mat[,"total"], type="o", lwd=1,
          col=ifelse(mat[,"alt"]/mat[,"total"] > uppers,"black",ifelse(mat[,"alt"]/mat[,"total"] < lowers,"grey","indianred")),
          pch=ifelse(mat[,"alt"]/mat[,"total"] > uppers,22,ifelse(mat[,"alt"]/mat[,"total"] < lowers,22,8)))
    legend("topright", c("expected germline", "expected somatic", "observed VAF",paste0((1-alpha)*100, "% CI"),
                         "high VAF", "low VAF"),
           col=c("black", "cadetblue", "indianred", "grey","black","grey"), lty=c(1,1,1,1,1,1), bty="n", pch=c(1,1,8,NA,22,22), lwd=2)

    mat$ccf <- mat.out$ccf
    mat$expected <- mat.out$vafs
    mat$lowers <- lowers
    mat$uppers <- uppers
    mat$clonal_status <- ifelse(mat[,"alt"]/mat[,"total"] > uppers,
                                "Artifact",
                                ifelse(mat[,"alt"]/mat[,"total"] < lowers,
                                       "Sub-clonal",
                                       "Neutral"))
  }

  else{
    n = nrow(mat.out)
    uppers = as.numeric(sapply(mat.out$vafsCI, function(x) get_confBound(x)$xx.u))
    lowers = as.numeric(sapply(mat.out$vafsCI, function(x) get_confBound(x)$xx.l))
    upperg = as.numeric(sapply(mat.out$vafgCI, function(x) get_confBound(x)$xx.u))
    lowerg = as.numeric(sapply(mat.out$vafgCI, function(x) get_confBound(x)$xx.l))

    mat$ccf <- mat.out$ccf
    mat$expected <- mat.out$vafs
    mat$lowers <- lowers
    mat$uppers <- uppers
    mat$obs <- mat[,"alt"]/mat[,"total"]
    mat$clonal_status <- ifelse(mat[,"alt"]/mat[,"total"] > uppers,
                                "Artifact",
                                ifelse(mat[,"alt"]/mat[,"total"] < lowers,
                                       "Sub-clonal",
                                       "Neutral"))
  }
  # return(mat.out)
  return(mat)


}
