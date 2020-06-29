###################################################
############### help functions ####################
###################################################

#' Compute CCF and multiplicity issue in major and minor CN.
#' @param vaf
#' @param tt
#' @param minor
#' @param purity
#' @param multiplicity
#' @return ccf, rawccf, out, multiplicity
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

computeCCF=function(vaf, tt, minor, purity, multiplicity=NULL){

  bs=vaf
  cg=tt
  major=tt-minor
  vlen=length(major)

  out=diff1=rep(list(NA),vlen)
  for(i in 1:vlen){
    #multiplicity. e.g, major-minor:3-1 multiplicity include 1, 2, 3
    majori=major[i]
    if(!is.na(majori)){
      mm=1:major[i]
      out[[i]]=bs[i]*(2*(1-purity[i])+cg[i]*purity[i])/mm/purity[i]
      diff1[[i]]=abs(out[[i]]-1) #distance to 1
    }
  }


  if(is.null(multiplicity)){
    column.list=lapply(diff1,which.min)
    #catch integer(0) to maintain list length, NA allowed when unlist
    column.list[sapply(column.list, function(x)length(x)==0)] = NA
    column.list=unlist(column.list)
  }else{
    column.list = multiplicity
  }

  #select the multiplicity that gives ccf closest to 1
  #ccf=out[cbind(1:nrow(out), column.list)]

  rawccf=unlist(lapply(1:length(out),function(x)out[[x]][column.list[x]]))

  ccf=rawccf
  ccf[major==0]=NA
  ccf[ccf>1]=1
  ccf[is.na(minor)]=NA

  return(list(ccf=ccf, rawccf=rawccf, out=out, multiplicity=column.list))
}

##################################################################

#' Confidence interval for CCF at a default alpha of 0.05 or 95% CI
#' @param alt
#' @param ref
#' @param tt
#' @param minor
#' @param purity
#' @param multiplicity
#' @param alpha
#' @return lower, upper
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

confCCF=function(alt, ref, tt, minor, purity, multiplicity, alpha=0.05){

  depth=alt+ref
  ci=binconf(alt, depth, alpha, method=c("exact"))
  lowerci=ci[,2]
  upperci=ci[,3]

  lower=computeCCF(vaf=lowerci, tt, minor, purity, multiplicity)$ccf
  upper=computeCCF(vaf=upperci, tt, minor, purity, multiplicity)$ccf

  return(list(lower=lower, upper=upper))

}


###################################################

#' Wrapper function to get necessary things from CCF
#' @param alt
#' @param total
#' @param tcn
#' @param lcn
#' @param purity
#' @return ccf, multiplicity, clonal_status, ccfCI
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

get_ccf<-function(alt, total, tcn,lcn, purity){

  fit=computeCCF(vaf=alt/total, tt=tcn, minor=lcn, purity=purity)
  #cancer cell fraction for the mutation
  ccf=round(fit$ccf,2) #assumes mut on minor copy
  #mutant copy number (this is the v you want)
  multiplicity=fit$multiplicity

  #calculate confidence interval for CCF
  conf=confCCF(alt=alt, ref=total-alt, tt=tcn, minor=lcn, purity=purity, multiplicity=multiplicity)
  lower.ccf=round(conf$lower,2)
  upper.ccf=round(conf$upper,2)
  #assign clonal status for each mutation: if the lower bound of 95% CI for the observed CCF is >= 75%, call it clonal mutation.
  clonal_status=ifelse(lower.ccf >= 0.75,"clonal","subclonal")
  clonal_status[lower.ccf < 0.75 & ccf >= 0.80]="likely clonal"
  ccfCI = paste0("[",lower.ccf, "-", upper.ccf,"]")
  return(list(ccf=ccf, multiplicity=multiplicity, clonal_status=clonal_status, ccfCI = ccfCI))

}


###################################################

#' Calculate expected germline and somatic VAFs along with their CI given tcn and multiplicity
#' @param alt
#' @param total
#' @param tcn
#' @param lcn
#' @param purity
#' @return ccf, multiplicity, clonal_status, ccfCI
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

get_expectedVAF<-function(multiplicity, tcn, purity, alt, total, alpha = 0.05) {

  vafg = ((purity*multiplicity)+ (1-purity))/((purity*tcn)+2*(1-purity))
  vafs = (purity*multiplicity)/((purity*tcn)+2*(1-purity))

  z = abs(qnorm(alpha/2))

  sigmag=sqrt(vafg*(1-vafg)/total)
  upperg=round(vafg + (z*sigmag),2); lowerg=round(vafg - (z*sigmag),2)

  sigmas=sqrt(vafs*(1-vafs)/total)
  uppers=round(vafs + (z*sigmas),2); lowers= round(vafs - (z*sigmas),2)
  flag.germline = ifelse( ((lowerg <= (alt/total)) &((alt/total) <= upperg)), TRUE, FALSE)

  return(list(vafg = round(vafg,2), vafgCI = paste0("[", lowerg, "-", upperg, "]"), vafs = round(vafs,2), vafsCI = paste0("[", lowers, "-", uppers, "]"), flag.germline=flag.germline))
}

##################################

#' Split confCI character
#' @param x
#' @return xx.l, xx.u
#' @export
#' @examples
#' @import
#' Hmisc
#' plyr

get_confBound<-function(x){

  xx = gsub("\\[|\\]","", x)
  xx.l = unlist(lapply(strsplit(xx, "-"), function(y) y[1]))
  xx.u = unlist(lapply(strsplit(xx,"-"), function(y) y[2]))
  return(list(xx.l=xx.l, xx.u=xx.u))
}


####################################

#' Split confCI character
#' @param tree
#' @param pdf
#' @param pdf.name
#' @param txt
#' @param txt.name
#' @param save
#' @param rdata.name
#' @return tree, mutation list
#' @export
#' @examples
#' @import
#' Canopy

canopy.plottree.mod <- function (tree, pdf = NULL, pdf.name = NULL, txt = NULL, txt.name = NULL,save = FALSE, rdata.name = NULL)
{
  if (is.null(pdf)) {
    pdf = FALSE
  }
  if (is.null(txt)) {
    txt = FALSE
  }
  if (pdf & is.null(pdf.name)) {
    stop("pdf.name has to be provided if pdf = TRUE!")
  }
  if (txt & is.null(txt.name)) {
    stop("txt.name has to be provided if txt = TRUE")
  }
  if (!is.null(pdf.name)) {
    pdf.split = strsplit(pdf.name, "\\.")[[1]]
    if (length(pdf.split) < 2 | pdf.split[2] != "pdf") {
      stop("pdf.name has to end with .pdf!")
    }
  }
  if (pdf) {
    pdf(file = pdf.name, height = 6, width = 6)
  }
  nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = c(3,
                                                                  3, 3), heights = c(1.3, 1, 1), respect = TRUE)
  par(mar = c(1, 7, 1, 10))
  K = ncol(tree$Z)
  plot(tree, label.offset = 0.1, type = "cladogram", direction = "d",
       show.tip.label = FALSE)
  nodelabels()
  tiplabels()
  snaedge = rep(NA, nrow(tree$sna))
  for (k in 1:nrow(tree$sna)) {
    snaedge[k] = intersect(which(tree$edge[, 1] == tree$sna[k,
                                                            2]), which(tree$edge[, 2] == tree$sna[k, 3]))
  }
  if (!is.null(tree$cna)) {
    cnaedge = rep(NA, nrow(tree$cna))
    for (k in 1:nrow(tree$cna)) {
      cnaedge[k] = intersect(which(tree$edge[, 1] == tree$cna[k,
                                                              2]), which(tree$edge[, 2] == tree$cna[k, 3]))
    }
  }
  else {
    cnaedge = NULL
  }
  edge.label = sort(unique(c(snaedge, cnaedge)))
  edgelabels(paste("mut", 1:length(edge.label), sep = ""),
             edge.label, frame = "n", col = 2, cex = 1.2)
  tiplabels("Normal", 1, adj = c(0.2, 1.5), frame = "n", cex = 1.2,
            col = 4)
  tiplabels(paste("Clone", 1:(K - 2), sep = ""), 2:(K - 1),
            adj = c(0.5, 1.5), frame = "n", cex = 1.2, col = 4)
  tiplabels(paste("Clone", (K - 1), sep = ""), K, adj = c(0.8,
                                                          1.5), frame = "n", cex = 1.2, col = 4)
  par(mar = c(1, 7, 0.5, 9.5))
  P = tree$P
  image(1:nrow(P), 1:ncol(P), axes = FALSE, ylab = "", xlab = "",
        P, breaks = 0:100/100, col = tim.colors(100))
  axis(4, at = 1:ncol(P), colnames(P), cex.axis = 1.2, las = 1,
       tick = FALSE)
  abline(h = seq(0.5, ncol(P) + 0.5, 1), v = seq(0.5, nrow(P) +
                                                   0.5, 1), col = "grey")
  for (i in 1:nrow(P)) {
    for (j in 1:ncol(P)) {
      txt.temp <- sprintf("%0.3f", P[i, j])
      if (P[i, j] <= 0.05 | P[i, j] >= 0.95) {
        text(i, j, txt.temp, cex = 0.7, col = "white")
      }
      else {
        text(i, j, txt.temp, cex = 0.7)
      }
    }
  }
  sna.name = rownames(tree$sna)
  cna.name = rownames(tree$cna)
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n",
       xaxt = "n", yaxt = "n")
  txt.output = matrix(nrow = length(edge.label), ncol = 1)
  for (i in 1:length(edge.label)) {
    txt.temp = paste("mut", i, ": ", paste(c(sna.name[which(snaedge ==
                                                              edge.label[i])], cna.name[which(cnaedge == edge.label[i])]),
                                           collapse = ", "), sep = "")
    text(x = 0, y = 0.95 - 0.1 * (i - 1), txt.temp, pos = 4,
         cex = 1.2)
    txt.output[i, 1] = txt.temp
  }
  if (txt) {
    write.table(txt.output, file = txt.name, col.names = FALSE,
                row.names = FALSE, quote = FALSE, sep = "\t")
  }

  if(save){

    mut.list = list()
    for (i in 1:length(edge.label)) {
      txt.temp = paste( paste(c(sna.name[which(snaedge ==
                                                 edge.label[i])], cna.name[which(cnaedge == edge.label[i])]),
                              collapse = ", "), sep = "")
      text(x = 0, y = 0.95 - 0.1 * (i - 1), txt.temp, pos = 4,
           cex = 1.2)
      # txt.output[i, 1] = txt.temp
      mut.list[[i]] <- txt.temp
    }
    save(mut.list, file = rdata.name)
  }

  if (!is.null(pdf.name)) {
    text(x = 0.5, y = 0.1, pdf.split[1], font = 2, cex = 1.2)
  }
  if (pdf) {
    dev.off()
  }
  par(mfrow = c(1, 1))
}



