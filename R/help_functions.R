###################################################
############### help functions ####################
###################################################

#' computeCCF
#'
#' Compute CCF and multiplicity issue in major and minor CN.
#'
#' @param vaf vaf
#' @param tt tt
#' @param minor minor
#' @param purity purity
#' @param multiplicity multiplicity
#' @return ccf, rawccf, out, multiplicity
#' @export
#'
#' @examples
#' \dontrun{
#' computeCCF(vaf, tt, minor, purity, multiplicity=NULL)
#' }
#' @import
#' Hmisc
#' dplyr

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

#' confCCF
#'
#' Confidence interval for CCF at a default alpha of 0.05 or 95% CI
#' @param alt alt
#' @param ref ref
#' @param tt tt
#' @param minor minor
#' @param purity purity
#' @param multiplicity multiplicity
#' @param alpha alpha
#' @return lower, upper
#' @export
#'
#' @examples
#' \dontrun{
#' computeCCF(vaf, tt, minor, purity, multiplicity=NULL)
#' }
#' @import
#' Hmisc
#' dplyr

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
#' @param alt alt
#' @param total total
#' @param tcn tcn
#' @param lcn lcn
#' @param purity purity
#' @return ccf, multiplicity, clonal_status, ccfCI
#' @export
#'
#' @examples
#' \dontrun{
#' computeCCF(vaf, tt, minor, purity, multiplicity=NULL)
#' }
#' @import
#' Hmisc
#' dplyr

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
#' \dontrun{
#' computeCCF(vaf, tt, minor, purity, multiplicity=NULL)
#' }
#' @import
#' Hmisc
#' dplyr

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
#' \dontrun{
#' computeCCF(vaf, tt, minor, purity, multiplicity=NULL)
#' }
#' @import
#' Hmisc
#' dplyr

get_confBound<-function(x){

  xx = gsub("\\[|\\]","", x)
  xx.l = unlist(lapply(strsplit(xx, "-"), function(y) y[1]))
  xx.u = unlist(lapply(strsplit(xx,"-"), function(y) y[2]))
  return(list(xx.l=xx.l, xx.u=xx.u))
}


####################################



