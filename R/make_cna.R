#' Create a Canopy ready CNA matrix and their corresponding variance matrices
#' @param cna.files Names of the files to be read in and processed. Default will be all the files found in the working directory.
#' @param path relative path to the folder containing the samples to be processed. Default is current directory.
#' @param sample.names Custom names to be given to the samples. Default is NULL.
#' @param cval Default is 100.
#' @param snp.nbhd Default is 250.
#' @param epsilon Default is 0.005.
#' @return WM, Wm,epsM, epsm
#' @export
#' @examples
#' @import
#' facets
#' pctGCdata
#' dplyr
#' dtplyr
#' gnomeR

make_cna_mat <- function(cna.files = list.files(), path = ".", sample.names = NULL,
                         cval = 100, snp.nbhd = 250,epsilon = 0.005){

  cnas <- lapply(cna.files, function(x){
    rcmat <- readSnpMatrix(filename=paste0(path,"/",x))
    nor.dp <- rcmat$NOR.DP
    tum.dp <- rcmat$TUM.DP

    set.seed(13692)
    xx=preProcSample(rcmat,snp.nbhd = snp.nbhd)
    oo=procSample(xx,cval=cval)

    M = 2^(oo$out$cnlr.median+1)/(1 + 1/(sqrt(exp(oo$out$mafR))))
    m = 2^(oo$out$cnlr.median+1) - M

    fit=emcncf(oo)
    cncf=data.frame(ID=rep(x, nrow(fit$cncf)), chrom=fit$cncf$chrom,
                    loc.start=fit$cncf$start, loc.end=fit$cncf$end,
                    nhet=fit$cncf$nhet, num.mark=fit$cncf$num.mark, fit$cncf[, 5:14],
                    seg.mean = log2(fit$cncf$tcn/fit$ploidy + 1 * 10^(-6)))
    return(cncf)
  })

  cncfs <- as.data.frame(do.call('rbind',cnas)) %>%
    dplyr::rename(sample = ID) %>%
    mutate(chrom = as.numeric(as.character(chrom)),
           start = as.numeric(as.character(start)),
           end = as.numeric(as.character(end)),
           num.mark = as.numeric(as.character(num.mark)),
           seg.mean = as.numeric(as.character(seg.mean))) %>%
    select(sample, chrom, start, end, num.mark, seg.mean)
  out <- CNregions.mod(seg = cncfs, epsilon = 0)

  temp <- as.data.frame(do.call('rbind',cnas)) %>%
    select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)

  outcome = cna.files
  names(outcome) <- cna.files
  facets.heatmap(seg = temp,epsilon = 0,outcome = outcome,patients=cna.files)

  info <- lapply(cna.files, function(x){
    rcmat <- readSnpMatrix(filename=paste0(path,"/",x))
    nor.dp <- rcmat$NOR.DP
    tum.dp <- rcmat$TUM.DP

    set.seed(13692)
    xx=preProcSample(rcmat,snp.nbhd = snp.nbhd)

    adj.count <- sum(xx$jointseg$rCountN)/sum(xx$jointseg$rCountT)

    mm <- lapply(colnames(out),function(y){

      temp <- strsplit(y, split = "\\.|-")
      chrm <- as.numeric(gsub("chr","",temp[[1]][1]))
      start <- as.numeric(temp[[1]][2])
      end <- as.numeric(temp[[1]][3])

      temp <- xx$jointseg %>%
        filter(chrom == chrm,
               maploc >= start,
               maploc <= end,
               het == 1) %>%
        mutate(
          M_T = ifelse(vafT > 0.5, rCountT * vafT, rCountT - rCountT * vafT),
          M_N = ifelse(vafN > 0.5, rCountN * vafN, rCountN - rCountN * vafN),
          m_T = rCountT - M_T,
          m_N = rCountN - M_N
        )

      if(nrow(temp) > 4){

        Nt = nrow(temp)

        W_M = 1/Nt * sum(temp$M_T/temp$M_N) * adj.count
        eps_M = sqrt(( sum((temp$M_T/temp$M_N)^2) - Nt*W_M^2 )/(Nt *(Nt-1)))

        W_m = 1/Nt * sum(temp$m_T/temp$m_N) * adj.count
        eps_m = sqrt(( sum((temp$m_T/temp$m_N)^2) - Nt*W_m^2 )/(Nt *(Nt-1)))
      }
      else{
        W_M  = W_m = eps_M = eps_m = NA
      }
      return(list("W_M" = W_M, "W_m" = W_m, "eps_M" = eps_M, "eps_m" = eps_m))
    })
    W_M <- sapply(mm,"[[","W_M")
    W_m <- sapply(mm,"[[","W_m")
    eps_M <- sapply(mm,"[[","eps_M")
    eps_m <- sapply(mm,"[[","eps_m")
    return(list("W_M" = W_M, "W_m" = W_m, "eps_M" = eps_M, "eps_m" = eps_m))
  })

  WM <- sapply(info,"[[","W_M")
  Wm <- sapply(info,"[[","W_m")
  epsM <- sapply(info,"[[","eps_M")
  epsm <- sapply(info,"[[","eps_m")

  if(is.null(sample.names))
    colnames(WM) <- colnames(Wm) <- colnames(epsM) <- colnames(epsm) <- cna.files
  else
    colnames(WM) <- colnames(Wm) <- colnames(epsM) <- colnames(epsm) <- sample.names
  rownames(WM) <- rownames(Wm) <- rownames(epsM) <- rownames(epsm) <- colnames(out)

  to.rm <- apply(WM, 1, anyNA)
  WM <- WM[-which(to.rm), ]
  Wm <- Wm[-which(to.rm), ]
  epsM <- epsM[-which(to.rm), ]
  epsm <- epsm[-which(to.rm), ]

  return(list("WM" = WM, "Wm" = Wm, "epsM" = epsM, "epsm" = epsm))
}
