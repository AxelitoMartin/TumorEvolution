#' Plot CNA heatmap
#' @param
#' @return facets heatmap
#' @export
#' @examples
#' @import
#' facets
#' pctGCdata
#' dplyr
#' dtplyr
#' gnomeR

plot.cna <- function(dat, patients = NULL, ordered = NULL, outcome = NULL){

  if(is.null(patients))
    patients <- rownames(dat)
  if (!is.null(outcome))
    names(outcome) <- patients
  if (!is.null(ordered))
    names(ordered) <- patients

  reducedM <- dat
  patients <- patients[match(rownames(reducedM), patients)]
  if (!is.null(outcome))
    outcome <- outcome[match(rownames(reducedM), names(outcome))]
  if (!is.null(ordered) && !is.null(outcome))
    ordered <- order(outcome)
  rownames(reducedM) <- abbreviate(rownames(reducedM),
                                   minlength = 10)
  imagedata = reducedM
  imagedata[imagedata > 1.5] = 1.5
  imagedata[imagedata < -1.5] = -1.5
  if (is.null(ordered)) {
    cl = stats::hclust(stats::dist(imagedata), method = "ward")
    imagedata.ordered = imagedata[cl$order, ]
    imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
  }
  if (!is.null(ordered)) {
    imagedata.ordered = imagedata[ordered, ]
    imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
  }
  chr = strsplit(colnames(imagedata), "\\.")
  chr = unlist(lapply(1:length(chr), function(x) chr[[x]][1]))
  chr = gsub("chr", "", chr)
  chr = as.numeric(chr)
  len = length(chr)
  chrom.ends <- rep(NA, length(table(chr)))
  d = 1
  for (r in unique(chr)) {
    chrom.ends[d] <- max(which(chr == r))
    d = d + 1
  }
  chrom.starts <- c(1, chrom.ends[-length(table(chr))] +
                      1)
  chrom.mids <- (chrom.starts + chrom.ends)/2
  bw = colorpanel(2, low = "white", high = "cadetblue4")
  colorkey = list(space = "right", height = 0.3, tick.number = 5)
  n <- nrow(reducedM)
  if (!is.null(outcome))
    x.lab <- outcome
  if (is.null(outcome))
    x.lab <- rep(" ", n)
  if (is.null(ordered))
    x.lab <- as.character(x.lab[cl$order])
  if (!is.null(ordered))
    x.lab <- as.character(x.lab[ordered])
  if (any(grepl("tcn", colnames(reducedM))) && any(grepl("ploidy",
                                                         colnames(reducedM))))
    scales = list(x = list(at = 1:n, labels = ploidy[cl$order],
                           rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                  z = list(at = n:1, labels = purity[cl$order],
                           rot = 90))
  else scales = list(x = list(at = 1:n, labels = x.lab,
                              rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                     z = list(at = n:1, labels = rep(1, n), rot = 90))
  my.panel.levelplot.2 <- function(...) {
    panel.levelplot(...)
    panel.abline(h = len - chrom.starts[-1], col = "gray",
                 lwd = 1)
    panel.scales = list(x = list(at = 1:n), y = list(at = len -
                                                       chrom.mids), z = list())
  }
  my.panel = my.panel.levelplot.2
  p = levelplot(imagedata.ordered, panel = my.panel, scales = scales,
                aspect = "fill", col.regions = bluered(256), xlab = "",
                ylab = "", colorkey = colorkey)
  return(list(p = p, out.cn = as.data.frame(dat$out.cn)))
}
