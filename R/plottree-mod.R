#' canopy.plottree.mod
#'
#' Split confCI character
#' @param tree
#' @param pdf
#' @param pdf.name
#' @param txt
#' @param txt.name
#' @param save
#' @param rdata.name
#' @param cex_names
#' @return tree, mutation list
#' @export
#' @examples
#' \dontrun{
#' canopy.plottree.mod(tree, pdf = NULL, pdf.name = NULL, txt = NULL, txt.name = NULL,save = FALSE, rdata.name = NULL)
#' }
#' @import
#' Canopy

canopy_plottree_mod <- function (tree, pdf = NULL, pdf.name = NULL,
                                 txt = NULL, txt.name = NULL,save = FALSE,
                                 rdata.name = NULL, cex_names = 1)
{
  # reformat tree mat #
  cl = stats::hclust(stats::dist(t(tree$P)), method = "ward")
  tree.ordered = t(t(tree$P)[cl$order, ])
  test_tree <- tree
  test_tree$P <- tree.ordered
  tree <- test_tree
  #####################

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
  ape::nodelabels()
  ape::tiplabels()
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
  ape::edgelabels(paste("mut", 1:length(edge.label), sep = ""),
             edge.label, frame = "n", col = 2, cex = 1.2)
  ape::tiplabels("Normal", 1, adj = c(0.2, 1.5), frame = "n", cex = 1.2,
            col = 4)
  ape::tiplabels(paste("Clone", 1:(K - 2), sep = ""), 2:(K - 1),
            adj = c(0.5, 1.5), frame = "n", cex = 1.2, col = 4)
  ape::tiplabels(paste("Clone", (K - 1), sep = ""), K, adj = c(0.8,
                                                          1.5), frame = "n", cex = 1.2, col = 4)
  par(mar = c(1, 7, 0.5, 9.5))
  P = tree$P
  image(1:nrow(P), 1:ncol(P), axes = FALSE, ylab = "", xlab = "",
        P, breaks = 0:100/100, col = fields::tim.colors(100))
  axis(4, at = 1:ncol(P), colnames(P), cex.axis = cex_names, las = 0.1,
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
