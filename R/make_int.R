#' Create a Canopy ready Y sna-cna interaction matrix
#' @param sna_obj An output object from the make_sna_mat() function.
#' @param cna_obj An output object from the make_cna_mat() function.
#' @return Y
#' @export
#' @examples
#' @import
#' dplyr
#' dtplyr


make_Y_mat <- function(sna_obj, cna_obj){
  X <- sna_obj$X
  WM <- cna_obj$WM

  ### Make overlap matrix Y ###
  Y <- as.data.frame(matrix(0L,nrow = nrow(X), ncol = nrow(WM)+1))
  rownames(Y) <- rownames(X)
  colnames(Y) <- c("Non-CNA",rownames(WM))

  all.cna <- as.data.frame(apply(gsub("chr","",do.call('rbind',strsplit(colnames(Y)[-1],"\\.|-"))),2,as.numeric))
  colnames(all.cna) <- c("chr","start","end")
  for(i in 1:nrow(Y)){
    temp <- as.numeric(strsplit(rownames(Y)[i],"_")[[1]])
    temp2 <- all.cna %>%
      filter(chr == temp[1],
             start <= temp[2],
             end >= temp[2])
    if(nrow(temp2)==0)
      Y[i,1] <- 1
    else{
      Y[i,match(paste0("chr",temp2[1],".",temp2[2],"-",temp2[3]),colnames(Y))] <- 1
    }
  }

  return(list("Y"=Y))
}

