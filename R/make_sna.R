#' Create a Canopy ready SNA matrices
#' @param mut.files Names of the files to be read in and processed. Default will be all the files found in the working directory.
#' @param path relative path to the folder containing the samples to be processed. Default is current directory.
#' @param sample.names Custom names to be given to the samples. Default is NULL.
#' @return X, R
#' @export
#' @examples
#' @import
#' dplyr
#' dtplyr


make_sna_mat <- function(mut.files = list.files(), path = ".", sample.names = NULL){

  for(i in 1:length(mut.files)){
    assign(paste0("s",i,"_sna"),read.delim(paste0(path,"/",mut.files[i])))
  }


  pos <- c()
  for(i in 1:length(mut.files)){
    pos <- c(pos,get(paste0("s",i,"_sna"))$pos)
  }
  pos <- unique(pos)
  R <- as.data.frame(matrix(ncol = length(mut.files), nrow = length(pos)))
  rownames(R) <- pos
  if(is.null(sample.names))
    colnames(R) <- mut.files
  else
    colnames(R) <- sample.names

  X <- as.data.frame(matrix(ncol = length(mut.files), nrow = length(pos)))
  rownames(X) <- pos
  colnames(X) <- colnames(R)

  for(i in 1:nrow(R)){
    j <- rownames(R)[i]
    for(f in 1:length(mut.files)){
      temp <- get(paste0("s",f,"_sna")) %>% filter(pos == j)
      if(nrow(temp)>0){
        R[i,f] <- as.numeric(temp$alt)
        X[i,f] <- as.numeric(temp$total)
      }
    }
  }
  X[is.na(X)] <- 0
  R[is.na(R)] <- 0

  pos <- c()
  for(i in 1:length(mut.files)){
    pos <- c(pos,paste0(get(paste0("s",i,"_sna"))$chr,"_",
             get(paste0("s",i,"_sna"))$pos))

  }
  pos <- unique(pos)

  rownames(X) <- pos
  rownames(R) <- pos

  return(list("X"=X, "R"=R))
}
