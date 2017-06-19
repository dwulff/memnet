#' Generalized topological overlap
#'
#' Computes the generalized topological overlap between pairs of vertices
#'
gtom = function(g, v1, v2, k = 2){

  # ----- extract adjacency matrix

  adj = as.matrix(as_adj(g))

  # ----- extract v1 and v2

  # v1
  if(is.character(v1)){
    if(v1 %in% rownames(adj)){
      v1 = which(v1 == rownames(adj))
      } else {
      stop('could not find v1')
      }
  } else if(is.numeric(v1)){
    if(v1 > nrow(adj)) stop('could not find v1')
  } else stop('could not find v1')

  # v2
  if(is.character(v2)){
    if(v2 %in% rownames(adj)){
      v2 = which(v2 == rownames(adj))
    } else {
      stop('could not find v2')
    }
  } else if(is.numeric(v2)){
    if(v2 > nrow(adj)) stop('could not find v2')
  } else stop('could not find v2')


  # ----- compute neighborhoods

  v1_hood = get_neighborhood(adj, v1, k)[,1]
  v2_hood = get_neighborhood(adj, v2, k)[,1]


  # ----- compute overlap

  # are v1 and v2 neighbors
  a = as.numeric(v1 %in% v2_hood)

  # remove v1 and v2
  v1_hood = v1_hood[-(which(v1 == v1_hood))]
  v2_hood = v2_hood[-(which(v2 == v2_hood))]

  # computer overlap
  num = length(intersect(v1_hood, v2_hood)) + a
  den = min(length(v1_hood), length(v2_hood)) + (1-a)
  tom = num / den

  # return
  return(tom)
  }

#' Generalized topological overlap
#'
#' Computes the generalized topological overlap between pairs of vertices
#'
mgtom = function(g, x, k = 2){

  # ----- extract adjacency matrix
  adj = as.matrix(as_adj(g))

  # store original x
  o = x


  # ----- extract pairs

  # v1
  if(!(is.matrix(x) | is.data.frame(x))) stop('x must be either matrix or data.frame')
  if(ncol(x) == 2){

    if(is.character(x[,1])){
      if(mean(x[,1] %in% rownames(adj)) != 1){
        stop('not all x found')
        } else {
        x[,1] = match_cn(x[,1],rownames(adj),1:nrow(adj))
        }
      } else if(is.numeric(x[,1])){
      if(max(x[,1]) > nrow(adj)) stop('some x are too large')
      } else stop('elements of x must be integer or character')
    if(is.character(x[,2])){
      if(mean(x[,2] %in% rownames(adj)) != 1){
        stop('not all x found')
        } else {
        x[,2] = match_cn(x[,2],rownames(adj),1:nrow(adj))
        }
      } else if(is.numeric(x[,2])){
      if(max(x[,2]) > nrow(adj)) stop('some x are too large')
      } else stop('elements of x must be integer or character')

    } else stop('x must be two column matrix or data.frame')

  # make numeric
  x = as.matrix(x)
  mode(x) = 'numeric'

  # ----- get gtom

  res = numeric(nrow(x))
  for(i in 1:nrow(x)){
    res[i] = gtom(g,x[i,1],x[i,2],k=k)
    }

  # name
  names(res) = paste(o[,1],o[,2],sep='_')

  return(res)
  }
