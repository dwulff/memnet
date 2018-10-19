#' Edge list to adjacency matrix
#'
#' Transforms an edge list as returned by \link{community_graph} into an
#' adjacency matrix or \code{adjlist}.
#'
#' @param edg character matrix with two columns containing the from and to
#'   nodes of the edges.
#' @param weight optional numeric vector specifiying the weights associated with
#'   each edge.
#' @param directed logical specifying whether edges should only be included
#'   for from-to or also fro to-from.
#' @param adjlist logical specifying whether to export the adjancy matrix as an
#'   \code{adjlist} as required by \link{fluency}.
#'
#' @return a numeric n x n adjacency matrix with n being the number of unique
#'   entries in \code{edg}.
#'
#' @export

edg_to_adj = function(edg, weight = NULL, directed = FALSE, adjlist = FALSE){

  # unique entries
  unis = unique(c(edg))

  # create matrix
  adj = matrix(0, ncol = length(unis), nrow = length(unis))
  rownames(adj) = unis
  colnames(adj) = unis

  # fill matrix
  if(!is.null(weight)){
    if(length(weight) == nrow(edg)) stop('dimensions of edg and weight do not match')
    adj[edg] = edg[,3]
    if(directed == FALSE){
      adj[edg[,2:1]] = edg[,3]
      }
    } else {
    adj[edg] = 1
    if(directed == FALSE){
      adj[edg[,2:1]] = 1
      }
    }

  # check if adjlist
  if(adjlist == TRUE) return(get_adjlist(adj))

  # out
  adj
  }

#' Edge list to adjlist
#'
#' Transforms an edge list as returned by \link{community_graph} into an
#' adjlist as required by, e.g., \link{fluency}.
#'
#' @param edg character matrix with two columns containing the from and to
#'   nodes of the edges.
#'
#' @return a numeric n x n adjacency matrix with n being the number of unique
#'   entries in \code{edg}.
#'
#' @export

edg_to_adjlist = function(edg){

  # out
  edg_to_adj(edg, adjlist = TRUE)
  }

#' Rename memnet objects
#'
#' Function replaces the index vectors created by \link{adjlist} or \link{fluency}
#' by their original names.
#'
#' @param dat lists of indices that correspond to positions in a name vector.
#' @param names character vector giving the names.
#'
#' @return lists of character vectors.
#'
#' @export
rename_memnet = function(dat, names){

  # test of list
  if(!is.list(dat)) stop('dat must be of type list.')

  # get inds
  inds = unlist(dat)
  if(!is.integer(inds)) stop('dat must contain integer vectors.')

  # test max index
  if(max(inds) > length(names)) stop('name vector is too short.')

  # rename vectors
  n = length(dat)
  for(i in 1:n){
    dat[[i]] = names[dat[[i]]]
  }

  # check for names
  list_names = names(dat)
  if(!is.null(names(dat))){
    names(dat) = as.integer(list_names)
    }

  # out
  dat
  }

