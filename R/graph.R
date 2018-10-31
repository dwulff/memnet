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
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency[1:3])
#'
#' # transfrom to adjacency matrix
#' edg_to_adj(edge_list)
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
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency[1:3])
#'
#' # transfrom to adjacency matrix
#' edg_to_adjlist(edge_list)
#'
#' @export
edg_to_adjlist = function(edg){

  # out
  edg_to_adj(edg, adjlist = TRUE)
  }

#' Restore names of memnet objects
#'
#' Function replaces the index vectors created by \link{get_adjlist}
#' \link{fluency}, or \link{search_rw} by their original names.
#'
#' @param dat lists of indices that correspond to positions in a name vector or
#'   matrix with two initial columns containing node indices.
#' @param names character vector giving the names.
#'
#' @return lists of character vectors.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency)
#'
#' # extract adjlist from community network
#' adjlist = edg_to_adjlist(edge_list)
#'
#  simulate fluency sequences
#' f = fluency(adjlist, c(10, 14, 16, 18))
#' restore_names(f, get_names(edge_list))
#'
#' @export
restore_names = function(dat, names){

  # test of list
  if(!is.list(dat) & !is.matrix(dat)){
    stop('dat must be of type list or matrix.')
    }

  if(is.list(dat)){

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

  } else {

    # round first two columns
    ind_cols = dat[,1:2]
    mode(ind_cols) = 'integer'

    # get inds
    inds = c(ind_cols)
    if(!is.integer(inds)) stop('first two columns of dat must contain integer vectors.')

    # test max index
    if(max(inds) > length(names)) stop('name vector is too short.')

    # rename vectors
    dat[, 1] = names[ind_cols[, 1]]
    dat[, 2] = names[ind_cols[, 2]]

  }

  # out
  dat
  }

#' Get node names of memnet objects
#'
#' Function extracts names from objects created from edgelists and adjacency
#' matrices.
#'
#' @param dat numeric or character matrix containing an edgelist or adjacency
#'   matrix.
#'
#' @return lists of character vectors.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency)
#'
#' # get names
#' get_names(edge_list)
#'
#' @export
get_names = function(dat){

  # check if matrix
  if(!is.matrix(dat)) stop('dat must be matrix.')

  # check if adj
  if(nrow(dat) == ncol(dat)){
    names = ncol(dat)
    } else {
    if(mode(dat) == 'character'){
      names = get_names_c(dat)
      } else {
      names = get_names_i(dat)
      }
    }

  # out
  names
  }



