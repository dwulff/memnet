
#g = graph_from_edgelist(cbind(c(1,1,2,3,2,3),c(1,2,2,1,3,4)),directed = F)
#as_edgelist(g)

#' Retrieve largest component
#'
#' Retrieves the largest component. In case of equally sized
#' components the first is component is retrieved.
#'
#' @param g graph object of class \code{igraph}.
#'
#' @return
#' A list containing the largest component as an \code{igraph}
#' object and the relative size of the largest component compared
#' to the entire graph.
#'
#' @export
l_comp = function(g){

  # retrieve components
  cmps = components(g)

  # determine max component size
  msiz = max(cmps$csize)

  # Select vertices of max component
  icmp = which.max(cmps$csize)
  ncmp = (1:length(cmps$membership))[cmps$membership == icmp]

  # return largest component and size of largest component
  out = list(induced_subgraph(g,ncmp), msiz / vcount(g))
  names(out) = c('largest component','proportion of graph')
  return(out)
  }


#' Maximum difference between cumulative degree distribution
#'
#' Determines the maximum difference between the cumulative
#' degree distributions of two graphs.
#'
#' @param g1 a graph object of class \code{igraph}.
#' @param g2 a graph object of class \code{igraph}.
#'
#' @return
#' Distance between cumulative degree distributions in terms
#' of the difference of cumulative proportions.
#'
#' @export
k_dist = function(g1, g2){
  da = degree.distribution(g1)
  dg = degree.distribution(g2)
  k_dist = trm(da,dg)
  names(k_dist) = 'k_dist'
  return(k_dist)
  }


#' Local clustering coefficient
#'
#' Computes the corrected average local clustering coffiecient
#' using \code{igraph}'s average \link[igraph]{transitivity}.
#'
#' Per default the corrected clustering coefficient ist computed
#' which is computed as the average transitivity divided
#' by the \code{average degree / n vertices}, i.e., the expected
#' average local transitivity for an Erdös-Renyi random graph.
#'
#' @param g a graph object of class \code{igraph}.
#'
#' @return the corrected local clustering coefficient and/or
#' the uncorrected clustering coefficient.
#'
#' @export
cc = function(g, type = 'corrected'){
  res = c()
  cc = transitivity(g,type='average')
  for(typ in type){
    if(typ == 'uncorrected'){
      res = c(res, cc)
      names(res)[length(res)] = 'cc'
      } else if(typ == 'corrected'){
      nv = vcount(g)
      ne = ecount(g)
      k  = ne / nv
      res = c(res, cc / (k / nv))
      names(res)[length(res)] = 'cc_c'
      }
    return(res)
    }
  }


#' Average shortest path length
#'
#' Computes the average shortest path length using \code{igraph}'s
#' \link[igraph]{average.path.length} function.
#'
#' Per default the corrected average shortest path length is computed
#' as the average shortest path divided by
#' \code{log(n vertices)/log(average degree)}, i.e., the expected
#' average shortest path length for an Erdös-Renyi random graph.
#'
#' @param g a graph object of class \code{igraph}.
#' @param directed a logical specifying whether the algorithm should
#' take edge weights into account \code{igraph}.
#'
#' @return Local clustering coefficient.
#'
#' @export
aspl = function(g, type = 'corrected',  directed = FALSE){
  res = c()
  aspl = average.path.length(g,directed = directed)
  for(typ in type){
    if(typ == 'uncorrected'){
      res = c(res, aspl)
      names(res)[length(res)] = 'aspl'
      } else if(typ == 'corrected'){
      n_v = vcount(g)
      n_e = ecount(g)
      k   = n_e / n_v
      res = c(res, aspl / (log(n_v) / log(k)))
      names(res)[length(res)] = 'aspl_c'
      }
    return(res)
    }
  }


#' Corrected average shortest path length
#'
#' Computes the corrected average shortest path length using
#' \code{igraph}'s \link[igraph]{average.path.length} function.
#' The correction is achieved by dividing average shortest path
#' length by the quotient of the average degree and the number of
#' vertices, i.e., the expected averafe transitivity for a
#' Erdös-Renyi random graph.

#'
#' @param g a graph object of class \code{igraph}.
#' @param directed a logical specifying whether the algorithm should
#' take edge weights into account \code{igraph}.
#'
#' @return Local clustering coefficient.
#'
#' @export
aspl_c = function(g, directed = FALSE){
  n_v = vcount(g)
  n_e = ecount(g)
  k   = n_e / n_v
  aspl = average.path.length(g,directed = directed)
  return(aspl / (log(n_v) / log(k)))
  }



#' Graph statistics
#'
#' Function computes various graph statistical measures
#'
#' @param g a graph object of class \code{igraph}.
#' @param gian logical specifying whether the largest component
#'   of \code{g} should be used instead of \code{g}.
#' @param directed logical specifying whether the graph should
#'   be interpreted as being directed. Applies for \code{aspl}
#'   and \code{assortivity}.
g_stat = function(g, giant = FALSE, directed = F){
  if(class(g) != 'igraph') stop('g must be of class igraph')
  comp   = l_comp(g)
  comp_p = comp[[2]]
  if(giant == TRUE){
    g = comp[[1]]
    }
  n_v = vcount(g)
  n_e = ecount(g)
  k   = n_e / n_v
  cc   = transitivity(g,type='average')
  cc_c = cc / (k / n_v)
  aspl   = average.path.length(g,directed = directed)
  aspl_c = aspl / (log(n_v) / log(n_e))
  smallw = cc_c / aspl_c

  assor  = assortativity_degree(g, directed = directed)
  res = c(n_v, n_e, k, cc, cc_c, aspl, aspl_c, smallw, assor, comp_p)
  names(res) = c('n_v', 'n_e', 'k', 'cc', 'cc_c', 'aspl', 'aspl_c', 'smallw', 'assor', 'comp_p')
  return(res)
  }





