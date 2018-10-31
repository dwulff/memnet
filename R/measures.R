
#' Retrieve largest component
#'
#' Retrieves the largest component. In case of equally sized
#' components the first is component is retrieved.
#'
#' @param adj numeric matrix representing the adjacency matrix.
#' @param weights numeric vector of edge weights. Optional.
#' @param mode character, either \code{"directed"} or \code{"undirected"},
#'   specifying whether the network should be interepeted as directed
#'   or undirected. Defaults to \code{"undirected"}.
#' @param igraph logical specifying whether the output should be of class
#'   \code{"igraph"}.
#'
#' @return
#' A list containing the, now, named adjacency matrix and a numeric value
#' indicating the size of the largest component relative to
#' to the entire graph.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency[1:10])
#'
#' # get adjacency matrix
#' adj = edg_to_adj(edge_list)
#'
#' # get largest component
#' l_comp(adj)
#'
#' @export
l_comp = function(adj, weights = NULL, mode = 'undirected', igraph = FALSE){

  # to igraph
  if(class(adj) != 'igraph'){
    g = igraph::graph_from_adjacency_matrix(adj, mode = mode)
    if(!is.null(weights)) igraph::E(g)$weight = weights
    } else {
    g = adj
    }

  # retrieve components
  cmps = igraph::components(g)

  # determine max component size
  msiz = max(cmps$csize)

  # Select nodes of max component
  icmp = which.max(cmps$csize)
  ncmp = (1:length(cmps$membership))[cmps$membership == icmp]

  # component
  comp = igraph::induced_subgraph(g,ncmp)

  # return igraph
  if(igraph == TRUE){
    out = list("adj" = comp,
               "rel_size" = msiz / igraph::vcount(g))
    return(out)
    }

  # return largest component and size of largest component
  out = list("adj" = igraph::get.adjacency(comp, sparse = FALSE),
             "rel_size" = msiz / igraph::vcount(g))
  out
  }


#' Maximum difference between cumulative degree distribution
#'
#' Determines the maximum difference between the cumulative
#' degree distributions of two graphs
#'
#' @inheritParams l_comp
#' @param adj_1 numeric matrix representing the adjacency matrix of graph 1.
#' @param adj_2 numeric matrix representing the adjacency matrix of graph 2.
#' @param weights_1,weights_2 numeric vector of edge weights for network 1 and
#'   2, respectively. Optional.
#' @param mode_1,mode_2 character, either \code{"directed"} or \code{"undirected"},
#'   specifying whether network 1 and 2 should be interepeted as directed
#'   or undirected, respectively. Defaults to \code{"undirected"}.
#'
#' @return
#' Distance between cumulative degree distributions.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list_1 = threshold_graph(animal_fluency[1:100])
#' edge_list_2 = threshold_graph(animal_fluency[101:200])
#'
#' # get adjacency matrices
#' adj_1 = edg_to_adj(edge_list_1)
#' adj_2 = edg_to_adj(edge_list_2)
#'
#' # get max degree distance
#' k_dist(adj_1, adj_2)
#'
#' @export
k_dist = function(adj_1, adj_2,
                  weights_1 = NULL, weights_2 = NULL,
                  mode_1 = 'undirected',
                  mode_2 = 'undirected'){

  # to igraph
  g1 = igraph::graph_from_adjacency_matrix(adj_1, mode = mode_1)
  if(!is.null(weights_1)) igraph::E(g1)$weight = weights_1
  g2 = igraph::graph_from_adjacency_matrix(adj_2, mode = mode_2)
  if(!is.null(weights_2)) igraph::E(g2)$weight = weights_2

  # get degree distance
  da = igraph::degree.distribution(g1)
  dg = igraph::degree.distribution(g2)
  k_dist = trm(da, dg)

  # out
  names(k_dist) = 'k_dist'
  k_dist
  }


#' Average local clustering (alc) coefficient
#'
#' Computes the uncorrected or corrected average local clustering coffiecient.
#'
#' The uncorrected clustering coefficent is computed according to Watts &
#' Strogatz (1998). The corrected clustering coefficient normalizes the
#' uncorrected one by the \code{average degree / n nodes}, i.e., the expected
#' average local clustering for an Erdös-Renyi random graph.
#'
#' @inheritParams l_comp
#' @param adj numeric matrix representing the adjacency matrix.
#' @param types character. Either \code{"uncorrected"} or \code{"corrected"}, or
#'   a vector containing both.
#'
#' @return the corrected local clustering coefficient and/or
#' the uncorrected clustering coefficient.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list = threshold_graph(animal_fluency)
#'
#' # get adjacency matrices
#' adj = edg_to_adj(edge_list)
#'
#' # get local average clustering coefficient
#' alc(adj)
#'
#' # get corrected local average clustering coefficient
#' alc(adj, types = 'corrected')
#'
#' @export
alc = function(adj, types = 'uncorrected', weights = NULL, mode = 'undirected'){

  # to igraph
  if(class(adj) != 'igraph'){
    g = igraph::graph_from_adjacency_matrix(adj, mode = mode)
    if(!is.null(weights)) igraph::E(g)$weight = weights
    } else {
    g = adj
    }

  res = c()
  cc = igraph::transitivity(g, type = 'localaverage')
  for(type in types){
    if(type == 'uncorrected'){
      res = c(res, cc)
      names(res)[length(res)] = 'cc'
      } else if(type == 'corrected'){
      nv = igraph::vcount(g)
      ne = igraph::ecount(g)
      k  = ne / nv
      res = c(res, cc / (k / nv))
      names(res)[length(res)] = 'cc_c'
      }
    }
  res
  }


#' Average shortest path length (aspl)
#'
#' Computes the average shortest path length (aspl) using \code{igraph}'s automatic
#' method.
#'
#' Per default the uncorrected apsl is computed. Otherwise, the uncorrected aspl
#' is normalized by \code{log(n nodes)/log(average degree)}, i.e., the expected
#' average shortest path length for an Erdös-Renyi random graph.
#'
#' @inheritParams l_comp
#' @param adj numeric matrix representing the adjacency matrix.
#' @param types character. Either \code{"uncorrected"} or \code{"corrected"}, or
#'   a vector containing both.
#'
#' @return Local clustering coefficient.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list = threshold_graph(animal_fluency)
#'
#' # get adjacency matrices
#' adj = edg_to_adj(edge_list)
#'
#' # get average shortest path length
#' aspl(adj)
#'
#' # get corrected average shortest path length
#' aspl(adj, types = 'corrected')
#'
#' @export
aspl = function(adj, types = 'uncorrected', weights = NULL, mode = 'undirected'){

  # to igraph
  if(class(adj) != 'igraph'){
    g = igraph::graph_from_adjacency_matrix(adj, mode = mode)
    if(!is.null(weights)) igraph::E(g)$weight = weights
    } else {
    g = adj
    }

  # get aspl
  res = c()
  aspl = igraph::mean_distance(g, directed = mode == 'directed')
  for(type in types){
    if(type == 'uncorrected'){
      res = c(res, aspl)
      names(res)[length(res)] = 'aspl'
      } else if(type == 'corrected'){
      n_v = igraph::vcount(g)
      n_e = igraph::ecount(g)
      k   = n_e / n_v
      res = c(res, aspl / (log(n_v) / log(k)))
      names(res)[length(res)] = 'aspl_c'
      }
    }
  res
  }

#' Network statistics
#'
#' Function computes various network measures including the size of the number of
#' nodes (|V|), the number of edges (|E|), the average degree (k),
#' uncorrected and corrected average local clustering (C, Cc), uncorrected and
#' corrected average shortest path lengths (L, Lc), the small-world index (S)
#' from Humphries & Gurney (2008), assortivity (A), and the size of largest
#' component (p) relative to the entire network.
#'
#' @inheritParams l_comp
#' @param adj numeric matrix representing the adjacency matrix.
#' @param giant logical specifying whether the should be computed for the
#'   largest component.
#'
#' @return A named numeric vector containing the computed network statistics.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list = threshold_graph(animal_fluency)
#'
#' # get adjacency matrices
#' adj = edg_to_adj(edge_list)
#'
#' # get structural overview
#' network_stats(adj)
#'
#' @export
network_stats = function(adj, giant = FALSE, weights = NULL, mode = 'undirected'){

  # to igraph
  if(class(adj) != 'igraph'){
    g = igraph::graph_from_adjacency_matrix(adj, mode = mode)
    if(!is.null(weights)) igraph::E(g)$weight = weights
    } else {
    g = adj
    }

  # get largest component
  comp   = l_comp(g, igraph = TRUE)
  comp_p = comp[[2]]
  if(giant == TRUE){
    g = comp[[1]]
  }

  # compute measures
  n_v = igraph::vcount(g)
  n_e = igraph::ecount(g)
  k   = n_e / n_v
  cc   = igraph::transitivity(g, type = 'localaverage')
  cc_c = cc / (k / n_v)
  aspl   = igraph::mean_distance(g, directed = mode == 'directed')
  aspl_c = aspl / (log(n_v) / log(n_e))
  smallw = cc_c / aspl_c
  assor  = igraph::assortativity_degree(g, directed = mode == 'directed')
  res = c(n_v, n_e, k, cc, cc_c, aspl, aspl_c, smallw, assor, comp_p)
  names(res) = c('|V|', '|E|', 'k',
                 'C', 'Cc',
                 'L', 'Lc',
                 'S', 'A',
                 'p')
  res
  }



#' Get common subgraph
#'
#' Function identifies and returns the common subgraphs of two networks.
#'
#' @inheritParams l_comp
#' @inheritParams k_dist
#'
#' @return List containing the two subgraphs.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list_1 = threshold_graph(animal_fluency[1:100])
#' edge_list_2 = threshold_graph(animal_fluency[101:200])
#'
#' # get adjacency matrices
#' adj_1 = edg_to_adj(edge_list_1)
#' adj_2 = edg_to_adj(edge_list_2)
#'
#' # get common subgraph
#' common_subgraphs(adj_1, adj_2)
#'
#' @export
common_subgraphs = function(adj_1, adj_2,
                            weights_1 = NULL, weights_2 = NULL,
                            mode_1 = 'undirected',
                            mode_2 = 'undirected'){

  # to igraph 1
  if(class(adj_1) != 'igraph'){
    g1 = igraph::graph_from_adjacency_matrix(adj_1, mode = mode_1)
    if(!is.null(weights_1)) igraph::E(g1)$weight = weights_1
    } else {
    g1 = adj_1
    }

  # to igraph 2
  if(class(adj_2) != 'igraph'){
    g2 = igraph::graph_from_adjacency_matrix(adj_2, mode = mode_2)
    if(!is.null(weights_2)) igraph::E(g2)$weight = weights_2
    } else {
    g2 = adj_2
    }

  # assign names
  if(is.null(igraph::V(g1)$name)) igraph::V(g1)$name = 1:nrow(adj_1)
  if(is.null(igraph::V(g2)$name)) igraph::V(g2)$name = 1:nrow(adj_2)

  # get common nodes
  common_nodes = intersect(igraph::V(g1)$name,
                           igraph::V(g2)$name)

  # common subgraph
  g1_com = igraph::induced_subgraph(g1, common_nodes)
  g2_com = igraph::induced_subgraph(g2, common_nodes)

  # out
  out = list(igraph::get.adjacency(g1_com, sparse = FALSE),
             igraph::get.adjacency(g2_com, sparse = FALSE))
  out
  }


#' Common subgraph statistics
#'
#' Function identifies the common subgraphs of two networks and returns the
#' networks statistics using \link{network_stats}.
#'
#' @inheritParams l_comp
#' @inheritParams k_dist
#' @inheritParams network_stats
#'
#' @return List containing the two subgraphs.
#'
#' @examples
#'
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge lists of fluency graphs
#' edge_list_1 = threshold_graph(animal_fluency[1:100])
#' edge_list_2 = threshold_graph(animal_fluency[101:200])
#'
#' # get adjacency matrices
#' adj_1 = edg_to_adj(edge_list_1)
#' adj_2 = edg_to_adj(edge_list_2)
#'
#' # get structural overview of both networks
#' common_subgraph_stats(adj_1, adj_2)
#'
#' @export
common_subgraph_stats = function(adj_1, adj_2,
                                 giant = FALSE,
                                 weights_1 = NULL, weights_2 = NULL,
                                 mode_1 = 'undirected',
                                 mode_2 = 'undirected'){

  # to igraph 1
  if(class(adj_1) != 'igraph'){
    g1 = igraph::graph_from_adjacency_matrix(adj_1, mode = mode_1)
    if(!is.null(weights_1)) igraph::E(g1)$weight = weights_1
    } else {
    g1 = adj_1
    }

  # to igraph 2
  if(class(adj_2) != 'igraph'){
    g2 = igraph::graph_from_adjacency_matrix(adj_2, mode = mode_2)
    if(!is.null(weights_2)) igraph::E(g2)$weight = weights_2
    } else {
    g2 = adj_2
    }

  # assign names
  if(is.null(igraph::V(g1)$name)) igraph::V(g1)$name = 1:nrow(adj_1)
  if(is.null(igraph::V(g2)$name)) igraph::V(g2)$name = 1:nrow(adj_2)

  # get common nodes
  common_nodes = intersect(igraph::V(g1)$name,
                           igraph::V(g2)$name)

  # common subgraph
  g1_com = igraph::induced_subgraph(g1, common_nodes)
  g2_com = igraph::induced_subgraph(g2, common_nodes)

  # compute stats
  res_1 = network_stats(g1_com)
  res_2 = network_stats(g2_com)

  # out
  out = rbind(res_1, res_2, res_1 - res_2)
  rownames(out) = c('adj_1','adj_2','adj_1-adj_2')
  out
  }





