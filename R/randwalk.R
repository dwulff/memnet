#' Random walk
#'
#' Random walk simulates the average number of steps between two vertices
#' for each pair of vertices in the graph.
#'
#' @param

rand_walk = function(g, agg = 'median', nrep = 100, nmax = 1000, pjump = 0, verbose = FALSE){

  # ----- extract graph
  if(class(g) != 'igraph'){
    if(is.matrix(g)){
      if(ncol(g) == 2) {
        g = graph_from_edgelist(g[,1:2])
        } else if(ncol(g) == nrow(g)){
        g = graph_from_adjacency_matrix(g)
        } else stop('g has inappropriate format')
      } else if(is.data.frame(g)) {
      g = graph_from_edgelist(as.matrix(g[,1:2]))
      } else {
      stop('g has inappropriate format')
      }
    }

  # ----- get adjacency matrix
  adj = as.matrix(as_adj(g))
  edglist = get_adjlist(adj)

  # ----- random walk
  steps = lapply(1:length(edglist),function(x){
    if(verbose == TRUE) cat('Processing vertex ',x,'\n')
    mrwalk(edglist,x,observe = (x+1):length(edglist),
           nmax = nmax, nrep = nrep, pjump = pjump)
    })
  steps = do.call(rbind, steps)

  # ----- aggregate
  df = as_data_frame(steps)
  if(agg == 'median') df = df %>% group_by(V1,V2) %>% summarise('RT' = median(V3))
  if(agg == 'mean')   df = df %>% group_by(V1,V2) %>% summarise('RT' = mean(V3))
  df = as.data.frame(df)

  # ----- rename

  if('name' %in% names(vertex_attr(g))){
    vnames = vertex_attr(g)$name
    }

  df[,1] = vnames[df[,1]]
  df[,2] = vnames[df[,2]]
  df$id  = ifelse(df[,1] < df[,2], paste(df[,1],df[,2],sep='_'),paste(df[,2],df[,1],sep='_'))

  # ----- return
  return(df)
  }
