

# -------------------------------------------------------------------------------------------------------------------- #
#
#          PLOTS
#
# -------------------------------------------------------------------------------------------------------------------- #

#' Plot graph
#'
#' Custom graph plot using \code{igraph}'s layout functions.
#'
#'
#' @param g graph object of class \code{igraph}. Alternatively, a
#' matrix or data.frame with the first two columns specifying the
#' from- and to-vertices.
#'

graph_plot = function(g,
                      vnames = NULL,
                      layout = layout.fruchterman.reingold,
                      grd_res = 48, horiz_lim = 3, vert_lim = 3,
                      ver_col = 'steelblue', edg_col = 'grey50', txt_col = 'black', shadow = T,
                      pt.cex = 2.5, cex = 1, lwd = 2){

  # ------ get graph

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

  # ------ extract edges and vertex names

  edg = as_edgelist(g)
  if(is.null(vnames)){
    if('name' %in% names(vertex_attr(g))){
      vnames = vertex_attr(g)$name
    } else stop('could not find vnames')
  } else if(length(V(g)) != length(vnames)) stop('length of vnames must match number of vertices')
  vnames = vertex_attr(g)[[1]]


  # ------ create grid

  grd    = expand.grid('x'=seq(-.05,1.05,length=grd_res),
                       'y'=seq(-.05,1.05,length=grd_res))
  grd[,c('row','col')] = expand.grid('row'=1:grd_res,'col'=1:grd_res)
  grd$id   = 1:nrow(grd)
  grd$pos  = apply(expand.grid(1:grd_res,1:grd_res),1,paste0,collapse='_')
  grd$free = rep(T,nrow(grd))


  # ------ get layout

  # apply layout function
  lyt = do.call(layout.fruchterman.reingold,list(g))

  # name and norm layout
  rownames(lyt) = vnames
  lyt[,1] = (lyt[,1]-min(lyt[,1]))/max(lyt[,1]-min(lyt[,1]))
  lyt[,2] = (lyt[,2]-min(lyt[,2]))/max(lyt[,2]-min(lyt[,2]))


  # ------ clear vertices space

  # limits of vertex block
  vert_lim = vert_lim
  horiz_lim = horiz_lim

  # iterate over vertices
  for(i in 1:nrow(lyt)){

    # find closest grid point
    pnt = lyt[i,]
    dif = sqrt((grd$x - pnt[1])**2 + (grd$y - pnt[2])**2)
    sel = which(dif==min(dif))[1]

    # clear rectangle around layout point
    pos = unlist(grd[sel,c('row','col')])
    rct = expand.grid(pos[1]+c(-horiz_lim : -1, 0, 1 : horiz_lim),
                      pos[2]+c(-vert_lim  : -1, 0, 1 : vert_lim))
    rem = apply(rct,1,paste0,collapse='_')
    grd$free[grd$pos %in% rem] <- F
    }


  # ------ place labels

  # init
  vert_lim = 1
  txt.pos = lyt

  # iterate over vertices
  for(i in 1:nrow(lyt)){

    # determine horizontal space around label
    wid = (.022*nchar(vnames[i])**.92)/2
    horiz_lim  = ceiling(wid*grd_res)

    # find closes non-taken grid position
    pnt = lyt[i,]
    tst = subset(grd,free == T)
    dif = sqrt((tst$x - pnt[1])**2 + (tst$y - pnt[2])**2)
    sel = which(dif==min(dif))[1]

    # clear space around label
    pos = unlist(tst[sel,c('row','col')])
    rct = expand.grid(pos[1]+c(-horiz_lim : -1, 0, 1 : horiz_lim),
                      pos[2]+c(-vert_lim  : -1, 0, 1 : vert_lim))
    rem = apply(rct,1,paste0,collapse='_')
    grd$free[ grd$pos %in% rem] = F

    # set text position
    id  = tst$id[sel]
    txt.pos[i,] = c(grd[id,1],grd[id,2])
    }

  # ------ Plot

  # set canvas
  plot.new()
  par(mar=c(.5,.5,1.5,.5))
  plot.window(c(-.06,1.06),c(-.06,1.06))

  # draw label lines
  apply(cbind(lyt,txt.pos),1,
        function(x) lines(c(x[1],x[3]),c(x[2],x[4]),lty=3,lwd=(lwd-.5)**.5,col='grey75'))
  points(txt.pos,pch=16,lwd=1,col='white',cex=1.5)

  # draw edges
  apply(edg,1,function(x){
    from=lyt[x[1],]
    to=lyt[x[2],];
    lines(c(from[1],to[1]),c(from[2],to[2]),
          lwd=lwd,col=edg_col)})

  # draw points
  if(length(ver_col) == 1) ver_col = rep(ver_col,nrow(lyt))
  sapply(1:nrow(lyt),function(i){
    if(shadow == TRUE) points(lyt[i,1]+.004, lyt[i,2]-.004,
                              pch=16,col=edg_col[1],cex=pt.cex)
    points(lyt[i,1],lyt[i,2],pch=16,lwd=1.3,col=ver_col[i],cex=pt.cex)
    NULL
    })

  # draw label background
  for(i in 1:nrow(txt.pos)){
    n = nchar(vnames[i])
    rect(txt.pos[i,1]-.007*n**.92,
         txt.pos[i,2]-.015,
         txt.pos[i,1]+.007*n**.92,
         txt.pos[i,2]+.015,
         col=rgb(1,1,1,alpha=.3),border=NA)
    }

  # draw
  text(txt.pos[,1],txt.pos[,2],vnames,font=1,cex=cex,col=txt_col)
  }


#' Neighborhood plot
#'
#' Plot k-Neighborhood of given vertex containing all vertices with distances
#' ≤ k from given vertex.
#'
#' @param g graph object of class \code{igraph}. Alternatively, a
#' matrix or data.frame with the first two columns specifying the
#' from- and to-vertices.
#' @param v integer specifying the row index (within the adjacency
#' matrix) of the vertex whose the neighborhood should be plotted.
#' Alternatively the vertex name, provided g allows for the extraction
#' of a name attribute.
#' @param k integer specifying the size of the neighborhood. Specifically,
#' the plot will contain all vertices that are \code{k} or fewer steps away
#' from \code{v}.
#' @param \dots arguments to be passed to \link{graph_plot}.
#'

neighborhood_plot = function(g, v, k, col = 'steelblue', ...){

  # ------ get graph

  if(class(g) != 'igraph'){
    if(is.matrix(g)){
      g = graph_from_edgelist(g[,1:2])
      } else if(is.data.frame(g)) {
      g = graph_from_edgelist(as.matrix(g[,1:2]))
      } else {
      stop('Non matching format for g')
      }
    }

  # ------ get neighborhood

  # get adjacency matrix
  adj = as.matrix(as_adj(g))

  # get index of v
  if(is.character(v)){
    if(!v %in% rownames(adj)){
      stop('vertex name not found')
      } else {
      v = which(v == rownames(adj))
      }
    } else if(is.numeric(v)){
    if(!v %in% 1:nrow(adj)) stop('vertex index is too large')
    }

  # get neighborhood
  hood = get_neighborhood(adj, v, k)
  cols = colorRampPalette(c('black',col,'white'))(max(hood[,2])+3)
  col  = cols[hood[,2] + 2]

  # ------ reduce graph

  adj = adj[hood[,1], hood[,1]]
  g   = graph_from_adjacency_matrix(adj)

  # ------ graph plot
  graph_plot(g, ver_col = col, txt_col = col, ...)
  }

