
# -------------------------------------------------------------------------------------------------------------------- #
#
#          HELPERS
#
# -------------------------------------------------------------------------------------------------------------------- #

# Circles --------

inPi = function(x) (x/360)*(pi*2)

x_circ  = function(n, deg, rad, orig, start){
  degs = seq(0,deg,length = n + 1)
  cbind(degs, cos(inPi(degs))*rad + orig[1])}

y_circP = function(x,rad,orig) orig[2] + sqrt(abs(rad^2-(x-orig[1])^2))
y_circM = function(x,rad,orig) orig[2] - sqrt(abs(rad^2-(x-orig[1])^2))

circle_raw<-function(n,deg,rad,orig){
  xs  = x_circ(n,-deg,rad,orig)
  ys1 = y_circP(xs[,2],rad,orig)
  ys2 = y_circM(xs[,2],rad,orig)
  sel = xs[,1]>180 & xs[,1]<360 | xs[,1]>-180 & xs[,1]<0 | xs[,1]>640 & xs[,1]<720
  ys = ys1
  ys[sel] = ys2[sel]
  cbind(xs[,2],ys)
}


circle = function(rad, orig, ..., n = 100){

  graphics::polygon(circle_raw(n,360,rad,orig),...)

}

# Other stuff --------

#' Fast general purpose color mixer
#'
#' Mixes two colors or matching vectors of colors according to some relative weight
#' and exports the result either in rgb or hex format.
#'
#' @param col_1,col_2 character vector of length one or of matching length containing
#'   colors either as a color name (see \link{colors}), rgb format (see \link{rgb}), or
#'   hex format.
#' @param weight numeric between 0 and 1 specifying the relative mixing weight for color
#'   one. E.g., \code{weight = .8} means that final color is composed of 80 percent
#'   \code{col_2} and 20 percent \code{col_1}.
#' @param format character string specifying the output format. Either \code{"hex"} or
#'   \code{"rgb"}.
#'
#' @return A vector of length \code{max(length(col_1), length(col_2))} containing the
#' mixed colors in the specified format.
#'
#' @examples
#'
#' # mix blue and red with more weight on blue
#' cmix('blue', 'red', .2)
#'
#' # mix blue and red with more weight on red
#' cmix('blue', 'red', .8)
#'
#' # mix blue and red and return as rgb
#' cmix('blue', 'red', .8, format = 'rgb')
#'
#' @export
cmix = function(col_1, col_2, weight, format = 'hex'){

  # get lens
  lens = c(ifelse(is.matrix(col_1),nrow(col_1),length(col_1)),
           ifelse(is.matrix(col_2),nrow(col_2),length(col_2)),
           length(weight))

  # test if lens fit
  if(max(lens[1:2]) %% min(lens[1:2]) != 0 | max(lens[1:2]) %% min(lens[3])){
    stop("Lengths of col_1, col_2, and weight are not equal or multiples of each other")
  }

  # test format
  if(!format %in% c('rgb','hex')) stop('Format can only be "hex" or "rgb"')

  # get target length
  target_length = max(lens)


  # col_1 to target length matrix
  if(is.character(col_1)){
    if(length(col_1) == 1){
      col_1 = grDevices::col2rgb(col_1)
      col_1 = matrix(c(col_1), ncol = 3, nrow = target_length)
      } else {
      n_rep = target_length / length(col_1)
      col_1 = t(sapply(rep(col_1, n_rep), grDevices::col2rgb))
      }
    } else{
    col_1 = matrix(c(col_1), ncol=3, nrow = target_length, byrow = FALSE)
    }

  # col_2 to target length matrix
  if(is.character(col_2)){
    if(length(col_2) == 1){
      col_2 = grDevices::col2rgb(col_2)
      col_2 = matrix(c(col_2), ncol = 3, nrow = target_length)
      } else {
      n_rep = target_length / length(col_2)
      col_2 = t(sapply(rep(col_2, n_rep), grDevices::col2rgb))
      }
    } else{
    col_2 = matrix(c(col_2), ncol=3, nrow = target_length, byrow = FALSE)
    }

  # expand weight
  weight = rep(weight, target_length / length(weight))


  # average colors
  col = col_1 * (1-weight) + col_2 * weight

  # out
  if(format == 'rgb') return(col / 255)
  grDevices::rgb(data.frame(col), maxColorValue = 255)

  }


# round even
round_even = function(x){
  rx = round(x)
  test = rx %% 2
  ifelse(test == 0,rx,
         ifelse( x < 0,
                 ifelse(x <  rx,floor(x),ceiling(x)),
                 ifelse(x >= rx,ceiling(x),floor(x))))
}


# find complementary color
compCol = function(col){

  # ----- extract color

  if(is.character(col)){
    col = grDevices::col2rgb(col)
  } else if(is.numeric(col)){
    if(length(col) == 3){
      col = matrix(col,nrow=3)
    } else stop('if numeric, col must be length 3')
  } else stop('col must be character or numeric')

  # ----- transform to complementary color

  compCol    = grDevices::rgb2hsv(col)
  compCol[1] = ifelse(compCol[1] > .5, compCol[1] - .5, compCol[1] + .5)
  compCol    = grDevices::hsv(compCol[1],compCol[2],compCol[3])

  # return
  return(compCol)
}

# get saturation
get_saturation = function(col){

  # ----- extract color

  if(is.character(col)){
    col = grDevices::col2rgb(col)
    } else if(is.numeric(col)){
    if(length(col) == 3){
      col = matrix(col, nrow=3)
      } else stop('if numeric, col must be length 3')
    } else stop('col must be character or numeric')

  # ----- transform to complementary color

  col    = grDevices::rgb2hsv(col)
  if(is.matrix(col)){
    return(col[2, ])
    } else {
    return(col[2])
    }
  }


# find complementary color
set_saturation = function(col, saturation = .5){

  # ----- extract color

  if(is.character(col)){
    col = grDevices::col2rgb(col)
  } else if(is.numeric(col)){
    if(length(col) == 3){
      col = matrix(col, nrow=3)
    } else stop('if numeric, col must be length 3')
  } else stop('col must be character or numeric')

  # check consistency
  if(!(ncol(col) %% length(saturation) == 0 |
    length(saturation) %% ncol(col) == 0)) stop('Non matching lengths.')
  if(any(saturation > 1)) stop('Saturation must be within 0 and 1.')

  # adjust length
  if(ncol(col) != length(saturation)){
    if(ncol(col) > length(saturation)) {
      saturation = rep(saturation, length(col) / length(saturation))
      }
    if(ncol(col) < length(saturation)) {
      col = matrix(c(col), nrow = 3, ncol = length(saturation) / ncol(col))
      }
    }

  # ----- transform to complementary color

  col    = grDevices::rgb2hsv(col)
  if(is.matrix(col)){
    col[2, ] = saturation
    col = apply(col, 2, function(x) grDevices::hsv(x[1], x[2], x[3]))
    } else {
    col[2] = saturation
    col    = grDevices::hsv(col[1],col[2],col[3])
    }

  # return
  col
  }




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
#' @param adj numeric matrix representing the adjacency matrix. Can also be an
#'   object of class \code{"igraph"} or an edge list, i.e., a two-column
#'   \code{matrix} or \code{data.frame} containing specifying the edges start
#'   and end points.
#' @param names optional character vector specifying the node names. Must be of
#'   appropriate length.
#' @param layout layout function from the \code{igraph} package. Default is
#'   \link[igraph]{layout.fruchterman.reingold}.
#' @param nod_col character vector of length 1 or length |V| specifying the
#'   node colors.
#' @param nod_cex numeric speciying the size of the node circles.
#' @param nod_shadow logical specifiying whether nodes should have shadows. Node
#'   shodow color is created from darkening \code{node_col}.
#' @param edg_col character vector of length 1 of length |E| specifying the edge
#'   line colors.
#' @param edg_lwd numeric vector of length 1 of length |E| specifying the edge
#'   line widths.
#' @param lab_col character vector of length 1 of length |V| specifying the
#'   text label colors.
#' @param lab_cex numeric vector of length 1 of length |V| specifying the
#'   text label sizes.
#' @param lab_lwd numeric vector of length 1 of length |V| specifying the
#'   width of the lines connecting the text label to the nodes.
#' @param lab_lcol character vector of length 1 of length |E| specifying
#'   specifying the color of the lines connecting the text label to the nodes.
#' @param lab_grid_size integer specifying the grid size used to place the node
#'   labels. Canvas is split in \code{lab_grid_size} and labels are placed into
#'   cells closest to the associated node.
#' @param lab_padding numeric vector of length 2 specifying the spacing among
#'   labels and between labels and nodes on the x and y dimension.
#'
#' @return nothing. A plot is created in \link{dev.cur}.
#'
#' @examples
#'
#' \dontrun{
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency[1:20])
#'
#' # get adjacency matrix
#' adj = edg_to_adj(edge_list)
#'
#' # plot
#' network_plot(adj)
#' }
#'
#' @export
network_plot = function(adj,
                        names = NULL,
                        layout = igraph::layout.fruchterman.reingold,
                        nod_col = "#E0EED4",
                        nod_cex = 3,
                        nod_shadow = T,
                        edg_col = 'grey25',
                        edg_lwd = 1.5,
                        lab_col = 'black',
                        lab_cex = 1,
                        lab_lwd = 1,
                        lab_lcol = 'grey25',
                        lab_grid_size = 48,
                        lab_padding = c(3, 3)){



  # ------ get graph

  if(class(adj) != 'igraph'){
    if(is.matrix(adj)){
      if(ncol(adj) == 2) {
        g = igraph::graph_from_edgelist(adj)
        } else if(ncol(adj) == nrow(adj)){
        g = igraph::graph_from_adjacency_matrix(adj)
        } else stop('adj has inappropriate format')
      } else if(is.data.frame(g)) {
        if(ncol(adj) == 2) {
        g = igraph::graph_from_edgelist(as.matrix(adj))
        } else stop('adj has inappropriate format')
      } else {
      stop('adj has inappropriate format')
      }
    } else {
    g = adj
    }

  # ------ extract edges and node names

  edg = igraph::as_edgelist(g)
  if(is.null(names)){
    if('name' %in% names(igraph::vertex_attr(g))){
      names = igraph::vertex_attr(g)$name
    } else {
      names = 1:igraph::vcount(g)
      message('labeled nodes 1:network_size')
    }
  } else if(igraph::vcount(g) != length(names)) stop('length of names must match number of nodes')



  # ------ create grid

  grd    = expand.grid('x'=seq(-.05,1.05,length=lab_grid_size),
                       'y'=seq(-.05,1.05,length=lab_grid_size))
  grd[,c('row','col')] = expand.grid('row'=1:lab_grid_size,'col'=1:lab_grid_size)
  grd$id   = 1:nrow(grd)
  grd$pos  = apply(expand.grid(1:lab_grid_size,1:lab_grid_size),1,paste0,collapse='_')
  grd$free = rep(T,nrow(grd))


  # ------ get layout

  # apply layout function
  lyt = do.call(layout, list(g))

  # name and norm layout
  rownames(lyt) = names
  lyt[,1] = (lyt[,1]-min(lyt[,1]))/max(lyt[,1]-min(lyt[,1]))
  lyt[,2] = (lyt[,2]-min(lyt[,2]))/max(lyt[,2]-min(lyt[,2]))


  # ------ clear nodes space

  # limits of node block
  vert_lim = lab_padding[1]
  horiz_lim = lab_padding[2]

  # iterate over nodes
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

  # iterate over nodes
  for(i in 1:nrow(lyt)){

    # determine horizontal space around label
    wid = (.022*nchar(names[i])**.92)/2
    horiz_lim  = ceiling(wid*lab_grid_size)

    # find closest non-taken grid position
    pnt = lyt[i,]
    tst = subset(grd, grd$free == T)
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
  graphics::plot.new()
  graphics::par(mar=c(.5,.5,1.5,.5))
  graphics::plot.window(c(-.06,1.06),c(-.06,1.06))

  # draw label lines
  if(length(lab_lwd) == 1) lab_lwd = rep(lab_lwd[1], nrow(lyt))
  if(length(lab_lcol) == 1) lab_lcol = rep(lab_lcol[1], nrow(lyt))
  lab_lcol
  for(i in 1:nrow(lyt)){
    graphics::lines(c(lyt[i, 1], txt.pos[i, 1]),
                    c(lyt[i, 2], txt.pos[i, 2]),
                    lty = 3,lwd = (lab_lwd[i] - .5) ** .5, col = lab_lcol[i])
    }
  graphics::points(txt.pos,pch=16,lwd=lab_lwd,col='white',cex = lab_cex)

  # draw edges
  if(length(edg_col) == 1) edg_col = rep(edg_col[1], nrow(edg))
  if(length(edg_lwd) == 1) edg_lwd = rep(edg_lwd[1], nrow(edg))
  for(i in 1:nrow(edg)){
    from=lyt[edg[i, 1],]
    to=lyt[edg[i, 2],];
    graphics::lines(c(from[1], to[1]), c(from[2], to[2]),
          lwd = edg_lwd[i], col = edg_col[i])
    }


  # draw points
  if(length(nod_cex) == 1) nod_cex = rep(nod_cex[1], nrow(lyt))
  if(length(nod_col) == 1) nod_col = rep(nod_col[1], nrow(lyt))
  for(i in 1:nrow(lyt)){
    if(nod_shadow == TRUE) {
      graphics::points(lyt[i,1]-.004, lyt[i,2]-.004,
      pch = 16, col = cmix(nod_col[i], 'black', .8), cex = nod_cex[i])
      }
    graphics::points(lyt[i, 1], lyt[i, 2], pch = 16,
           col = nod_col[i], cex = nod_cex[i])
  }


  # draw label background
  for(i in 1:nrow(txt.pos)){
    n = nchar(names[i])
    graphics::rect(txt.pos[i,1]-.007*n**.92,
                   txt.pos[i,2]-.015,
                   txt.pos[i,1]+.007*n**.92,
                   txt.pos[i,2]+.015,
                   col=grDevices::rgb(1,1,1,alpha=.3),border=NA)
    }

  # draw labels
  if(length(lab_cex) == 1) lab_cex = rep(lab_cex[1], nrow(lyt))
  if(length(lab_col) == 1) lab_col = rep(lab_col[1], nrow(lyt))
  graphics::text(txt.pos[, 1], txt.pos[, 2], labels = names,
                 font = 1, cex = lab_cex, col = lab_col)
  }



#' Neighborhood plot
#'
#' Plot k-Neighborhood of given node containing all nodes with distances
#' larger than k from given node.
#'
#' @param adj numeric matrix representing the adjacency matrix. Can also be an
#'   object of class \code{"igraph"} or an edge list, i.e., a two-column
#'   \code{matrix} or \code{data.frame} containing specifying the edges start
#'   and end points.
#' @param names optional character vector specifiying the node names.
#' @param node integer specifying the row index (within the adjacency
#'   matrix) of the node whose the neighborhood should be plotted.
#'   Alternatively the node name.
#' @param k integer specifying the size of the neighborhood. Specifically,
#'   the plot will contain all nodes that are \code{k} or fewer steps away
#'   from \code{v}.
#' @param nod_col character vector of length 1 specifying the node colors.
#' @param nod_shading logical specifying whether the node colors should be shaded
#'   as a function of the distance to \code{node}.
#' @param \dots arguments to be passed to \link{network_plot}.
#'
#' @return nothing. A plot is created in \link{dev.cur}.
#'
#'@examples
#'
#' \dontrun{
#' # get fluency data
#' data(animal_fluency)
#'
#' # edge list of fluency graph
#' edge_list = threshold_graph(animal_fluency[1:40])
#'
#' # get adjacency matrix
#' adj = edg_to_adj(edge_list)
#'
#' # plot
#' neighborhood_plot(adj, node = 'dog', k = 2)
#' }
#'
#' @export
neighborhood_plot = function(adj, names = NULL, node, k = 2, nod_col = "#E0EED4", nod_shading = TRUE, ...){

  # names function
  get_graph_names = function(g, names){
    if(is.null(names)){
      if('name' %in% names(igraph::vertex_attr(g))){
        names = igraph::vertex_attr(g)$name
      } else {
        names = 1:igraph::vcount(g)
        message('labeled nodes 1:network_size')
      }
    } else if(igraph::vcount(g) != length(names)) stop('length of vnames must match number of nodes')

    # out
    names
    }

  # ------ get graph

  if(class(adj) != 'igraph'){
    if(is.matrix(adj)){
      if(ncol(adj) == 2) {
        g = igraph::graph_from_edgelist(adj)
        adj = igraph::get.adjacency(g, sparse = FALSE)
        names = get_graph_names(g, names)
        } else if(nrow(adj) != ncol(adj)){
        stop('adj has inappropriate format')
        } else {
        if(is.null(rownames(adj))){
          names = 1:nrow(adj)
          message('labeled nodes 1:network_size')
          } else {
          names = rownames(adj)
          }
        }
      } else if(is.data.frame(g)) {
      if(ncol(adj) == 2) {
        g = igraph::graph_from_edgelist(as.matrix(adj))
        adj = igraph::get.adjacency(g, sparse = FALSE)
        names = get_graph_names(g, names)
        } else stop('adj has inappropriate format')
      } else {
      stop('adj has inappropriate format')
      }
    } else {
    names = get_graph_names(g, names)
    adj = igraph::get.adjacency(g, sparse = FALSE)
    }

  # ------ get neighborhood

  # get index of node
  if(is.character(node)){
    if(!node %in% rownames(adj)){
      stop('node name not found')
      } else {
      node = which(node == rownames(adj))
      }
    } else if(is.numeric(node)){
    if(!node %in% 1:nrow(adj)) stop('node index is too large')
    }

  # get neighborhood
  hood = get_neighborhood(adj, node, k)
  if(nod_shading == TRUE){
    hood[,2] = abs(hood[,2] - max(hood[,2]))
    cols = set_saturation(nod_col, (hood[,2] + .3)/max((hood[,2] + 1) + .6))
    } else {
    cols = rep(nod_col, nrow(hood))
    }

  # ------ reduce graph

  colnames(adj) = names
  adj = adj[hood[,1], hood[,1]]

  # ------ graph plot
  network_plot(adj, nod_col = cols, ...)
  }


#'
#' #' Circular network plot
#' #'
#' #' Plot the network as a circular graph
#' #'
#' #' @param g graph object of class \code{igraph}. Alternatively, a
#' #' matrix or data.frame with the first two columns specifying the
#' #' from- and to-nodes.
#' #' @param v integer specifying the row index (within the adjacency
#' #' matrix) of the node whose the neighborhood should be plotted.
#' #' Alternatively the node name, provided g allows for the extraction
#' #' of a name attribute.
#' #' @param k integer specifying the size of the neighborhood. Specifically,
#' #' the plot will contain all nodes that are \code{k} or fewer steps away
#' #' from \code{v}.
#' #' @param \dots arguments to be passed to \link{graph_plot}.
#' #'
#'
#' # load data
#' sim = readRDS('2 Clean Data/Tablet_SimRatings.RDS')
#'
#' sim_cats = readr::read_delim('0 Material/sim_words_translation.txt',
#'                              delim=' ',
#'                              col_names=c('word','translation','category'))
#'
#' # sim$in_cat = f.utils::replace_cc(sim$left_word,
#' #                                  sim_cats$word,
#' #                                  sim_cats$category) ==
#' #   f.utils::replace_cc(sim$right_word,
#' #                       sim_cats$word,
#' #                       sim_cats$category)
#' # print(sim %>% group_by(group, in_cat, pair_label) %>%
#' #   summarize(mr = mean(norm_rating)) %>%
#' #   filter(!in_cat & mr > .6),n=100)
#' #
#'
#' # pair label
#' sim$pair_label = ifelse(sim$left_word>sim$right_word,
#'                         paste(sim$left_word,sim$right_word,sep='_'),
#'                         paste(sim$right_word,sim$left_word,sep='_'))
#'
#' # extract test
#' sim_test = subset(sim, part == 'test')
#'
#' # aggregate
#' sim_agg = sim_test %>%
#'   group_by(group, pair_label) %>%
#'   summarize(weight = mean(norm_rating)) %>%
#'   mutate('from' = sapply(stringr::str_split(pair_label,'_'),`[`,1),
#'          'to'   = sapply(stringr::str_split(pair_label,'_'),`[`,2)) %>%
#'   ungroup()
#'
#'
#'
#' #
#' nodes = unique(c(sim_agg$from, sim_agg$to))
#'
#' l = memnet:::circle_raw(length(nodes),360,1,c(0,0))[-(length(nodes)+1),]
#' l_t = memnet:::circle_raw(length(nodes),360,1.07,c(0,0))[-(length(nodes)+1),]
#' l_t2 = memnet:::circle_raw(length(nodes),360,1.14,c(0,0))[-(length(nodes)+1),]
#'
#' cats = read.table('~/Dropbox (2.0)/Work/Projects/Memory/-- AgingLexicon/0 Material/Tablet_categories.txt',sep='\t')
#' n_cats = memnet:::match_cc(nodes,cats[,1],cats[,2])
#'
#' # order nodes
#' nodes = nodes[order(n_cats)]
#' n_cats = n_cats[order(n_cats)]
#'
#' sim_agg$from_ind = f.utils::replace_cn(sim_agg$from, nodes, 1:length(nodes))
#' sim_agg$to_ind = f.utils::replace_cn(sim_agg$to, nodes, 1:length(nodes))
#'
#' # translate
#' trans = read.table('~/Dropbox (2.0)/Work/Projects/Memory/-- AgingLexicon/0 Material/sim_words_translation.txt')
#' node_trans = f.utils::replace_cc(nodes,trans[,1],trans[,2])
#'
#' # split and process
#' sim_old = sim_agg %>% filter(group == 'old') %>% select(from_ind, to_ind, weight)
#' sim_young = sim_agg %>% filter(group != 'old') %>% select(from_ind, to_ind, weight)
#'
#'
#' # multiply colors
#' node_col = rgb(224,238,212,maxColorValue = 255)
#' #if(length(node_col) != nrow(l)) node_col = rep(node_col[1], nrow(l))
#'
#'
#' cmix = function(x,y,r=.5) {
#'   rgb(data.frame(matrix(rowSums(cbind(r*col2rgb(x),(1-r)*col2rgb(y))),nrow=1)),
#'       maxColorValue=255)
#' }
#'
#' pal = c(cmix(node_col,'black',.7),cmix(node_col, 'white', .7),cmix(node_col,'black',.7),cmix(node_col,'white',.7))
#' c_cols = colorRampPalette(pal)(length(unique(n_cats)))
#' n_cols = memnet:::match_cc(n_cats,unique(n_cats),c_cols)
#' node_col = n_cols
#'
#' # older adults
#'
#' png('5 Figures/older_sim_network.tiff',width=1000,height=1000,bg='transparent')
#'
#' plot.new();par(mar=c(0,0,0,0));plot.window(xlim=range(l_t[,1])*1.6,ylim=range(l_t[,2])*1.6)
#' points(l_t, col = node_col, pch = 16, cex = 5)
#'
#' for(i in 1:nrow(sim_old)) {
#'   if(sim_old$weight[i] > .4){
#'     from = l[sim_old$from_ind[i],]
#'     to = l[sim_old$to_ind[i],]
#'     lines(c(from[1],to[1]),c(from[2],to[2]), lwd = (sim_old$weight[i] + .5)**5,col=cmix(rgb(224,238,212,maxColorValue = 255),'black',.25))
#'   }
#' }
#'
#' for(i in 1:length(nodes)) text(l_t2[i,1],l_t2[i,2],labels=node_trans[i],srt=seq(360,0,length.out = length(nodes))[i],adj=0, cex= 2)
#'
#' dev.off()
#'
#'
