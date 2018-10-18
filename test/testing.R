
require(Rcpp)
require(memnet)

adj = grow_ws(n = 1000,k = 4, p = .5)
adj_2 = grow_ws(n = 1000,k = 4, p = .2)

network_plot(adj,nod_col = grey(seq(.2,.8, length=10)))

k_dist(net1,net2)

network_stats(net1)

common_subgraphs(net1, net2)
common_subgraph_stats(net1, net2)

aspl(adj_1, c('uncorrected','corrected'))
alc(adj_1, c('uncorrected','corrected'))


memnet::network_plot(adj)

neighborhood_plot(adj, 5, k = 2, lab_lcol = 'grey25')

