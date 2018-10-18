
require(Rcpp)
require(memnet)

adj_1 = grow_ws(n = 1000,k = 4, p = .2)
adj_2 = grow_ws(n = 1000,k = 4, p = .2)

k_dist(net1,net2)

network_stats(net1)

common_subgraphs(net1, net2)
common_subgraph_stats(net1, net2)

aspl(adj_1, c('uncorrected','corrected'))
alc(adj_1, c('uncorrected','corrected'))
