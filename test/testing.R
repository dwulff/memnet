
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

data("animal_fluency")


# get data
data("animal_fluency")
age = as.numeric(names(animal_fluency))

# infer networks for age > 70
net_comunity = community_graph(animal_fluency[age > 70])
net_threshold = threshold_graph(animal_fluency[age > 70])
net_rw = rw_graph(animal_fluency[age > 70])

# show stats
network_stats(edg_to_adj(net_comunity))
network_stats(edg_to_adj(net_threshold))
network_stats(edg_to_adj(net_rw))

# plot
network_plot(edg_to_adj(net_comunity), nod_cex = 2, lab_cex = 1)
network_plot(edg_to_adj(net_threshold), nod_cex = 2, lab_cex = 1)
network_plot(edg_to_adj(net_rw), nod_cex = 2, lab_cex = 1, lab_grid_size = 70)

# inspect words
neighborhood_plot(edg_to_adj(net_comunity), k = 3, node = 'ant', nod_cex = 3, lab_cex = 1)
adj = edg_to_adj(net_comunity)
node = 'dog'


n = grow_ws(100, 5)
min(unlist((get_adjlist(n))))
network_plot(n)
get_adjlist(n)

f = fluency(get_adjlist(n), rep(10,100))
min(unlist(f))


# get data
data("animal_fluency")
age = as.numeric(names(animal_fluency))

# infer networks for age > 70
net_comunity = community_graph(animal_fluency[age > 70])

adj = edg_to_adj(net_comunity)
adjlist = get_adjlist(adj)

min(unlist(adjlist))

f = fluency(adjlist, rep(10,100))

rename_memnet(f, rownames(adj))

min(unlist(f))

