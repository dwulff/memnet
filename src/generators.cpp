#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////
//
//    HELPERS
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
std::vector<int> shuffle(std::vector<int> & v){
  int n = v.size();

  for (int i = 0; i < n - 1; ++i){
    // generate a random number j such that i <= j < n and
    // swap the element present at index j with the element
    // present at current index i
    int j = i + Rcpp::sample(100000000,1)[0] % (  n - i);
    std::swap(v[i], v[j]);
  }
  return v;
}


////////////////////////////////////////////////////////////////////////////////
//
//    STEYVERS & TENENBAUM
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericMatrix seed(int n, int m){
  int i,j;
  NumericMatrix adj(n,n);
  for(i = 0; i < (m - 1); i++){
    for(j = (i + 1); j < m; j++){
      if(j != i){
        adj(i,j) = adj(j,i) = 1;
        }
      }
    }
  return adj;
  }



// [[Rcpp::export]]
int sm(NumericVector x){
  int i,s = 0,n = x.size();
  for(i = 0; i < n; i++){
    s += x[i];
    }
  return s;
  }

// [[Rcpp::export]]
std::vector<int> getdegrees(NumericMatrix adj, int pos){
  int i;
  std::vector<int> degrees;
  for(i = 0; i < pos; i++){
    degrees.push_back(sm(adj(_,i)));
    }
  return degrees;
  }


std::vector<int> getneighbors(NumericMatrix adj,int node, int pos){
  int i;
  NumericVector vec = adj(_,node);
  std::vector<int> neighbors;
  for(i = 0; i < pos; i++){
    if(vec[i] == 1){
      neighbors.push_back(i);
      }
    }
  return neighbors;
  }

// [[Rcpp::export]]
std::vector<int> getnonneighbors(NumericMatrix adj,int node){
  int i, n = adj.nrow();
  NumericVector vec = adj(_,node);
  std::vector<int> neighbors;
  for(i = 0; i < n; i++){
    if(vec[i] == 0){
      neighbors.push_back(i);
    }
  }
  return neighbors;
}

// // [[Rcpp::export]]
// int randint(int n){
//   return std::rand() % n;
//   }

// [[Rcpp::export]]
int randint(int n){
  return (Rcpp::sample(n, 1))[0] - 1;
}

// [[Rcpp::export]]
int selectnode(std::vector<int> ps){
  int k, sum = 0, i = 0, n = ps.size();
  double v, r = double((Rcpp::runif(1,0,1))[0]);
  for(k = 0; k < n; k++){
    sum += ps[k];
    }
  v = ps[0] / double(sum);
  if(sum != 0){
    while(i < n && v <= r){
      i++;
      v += double(ps[i]) / double(sum);
      }
    } else {
    i = randint(n);
    }
  return i;
  }

// [[Rcpp::export]]
int selectnode_power(std::vector<int> ps, double power){
  int j, sum = 0, i = 0, n = ps.size();
  double v, r = double((Rcpp::runif(1,0,1))[0]);
  for(j = 0; j < n; j++){
    ps[j] = std::pow(ps[j], power);
    sum += ps[j];
  }
  v = ps[0] / double(sum);
  if(sum != 0){
    while(i < n && v <= r){
      i++;
      v += double(ps[i]) / double(sum);
    }
  } else {
    i = randint(n);
  }
  return i;
}


//' Steyvers and Tenenbaum (2004) network growth model
//'
//' Grow networks using Steyvers and Tenenbaum (2004) model, which combines
//' preferential attachment with a triad formation.
//'
//' @param n Integer. Number of nodes in the network.
//' @param m Integer. Number of edges added for each incoming node.
//'
//' @return n x n adjacency matrix.
//'
//' @references
//' Steyvers, M., & Tenenbaum, J. B. (2005). The large‐scale structure of
//' semantic networks: Statistical analyses and a model of semantic growth.
//' Cognitive science, 29(1), 41-78.
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' @examples
//' # generate small graph
//' grow_st(n = 6, m = 2)
//'
//' \dontrun{
//' # generate large graph
//' grow_st(n = 100, m = 10)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericMatrix grow_st(int n = 100, int m = 5){
  if(n < m + 1) return 0;
  int i,j, node, connect;
  NumericMatrix adj = seed(n, m + 1);
  std::vector<int> degrees(m+1,m);
  std::vector<int> neighbors;
  for(i = m + 1; i < n; i++){
    node      = selectnode(degrees);
    neighbors = getneighbors(adj, node, i);
    neighbors = shuffle(neighbors);
    if(neighbors.size() < m) return 0;
    for(j = 0; j < m; j++){
      connect = neighbors[j];
      adj(i, connect) = adj(connect, i) = 1;
      degrees[connect]++;
      }
    degrees.push_back(m);
    }
  return adj;
  }


////////////////////////////////////////////////////////////////////////////////
//
//    HOLME & KIM (1999) SCALE FREE WITH TUNABLE CLUSTERING
//
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
NumericMatrix emptyseed(int n){
  NumericMatrix adj(n,n);
  return adj;
  }

// [[Rcpp::export]]
double puni(){
  return double((Rcpp::runif(1,0,1))[0]);
  }

// [[Rcpp::export]]
int unconnectedneighbor(NumericMatrix adj, int from, int to){
  int i, test;
  std::vector<int> neighbors = getneighbors(adj, to, from);
  neighbors = shuffle(neighbors);
  //std::cout << neighbors.size() << '\n';
  for(i = 0; i < neighbors.size(); i++){
    test = adj(neighbors[i],from);
    if(test == 0) return neighbors[i];
    }
  return -1;
  }


// // [[Rcpp::export]]
// void test(int n = 100, int m = 5){
//   std::vector<int> neighbors, degrees(m,0.0);
//   for(int i = 0; i < n; ++i) std::cout << selectnode(degrees) << '\n';
//   }

//' Holme and Kim (2002) network growth model
//'
//' Grow networks using Holme & Kim's (2002) model, which combines preferential
//' attachment with tunable triad formation flexibly controlling the amount of
//' clustering in the network via \code{p}.
//'
//' @param n Integer. Number of nodes in the network.
//' @param m Integer. Number of edges added for each incoming node.
//' @param p Numeric. Proability that a triad formation step follows a preferential
//' attachment step.
//'
//' @return n x n adjacency matrix.
//'
//' @references
//' Holme, P., & Kim, B. J. (2002). Growing scale-free networks with tunable
//' clustering. Physical review E, 65(2), 026107.
//'
//' @examples
//' # generate small graph
//' grow_hk(n = 6, m = 2, p = .1)
//'
//' \dontrun{
//' # generate large graph, low clustering
//' grow_hk(n = 100, m = 10, p = .1)
//'
//' # generate large graph, high clustering
//' grow_hk(n = 100, m = 10, p = .9)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericMatrix grow_hk(int n = 100, int m = 5, double p = 5){
  int i, f, node;
  NumericMatrix adj = seed(n,m);
  std::vector<int> neighbors, degrees(m,0.0);
  for(i = 0; i < m; i++) degrees[i] = m - 1;
  for(f = m; f < n; f++){
    int im = 0;
    while(im < m){
      std::vector<int> degrees_cur(degrees);
      node      = selectnode(degrees_cur);
      //std::cout << f << '\t' << node << '\t' << adj(f, node) << '\n';
      if(adj(node, f) != 1){
        im++;
        adj(f, node)  = adj(node, f) = 1;
        degrees[node]++;
        degrees_cur[node] = 0;
        }
      if(runi() < p && im < m){
        neighbors = getneighbors(adj, node, f);
        neighbors = shuffle(neighbors);
        int nneigh = neighbors.size();
        //for(int j = 0; j < nneigh; ++j) std::cout << neighbors[j] << ' ';
        //std::cout << '\n';
        for(int i = 0; i < nneigh; ++i){
          node = neighbors[i];
          if(adj(node, f) != 1){
            im++;
            adj(f, node)  = adj(node, f) = 1;
            degrees[node]++;
            degrees_cur[node] = 0;
            break;
            }
          }
        }
      }
    degrees.push_back(m);
    }
  return adj;
  }

////////////////////////////////////////////////////////////////////////////////
//
//    BARABÁSI & ALBERT
//
////////////////////////////////////////////////////////////////////////////////


//' Barabási & Albert (2002) network growth model
//'
//' Grow networks using Barabási & Alberts's (2002) preferential
//' attachment model.
//'
//' @param n Integer. Number of nodes in the network.
//' @param m Integer. Number of edges added for each incoming node.
//' @param power Numeric. Controls the selection of nodes by raising the degree
//' to this power.
//'
//' @return n x n adjacency matrix.
//'
//' @references
//' Barabási, A. L., & Albert, R. (1999). Emergence of scaling in random
//' networks. Science, 286(5439), 509-512.
//'
//' @examples
//' # generate small graph
//' grow_ba(n = 6, m = 2)
//'
//' \dontrun{
//' # generate large graph, flat degree distribution
//' grow_ba(n = 100, m = 10, p = .1)
//'
//' # generate large graph, steep degree distribution
//' grow_ba(n = 100, m = 10, p = 10)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericMatrix grow_ba(int n = 100, int m = 5, double power = 1){
  int i, f, node;
  NumericMatrix adj = seed(n,m);
  std::vector<int> neighbors, degrees(m,0.0);
  for(i = 0; i < m; i++) degrees[i] = m - 1;
  for(f = m; f < n; f++){
    int im = 0;
    while(im < m){
      std::vector<int> degrees_cur(degrees);
      if(power == 1){
        node = selectnode(degrees_cur);
        } else {
        node = selectnode_power(degrees_cur, power);
        }
      //std::cout << f << '\t' << node << '\t' << adj(f, node) << '\n';
      if(adj(node, f) != 1){
        im++;
        adj(f, node)  = adj(node, f) = 1;
        degrees[node]++;
        degrees_cur[node] = 0;
      }
    }
    degrees.push_back(m);
  }
  return adj;
}



////////////////////////////////////////////////////////////////////////////////
//
//    Watts & Strogatz
//
////////////////////////////////////////////////////////////////////////////////

// container
struct container{
  int first, second, third;
};

// comparer
struct comparer{
  inline bool operator()(const container& one, const container& two){
    if(one.first == two.first){
      if(one.second == two.second){
        return one.third < two.third;
      }
      return one.second < two.second;
    }
    return one.first < two.first;
  }
};


void sort_3(IntegerMatrix& mat, IntegerVector by){

  // vectors to struct
  std::vector<container> vec(mat.nrow());
  for(int i = 0; i < vec.size(); ++i){
    vec[i].first = mat(i, by[0]);
    vec[i].second = mat(i, by[1]);
    vec[i].third = mat(i, by[2]);
  }

  // sort
  sort(vec.begin(),vec.end(),comparer());

  // back to matrix
  for(int i = 0; i < vec.size(); ++i){
    mat(i, by[0]) = vec[i].first;
    mat(i, by[1]) = vec[i].second;
    mat(i, by[2]) = vec[i].third;
  }

}

//' Watts & Strogatz (2002) network growth model
//'
//' Grow networks using Watts & Strogatz (1999) growth model, which constructs
//' in-between regular lattices and random networks by re-wiring edges of a
//' regular lattice with probability \code{p}.
//'
//' @param n Integer. Number of nodes in the network.
//' @param k Integer. Number of edges added for each incoming node. Can only be
//'   even.
//' @param p Numeric. Proability that an edge e_ij is rewired to e_ik with k being
//'   randomly drawn from the set of nodes.
//'
//' @return n x n adjacency matrix.
//'
//' @references
//' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of ‘small-world’
//' networks. Nature, 393(6684), 440-442.
//'
//' @examples
//' # generate small, mildly random graph
//' grow_ws(n = 6, k = 2, p = .2)
//'
//' \dontrun{
//' # generate large, mildly random graph
//' grow_ws(n = 100, k = 10, p = .1)
//'
//' # generate large, highly random graph
//' grow_ws(n = 100, k = 10, p = 10)
//' }
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix grow_ws(int n = 100, int k = 10, double p = .2){

  // extract edges
  int ind = 0, max_dist = k / 2;
  int n_edges = n * max_dist;
  IntegerMatrix edg(n_edges, 3);
  for(int i = 0; i < n - 1; ++i){
    for(int j = i + 1; j < n; ++j){
      int dist = std::min(std::abs(j - i),std::abs((j - n) - i));
      if(dist < max_dist + 1){
        edg(ind, 0) = i;
        edg(ind, 1) = j;
        edg(ind, 2) = dist;
        ind++;
      }
    }
  }

  // sort
  sort_3(edg, IntegerVector::create(2, 0, 1));

  // fill matrix
  IntegerMatrix adj(n, n);
  for(int i = 0; i < n_edges; ++i){
    adj(edg(i,0),edg(i,1)) = adj(edg(i,1), edg(i,0)) = 1;
  }

  // rewire
  for(int i = 0; i < n_edges; ++i){
    double p_cur = double((Rcpp::runif(1,0,1))[0]);
    if(p > p_cur){
      int new_j = randint(n);
      if(new_j != edg(i, 0) && adj(edg(i, 0), new_j) == 0){
        adj(edg(i,0), edg(i,1)) = 0;
        adj(edg(i,1), edg(i,0)) = 0;
        adj(edg(i,0), new_j) = 1;
        adj(new_j, edg(i,0)) = 1;
      }
    }
  }

  return adj;
}


////////////////////////////////////////////////////////////////////////////////
//
//    Regular lattice
//
////////////////////////////////////////////////////////////////////////////////

//' Regular lattice network model
//'
//' Grow regular lattice networks, in which every node is connected to m neighbors.
//'
//' @inheritParams grow_ws
//'
//' @return n x n adjacency matrix.
//'
//' @references
//' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of ‘small-world’
//' networks. Nature, 393(6684), 440-442.
//'
//' @examples
//' # generate small lattice
//' grow_lattice(n = 6, k = 2)
//'
//' \dontrun{
//' # generate large lattice
//' grow_lattice(n = 100, k = 10, p = .1)
//' }
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix grow_lattice(int n = 100, int k = 10){

  // extract edges
  int ind = 0, max_dist = k / 2;
  int n_edges = n * max_dist;
  IntegerMatrix edg(n_edges, 3);
  for(int i = 0; i < n - 1; ++i){
    for(int j = i + 1; j < n; ++j){
      int dist = std::min(std::abs(j - i),std::abs((j - n) - i));
      if(dist < max_dist + 1){
        edg(ind, 0) = i;
        edg(ind, 1) = j;
        edg(ind, 2) = dist;
        ind++;
      }
    }
  }

  // fill matrix
  IntegerMatrix adj(n, n);
  for(int i = 0; i < n_edges; ++i){
    adj(edg(i,0),edg(i,1)) = adj(edg(i,1), edg(i,0)) = 1;
  }


  return adj;
}


