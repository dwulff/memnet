#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          HELPERS
//
////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
std::vector<int> unique_int(std::vector<int> v){
  std::sort(v.begin(), v.end());
  std::vector<int>::iterator last;
  last = std::unique(v.begin(), v.end());
  v.erase(last, v.end());
  return v;
  }

// [[Rcpp::export]]
int rint(int n){
  return rand() % n;
  }

// [[Rcpp::export]]
double runi(){
  double t,p;
  t = std::rand() % 100001;
  p = t/100000;
  return p;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          GRAPH FUNCTIONS
//
////////////////////////////////////////////////////////////////////////////////////////////////////

//' Get adjacency list
//'
//' Get list containing adjacent, i.e., neighboring, nodes for each node in the
//' graph. Nodes are returned as their row indices of the adjacency matrix.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//'
//' @return A list of vectors containing the indices of adjacent nodes.
//'
//' @export
// [[Rcpp::export]]
GenericVector get_adjlist(NumericMatrix adj){
  int i, j, n = adj.nrow();
  GenericVector adjlist(n);
  std::vector<int> neighbors;
  NumericVector col(n), row(n);
  for(i = 0; i < n; i++){
    col = adj(_,i);
    row = adj(i,_);
    for(j = 0; j < n; j++){
      if(col[j] == 1 || row[j] == 1){
        neighbors.push_back(j + 1);
        }
      }
    adjlist[i] = neighbors;
    neighbors.clear();
    }
  return adjlist;
  }


//' Get neighbors \code{k} or fewer steps away
//'
//' Function iterates over graph to identify for a given node all nodes that are
//' no more than \code{k} steps apart.
//'
//' \code{k < 1} will be set to \code{k = 1}.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//' @param start integer specifying the row of the start node in the adjacency matrix.
//' @param k integer specifying the maximum distance to the start node.
//'
//' @return A character vector containing nodes \code{k} or fewer steps away
//' from \code{start}.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix get_neighborhood(NumericMatrix adj, int start, int k){
  int ik = 1, n = adj.nrow();
  std::vector<int> todos;
  std::vector<int> new_todos;
  std::vector<int> neighbors;
  std::map<int, int> nodes;
  NumericVector col(n), row(n);

  // init nodes
  nodes[start] = 0;

  // init todos
  col = adj(_,start - 1);
  row = adj(start - 1,_);
  for(int j = 0; j < n; j++){
    if(col[j] == 1 || row[j] == 1){
      todos.push_back(j);
      if(!inmap_2int(j+1, nodes)) nodes[j + 1] =  ik;
      }
    }

  // iterate
  while(ik < k){
    ik++;
    int n_todo = todos.size();
    for(int i = 0; i < n_todo; ++i){
      col = adj(_,todos[i]); // consider for allowing directions
      row = adj(todos[i],_); // (i.e., including only one of these terms)
      for(int j = 0; j < n; j++){
        if(col[j] == 1 || row[j] == 1){
          if(!inmap_2int(j+1, nodes)){
            new_todos.push_back(j);
            nodes[j + 1] =  ik;
            }
          }
        }
      }
    todos = new_todos;
    }

  NumericMatrix res(nodes.size(),2);
  std::map<int, int>::const_iterator it;
  int i = 0;
  for(it = nodes.begin(); it != nodes.end(); ++it, ++i){
    res(i,0) = it->first;
    res(i,1) = it->second;
    }

  return res;
  }

//' Get vector of neighbors exactly k steps away
//'
//' Function iterates over graph to identify for a given node all nodes that are
//' exactly k steps apart.
//'
//' \code{k < 1} will be set to \code{k = 1}.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//' @param start integer specifying the row index of the start node in the
//'   adjacency matrix.
//' @param k integer specifying the exact distance to the start node.
//'
//' @return A character vector containing nodes \code{k} or fewer steps away
//' from \code{start}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> get_kneighbors(NumericMatrix adj, int start, int k){
  int ik = 1, n = adj.nrow();
  std::vector<int> todos;
  std::vector<int> new_todos;
  std::vector<int> neighbors;
  std::vector<int> nodes;
  std::vector<int> kneighbors;
  NumericVector col(n), row(n);

  // init nodes
  nodes.push_back(start);

  // init todos
  col = adj(_,start - 1);
  row = adj(start - 1,_);
  for(int j = 0; j < n; j++){
    if(col[j] == 1 || row[j] == 1){
      todos.push_back(j);
      nodes.push_back(j + 1);
      if(ik == k) kneighbors.push_back(j + 1);
      }
    }

  // iterate
  while(ik < k){
    ik++;
    int n_todo = todos.size();
    for(int i = 0; i < n_todo; ++i){
      col = adj(_,todos[i]); // consider for allowing directions
      row = adj(todos[i],_); // (i.e., including only one of these terms)
      for(int j = 0; j < n; j++){
        if(col[j] == 1 || row[j] == 1){
          if(!inset(j+1,nodes)) new_todos.push_back(j);
          if(ik == k && !inset(j+1,nodes)) kneighbors.push_back(j + 1);
          nodes.push_back(j + 1);
          }
        }
      }
    todos = new_todos;
    }
  return unique_int(kneighbors);
  }


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          DEGREE DISTRIBUTION
//
////////////////////////////////////////////////////////////////////////////////////////////////////


// get position in cumulative distribution
// [[Rcpp::export]]
double prbs(std::vector<double> x, double p){
  int i = 0;
  double rp = x[0], n = double(x.size());
  while(rp <= p){
    rp += double(x[i]);
    i++;
    }
  return double(i) / n;
  }

// determine the largest difference between 2 cumulative distributions
// [[Rcpp::export]]
double trm(std::vector<double> x, std::vector<double> y){
  int i,n,nx = x.size(),ny = y.size();
  double p, d, dmax = 0;
  std::vector<double> cumx = csum(x), cumy = csum(y);
  n = cumx.size();
  for(i = 1; i < n; i++){
    p = cumx[i];
    d = std::abs(prbs(cumx,p) - prbs(cumy,p));
    if(d > dmax) dmax = d;
    }
  n = cumy.size();
  for(i = 1; i < n; i++){
    p = cumy[i];
    d = std::abs(prbs(cumx,p) - prbs(cumy,p));
    if(d > dmax) dmax = d;
    }
  return dmax;
  }



