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
//' Get list containing adjacent vertices for each vertex in the graph.
//'
//' Adjacent vertices are returned in terms of their row index of the adjacency matrix.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//'
//' @return A list of vectors containing the indices of adjacent vertices.
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
        neighbors.push_back(j);
        }
      }
    adjlist[i] = neighbors;
    neighbors.clear();
    }
  return adjlist;
  }


//' Get neighbors k or fewer steps away
//'
//' Function iterates over graph to identify for a given vertex all vertices that are
//' no more than k steps apart.
//'
//' k = 0 will return be interpreted as k = 1.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//' @param start integer specifying the row of the start vertex in the adjacency matrix.
//' @param k integer specifying the maximum distance to the start vertex.
//'
//' @return A character vector containing vertices \code{k} or fewer steps away
//' from \code{start}.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix get_neighborhood(NumericMatrix adj, int start, int k){
  int ik = 1, n = adj.nrow();
  std::vector<int> todos;
  std::vector<int> new_todos;
  std::vector<int> neighbors;
  std::map<int, int> vertices;
  NumericVector col(n), row(n);

  // init vertices
  vertices[start] = 0;

  // init todos
  col = adj(_,start - 1);
  row = adj(start - 1,_);
  for(int j = 0; j < n; j++){
    if(col[j] == 1 || row[j] == 1){
      todos.push_back(j);
      if(!inmap_2int(j+1, vertices)) vertices[j + 1] =  ik;
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
          if(!inmap_2int(j+1, vertices)){
            new_todos.push_back(j);
            vertices[j + 1] =  ik;
            }
          }
        }
      }
    todos = new_todos;
    }

  NumericMatrix res(vertices.size(),2);
  std::map<int, int>::const_iterator it;
  int i = 0;
  for(it = vertices.begin(); it != vertices.end(); ++it, ++i){
    res(i,0) = it->first;
    res(i,1) = it->second;
    }

  return res;
  }

//' Get vector of neighbors exactly k steps away
//'
//' Function iterates over graph to identify for a given vertex all vertices that are
//' exactly k steps apart.
//'
//' k = 0 will return be interpreted as k = 1.
//'
//' @param adj numeric matrix specifying the adjacency matrix.
//' @param start integer specifying the row of the start vertex in the adjacency matrix.
//' @param k integer specifying the exact distance to the start vertex.
//'
//' @return A character vector containing vertices \code{k} or fewer steps away
//' from \code{start}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> get_kneighbor(NumericMatrix adj, int start, int k){
  int ik = 1, n = adj.nrow();
  std::vector<int> todos;
  std::vector<int> new_todos;
  std::vector<int> neighbors;
  std::vector<int> vertices;
  std::vector<int> kneighbors;
  NumericVector col(n), row(n);

  // init vertices
  vertices.push_back(start);

  // init todos
  col = adj(_,start - 1);
  row = adj(start - 1,_);
  for(int j = 0; j < n; j++){
    if(col[j] == 1 || row[j] == 1){
      todos.push_back(j);
      vertices.push_back(j + 1);
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
          if(!inset(j+1,vertices)) new_todos.push_back(j);
          if(ik == k && !inset(j+1,vertices)) kneighbors.push_back(j + 1);
          vertices.push_back(j + 1);
          }
        }
      }
    todos = new_todos;
    }
  return unique_int(kneighbors);
  }

