#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


// STEYVERS & TENENBAUM

// [[Rcpp::export]]
NumericMatrix seed(int n, int m){
  int i,j;
  NumericMatrix adj(n,n);
  for(i = 0; i < m; i++){
    for(j = 0; j < m; j++){
      if(j != i){
        adj(i,j) = 1;
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

int randint(int n){
  return rand() % n;
  }

// [[Rcpp::export]]
int selectnode(std::vector<int> ps){
  int k, sum = 0, i = 0, n = ps.size();
  double v, r = double(std::rand()) / RAND_MAX;
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
NumericMatrix stgame(int n, int m){
  if(n < m + 1) return 0;
  int i,j, node, connect;
  NumericMatrix adj = seed(n,m + 1);
  std::vector<int> degrees(m+1,m);
  std::vector<int> neighbors;
  for(i = m + 1; i < n; i++){
    node      = selectnode(degrees);
    neighbors = getneighbors(adj, node, i);
    std::random_shuffle(neighbors.begin(),neighbors.end());
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


// // HOLME & KIM (1999) SCALE FREE WITH TUNABLE CLUSTERING


// [[Rcpp::export]]
NumericMatrix emptyseed(int n){
  NumericMatrix adj(n,n);
  return adj;
  }

// [[Rcpp::export]]
double puni(){
  return double(std::rand()) / RAND_MAX;
  }

// [[Rcpp::export]]
int unconnectedneighbor(NumericMatrix adj, int from, int to){
  int i, test;
  std::vector<int> neighbors = getneighbors(adj, to, from);
  std::random_shuffle(neighbors.begin(),neighbors.end());
  //std::cout << neighbors.size() << '\n';
  for(i = 0; i < neighbors.size(); i++){
    test = adj(neighbors[i],from);
    if(test == 0) return neighbors[i];
    }
  return -1;
  }


// [[Rcpp::export]]
void test(int n = 100, int m = 5){
  std::vector<int> neighbors, degrees(m,0.0);
  for(int i = 0; i < n; ++i) std::cout << selectnode(degrees) << '\n';
  }

// [[Rcpp::export]]
NumericMatrix hkgame(int n, int m, double p){
  int i, f, node;
  NumericMatrix adj = seed(n,m);
  std::vector<int> neighbors, degrees(m,0.0);
  for(i = 0; i < m; i++) degrees[i] = m - 1;
  for(f = m; f < n; f++){
    int im = 0, t = 0;
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
        std::random_shuffle(neighbors.begin(),neighbors.end());
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




// // [[Rcpp::export]]
// NumericMatrix hkgame(int n, int m, double p, int pw = 1){
//   int i, f, j, pos, to, ncon, uncneighbor;
//   NumericMatrix adj = emptyseed(n);
//   std::vector<int> neighbors, connect, activen, degrees(n,0);
//   // for(i = 0; i < m; i++) degrees[i] = m - 1;
//   // for(i = 0; i < m; i++) degrees[i] = 0;
//   for(f = m; f < n; f++){
//     activen      = degrees;
//     connect      = getneighbors(adj,f,f);
//     ncon         = connect.size();
//     for(i = 0; i < ncon; i++){
//       pos = connect[i];
//       activen[pos] = 0;
//       }
//     to          = selectnode(activen);
//     adj(f, to)  = adj(to, f) = 1;
//     degrees[to] ++;
//     activen[to] = 0;
//     for(j = 1; j < m; j++){
//       if(p > puni()){
//         uncneighbor = unconnectedneighbor(adj, f, to);
//         if(uncneighbor >= 0){
//           adj(uncneighbor, f) = adj(uncneighbor, f) = 1;
//           degrees[uncneighbor] ++;
//           activen[uncneighbor] = 0;
//           } else {
//           j -= 1;
//           }
//         } else {
//         to            = selectnode(activen);
//         adj(f, to)    = adj(to, f) = 1;
//         degrees[to]   ++;
//         activen[to]   = 0;
//         }
//       }
//     }
//   return adj;
//   }






//
//
//
//
// //std::chrono::high_resolution_clock::time_point t_start, t_end;
// //t_start = std::chrono::high_resolution_clock::now();
//
// //t_end = std::chrono::high_resolution_clock::now();
// //std::cout << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms\n";
//
// //
// // if(p > .9 && p < 1){
// //   std::cout << "Warning: code may run forever\n";
// // }
// // if(p == 1){
// //   std::cout << "Can not work. Reduce p!\n";
// //   return 0;
// // }
