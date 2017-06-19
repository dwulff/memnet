#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          SEARCH HELPERS
//
////////////////////////////////////////////////////////////////////////////////////////////////////

// get list of neighbors of a vertex
// [[Rcpp::export]]
std::vector<int> getneighbors(GenericVector adjlist, int pos){
  return adjlist[pos];
  }

// draw neighbor at random
// [[Rcpp::export]]
int getnext(std::vector<int> neighbors){
  int nn, next;
  nn   = neighbors.size();
  next = rint(nn);
  return neighbors[next];
  }


// unique and cut
// return unique entries with length n or length(unique) if length(unique) < n
// used for fast search
// [[Rcpp::export]]
std::vector<int> unicut(std::vector<int> vs, int n){
  //std::sort(v.begin(), v.end());
  std::vector<int> nv;
  int vn = vs.size();
  for(int i = 0; i < vn; ++i){
    if(nv.size() >= n) break;
    int v = vs[i];
    if(inset(v,nv) == false){
      nv.push_back(v);
      //std::cout << v << '\n';
      }
    }
  return nv;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          FLUENCY PRODUCTIONS
//
////////////////////////////////////////////////////////////////////////////////////////////////////

//' Verbal fluency generator
//'
//' Generates verbal fluency data using a switcher-random walk process.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting a neighbor with
//' probability 1-\code{pjump} or jumps to a random place in the network
//' with probability \code{pjump}. Where the random walk process enters
//' the network and where it jumps to is additionally controlled
//' by \code{type}. Neighbors are always selected uniformly.
//'
//' @param adjlist a list containing row indices for adjacent vertices as created
//'   by \link{get_adjlist}.
//' @param n integer specifying the number of productions.
//' @param pjump numeric specifying the probability of a jump.
//' @param type integer controlling network start and jump vertices.
//'   For \code{type = 0} the process selects the start vertex and any jump
//'   vertices proportional to their degree. For \code{type = 1} the process
//'   selects a random vertex to serve as the start vertex and the jump vertex.
//'   For \code{type = 2} the process selects the start and any jump vertices
//'   uniformly.
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> fluency(GenericVector adjlist, int n, double pjump, int type){
  int i, j, start, npos, reset, nitem = adjlist.size();
  bool jump = false;
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;
  if(type == 1 || type == 2){
    start = rint(nitem);
    items.push_back(start);
    } else {
    for(i = 0; i < nitem; i++){
      neighbors = getneighbors(adjlist,i);
      degrees.push_back(neighbors.size());
      }
    start = smpl(degrees);
    items.push_back(start);
    }
  for(j = 1; j < n; j++){
    i = 0;
    if(runi() > pjump || jump == true){
      jump = false;
      if(type == 1) {
        reset = start;
        }
      else if (type == 2) {
        reset = rint(nitem);
        }
      else {
        reset = smpl(degrees);
        }
      neighbors = getneighbors(adjlist,reset);
      npos      = getnext(neighbors);
      while(inset(npos,items)){
        i++;
        neighbors = getneighbors(adjlist,npos);
        npos      = getnext(neighbors);
        if(i == 10) {
          jump = true;
          break;
          }
        }
      } else {
      neighbors = getneighbors(adjlist,items[j-1]);
      npos      = getnext(neighbors);
      while(inset(npos,items)){
        i++;
        neighbors = getneighbors(adjlist,npos);
        npos      = getnext(neighbors);
        if(i == 10) {
          jump = true;
          break;
        }
      }
    }
    if(jump == false) {
      items.push_back(npos);
      } else {
      j = j - 1;
      }
    }
  return items;
  }


//' Verbal fluency generator wrapper
//'
//' Generates multiple verbal fluency sequences using \code{fluency}.
//'
//' For details see \link{fluency}.
//'
//' @inheritParams fluency
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
GenericVector mfluency(GenericVector adjlist, NumericVector n, double pjump = 0, int type = 1){
  int i, nsub = n.size();
  GenericVector data(nsub);
  for(i = 0; i < nsub; i++){
    data[i] = tostring(fluency(adjlist,n[i],pjump,type));
    }
  return data;
  }


//' Fast verbal fluency generator
//'
//' Generates verbal fluency data using a switcher-random walk process.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting neighbors with
//' probability 1-\code{pjump} or jumps to a random place in the network
//' with probability \code{pjump}. How the random walk process enters
//' the network and how it jumps to is additionally controlled
//' by \code{random}. Neighbors are always selected uniformly.
//'
//' In contrast to \link{fluency}, does not check at every step whether
//' the sampled neighbor is already in the list of productions. Instead
//' \code{ffluency} simply returns the list of unique productions. This means
//' that if repetitions occur \code{ffluency} will produce sequences of length
//' \code{min(n*3 - k, n)} where k is the number of repeptitions.
//'
//' @inheritParams fluency
//' @param n integer specifying the number of productions.
//' @param random bool controlling jump vertices.
//'   For \code{random = TRUE} the process selects jump at random. For
//'   \code{random = FALSE} the process always jumps back to the start vertex. The
//'   start vertex is always selected at random.
//'
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> ffluency(GenericVector adjlist, int n, double pjump, bool random, bool pref_start = false){
  int j, start, npos, reset, nitem = adjlist.size();
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;
  if(pref_start == false){
    start = rint(nitem);
    items.push_back(start);
    } else {
    for(int i = 0; i < nitem; i++){
      neighbors = getneighbors(adjlist,i);
      degrees.push_back(neighbors.size());
      }
    start = smpl(degrees);
    items.push_back(start);
    }
  for(j = 1; j < n * 3; j++){
    if(runi() < pjump){
      if(random == false) {
        reset = start;
        } else {
        reset = rint(nitem);
        }
      neighbors = getneighbors(adjlist,reset);
      npos      = getnext(neighbors);
      items.push_back(npos);
      } else {
      neighbors = getneighbors(adjlist,items[j-1]);
      npos      = getnext(neighbors);
      items.push_back(npos);
      }
    }
  return unicut(items,n);
  }


//' Fast verbal fluency generator wrapper
//'
//' Generates multiple verbal fluency sequences using \code{fluency}.
//'
//' For details see \link{ffluency}.
//'
//' @inheritParams ffluency
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
GenericVector mfsearch(GenericVector adjlist, NumericVector n, double pjump = 0, bool random = false, bool pref_start = false){
  int i, nsub = n.size();
  GenericVector data(nsub);
  for(i = 0; i < nsub; i++){
    data[i] = tostring(ffluency(adjlist,n[i],pjump,random,pref_start));
    }
  return data;
  }


//' Exhaustive verbal fluency generator
//'
//' Generates verbal fluency data using a switcher-random walk process.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting neighbors with
//' probability 1-\code{pjump} or jumps to a random place in the network
//' with probability \code{pjump}. How the random walk process enters
//' the network and how it jumps to is additionally controlled
//' by \code{random}. Neighbors are always selected uniformly.
//'
//' In contrast to \link{fluency}, does not check at every step whether
//' the sampled neighbor is already in the list of productions. Instead
//' \code{ffluency} simply returns the list of unique productions. This means
//' that if repetitions occur \code{ffluency} will produce sequences of length
//' \code{min(n*3 - k, n)} where k is the number of repeptitions.
//'
//' @inheritParams fluency
//' @param n integer specifying the number of productions.
//' @param random bool controlling jump vertices.
//'   For \code{random = TRUE} the process selects jump at random. For
//'   \code{random = FALSE} the process always jumps back to the start vertex. The
//'   start vertex is always selected at random.
//'
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> efluency(GenericVector adjlist, int n, double pjump, bool random, bool pref_start = false){
  int j, start, npos, reset, nitem = adjlist.size();
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;
  if(pref_start == false){
    start = rint(nitem);
    items.push_back(start);
  } else {
    for(int i = 0; i < nitem; i++){
      neighbors = getneighbors(adjlist,i);
      degrees.push_back(neighbors.size());
    }
    start = smpl(degrees);
    items.push_back(start);
  }
  int cnt = 1;
  while(true){
    if(runi() < pjump){
      if(random == false) {
        reset = start;
        } else {
        reset = rint(nitem);
        }
      neighbors = getneighbors(adjlist,reset);
      npos      = getnext(neighbors);
      if(inset(npos,items) == false) cnt++;
      items.push_back(npos);
      } else {
      neighbors = getneighbors(adjlist,items.back());
      npos      = getnext(neighbors);
      if(inset(npos,items) == false) cnt++;
      items.push_back(npos);
      }
    //std::cout << cnt << '\n';
    if(cnt == n) break;
    }
  return unicut(items,n);
  }


//' Verbal fluency generator wrapper
//'
//' Generates multiple verbal fluency sequences using \code{efluency}.
//'
//' For details see \link{efluency}.
//'
//' @inheritParams efluency
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
GenericVector mesearch(GenericVector adjlist, NumericVector n, double pjump = 0, bool random = false, bool pref_start = false){
  int i, nsub = n.size();
  GenericVector data(nsub);
  for(i = 0; i < nsub; i++){
    data[i] = tostring(efluency(adjlist,n[i],pjump,random,pref_start));
  }
  return data;
}



//' Verbal fluency step counter
//'
//' Generates verbal fluency data using a switcher-random walk process.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting neighbors with
//' probability 1-\code{pjump} or jumps to a random place in the network
//' with probability \code{pjump}. How the random walk process enters
//' the network and how it jumps to is additionally controlled
//' by \code{random}. Neighbors are always selected uniformly.
//'
//' In contrast to \link{fluency}, does not check at every step whether
//' the sampled neighbor is already in the list of productions. Instead
//' \code{ffluency} simply returns the list of unique productions. This means
//' that if repetitions occur \code{ffluency} will produce sequences of length
//' \code{min(n*3 - k, n)} where k is the number of repeptitions.
//'
//' @inheritParams fluency
//' @param n integer specifying the number of productions.
//' @param random bool controlling jump vertices.
//'   For \code{random = TRUE} the process selects jump at random. For
//'   \code{random = FALSE} the process always jumps back to the start vertex. The
//'   start vertex is always selected at random.
//'
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
int sfluency(GenericVector adjlist, int n, double pjump, bool random, bool pref_start = false){
  int j, start, npos, reset, nitem = adjlist.size();
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;
  if(pref_start == false){
    start = rint(nitem);
    items.push_back(start);
  } else {
    for(int i = 0; i < nitem; i++){
      neighbors = getneighbors(adjlist,i);
      degrees.push_back(neighbors.size());
    }
    start = smpl(degrees);
    items.push_back(start);
  }
  int cnt = 1;
  while(true){
    if(runi() < pjump){
      if(random == false) {
        reset = start;
      } else {
        reset = rint(nitem);
      }
      neighbors = getneighbors(adjlist,reset);
      npos      = getnext(neighbors);
      if(inset(npos,items) == false) cnt++;
      items.push_back(npos);
    } else {
      neighbors = getneighbors(adjlist,items.back());
      npos      = getnext(neighbors);
      if(inset(npos,items) == false) cnt++;
      items.push_back(npos);
    }
    //std::cout << cnt << '\n';
    if(cnt == n) break;
  }

  return items.size() - n;
}


//' Exhaustive verbal fluency generator wrapper
//'
//' Generates multiple verbal fluency sequences using \code{efluency}.
//'
//' For details see \link{efluency}.
//'
//' @inheritParams efluency
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> mssearch(GenericVector adjlist, NumericVector n, double pjump = 0, bool random = false, bool pref_start = false){
  int i, nsub = n.size();
  std::vector<int> data(nsub);
  for(i = 0; i < nsub; i++){
    data[i] =sfluency(adjlist,n[i],pjump,random,pref_start);
  }
  return data;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          RANDOM WALKS
//
////////////////////////////////////////////////////////////////////////////////////////////////////


//' Random walk
//'
//' Traverses a network using a switcher-random walk process and records the earliest
//' visit to vertices of interest.
//'
//' Beginning at a given start vertext, function traverses a network using switcher
//' random walk and records for each of a list of vertices of interest the index at which
//' the respective vertices have been visited first.
//'
//' If a vertex specified in \code{observe} has never been visited then the function
//' returns \code{nmax} for that vertex.
//'
//' @inheritParams fluency
//' @param start index of the start vertix.
//' @param observe integer vector specifying the vertices whose first visits should be recorded.
//' @param nmax integer specifying the maximum number of steps.
//'
//' @return Numeric, 3 column matrix containing in each row the start vertex, the end vertex, and
//' the (minimum) number of steps it took to reach the end vertext from the start vertex.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix rwalk(GenericVector adjlist, int start, std::vector<int> observe, int nmax = 1000, double pjump = 0){

  // scalars
  int i, npos = start - 1, ind = 0, nitem = adjlist.size(), nobs = observe.size();

  // containers
  std::vector<int> neighbors;
  std::map<std::pair<int, int>, int > m;
  std::pair<int, int> pair;

  // adjust arguments
  start = start - 1;
  for(i = 0; i < nobs; ++i) observe[i]--;

  // fill map with nmax
  for(i = 0; i < nobs; ++i){
    pair = std::make_pair(start, observe[i]);
    m[pair] = nmax;
    }

  // iterate until nmax
  for(i = 0; i < nmax; ++i){
    ind++;

    // jump
    if(runi() < pjump){
      neighbors = getneighbors(adjlist,rint(nitem));
      npos      = getnext(neighbors) - 1;
      if(inset(npos,observe)){
        pair = std::make_pair(start, npos);
        if(m[pair] == nmax) m[pair] = ind;
        }

      // no jump
      } else {
      neighbors = getneighbors(adjlist,npos);
      //std::cout<< neighbors[1] << '_' << neighbors.size() << '\n';
      npos      = getnext(neighbors) - 1;
      if(inset(npos,observe)){
        pair = std::make_pair(start, npos);
        if(m[pair] == nmax) m[pair] = ind;
        }
      }
    }

  // extract and store results
  NumericMatrix res(nobs,3);
  std::map<std::pair<int, int>, int>::const_iterator it;
  for(it = m.begin(), i = 0; it != m.end(); ++it, ++i){
    std::pair<int, int> value = it->first;
    res(i,0) = value.first + 1;
    res(i,1) = value.second + 1;
    res(i,2) = it->second;
    }

  return res;
  }

//' Random walk
//'
//' Traverses a network using a switcher-random walk process and records the earliest
//' visit to vertices of interest.
//'
//' Beginning at a given start vertext, function traverses a network using switcher
//' random walk and records for each of a list of vertices of interest the index at which
//' the respective vertices have been visited first.
//'
//' If a vertex specified in \code{observe} has never been visited then the function
//' returns \code{nmax} for that vertex.
//'
//' @inheritParams fluency
//' @param start index of the start vertix.
//' @param observe integer vector specifying the vertices whose first visits should be recorded.
//' @param nmax integer specifying the maximum number of steps.
//'
//' @return Numeric, 3 column matrix containing in each row the start vertex, the end vertex, and
//' the (minimum) number of steps it took to reach the end vertext from the start vertex.
//'
//' @export
// [[Rcpp::export]]

NumericMatrix mrwalk(GenericVector adjlist,
                     int start, std::vector<int> observe,
                     int nrep = 100, bool aggregate = true,
                     int nmax = 1000, double pjump = 0){

  // scalars
  int nobs = observe.size();

  // containers
  if(aggregate == true){
    NumericMatrix agg = rwalk(adjlist,start,observe,nmax,pjump);
    for(int i = 1; i < nrep; ++i){
      NumericMatrix tmp = rwalk(adjlist,start,observe,nmax,pjump);
      for(int j = 0; j < nobs; ++j){
        agg(j,2) +=  tmp(j,2);
        }
      }
    for(int j = 0; j < nobs; ++j) agg(j,2) /= nrep;
    return agg;
    } else {
    NumericMatrix res(nrep * nobs,3);
    int ind = 0;
    for(int i = 1; i < nrep; ++i){
      NumericMatrix tmp = rwalk(adjlist,start,observe,nmax,pjump);
      for(int j = 0; j < nobs; ++j){
        res(ind,_) = tmp(j,_);
        ind++;
        }
      }
    return res;
    }
  return 0;
  }


