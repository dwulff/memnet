// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/lexical_cast.hpp>
#include "helpers.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          SEARCH HELPERS
//
////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int smpl(std::vector<double> ps){
  int k, i = 0, n = ps.size();
  double v, sum = 0, r = double((Rcpp::runif(1,0,1))[0]);
  for(k = 0; k < n; k++){
    sum += ps[k];
  }
  v = ps[0] / double(sum);
  while(i < n && v <= r){
    i++;
    v += ps[i] / double(sum);
  }
  return i;
}

// [[Rcpp::export]]
std::string my_to_string(int value) {
  std::string stri = boost::lexical_cast<std::string>(value);
  return stri;
  }

// convert integer vector into Character vector
// [[Rcpp::export]]
Rcpp::CharacterVector to_str(std::vector<int> items){
  int n = items.size();
  Rcpp::CharacterVector res(n);
  for(int i = 0; i < n; i++) {
    //res[i] = std::to_string(items[i]);
    res[i] = my_to_string(items[i]);
  }
  return res;
}

// get list of neighbors of a node
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


// [[Rcpp::export]]
GenericVector adjlist_minus1(GenericVector &adjlist){
  int n = adjlist.size();
  GenericVector new_adjlist(n);
  for(int i  = 0; i < n; ++i){
    std::vector<int> neighbors = adjlist[i];
    int n_neighbors = neighbors.size();
    for(int j = 0; j < n_neighbors; ++j){
      neighbors[j]--;
      }
    new_adjlist[i] = neighbors;
    }
  return new_adjlist;
  }

// [[Rcpp::export]]
void add_1(std::vector<int> &items){
  int n = items.size();
  for(int i  = 0; i < n; ++i) items[i]++;
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
//' probability \code{1-pjump} or jumps to a random place in the network
//' with probability \code{pjump}. Where the random walk process enters
//' the network and where it jumps to is further controlled
//' by \code{type}. Neighbors are always selected uniformly.
//'
//' @param adj_list a list containing row indices of nodes adjacent node to the ith
//'   node as created by \link{get_adjlist}.
//' @param n integer specifying the number of unique productions.
//' @param pjump numeric specifying the probability of a jump.
//' @param type integer controlling network start and jump nodes.
//'   For \code{type = 0} the process selects the start node and any jump
//'   nodes proportional to their degree. For \code{type = 1} the process
//'   selects a random node to serve both as the start node and the jump node.
//'   For \code{type = 2} the process selects the start and any jump nodes
//'   uniformly at random.
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # create verbal fluency sequence
//' one_fluency(get_adjlist(network), 10)
//'
//' # create verbal fluency sequence
//' # with high jump probability
//' one_fluency(get_adjlist(network), 10, pjump = .5)
//'
//' @export
// [[Rcpp::export]]
std::vector<int> one_fluency(GenericVector adj_list,
                             int n,
                             double pjump = 0,
                             int type = 0){
  int i, j, start, npos, reset, nitem = adj_list.size();
  bool jump = false;
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;

  // reset adjlist
  GenericVector adjlist = adjlist_minus1(adj_list);

  // enter process
  if(type == 1 || type == 2){
    start = rint(nitem);
    items.push_back(start);
    } else {
    for(i = 0; i < nitem; i++){
      neighbors = getneighbors(adjlist, i);
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
      neighbors = getneighbors(adjlist, reset);
      npos      = getnext(neighbors);
      while(inset(npos,items)){
        i++;
        neighbors = getneighbors(adjlist, npos);
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

  //add one
  add_1(items);

  //out
  return items;
  }


//' Repeated verbal fluency generator.
//'
//' Generates multiple verbal fluency sequences using \link{one_fluency}.
//'
//' For details see \link{one_fluency}.
//'
//' @inheritParams one_fluency
//' @param adjlist a list containing row indices of nodes adjacent node to the ith
//'   node as created by \link{get_adjlist}.
//' @param n integer vector specifying for each sequence the number of
//'   unique productions.
//' @param string logical specifying whether the output should be of mode character.
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 3)
//'
//' # create verbal fluency sequences
//' fluency(get_adjlist(network), c(10, 10))
//'
//' # create verbal fluency sequence
//' # with high jump probability
//' fluency(get_adjlist(network), c(10, 10), pjump = .5)
//'
//' @export
// [[Rcpp::export]]
GenericVector fluency(GenericVector adjlist,
                      NumericVector n,
                      double pjump = 0,
                      int type = 0,
                      bool string = false){
  int i, nsub = n.size();
  GenericVector data(nsub);
  for(i = 0; i < nsub; i++){
    if(string == true){
    data[i] = to_str(one_fluency(adjlist,n[i],pjump,type));
    } else {
    data[i] = one_fluency(adjlist,n[i],pjump,type);
    }
  }
  return data;
  }


//' Fast verbal fluency generator
//'
//' Generates verbal fluency data using a switcher-random walk process.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting a neighbor with
//' probability \code{1-pjump} or jumps to a random place in the network
//' with probability \code{pjump}. Where the random walk process enters
//' the network and where it jumps to is further controlled
//' by \code{type}. Neighbors are always selected uniformly.
//'
//' In contrast to \link{fluency}, this function does not check at every step
//' whether the sampled neighbor is already in the list of productions. Instead,
//' \code{ffluency} simply returns the list of unique productions. This means
//' that if repetitions occur \code{ffluency} will produce sequences of length
//' \code{min(n * 3 - k, n)} where k is the number of repeptitions.
//'
//' @inheritParams one_fluency
//' @param n integer specifying the maximum number of productions. Function may
//' return fewer than \code{n}.
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # create verbal fluency sequences
//' one_ffluency(get_adjlist(network),  10)
//'
//' # create verbal fluency sequence
//' # with high jump probability
//' one_ffluency(get_adjlist(network), 10, pjump = .5)
//'
//' @export
// [[Rcpp::export]]
std::vector<int> one_ffluency(GenericVector adj_list,
                              int n,
                              double pjump = 0,
                              int type = 0){
  int i, j, start, npos, reset, nitem = adj_list.size();
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;

  // reset adjlist
  GenericVector adjlist = adjlist_minus1(adj_list);

  // enter process
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
  for(j = 1; j < n * 3; j++){
    if(runi() < pjump){
      if(type == 1) {
        reset = start;
        }
      else if (type == 2) {
        reset = rint(nitem);
        }
      else {
        reset = smpl(degrees);
        }
      neighbors = getneighbors(adjlist, reset);
      npos      = getnext(neighbors);
      items.push_back(npos);
      } else {
      neighbors = getneighbors(adjlist,items[j-1]);
      npos      = getnext(neighbors);
      items.push_back(npos);
      }
    }

  // reduce to unique
  items = unicut(items,n);

  // add 1
  add_1(items);

  // out
  return items;
  }


//' Fast verbal fluency generator
//'
//' Generates multiple verbal fluency sequences using \link{one_ffluency}.
//'
//' For details see \link{one_ffluency}.
//'
//' @inheritParams one_ffluency
//' @inheritParams fluency
//' @param n integer vector specifying for each sequence the maximum numbers of
//' productions. Function may return fewer than \code{n}.
//' @param string logical specifying whether the output should be of mode character.
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # create verbal fluency sequences
//' ffluency(get_adjlist(network), c(10, 10))
//'
//' # create verbal fluency sequence
//' # with high jump probability
//' ffluency(get_adjlist(network), c(10, 10), pjump = .5)
//'
//' @export
// [[Rcpp::export]]
GenericVector ffluency(GenericVector adjlist,
                       NumericVector n,
                       double pjump = 0,
                       int type = 0,
                       bool string = false){
  int i, nsub = n.size();
  GenericVector data(nsub);
  for(i = 0; i < nsub; i++){
    if(string == true){
      data[i] = to_str(one_ffluency(adjlist,n[i],pjump,type));
      } else {
      data[i] = one_ffluency(adjlist,n[i],pjump,type);
      }
    }
  return data;
  }


// //' Exhaustive verbal fluency generator
// //'
// //' Generates verbal fluency data using a switcher-random walk process.
// //'
// //' Function produces verbal fluency data via a switcher random walk
// //' process that traverses the network by selecting a neighbor with
// //' probability \code{1-pjump} or jumps to a random place in the network
// //' with probability \code{pjump}. Where the random walk process enters
// //' the network and where it jumps to is further controlled
// //' by \code{type}. Neighbors are always selected uniformly.
// //'
// //' In contrast to \link{fluency}, does not check at every step whether
// //' the sampled neighbor is already in the list of productions. In contrast to
// //' \link{ffluency}, this function will always produce sequences of length
// //' \code{n}, by continuing the search until \code{n} nodes have been visited.
// //'
// //' @inheritParams fluency
// //'
// //' @return Integer vector containing the indices of the fluency productions.
// //'   Indices refer to the row of the item in the original adjacency matrix. See
// //'   \link{get_adjlist}.
// //'
// //' @export
// // [[Rcpp::export]]
// std::vector<int> efluency(GenericVector adjlist, int n, double pjump, int type){
//   int i, start, npos, reset, nitem = adjlist.size();
//   std::vector<int> items;
//   std::vector<int> neighbors;
//   std::vector<double> degrees;
//   if(type == 1 || type == 2){
//     start = rint(nitem);
//     items.push_back(start);
//   } else {
//     for(i = 0; i < nitem; i++){
//       neighbors = getneighbors(adjlist,i);
//       degrees.push_back(neighbors.size());
//     }
//     start = smpl(degrees);
//     items.push_back(start);
//   }
//   int cnt = 1;
//   while(true){
//     if(runi() < pjump){
//       if(type == 1) {
//         reset = start;
//         }
//       else if (type == 2) {
//         reset = rint(nitem);
//         }
//       else {
//         reset = smpl(degrees);
//         }
//       neighbors = getneighbors(adjlist,reset);
//       npos      = getnext(neighbors);
//       if(inset(npos,items) == false) cnt++;
//       items.push_back(npos);
//       } else {
//       neighbors = getneighbors(adjlist,items.back());
//       npos      = getnext(neighbors);
//       if(inset(npos,items) == false) cnt++;
//       items.push_back(npos);
//       }
//     //std::cout << cnt << '\n';
//     if(cnt == n) break;
//     }
//   return unicut(items,n);
//   }
//
//
// //' Verbal fluency generator wrapper
// //'
// //' Generates multiple verbal fluency sequences using \code{efluency}.
// //'
// //' For details see \link{efluency}.
// //'
// //' @inheritParams efluency
// //'
// //' @return List of character vectors containing the indices of the fluency productions.
// //'   Indices refer to the row of the item in the original adjacency matrix. See
// //'   \link{get_adjlist}.
// //'
// //' @export
// // [[Rcpp::export]]
// GenericVector n_esearch(GenericVector adjlist,
//                         NumericVector n,
//                         double pjump = 0,
//                         int type = 1){
//   int i, nsub = n.size();
//   GenericVector data(nsub);
//   for(i = 0; i < nsub; i++){
//     data[i] = to_str(efluency(adjlist,n[i], pjump, type));
//   }
//   return data;
// }



//' Verbal fluency step counter
//'
//' Generates verbal fluency data using a switcher-random walk process and counts
//' the number of steps required to produce \code{n} unique responses.
//'
//' Function produces verbal fluency data via a switcher random walk
//' process that traverses the network by selecting a neighbor with
//' probability \code{1-pjump} or jumps to a random place in the network
//' with probability \code{pjump}. Where the random walk process enters
//' the network and where it jumps to is further controlled
//' by \code{type}. Neighbors are always selected uniformly.
//'
//' In contrast to \link{fluency} and \code{ffluency}, returns the number of steps
//' required to produce a sequence of unique productions, rather than the
//' productions itself.
//'
//' @inheritParams one_fluency
//' @param n integer specifying the number of productions.
//'
//'
//' @return Integer vector containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # count number of steps needed to create sequence
//' one_fluency_steps(get_adjlist(network), 10)
//'
//' # count number of steps needed to create sequence
//' # with high jump probability
//' one_fluency_steps(get_adjlist(network), 10, pjump = .5)
//'
//' @export
// [[Rcpp::export]]
int one_fluency_steps(GenericVector adj_list, int n, double pjump = 0 , int type = 0){
  int i, start, npos, reset, nitem = adj_list.size();
  std::vector<int> items;
  std::vector<int> neighbors;
  std::vector<double> degrees;

  // reset adjlist
  GenericVector adjlist = adjlist_minus1(adj_list);

  // enter process
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
  int cnt = 1;
  while(true){
    if(runi() < pjump){
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
      if(inset(npos, items) == false) cnt++;
      items.push_back(npos);
    } else {
      neighbors = getneighbors(adjlist,items.back());
      npos      = getnext(neighbors);
      if(inset(npos, items) == false) cnt++;
      items.push_back(npos);
    }
    //std::cout << cnt << '\n';
    if(cnt == n) break;
  }

  return items.size() - n;
}


//' Verbal fluency step counter
//'
//' Repeatedly generates verbal fluency data using \link{one_fluency_steps} and
//' counts the number of steps required to produce \code{n} unique responses.
//'
//' For details see \link{one_fluency_steps}.
//'
//' @inheritParams one_fluency_steps
//' @inheritParams fluency
//' @param n integer vector specifying the numbers of production.
//'
//' @return List of character vectors containing the indices of the fluency productions.
//'   Indices refer to the row of the item in the original adjacency matrix. See
//'   \link{get_adjlist}.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # count number of steps needed to create sequence
//' fluency_steps(get_adjlist(network), c(10, 10))
//'
//' # count number of steps needed to create sequence
//' # with high jump probability
//' fluency_steps(get_adjlist(network), c(10, 10), pjump = .5)
//'
//' @export
// [[Rcpp::export]]
std::vector<int> fluency_steps(GenericVector adjlist,
                               NumericVector n,
                               double pjump = 0,
                               int type = 0){
  int i, nsub = n.size();
  std::vector<int> data(nsub);
  for(i = 0; i < nsub; i++){
    data[i] = one_fluency_steps(adjlist,n[i],pjump,type);
  }
  return data;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          RANDOM WALKS
//
////////////////////////////////////////////////////////////////////////////////////////////////////


//' Search network using switcher-random walk process
//'
//' Traverses a network using a switcher-random walk process and records the
//' number of steps required to get from node \code{start} to node \code{observe}.
//'
//' If a node specified in \code{observe} has never been visited then the function
//' returns \code{nmax} for that node.
//'
//' @inheritParams one_fluency
//' @param start integer vector of length 1 specifying the index of
//'   the start node.
//' @param observe integer vector of length 1 or larger specifying the target
//'   end nodes.
//' @param nmax integer specifying the maximum number of steps before search
//'   terminates.
//'
//' @return Numeric, 3 column matrix containing in each row the start node, the
//' end node, and the (minimum) number of steps it took to reach the end node
//' from the start node.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # observe number of steps from node 2
//' # to nodes 3, 4, and 5
//' one_search(get_adjlist(network), 2, c(3, 4, 5))
//'
//' # observe number of steps from node 2 to nodes 3, 4, and 5
//' # with high jump probability
//' one_search(get_adjlist(network), start = 2, observe = c(3, 4, 5), pjump = .5)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix one_search(GenericVector adj_list,
                         int start,
                         std::vector<int> observe,
                         int nmax = 1000,
                         double pjump = 0,
                         int type = 0){

  // scalars
  int i, reset, npos, ind = 0, nitem = adj_list.size(), nobs = observe.size();

  // containers
  std::vector<int> neighbors;
  std::map<std::pair<int, int>, int > m;
  std::pair<int, int> pair;

  // reset adjlist
  GenericVector adjlist = adjlist_minus1(adj_list);

  // adjust arguments
  start = start - 1;
  for(i = 0; i < nobs; ++i) observe[i]--;
  npos = start;

  // fill map with nmax
  for(i = 0; i < nobs; ++i){
    pair = std::make_pair(start, observe[i]);
    m[pair] = nmax;
    }

  // find degrees if necessary
  std::vector<double> degrees;
  for(i = 0; i < nitem; ++i){
    neighbors = getneighbors(adjlist, i);
    degrees.push_back(neighbors.size());
    }

  // iterate until nmax
  for(i = 0; i < nmax; ++i){
    ind++;

    // jump
    if(runi() < pjump){
      if(type == 1) {
        reset = start;
        }
      else if (type == 2) {
        reset = rint(nitem);
        }
      else {
        reset = smpl(degrees);
        }
      neighbors = getneighbors(adjlist, reset);
      npos      = getnext(neighbors);
      if(inset(npos,observe)){
        pair = std::make_pair(start, npos);
        if(m[pair] == nmax) m[pair] = ind;
        }

      // no jump
      } else {
      neighbors = getneighbors(adjlist,npos);
      npos      = getnext(neighbors);
      if(inset(npos,observe)){
        pair = std::make_pair(start, npos);
        if(m[pair] == nmax) m[pair] = ind;
        }
      }
    }

  // extract and store results
  NumericMatrix res(nobs, 3);
  std::map<std::pair<int, int>, int>::const_iterator it;
  for(it = m.begin(), i = 0; it != m.end(); ++it, ++i){
    std::pair<int, int> value = it->first;
    res(i,0) = value.first + 1;
    res(i,1) = value.second + 1;
    res(i,2) = it->second;
    }

  return res;
  }

//' Search network using switcher-random walk process
//'
//' Traverses a network using a switcher-random walk process and records the
//' number of steps required to get from node \code{start} to node \code{observe}.
//'
//' If a node specified in \code{observe} has never been visited then the function
//' returns \code{nmax} for that node.
//'
//' @inheritParams one_search
//' @inheritParams fluency
//' @param start integer vector of length 1 or larger specifying the index of
//'   the start node.
//' @param observe integer vector of length 1 or larger specifying the target
//'   end nodes.
//' @param nmax integer specifying the maximum number of steps before search
//'   terminates.
//'
//' @return Numeric, 3 column matrix containing in each row the start node, the
//' end node, and the (minimum) number of steps it took to reach the end node
//' from the start node.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # observe number of steps from node 2 and 6
//' # to nodes 3, 4, and 5
//' search_rw(get_adjlist(network), c(2, 6), c(3, 4, 5))
//'
//' # observe number of steps from node 2 and 6 to nodes 3, 4, and 5
//' # with high jump probability
//' search_rw(get_adjlist(network), start = c(2, 6), observe = c(3, 4, 5), pjump = .5)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix search_rw(GenericVector adjlist,
                        std::vector<int> start,
                        std::vector<int> observe,
                        int nmax = 1000,
                        double pjump = 0,
                        int type = 0){

  // scalars
  int i, j, ind = 0, n_start = start.size(), n_observe = observe.size();

  // create container for all
  NumericMatrix res(n_start * n_observe, 3);

  // iterate through starts
  for(i = 0; i < n_start; ++i){
    NumericMatrix one = one_search(adjlist,start[i],observe,nmax,pjump,type);

    // iterate through ones
    for(j = 0; j < n_observe; ++j){
      res(ind, _) = one(j, _);
      ind++;
      }
    }

  return res;
  }


//' Search network repeatedly using switcher-random walk process
//'
//' Traverses a network using a switcher-random walk process repeatedly, records
//' the earliest visit to nodes of interest and averages the result.
//'
//' Beginning at a given start node, function traverses a network using switcher
//' random walk and records for each of a list of nodes of interest the index at which
//' the respective nodes have been visited first.
//'
//' If a node specified in \code{observe} has never been visited then the function
//' returns \code{nmax} for that node.
//'
//' @inheritParams fluency
//' @inheritParams search_rw
//' @param observe integer vector specifying the nodes whose first visits should be recorded.
//' @param nmax integer specifying the maximum number of steps.
//' @param nrep integer specifying the number of iterations across which
//'   aggregates are computed.
//'
//' @return Numeric, 3 column matrix containing in each row the start node, the end node, and
//' the (minimum) number of steps it took to reach the end node from the start node.
//'
//' @references
//'
//' Wulff, D. U., Hills, T., & Mata, R. (2018, October 29). Structural
//' differences in the semantic networks of younger and older adults.
//' https://doi.org/10.31234/osf.io/s73dp
//'
//' Goñi, J., Martincorena, I., Corominas-Murtra, B., Arrondo, G., Ardanza-
//' Trevijano, S., & Villoslada, P. (2010). Switcher-random-walks: A cognitive-
//' inspired mechanism for network exploration. International Journal of
//' Bifurcation and Chaos, 20(03), 913-922.
//'
//' @examples
//' # generate watts strogatz graph
//' network = grow_ws(n = 100, k = 10)
//'
//' # determine mean number of steps from node 2 and 6
//' # to nodes 3, 4, and 5
//' search_rw_mean(get_adjlist(network), c(2, 6), c(3, 4, 5))
//'
//' # determine mean number of steps from node 2 and 6 to nodes 3, 4, and 5
//' # with high jump probability
//' search_rw_mean(get_adjlist(network), start = c(2, 6), observe = c(3, 4, 5), pjump = .5)
//'
//' @export
// [[Rcpp::export]]

NumericMatrix search_rw_mean(GenericVector adjlist,
                             std::vector<int> start,
                             std::vector<int> observe,
                             int nmax = 1000,
                             double pjump = 0,
                             int type = 0,
                             int nrep = 100){

  // first time
  NumericMatrix res = search_rw(adjlist, start, observe, nmax, pjump, type);

  // iterate
  for(int i = 0; i < nrep; i++){
    NumericMatrix res_i = search_rw(adjlist, start, observe, nmax, pjump, type);
    res(_,2) = res(_,2) + res_i(_,2);
    }

  // normalize
  res(_,2) = res(_,2) / nrep;

  return res;
  }


