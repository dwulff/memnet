#include <Rcpp.h>
#include <time.h>
#include "helpers.h"
using namespace Rcpp;


// identifies unique set of string vector (productions of one individual)
// and returns ordered vector
// [[Rcpp::export]]
std::vector<std::string> set(std::vector<std::string> v){
  std::sort(v.begin(), v.end());
  std::vector<std::string>::iterator last;
  last = std::unique(v.begin(), v.end());
  v.erase(last, v.end());
  return v;
  }

// get unique set from list of string vectors (productions of all individuals)
// [[Rcpp::export]]
std::vector<std::string> mset(GenericVector dat){
  int i, ndat = dat.size();
  std::vector<std::string> allitems;
  std::vector<std::string> subdat;
  for(i = 0; i < ndat; i++){
    subdat = dat[i];
    allitems.insert(allitems.end(), subdat.begin(), subdat.end());
    }
  return set(allitems);
  }

// find index of entry string vector
// [[Rcpp::export]]
int indx(std::string s, std::vector<std::string> set){
  int pos;
  std::vector<std::string>::iterator iter;
  iter = std::find(set.begin(), set.end(), s);
  pos  = std::distance(set.begin(), iter);
  return pos;
  }


// calucalte all distances between productions
// [[Rcpp::export]]
GenericVector lags(GenericVector dat, int  l){
  int i, j, k, n, lim, ndat = dat.size();
  CharacterVector subdat;
  GenericVector res(2);
  std::vector<std::string> ids;
  std::vector<int> dists;
  std::string itemj, itemk;
  for(i = 0; i < ndat; i++){
    subdat = dat[i];
    n = subdat.size();
    for(j = 0; j < n; j++){
      lim = j + l + 1;
      if(lim > n){
        lim = n;
        }
      for(k = j + 1; k < lim; k++){
        itemj = std::string(subdat[j]);
        itemk = std::string(subdat[k]);
        if(itemj < itemk){
          ids.push_back(itemj + '-' + itemk);
          } else {
          ids.push_back(itemk + '-' + itemj);
          }
        dists.push_back(k - j);
        }
      }
    }
  res[0] = ids;
  res[1] = dists;
  return res;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          STRING HANDLING
//
////////////////////////////////////////////////////////////////////////////////////////////////////

// split string by delimiter
// necessary to split ids from lags function
// [[Rcpp::export]]
std::vector<std::string> strsplit(const std::string & s, const std::string & delim) {
  int end = 0, start = 0;
  std::vector<std::string> parts;
  while( end != std::string::npos ) {
    end   = s.find( delim, start );
    parts.push_back(s.substr(start,end));
    start = int(end) + delim.size();
    }
  return parts;
  }

// translate word productions into indices
// based on vector of unique productions producted
// by set
// [[Rcpp::export]]
NumericMatrix getinds(std::vector<std::string> pairs,
                      std::vector<std::string> unis){
  int i, n = pairs.size();
  NumericMatrix inds(n,2);
  std::vector<std::string> parts;
  std::string pair, start, end;
  for(i = 0; i < n; i++){
    pair  = pairs[i];
    parts = strsplit(pair,"-");
    start = parts[0];
    end   = parts[1];
    inds(i,0) = indx(start,unis);
    inds(i,1) = indx(end,unis);
    }
  return inds;
  }

// split word pairs into separate strings
// [[Rcpp::export]]
CharacterMatrix getpairs(std::vector<std::string> spairs, std::string del){
  int i, n = spairs.size();
  CharacterMatrix res(n,2);
  std::vector<std::string> parts;
  std::string pair;
  for(i = 0; i < n; i++){
    pair = std::string(spairs[i]);
    parts = strsplit(pair, del);
    res(i,0) = parts[0];
    res(i,1) = parts[1];
  }
  return res;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          COUNTS AND PROBABILITIES
//
////////////////////////////////////////////////////////////////////////////////////////////////////


// count occurences in string vector (productions of one individual)
// returns counts in alphabetic order
// [[Rcpp::export]]
std::vector<int> count(std::vector<std::string> v){
  std::vector<int> cnts;
  std::map<std::string, int> counter;
  for(std::vector<std::string>::iterator sit = v.begin(); sit != v.end(); ++sit){
    counter[*sit]++;
    }
  for(std::map<std::string,int>::iterator mit = counter.begin(); mit != counter.end(); ++mit) {
    cnts.push_back( mit -> second );
    }
  return cnts;
  }

// count occurences across lists of string vectors (productions of
// all individuals)
// returns counts in alphabetic order
// [[Rcpp::export]]
std::vector<int> mcount(GenericVector dat){
  int i, ndat = dat.size();
  std::vector<std::string> allitems;
  std::vector<std::string> subdat;
  for(i = 0; i < ndat; i++){
    subdat = dat[i];
    allitems.insert(allitems.end(), subdat.begin(), subdat.end());
    }
  return count(allitems);
  }


// turns counts in relative frequencies
// [[Rcpp::export]]
NumericVector getprob(std::vector<int> counts, double N){
  int i, n = counts.size();
  NumericVector probs(n);
  for(i = 0; i < n; i++){
    probs[i] = double(counts[i]) / N;
  }
  return probs;
}

// calculates probability to be in window
// see Goni et al. (2009)
// [[Rcpp::export]]
double pinwin(double n, double l){
  return (2 / ( n * ( n - 1 ))) * (l * n - (l * (l + 1)) / 2);
}

// determine average probability to be in window
// [[Rcpp::export]]
double mpinwin(NumericVector ns, double l){
  int  i, nns = ns.size();
  double n, p = 0;
  for(i = 0; i < nns; i++){
    n = double(ns[i]);
    p += pinwin(n,l);
  }
  return p / double(nns);
}

// determine the number of productions across lists of
// prodcutions (all individuals)
// [[Rcpp::export]]
NumericVector lens(GenericVector dat){
  int i, n = dat.size();
  NumericVector lengths(n);
  CharacterVector subdat;
  for(i = 0; i < n; i++){
    subdat = dat[i];
    lengths[i] = subdat.size();
  }
  return lengths;
}

// determine mean number of productions
// [[Rcpp::export]]
double mlength(GenericVector dat){
  int i, n = 0, ndat = dat.size();
  CharacterVector subdat;
  for(i = 0; i < ndat; i++){
    subdat = dat[i];
    n += subdat.size();
  }
  return n/double(ndat);
}

// calculate probability of a link as the product
// of the base probability of a pair's first production,
// the probility of the pair's second production, and the
// probability that they co-occur within the window size.
// [[Rcpp::export]]
NumericVector getplink(NumericMatrix inds,
                       NumericVector probs,
                       double pinwin){
  int i, indstart, indend, n = inds.nrow();
  double pstart, pend;
  NumericVector plinked(n);
  std::vector<std::string> parts;
  std::string pair, start, end;
  for(i = 0; i < n; i++){
    indstart = inds(i,0);
    indend   = inds(i,1);
    pstart   = probs[indstart];
    pend     = probs[indend];
    plinked[i] = pstart * pend * pinwin;
    }
  return plinked;
  }

// point probability of the binomial distribution
// [[Rcpp::export]]
double dbinom(int k, int n, double p){
  double coef = noverk(n,k);
  double lik  = pow(p, k) * pow(1.0-p, n-k);
  return coef * lik;
  }

// point probability of the cumulative binomial distribution
// [[Rcpp::export]]
double pbinom(int k, int n, double p){
  int i;
  double prob = 0;
  if(k < n - k){
    for(i = 0; i <= k; i++){
      prob += dbinom(double(i),n,p);
      }
    prob = 1.0 - prob;
    } else {
    for(i = 0; i <= n-k; i++){
      prob += dbinom(double(i),n,p);
      }
    }
  return prob;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          GONI GRAPH
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//' Create Goni graph
//'
//' Creare Goni graph from verbal fluency data. Function creates a graph
//' by adding edges for words that occur within a window size \code{l}
//' more frequently than \code{min_cooc} and \code{100(1-crit)%} of
//' chance productions.
//'
//'
//' @param dat
//' @param l an integer specifying the window size. The internal upper limit
//'   of \code{l} is the number of productions.
//' @param crit a numeric within \code{[0,1]} specifiying the type-1 error
//'   rate of including an edge between unconnected words.
//' @param min_cooc integer specifying the minimum number of times two words
//'   have to coocur within a window size of \code{l} to consider including
//'   an edge between them.
//'
//' @return
//' A matrix
//'
// [[Rcpp::export]]
CharacterMatrix goni_graph(
  GenericVector dat,
  int l        = 3,
  int min_cooc = 1,
  double crit  = .05
  ){
  //std::chrono::high_resolution_clock::time_point t_start, t_end;
  int i, nuni, ndat, npairs, cnt, igr = 0;
  double probinwin, p, ptest;
  std::vector<std::string> unis, pairsinwin, unipairsinwin;
  NumericVector freqs, probs, cntpairsinwin, plinked;
  CharacterMatrix pairs;
  NumericMatrix pairinds;
  std::vector<std::string> starts;
  std::vector<std::string> ends;
  std::vector<int> cnts;

  //t_start = std::chrono::high_resolution_clock::now();

  //number of tests
  ndat = dat.size();
  unis = mset(dat);   //in alphabetic order
  nuni = unis.size();

  //determine frequency of pairs
  pairsinwin    = lags(dat, l)[0];
  unipairsinwin = set(pairsinwin);
  cntpairsinwin = count(pairsinwin);

  //determine probabilities
  cnts      = mcount(dat);
  probs     = getprob(cnts,ndat);
  probinwin = mpinwin(lens(dat), l);
  pairinds  = getinds(unipairsinwin,unis);
  plinked   = getplink(pairinds,probs,probinwin);

  //determine probabilities
  pairs     = getpairs(unipairsinwin, "-");
  npairs    = plinked.size();

  //loop over pairs
  for(i = 0; i < npairs; i++){
    cnt = cntpairsinwin[i];
    p   = plinked[i];
    //higher count than min?
    if(cnt > min_cooc){
      if(cnt > ndat){
        cnt = ndat;
        }
      ptest = pbinom(cnt,ndat,double(p));
      //std::cout << pairs;
      if(ptest < crit){
        starts.push_back(std::string(pairs(i,0)));
        ends.push_back(std::string(pairs(i,1)));
        igr ++;
        }
      }
    }

  //fill graph
  CharacterMatrix graph(igr, 2);
  for(i = 0; i < igr; i++){
    graph(i,0) = starts[i];
    graph(i,1)  = ends[i];
    }

  //t_end = std::chrono::high_resolution_clock::now();
  //std::cout << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms\n";

  return graph;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          RW GRAPH
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//' Create random walk graph
//'
//' Creare random walk graph from verbal fluency data. Function creates a graph
//' by adding edges for words that occur within a window size of 1.
//'
//' @param dat
//'
//' @return
//' A matrix
//'
// [[Rcpp::export]]
CharacterMatrix rw_graph(
    GenericVector dat
    ){
  //std::chrono::high_resolution_clock::time_point t_start, t_end;
  int i, npairs;
  std::vector<std::string> pairsinwin, unipairsinwin;
  CharacterMatrix pairs;

  //t_start = std::chrono::high_resolution_clock::now();

  //determine frequency of pairs
  pairsinwin    = lags(dat, 1)[0];
  unipairsinwin = set(pairsinwin);
  pairs     = getpairs(unipairsinwin, "-");
  npairs    = pairs.nrow();

  //fill graph
  CharacterMatrix graph(npairs, 2);
  for(i = 0; i < npairs; i++){
    graph(i,0) = std::string(pairs(i,0));
    graph(i,1)  = std::string(pairs(i,1));
    }

  //t_end = std::chrono::high_resolution_clock::now();
  //std::cout << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms\n";

  return graph;
  }

