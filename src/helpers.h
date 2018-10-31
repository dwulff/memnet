#ifndef __UTILITIES__
#define __UTILITIES__



//////////////////////////////////////////////////////////////////////////////
//
//          GENERIC UTILITIES
//
//////////////////////////////////////////////////////////////////////////////

// cumulative sum
inline std::vector<double> csum(std::vector<double> x){
  int i, j, n = x.size();
  double cp = 0;
  std::vector<double> cumx(n);
  for(i = 0; i < n; i++){
    for(j = 0; j <= i; j++){
      cp += double(x[j]);
      }
    cumx[i] = cp;
    cp = 0;
    }
  return cumx;
  }

// create diagonal matrix
// with two types of entries for diagnoal and non-diagnoal entries
inline Rcpp::NumericMatrix diagmat(int g, double pout, double pin){
  int i,j;
  Rcpp::NumericMatrix dm(g,g);
  for(i = 0; i < g; i++){
    for(j = 0; j < g; j++){
      if(i == j){
        dm(i,j) = pin;
        } else {
        dm(i,j) = pout;
        }
      }
    }
  return(dm);
  }


// get position in cumulative distribution
double prbs(std::vector<double> x, double p);

// determine the largest difference between 2 cumulative distributions
double trm(std::vector<double> x, std::vector<double> y);

// random integer
// inline int rndint(int n){
// return std::rand() % n;
//  }

inline int rndint(int n){
  return (Rcpp::sample(n, 1))[0] - 1;
}

// allocate
//inline Rcpp::NumericVector allct(int n, int g, int mx){
//  int i = 0, grs;
// int total = n - g;
// if(mx < 1) mx = 1;
//  Rcpp::NumericVector sizes(g);
//  while( i < (g-1) ){
//    if(total > 0){
//      if(total >= mx) {
//        grs = rndint(mx-1);
//        sizes[i] = grs + 1;
//        } else {
//        grs = rndint(total-1);
//        sizes[i] = grs + 1;
//        }
//      total = total - grs;
//      } else {
//      sizes[i] = 1;
//      }
//    i++;
//    }
//  if(total > 0){
//    sizes[i] = total + 1;
//    } else {
//    sizes[i] = 1;
//    }
//  return sizes;
//  }

// draw random integer between 0 and (n-1)
int rint(int n);

// draw random from uniform [0,1]
double runi();

// sampl according to vector of values (interpreted as proportional to probability)
int smpl(std::vector<double> ps);

// unique
std::vector<int> unique_int(std::vector<int> v);

// shuffle
std::vector<int> shuffle(std::vector<int> & v);

//////////////////////////////////////////////////////////////////////////////
//
//          GRAPH UTILITIES
//
//////////////////////////////////////////////////////////////////////////////


// join two integer vectors
inline std::vector<int> join_int(std::vector<int> a, std::vector<int> b){
  std::vector<int> ab;
  ab.reserve( a.size() + b.size() ); // preallocate memory
  ab.insert( ab.end(), a.begin(), a.end() );
  ab.insert( ab.end(), b.begin(), b.end() );
  return ab;
  }

// reduce integer vector to unique items
//std::vector<int> unique_int(std::vector<int> v)


//////////////////////////////////////////////////////////////////////////////
//
//          SEARCH UTILITIES
//
//////////////////////////////////////////////////////////////////////////////

// test if element is in set
inline bool inset(int el, std::vector<int> set){
  return std::find(set.begin(), set.end(), el) != set.end();
  }

// test if element is in set
inline bool inset_str(std::string el, std::vector<std::string> set){
  return std::find(set.begin(), set.end(), el) != set.end();
  }

// test if element is in set
inline bool inset_strset(std::string el, std::set<std::string> set){
  return set.find(el) != set.end();
  }

// test if edge is in pair
inline bool inmap_2int(int el, std::map<int, int> m){
  return m.find(el) != m.end();
  }

//////////////////////////////////////////////////////////////////////////////
//
//          GONI UTILITIES
//
//////////////////////////////////////////////////////////////////////////////

// n over k
double noverk(int n, int k);

#endif // __UTILITIES__
