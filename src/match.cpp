#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double m_nn(double ele, std::vector<double> set, std::vector<double> rep) {
  std::vector<double>::iterator it;
  it = std::find(set.begin(), set.end(), ele);
  int i = std::distance(set.begin(), it);
  if(i == rep.size()) return NA_REAL;
  return rep[i];
  }

// [[Rcpp::export]]
std::string m_nc(double ele, std::vector<double> set, std::vector<std::string> rep) {
  std::vector<double>::iterator it;
  it = std::find(set.begin(), set.end(), ele);
  int i = std::distance(set.begin(), it);
  if(i == rep.size()) return "NA";
  return rep[i];
}

// [[Rcpp::export]]
double m_cn(std::string ele, std::vector<std::string> set, std::vector<double> rep) {
  std::vector<std::string>::iterator it;
  it = std::find(set.begin(), set.end(), ele);
  int i = std::distance(set.begin(), it);
  if(i == rep.size()) return NA_REAL;
  return rep[i];
}

// [[Rcpp::export]]
std::string m_cc(std::string ele, std::vector<std::string> set, std::vector<std::string> rep) {
  std::vector<std::string>::iterator it;
  it = std::find(set.begin(), set.end(), ele);
  int i = std::distance(set.begin(), it);
  if(i == rep.size()) return "NA";
  return rep[i];
}

// [[Rcpp::export]]
std::vector<double> match_nn(std::vector<double> elems, std::vector<double> set, std::vector<double> rep) {
  int n = elems.size();
  std::vector<double> reps(n);
  for(int i = 0; i < n; ++i){
    reps[i] = m_nn(elems[i],set,rep);
    }
  return reps;
  }

// [[Rcpp::export]]
std::vector<std::string> match_nc(std::vector<double> elems, std::vector<double> set, std::vector<std::string> rep) {
  int n = elems.size();
  std::vector<std::string> reps(n);
  for(int i = 0; i < n; ++i){
    reps[i] = m_nc(elems[i],set,rep);
  }
  return reps;
}

// [[Rcpp::export]]
std::vector<double> match_cn(std::vector<std::string> elems, std::vector<std::string> set, std::vector<double> rep) {
  int n = elems.size();
  std::vector<double> reps(n);
  for(int i = 0; i < n; ++i){
    reps[i] = m_cn(elems[i],set,rep);
  }
  return reps;
}

// [[Rcpp::export]]
std::vector<std::string> match_cc(std::vector<std::string> elems, std::vector<std::string> set, std::vector<std::string> rep) {
  int n = elems.size();
  std::vector<std::string> reps(n);
  for(int i = 0; i < n; ++i){
    reps[i] = m_cc(elems[i],set,rep);
  }
  return reps;
}